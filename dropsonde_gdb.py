# dropsonde_gdb.py -- scheme-boundary argument dumper. Runs INSIDE gdb:
#
#   DROPSONDE_CONFIG=<config.json> gdb --batch -x dropsonde_gdb.py cesm.exe
#
# (normally launched by ./dropsonde). Must stay compatible with gdb's
# embedded Python (>= 3.6, e.g. Derecho's 3.6.15) and use only the stdlib.
#
# Strategy: a breakpoint at every <scheme>_run entry. At entry we walk the
# frame's dummy arguments, probe each array's base address and per-subscript
# byte strides empirically (element-address arithmetic, so we never depend
# on gdb's descriptor internals), dump the bytes, and plant a
# FinishBreakpoint that re-reads the same addresses at scheme exit.
# Argument storage is caller-owned, so entry-time addresses remain valid at
# exit. The FinishBreakpoint's stop() does memory reads and file writes
# only (safe inside stop callbacks) and never halts the run.

import json
import os
import re
import struct
import traceback

import gdb

with open(os.environ["DROPSONDE_CONFIG"]) as _cfg_f:
    CFG = json.load(_cfg_f)
ROLE = CFG["role"]            # "cam" | "sima"
OUT = CFG["out_dir"]
SCHEMES = CFG["schemes"]      # unique scheme names, suite order
KILL_AFTER = CFG.get("kill_after_steps", 0)  # 0 = run to natural exit
# {scheme: portable_symbol}: schemes whose breakpoint is planted on a shared
# portable subroutine (called by both the CAM and CAM-SIMA wrappers) instead
# of the SIMA-only <scheme>_run CCPP wrapper. See parse_portable_map in the
# dropsonde driver.
PORTABLE = CFG.get("portable", {}) or {}
# {derived_type_name_lower: [component paths]}: derived-type dummies whose
# type matches a key are expanded into per-component pseudo-args
# (e.g. elem%state%v) instead of being skipped. Authored as a JSON spec
# (--capture); needed for non-CCPP code like the SE dycore, where the model
# state travels in element_t arrays rather than plain-array arguments.
CAPTURE = CFG.get("capture", {}) or {}

# gdb (Derecho gfortran, calibrated 2026-06-12) reports Fortran array
# ranges outermost-first, i.e. REVERSED relative to subscript order, while
# our stride probes are per subscript position. Each array's extents/stride
# pairing is additionally cross-checked against gdb's type sizeof in
# array_plan() and flipped automatically when the check is conclusive;
# DIM_ORDER is only the fallback for ambiguous (square/strided) cases.
DIM_ORDER = "reversed"

MANIFEST = {
    "role": ROLE,
    "dim_order": DIM_ORDER,
    "breakpoints": {},        # scheme -> resolved spec | "missing"
    "constituents": None,     # ordered names (CAM short / CCPP standard)
    "hits": [],               # one record per scheme entry, execution order
    "notes": [],
}
FILE_IDX = [0]
HIT_COUNT = {}

SCALAR_FMT = {"f8": "d", "f4": "f", "i8": "q", "i4": "i", "i2": "h", "i1": "b"}

# gdb's prologue skip lands the scheme-entry breakpoint INSIDE gfortran's
# argument-descriptor setup, where assumed-shape dummies aren't readable yet
# ("Location address is not set") and single-stepping does not reliably escape
# (stepping bounces among the subroutine-declaration line numbers for the
# ~66-dummy park_macrophysics_run). We instead advance to the first executable
# statement of the body, where the dummies are live. CCPP _run/_init
# subroutines open with errmsg=''/errflg=0, so the body line is found by
# scanning the source from the subroutine's own definition line for that
# idiom. Cached per scheme. (Subroutines large enough that gfortran -O0 drops
# the DWARF locations of their stack-passed dummies altogether additionally
# need the ABI fallback in handle_entry; the body is a valid vantage for it.)
_ERR_INIT_RE = re.compile(r"^\s*(errmsg|errflg|errcode|errstat)\s*=", re.I)
_BODY_LINE = {}  # scheme -> (filename, line) | None


def note(msg):
    MANIFEST["notes"].append(msg)
    gdb.write("dropsonde[{}]: {}\n".format(ROLE, msg))


# --------------------------------------------------------------------------
# compiler detection
#
# Most of the capture is compiler-agnostic: gdb's by-name element addressing
# (array_plan) and scalar reads work for gfortran and Intel alike. A few
# memory-layout-specific paths -- the manual descriptor / deferred-length
# character decode, and the stack-arg ABI fallback -- need to know which
# compiler emitted the binary. We read it from the DWARF producer string
# (e.g. "GNU Fortran ..." vs "Intel(R) Fortran ... Version 19.1..."), which
# gdb exposes per compilation unit as symtab.producer. Requires a selected
# frame (a symtab to read it from), so it is resolved lazily at the first
# scheme hit and cached.
# --------------------------------------------------------------------------

_COMPILER = [None]

# Optimized-binary (production, e.g. gnu -O2) support. gfortran's DWARF
# producer string carries the compile flags, so optimization is detected
# pre-run from the CU of the first scheme symbol that resolves; None means
# the probe could not decide (treated as not optimized). At -O2 the normal
# capture path is unusable: the prologue-skip breakpoint can land where the
# DWARF locations of descriptor-backed dummies are invalid (bounds read as
# garbage), and Fortran subscript evaluation (parse_and_eval("a(1,1)"))
# aborts gdb 16.2 outright with an internal error on the dynamic stride
# properties gfortran emits at -O2. So optimized binaries are captured at
# the raw function entry instruction instead, where the System V AMD64 ABI
# guarantees every by-reference argument pointer's location; DWARF is used
# only for static facts (arg names, declaration order, element types,
# assumed- vs explicit-shape).
OPTIMIZED = [None]
_OPT_FLAG_RE = re.compile(r"(?:^|\s)-O([123sz])(?:\s|$)")


def compiler():
    """'gnu' | 'intel' | 'nag' | 'unknown' for the binary under debug."""
    if _COMPILER[0] is not None:
        return _COMPILER[0]
    prod = ""
    try:
        prod = gdb.newest_frame().find_sal().symtab.producer or ""
    except Exception:
        pass
    low = prod.lower()
    if "intel" in low:
        c = "intel"
    elif "nag" in low:
        c = "nag"
    elif "gnu" in low or "gcc" in low:
        c = "gnu"
    else:
        c = "unknown"
    _COMPILER[0] = c
    note("compiler: {} (producer: {})".format(c, prod[:70] or "?"))
    return c


# --------------------------------------------------------------------------
# type and memory helpers
# --------------------------------------------------------------------------

def dtype_of(t):
    if t.code == gdb.TYPE_CODE_FLT:
        return "f{}".format(t.sizeof)
    if t.code in (gdb.TYPE_CODE_INT, gdb.TYPE_CODE_BOOL):
        return "i{}".format(t.sizeof)
    return None


def is_character(t):
    return "character" in str(t).lower()


def array_dims(t):
    """Nested array ranges as [(lo, hi), ...] plus the element type."""
    dims = []
    t = t.strip_typedefs()
    while t.code == gdb.TYPE_CODE_ARRAY:
        lo, hi = t.range()
        dims.append((int(lo), int(hi)))
        t = t.target().strip_typedefs()
    return dims, t


def elem_addr(name, subs):
    v = gdb.parse_and_eval("{}({})".format(
        name, ",".join(str(s) for s in subs)))
    if v.address is None:
        raise gdb.error("no address for element of " + name)
    return int(v.address)


def array_plan(name, val):
    """Build a capture plan {addr, elsize, dtype, los, extents, strides}.

    Strides are bytes per increment of each subscript (source order),
    probed via element addresses so strided actual arguments (array
    slices) are handled correctly.
    """
    dims, base = array_dims(val.type)
    dt = dtype_of(base)
    if dt is None:
        return None, "unsupported element type: {}".format(base)
    if DIM_ORDER == "reversed":
        dims = dims[::-1]
    los = [d[0] for d in dims]
    extents = [d[1] - d[0] + 1 for d in dims]
    elsize = base.sizeof
    sizeof = val.type.sizeof
    addr0 = elem_addr(name, los)
    strides = []
    for k in range(len(dims)):
        if extents[k] <= 1:
            strides.append(0)
            continue
        subs = list(los)
        subs[k] = los[k] + 1
        strides.append(elem_addr(name, subs) - addr0)
    if any(s < 0 for s in strides):
        return None, "negative stride (reversed slice?): {}".format(strides)

    # Cross-check the extents<->stride pairing against gdb's logical array
    # size and flip when (and only when) the flipped pairing matches.
    # Strided dummies legitimately fail both checks; they keep the
    # DIM_ORDER pairing, which is still self-consistent across both runs.
    def span(exts):
        return sum(s * (e - 1) for s, e in zip(strides, exts)) + elsize

    if len(dims) > 1 and span(extents) != sizeof:
        flipped = extents[::-1]
        if span(flipped) == sizeof and los == los[::-1]:
            extents = flipped
    return {"addr": addr0, "elsize": elsize, "dtype": dt, "los": los,
            "extents": extents, "strides": strides}, None


# A synthesized derived-type plan can span the whole element array (GBs)
# for a component that is a small fraction of it; batch the outermost
# (largest-stride) dimension so no single gdb read exceeds this.
CHUNK_SPAN = 64 << 20


def read_array(plan):
    """Read the logical array contents as contiguous subscript-major bytes."""
    if plan.get("addrs") is not None:
        # per-element plans (allocatable derived-type component): same
        # inner shape at per-element heap addresses, outermost dim last
        inner = plan["inner"]
        return b"".join(read_array(dict(inner, addr=a))
                        for a in plan["addrs"])
    extents = plan["extents"]
    strides = plan["strides"]
    elsize = plan["elsize"]
    total = 1
    for e in extents:
        total *= e
    if total == 0:
        return b""
    span = sum(s * (e - 1) for s, e in zip(strides, extents)) + elsize
    if (span > CHUNK_SPAN and len(extents) > 1 and strides[-1] > 0 and
            strides[-1] == max(strides)):
        nb = max(1, CHUNK_SPAN // strides[-1])
        out = []
        for j0 in range(0, extents[-1], nb):
            cnt = min(nb, extents[-1] - j0)
            sub = dict(plan, addr=plan["addr"] + j0 * strides[-1],
                       extents=extents[:-1] + [cnt])
            out.append(read_array(sub))
        return b"".join(out)
    raw = bytes(gdb.selected_inferior().read_memory(plan["addr"], span))

    contig = True
    expect = elsize
    for s, e in zip(strides, extents):
        if e > 1 and s != expect:
            contig = False
            break
        expect *= e
    if contig:
        return raw[:total * elsize]

    out = bytearray(total * elsize)
    pos = 0
    for li in range(total):
        rem = li
        off = 0
        for s, e in zip(strides, extents):
            off += (rem % e) * s
            rem //= e
        out[pos:pos + elsize] = raw[off:off + elsize]
        pos += elsize
    return bytes(out)


def read_scalar(addr, dt):
    raw = bytes(gdb.selected_inferior().read_memory(addr, int(dt[1:])))
    return struct.unpack("=" + SCALAR_FMT[dt], raw)[0]


def write_blob(scheme, hit, arg, phase, data):
    FILE_IDX[0] += 1
    fn = "{:05d}_{}_h{}_{}_{}.bin".format(
        FILE_IDX[0], scheme, hit, arg.replace("%", "."), phase)
    with open(os.path.join(OUT, fn), "wb") as f:
        f.write(data)
    return fn


# --------------------------------------------------------------------------
# derived-type dummy expansion (--capture spec)
#
# Non-CCPP code like the SE dycore passes its state in derived types
# (element_t arrays: elem(:)%state%v ...) rather than plain-array arguments.
# A capture spec maps a derived-type NAME to the component paths worth
# comparing; matching dummies are expanded into per-component pseudo-args
# ("elem%state%v") whose records look exactly like plain array/scalar args,
# so the exit re-read and the differ are unchanged. Components at a fixed
# offset in the element (compile-time shape) become one strided plan over
# all elements; allocatable components live at per-element heap addresses
# and get a per-element address list (read_array's "addrs" path). Addresses
# are probed empirically per hit, like every other capture.
# --------------------------------------------------------------------------

def _capturable_struct(t):
    """Lowercase derived-type name if `t` is a struct, pointer-to-struct, or
    array-of-struct dummy; None otherwise."""
    t = t.strip_typedefs()
    if t.code == gdb.TYPE_CODE_PTR:
        t = t.target().strip_typedefs()
    if t.code == gdb.TYPE_CODE_ARRAY:
        t = _element_type(t)
    if t.code not in (gdb.TYPE_CODE_STRUCT, gdb.TYPE_CODE_UNION):
        return None
    name = t.name or t.tag or ""
    name = name.split("::")[-1].strip().lower()
    return name or None


def _expr_plan(expr):
    """('array', plan) | ('scalar', (addr, dtype)) | (None, None) plus an
    error string, for a component-path gdb expression like 'tl%n0' or
    'elem(3)%state%v'."""
    v = gdb.parse_and_eval(expr)
    t = v.type.strip_typedefs()
    if is_character(t):
        return None, None, "character component not captured"
    if t.code == gdb.TYPE_CODE_ARRAY:
        plan, err = array_plan(expr, v)
        if plan is None:
            return None, None, err
        return "array", plan, None
    if t.code in (gdb.TYPE_CODE_STRUCT, gdb.TYPE_CODE_UNION,
                  gdb.TYPE_CODE_PTR):
        return None, None, ("component is a derived type; extend the "
                            "capture path into its numeric members")
    dt = dtype_of(t)
    if dt is None:
        return None, None, "unsupported component type: {}".format(t)
    if v.address is None:
        return None, None, "component has no address"
    return "scalar", (int(v.address), dt), None


def _as_rank0_plan(payload):
    """A scalar component as a rank-0 plan, so array-of-struct expansion
    can treat scalar and array components uniformly (the element dimension
    is appended either way)."""
    addr, dt = payload
    return {"addr": addr, "elsize": int(dt[1:]), "dtype": dt,
            "los": [], "extents": [], "strides": []}


def _expand_struct_scalar(scheme, hit, name, paths, rec):
    """Pseudo-args for the components of a (scalar) derived-type dummy."""
    for path in paths:
        key = "{}%{}".format(name, path)
        info = {"kind": "skipped"}
        rec["args"][key] = info
        try:
            kind, payload, err = _expr_plan(key)
            if kind == "array":
                _record_array(info, scheme, hit, key, payload)
            elif kind == "scalar":
                addr, dt = payload
                info["kind"] = "scalar"
                info["dtype"] = dt
                info["addr"] = addr
                info["entry_value"] = read_scalar(addr, dt)
            else:
                info["why"] = err
        except Exception as exc:
            info["kind"] = "error"
            info["why"] = str(exc)


def _expand_struct_array(scheme, hit, name, val, paths, rec):
    """Pseudo-args for the components of an array-of-struct dummy (e.g.
    elem(:)). Each component is captured across ALL elements, the element
    index appended as the outermost (last, slowest-varying) dimension."""
    dims, _base = array_dims(val.type)
    if len(dims) != 1:
        raise gdb.error("rank-{} derived-type array not supported".format(
            len(dims)))
    lo, hi = dims[0]
    n = hi - lo + 1
    if n <= 0:
        for path in paths:
            rec["args"]["{}%{}".format(name, path)] = {
                "kind": "skipped", "why": "empty derived-type array"}
        return
    estride = None
    if n > 1:
        a0 = gdb.parse_and_eval("{}({})".format(name, lo)).address
        a1 = gdb.parse_and_eval("{}({})".format(name, lo + 1)).address
        if a0 is not None and a1 is not None:
            estride = int(a1) - int(a0)
    for path in paths:
        key = "{}%{}".format(name, path)
        info = {"kind": "skipped"}
        rec["args"][key] = info
        try:
            kind, payload, err = _expr_plan(
                "{}({})%{}".format(name, lo, path))
            if kind is None:
                info["why"] = err
                continue
            p0 = payload if kind == "array" else _as_rank0_plan(payload)
            if n > 1:
                # shape must be uniform across elements (same allocation
                # code ran for each); spot-check the last element
                kind_l, payload_l, err_l = _expr_plan(
                    "{}({})%{}".format(name, hi, path))
                if kind_l != kind:
                    raise gdb.error("component kind varies across "
                                    "elements: {}".format(err_l or kind_l))
                pl = (payload_l if kind == "array"
                      else _as_rank0_plan(payload_l))
                if (pl["extents"] != p0["extents"] or
                        pl["strides"] != p0["strides"]):
                    raise gdb.error(
                        "component shape varies across elements: "
                        "{}/{} vs {}/{}".format(
                            p0["extents"], p0["strides"],
                            pl["extents"], pl["strides"]))
            delta = None
            if n > 1:
                kind1, payload1, _err1 = _expr_plan(
                    "{}({})%{}".format(name, lo + 1, path))
                if kind1 == kind:
                    p1 = (payload1 if kind == "array"
                          else _as_rank0_plan(payload1))
                    delta = p1["addr"] - p0["addr"]
            if n == 1 or (delta is not None and delta == estride and
                          delta > 0):
                # fixed-offset component: one strided plan covers all
                # elements
                plan = {"addr": p0["addr"], "elsize": p0["elsize"],
                        "dtype": p0["dtype"], "los": p0["los"] + [lo],
                        "extents": p0["extents"] + [n],
                        "strides": p0["strides"] + [estride or 0]}
            else:
                # allocatable/pointer component: per-element heap
                # addresses; probe each element's base once
                addrs = [p0["addr"]]
                for i in range(lo + 1, hi + 1):
                    e = "{}({})%{}".format(name, i, path)
                    if kind == "array":
                        addrs.append(elem_addr(e, p0["los"]))
                    else:
                        addrs.append(
                            int(gdb.parse_and_eval(e).address))
                plan = {"elsize": p0["elsize"], "dtype": p0["dtype"],
                        "los": p0["los"] + [lo],
                        "extents": p0["extents"] + [n],
                        "strides": p0["strides"] + [0],
                        "inner": p0, "addrs": addrs}
            _record_array(info, scheme, hit, key, plan)
        except Exception as exc:
            info["kind"] = "error"
            info["why"] = str(exc)


# --------------------------------------------------------------------------
# entry / exit capture
# --------------------------------------------------------------------------

def arg_symbols(frame):
    blk = frame.block()
    while blk is not None and blk.function is None:
        blk = blk.superblock
    if blk is None:
        return []
    return [s for s in blk if s.is_argument]


_END_SUB_RE = re.compile(r"^\s*end\s+subroutine\b", re.I)

# Heuristic classification of a source line as executable, used only to
# sanity-check the errmsg-sentinel assumption below. Declarations carry '::'
# or open with a declaration keyword; executable statements open with a
# control keyword or look like an assignment.
_DECL_KEYWORD_RE = re.compile(
    r"^\s*(real|integer|logical|character|complex|double\s+precision|"
    r"type\s*[(,:]|class\s*[(,:]|use\b|implicit\b|save\b|parameter\b|"
    r"dimension\b|external\b|intrinsic\b|pointer\b|allocatable\b|target\b|"
    r"optional\b|public\b|private\b|protected\b|procedure\b|import\b|"
    r"interface\b|abstract\b|end\b|data\b|equivalence\b|common\b|"
    r"namelist\b|format\b|entry\b|include\b|sequence\b|contiguous\b|"
    r"volatile\b|value\b|enum\b|enumerator\b)", re.I)
_EXEC_LOOK_RE = re.compile(
    r"^\s*(call\b|if\s*\(|do\b|select\s+case|where\s*\(|associate\s*\(|"
    r"forall\s*\(|block\b|[a-z]\w*(\s*\([^)]*\))?\s*=[^=>])", re.I)


def _exec_lines_before(lines, start, stop):
    """1-based line numbers in [start, stop) that look like executable
    statements. The drive loop assumes the errmsg=''/errflg=0 sentinel is the
    FIRST executable statement of the body; anything executable before it
    means entry capture happens mid-body (intent(inout) inputs possibly
    already mutated) -- a convention, not a guarantee, for portable (non-CCPP
    -wrapper) subroutines, so it is checked and warned about rather than
    trusted. Skips blanks, comments, preprocessor lines, continuation lines,
    and declarations."""
    suspects = []
    cont = False
    for idx in range(start - 1, stop - 1):
        code = lines[idx].split("!", 1)[0].strip()
        was_cont = cont
        cont = code.endswith("&")
        if was_cont or not code or code.startswith("#"):
            continue
        if "::" in code or _DECL_KEYWORD_RE.match(code):
            continue
        if _EXEC_LOOK_RE.match(code):
            suspects.append(idx + 1)
    return suspects


def _first_body_line(scheme, frame):
    """(filename, line) of the first executable statement of the scheme's
    body, or None. Dummy arguments are only reliably readable once execution
    is in the body (past the compiler-generated descriptor setup), so the
    drive loop advances here before capturing. CCPP subroutines open with
    errmsg=''/errflg=0; we scan the source from the subroutine's OWN
    definition line (so a sibling subroutine earlier in the same file -- e.g.
    the _init above _run -- is not matched). Stop at 'end subroutine' so
    schemes that lack the errmsg idiom don't match a sibling's body line
    (which would cause advance to overshoot into the caller frame).
    Result is cached per scheme."""
    if scheme in _BODY_LINE:
        return _BODY_LINE[scheme]
    result = None
    try:
        start = frame.function().line
        fname = frame.find_sal().symtab.fullname()
        with open(fname) as fh:
            lines = fh.readlines()
        end = len(lines) + 1
        for idx in range(start - 1, len(lines)):
            if _END_SUB_RE.match(lines[idx]):
                end = idx + 1
                break
            if _ERR_INIT_RE.match(lines[idx]):
                result = (fname, idx + 1)
                break
        if result is not None:
            suspects = _exec_lines_before(lines, start, result[1])
            if suspects:
                note("{}: WARNING: {} executable-looking line(s) precede "
                     "the errmsg/errflg init at {}:{} (first: line {}); "
                     "entry capture happens there and may reflect mid-body "
                     "state".format(scheme, len(suspects),
                                    os.path.basename(fname), result[1],
                                    suspects[0]))
        else:
            # Non-CCPP subroutines (e.g. the SE dycore) have no errmsg
            # idiom; the first executable-looking statement is an equally
            # safe advance target -- 'advance' stops BEFORE the line runs,
            # so the dummies are live and still unmodified.
            suspects = _exec_lines_before(lines, start, end)
            if suspects:
                result = (fname, suspects[0])
                note("{}: no errmsg idiom; using first executable "
                     "statement".format(scheme))
    except Exception as exc:
        note("{}: body-line scan failed: {}".format(scheme, exc))
    if result is not None:
        note("{}: capturing at body {}:{}".format(
            scheme, os.path.basename(result[0]), result[1]))
    _BODY_LINE[scheme] = result
    return result


# --------------------------------------------------------------------------
# ABI-level argument capture (gfortran / System V AMD64 fallback)
#
# A few CCPP schemes are large enough that gfortran -O0 emits an EMPTY
# DW_AT_location for every stack-passed dummy: `info address numliq` reports
# "optimized out" and gdb cannot read the dummy by name ANYWHERE in the
# routine -- not a prologue or stopping-point issue (park_macrophysics_run:
# 66 dummies, ~1500 lines, a ~30 KB frame from its (ncol,pver) automatic
# arrays). The data is still on the stack where the ABI put it, so we recover
# it directly. Every Fortran dummy is passed by reference (one 8-byte
# pointer); System V AMD64 puts the first 6 in registers and the rest as
# consecutive 8-byte slots from rbp+16. So dummy i (1-based, declaration
# order) with i > 6 is the pointer at [rbp + 16 + 8*(i-7)] -- a descriptor
# pointer for assumed-shape arrays, a value pointer for scalars. Confirmed on
# gcc 14.3.0. Assumes all-by-reference dummies (the CCPP norm); by-value float
# arguments would shift the integer-register accounting and are not handled.
# --------------------------------------------------------------------------

import platform as _platform
_ARCH = _platform.machine()

N_INT_REG = 6          # System V AMD64 integer/pointer argument registers
STACK_ARG0 = 16        # first stack arg at rbp+16 (saved rbp, then return addr)

# gfortran array descriptor, LP64, gcc 9-14 (byte offsets):
_D_BASE = 0            # void *base_addr (points at the element at the lbounds)
_D_ELEM_LEN = 16       # size_t elem_len
_D_RANK = 28           # signed char rank
_D_DIM0 = 40           # dim[]: {index_type stride (ELEMENTS), lbound, ubound}
_D_DIM = 24            #   3 * 8 bytes per dimension


def _read_u64(addr):
    return struct.unpack("=Q", bytes(
        gdb.selected_inferior().read_memory(addr, 8)))[0]


def _read_i64(addr):
    return struct.unpack("=q", bytes(
        gdb.selected_inferior().read_memory(addr, 8)))[0]


def _frame_base(frame):
    return int(frame.read_register("rbp")) & 0xFFFFFFFFFFFFFFFF


def _element_type(t):
    """Strip array codes to the element type without touching dynamic ranges."""
    t = t.strip_typedefs()
    while t.code == gdb.TYPE_CODE_ARRAY:
        t = t.target().strip_typedefs()
    return t


def _descriptor_plan(desc_addr, elem_type):
    """An array_plan()-shaped dict built from a raw gfortran descriptor."""
    dt = dtype_of(elem_type)
    if dt is None:
        return None, "unsupported element type: {}".format(elem_type)
    rank = struct.unpack("=b", bytes(
        gdb.selected_inferior().read_memory(desc_addr + _D_RANK, 1)))[0]
    if rank < 1 or rank > 7:
        return None, "implausible descriptor rank: {}".format(rank)
    base = _read_u64(desc_addr + _D_BASE)
    if base == 0:
        return None, "null descriptor base_addr"
    elem_len = _read_i64(desc_addr + _D_ELEM_LEN)
    los, extents, strides = [], [], []
    for d in range(rank):
        b = desc_addr + _D_DIM0 + _D_DIM * d
        st = _read_i64(b)
        lb = _read_i64(b + 8)
        ub = _read_i64(b + 16)
        los.append(lb)
        extents.append(ub - lb + 1)
        strides.append(st * elem_len)
    if any(e < 0 for e in extents):
        return None, "negative extent in descriptor: {}".format(extents)
    return {"addr": base, "elsize": elem_type.sizeof, "dtype": dt,
            "los": los, "extents": extents, "strides": strides}, None


def _capture_abi(scheme, hit, name, sym, idx, frame, info):
    """Recover a stack-passed dummy whose DWARF location gfortran dropped,
    reading it straight from the ABI argument slot. Fills `info` exactly like
    handle_entry's normal path, so the exit re-read and differ are unchanged."""
    if _ARCH not in ("x86_64", "AMD64"):
        info["kind"] = "error"
        info["why"] = ("ABI fallback only supports x86-64 (running on {})"
                       .format(_ARCH))
        return
    # The descriptor decode below is the gfortran array-descriptor layout.
    # Intel uses a different dope-vector format and (so far) has not dropped
    # any dummy locations, so this path is gfortran-only; bail with a clear
    # reason rather than misread an Intel descriptor as a gfortran one.
    if compiler() not in ("gnu", "unknown"):
        info["kind"] = "error"
        info["why"] = "ABI fallback unsupported under {}".format(compiler())
        return
    if name.startswith("_"):              # gfortran hidden character-length arg
        info["why"] = "hidden character-length arg (by value)"
        return
    st = sym.type.strip_typedefs()
    if is_character(st):
        info["why"] = "character arg via ABI not captured"
        return
    slot = _frame_base(frame) + STACK_ARG0 + 8 * (idx - N_INT_REG - 1)
    ptr = _read_u64(slot)                 # by reference: descriptor/value pointer
    info["abi"] = True
    if st.code == gdb.TYPE_CODE_ARRAY:
        elem = _element_type(st)
        if elem.code in (gdb.TYPE_CODE_STRUCT, gdb.TYPE_CODE_UNION):
            info["why"] = "derived-type array via ABI not captured"
            return
        # Known gap: an EXPLICIT-shape stack-passed dummy (b(ncol,pver)) is a
        # raw data pointer, not a descriptor, and gdb's type info here does
        # not distinguish the two. Decoding data bytes as a descriptor is
        # expected to fail _descriptor_plan's plausibility gates (rank 1-7 at
        # +28, non-null base, non-negative extents) -- probabilistic, not
        # principled. If a capture from this path ever looks wrong, check
        # whether the dummy is explicit-shape.
        plan, err = _descriptor_plan(ptr, elem)
        if plan is None:
            info["kind"] = "error"
            info["why"] = "abi: " + err
            return
        info["kind"] = "array"
        info["dtype"] = plan["dtype"]
        info["los"] = plan["los"]
        info["extents"] = plan["extents"]
        info["strides"] = plan["strides"]
        info["plan"] = plan
        info["entry_file"] = write_blob(
            scheme, hit, name, "in", read_array(plan))
    elif st.code in (gdb.TYPE_CODE_STRUCT, gdb.TYPE_CODE_UNION):
        info["why"] = "derived-type arg via ABI not captured"
    else:
        dt = dtype_of(st)
        if dt is None:
            info["kind"] = "error"
            info["why"] = "abi: unsupported scalar type: {}".format(st)
            return
        info["kind"] = "scalar"
        info["dtype"] = dt
        info["addr"] = ptr
        info["entry_value"] = read_scalar(ptr, dt)


# --------------------------------------------------------------------------
# optimized-binary capture (raw-entry ABI reads; see OPTIMIZED above)
# --------------------------------------------------------------------------

ARG_REGS = ("rdi", "rsi", "rdx", "rcx", "r8", "r9")


def _cu_producer(spec):
    """DWARF producer string of the compilation unit defining `spec`, or
    None. Works pre-run off the static symbol table."""
    try:
        _unparsed, sals = gdb.decode_line(spec)
        return sals[0].symtab.producer or None
    except Exception:
        return None


def probe_optimized(specs_per_scheme):
    """Set OPTIMIZED[0] from the first scheme spec whose CU producer is
    readable. gfortran producers carry the compile flags (e.g. 'GNU
    Fortran2008 14.3.0 ... -O2 -ffp-contract=off -g')."""
    for specs in specs_per_scheme:
        for spec in specs:
            if callable(spec):
                continue        # _qualified sweep: too costly for a probe
            prod = _cu_producer(spec)
            if prod is None:
                continue
            m = _OPT_FLAG_RE.search(prod)
            OPTIMIZED[0] = bool(m)
            note("optimization probe: {} (producer: {})".format(
                "-O" + m.group(1) if m else "not optimized", prod[:90]))
            return
    note("optimization probe: no scheme CU producer readable; assuming "
         "not optimized")


def _entry_addr(spec):
    """Raw entry-point address of the function `spec` resolves to (the
    function block's start), or None. 'break <spec>' would stop at the
    post-prologue-skip address, where -O2 DWARF locations may not be valid
    yet; the ABI argument slots are only guaranteed at the first
    instruction."""
    try:
        _unparsed, sals = gdb.decode_line(spec)
    except gdb.error:
        return None
    if not sals:
        return None
    if len(sals) > 1:
        note("WARNING: {} resolves to {} locations; using the first".format(
            spec, len(sals)))
    try:
        blk = gdb.block_for_pc(sals[0].pc)
    except Exception:
        return None
    while blk is not None and blk.function is None:
        blk = blk.superblock
    if blk is None:
        return None
    return int(blk.start)


def _make_entry_bp_optimized(specs, scheme, target):
    """EntryBP at the raw entry address of the first spec that resolves."""
    for spec in specs:
        if callable(spec):
            spec = spec()
        if not spec:
            continue
        addr = _entry_addr(spec)
        if addr is None:
            continue
        return EntryBP("*{:#x}".format(addr), scheme, target)
    return None


def _abi_arg_ptr(frame, idx):
    """Pointer argument `idx` (1-based, declaration order) read at the
    function's FIRST instruction, where the System V AMD64 locations are
    unconditionally valid: args 1-6 in registers, 7+ in consecutive 8-byte
    stack slots above the return address at [rsp]."""
    if idx <= N_INT_REG:
        return int(frame.read_register(ARG_REGS[idx - 1])) \
            & 0xFFFFFFFFFFFFFFFF
    rsp = int(frame.read_register("rsp")) & 0xFFFFFFFFFFFFFFFF
    return _read_u64(rsp + 8 * (idx - N_INT_REG))


def _record_array(info, scheme, hit, name, plan):
    """Fill an arg record from a capture plan and dump the entry bytes."""
    info["kind"] = "array"
    info["dtype"] = plan["dtype"]
    info["los"] = plan["los"]
    info["extents"] = plan["extents"]
    info["strides"] = plan["strides"]
    info["plan"] = plan
    info["entry_file"] = write_blob(scheme, hit, name, "in", read_array(plan))


def _static_shape_plan(ptr, st):
    """Capture plan for a compile-time-constant explicit-shape dummy: a raw
    contiguous data pointer, shape fully known from the static DWARF type."""
    dims, base = array_dims(st)
    dt = dtype_of(base)
    if dt is None:
        return None, "unsupported element type: {}".format(base)
    if DIM_ORDER == "reversed":
        dims = dims[::-1]
    los = [d[0] for d in dims]
    extents = [d[1] - d[0] + 1 for d in dims]
    if any(e <= 0 for e in extents):
        return None, "implausible static extents: {}".format(extents)
    strides = []
    acc = base.sizeof
    for e in extents:
        strides.append(acc)
        acc *= e
    return {"addr": ptr, "elsize": base.sizeof, "dtype": dt, "los": los,
            "extents": extents, "strides": strides}, None


def handle_entry_optimized(scheme, frame):
    """Entry capture for an optimized binary, stopped at the function's
    first instruction. Mirrors handle_entry's record shape exactly (the
    exit re-read and the differ are unchanged), but never evaluates a
    dummy's value or subscripts through gdb -- every dynamic fact comes
    from the ABI argument pointers and the gfortran descriptors they
    reference."""
    hit = HIT_COUNT.get(scheme, 0)
    HIT_COUNT[scheme] = hit + 1
    try:
        caller = frame.older().name() or "?"
    except Exception:
        caller = "?"
    rec = {"scheme": scheme, "hit": hit, "step": CURRENT_STEP[0],
           "caller": caller, "args": {}}
    frame.select()
    for idx, sym in enumerate(arg_symbols(frame), start=1):
        name = sym.name
        info = {"kind": "skipped", "abi": True}
        rec["args"][name] = info
        try:
            if name.startswith("_"):
                info["why"] = "hidden character-length arg (by value)"
                continue
            st = sym.type.strip_typedefs()
            ptr = _abi_arg_ptr(frame, idx)
            if ptr == 0:
                info["why"] = "null argument pointer (absent optional?)"
                continue
            if is_character(st):
                info["why"] = "character arg not captured in optimized mode"
                continue
            if st.code == gdb.TYPE_CODE_ARRAY:
                elem = _element_type(st)
                if elem.code in (gdb.TYPE_CODE_STRUCT, gdb.TYPE_CODE_UNION):
                    info["why"] = "derived-type array not captured"
                    continue
                if ":" in str(sym.type):
                    # assumed shape: the argument is a descriptor pointer
                    plan, err = _descriptor_plan(ptr, elem)
                elif not getattr(st, "dynamic", False):
                    plan, err = _static_shape_plan(ptr, st)
                else:
                    info["why"] = ("explicit-shape dummy with runtime "
                                   "bounds not captured in optimized mode")
                    continue
                if plan is None:
                    info["kind"] = "error"
                    info["why"] = err
                    continue
                _record_array(info, scheme, hit, name, plan)
            elif st.code in (gdb.TYPE_CODE_STRUCT, gdb.TYPE_CODE_UNION):
                info["why"] = "derived-type arg not captured"
            else:
                dt = dtype_of(st)
                if dt is None:
                    info["why"] = "unsupported type: {}".format(st)
                    continue
                info["kind"] = "scalar"
                info["dtype"] = dt
                info["addr"] = ptr
                info["entry_value"] = read_scalar(ptr, dt)
        except Exception as exc:
            info["kind"] = "error"
            info["why"] = str(exc)
    MANIFEST["hits"].append(rec)
    return rec


def handle_entry(scheme, frame):
    hit = HIT_COUNT.get(scheme, 0)
    HIT_COUNT[scheme] = hit + 1
    try:
        caller = frame.older().name() or "?"
    except Exception:
        caller = "?"
    rec = {"scheme": scheme, "hit": hit, "step": CURRENT_STEP[0],
           "caller": caller, "args": {}}
    frame.select()
    for idx, sym in enumerate(arg_symbols(frame), start=1):
        name = sym.name
        info = {"kind": "skipped"}
        rec["args"][name] = info
        try:
            val = sym.value(frame)
            t = val.type.strip_typedefs()
            stname = _capturable_struct(t)
            if stname is not None:
                paths = CAPTURE.get(stname)
                if not paths:
                    info["why"] = ("derived type {} (no capture spec "
                                   "entry)".format(stname))
                elif t.code == gdb.TYPE_CODE_ARRAY:
                    info["why"] = ("derived-type array {}: expanded into "
                                   "{} component pseudo-args".format(
                                       stname, len(paths)))
                    _expand_struct_array(scheme, hit, name, val, paths, rec)
                else:
                    info["why"] = ("derived type {}: expanded into {} "
                                   "component pseudo-args".format(
                                       stname, len(paths)))
                    _expand_struct_scalar(scheme, hit, name, paths, rec)
            elif t.code == gdb.TYPE_CODE_ARRAY and not is_character(t):
                plan, err = array_plan(name, val)
                if plan is None:
                    info["why"] = err
                    continue
                data = read_array(plan)
                info["kind"] = "array"
                info["dtype"] = plan["dtype"]
                info["los"] = plan["los"]
                info["extents"] = plan["extents"]
                info["strides"] = plan["strides"]
                info["plan"] = plan
                info["entry_file"] = write_blob(scheme, hit, name, "in", data)
            elif is_character(t):
                info["kind"] = "char"
                try:
                    info["entry_value"] = val.string().rstrip()
                except Exception:
                    if val.address is not None and t.sizeof > 0:
                        raw = bytes(gdb.selected_inferior().read_memory(
                            int(val.address), t.sizeof))
                        info["entry_value"] = raw.decode("latin-1").rstrip()
                if val.address is not None and t.sizeof > 0:
                    info["addr"] = int(val.address)
                    info["len"] = t.sizeof
            elif dtype_of(t) is not None:
                dt = dtype_of(t)
                info["kind"] = "scalar"
                info["dtype"] = dt
                info["entry_value"] = (float(val) if dt[0] == "f"
                                       else int(val))
                if val.address is not None:
                    info["addr"] = int(val.address)
            else:
                info["why"] = "unsupported type: {}".format(t)
        except Exception as exc:
            # gfortran -O0 drops the DWARF location of stack-passed dummies of
            # very large subroutines; recover those (declaration index > 6)
            # straight from the ABI argument slots.
            if idx > N_INT_REG:
                try:
                    _capture_abi(scheme, hit, name, sym, idx, frame, info)
                    continue
                except Exception as exc2:
                    info["kind"] = "error"
                    info["why"] = "abi fallback failed: {}".format(exc2)
                    continue
            info["kind"] = "error"
            info["why"] = str(exc)
    MANIFEST["hits"].append(rec)
    return rec


class ExitBP(gdb.FinishBreakpoint):
    """Re-reads entry-time addresses when the scheme returns."""

    def __init__(self, frame, rec):
        super(ExitBP, self).__init__(frame, internal=True)
        self.rec = rec

    def stop(self):
        rec = self.rec
        for name, info in rec["args"].items():
            try:
                kind = info.get("kind")
                if kind == "array":
                    data = read_array(info["plan"])
                    info["exit_file"] = write_blob(
                        rec["scheme"], rec["hit"], name, "out", data)
                elif kind == "scalar" and "addr" in info:
                    info["exit_value"] = read_scalar(
                        info["addr"], info["dtype"])
                elif kind == "char" and "addr" in info:
                    raw = bytes(gdb.selected_inferior().read_memory(
                        info["addr"], info["len"]))
                    info["exit_value"] = raw.decode("latin-1").rstrip()
            except Exception as exc:
                info["exit_error"] = str(exc)
        rec["complete"] = True
        return False

    def out_of_scope(self):
        note("finish breakpoint out of scope: {} h{}".format(
            self.rec["scheme"], self.rec["hit"]))


# --------------------------------------------------------------------------
# constituent name capture (once, at first scheme hit, i.e. post-init)
# --------------------------------------------------------------------------

def _string_at(expr):
    v = gdb.parse_and_eval(expr)
    try:
        s = v.string().strip()
        if s:
            return s
    except Exception:
        pass
    # Let gdb's own printer resolve descriptors it understands.
    out = gdb.execute("print {}".format(expr), to_string=True)
    m = re.search(r"'(.*)'", out, re.S)
    if m:
        return m.group(1).strip()
    # Deferred-length allocatable component: gdb sees only the data
    # pointer ("PTR TO -> character*0").  gfortran stores the length in a
    # hidden sibling member _<name>_length appended after the type's own
    # visible components (named in DWARF, readable like class %_data).
    if v.type.code == gdb.TYPE_CODE_PTR:
        ptr = int(v)
        if ptr == 0:
            return ""
        head, sep, last = expr.rpartition("%")
        if sep:
            try:
                length = int(gdb.parse_and_eval(
                    "{}%_{}_length".format(head, last)))
            except gdb.error:
                length = -1
            if 0 < length <= 1024:
                raw = bytes(gdb.selected_inferior().read_memory(ptr, length))
                return raw.decode("latin-1").strip()
        # last resort: printable run at the data pointer; "?" marks the
        # name as heuristic (it will show as unmatched in the report)
        raw = bytes(gdb.selected_inferior().read_memory(ptr, 256))
        m = re.match(rb"[A-Za-z][A-Za-z0-9_]*", raw)
        if m:
            return m.group(0).decode("latin-1") + "?"
    raise gdb.error("could not read string {}: {}".format(
        expr, out.strip()[:120]))


def _field_offset(t, name):
    """(byte offset, field type) of a named component of a derived type."""
    for f in t.strip_typedefs().fields():
        if f.name == name:
            return f.bitpos // 8, f.type
    raise gdb.error("no component {!r} in {}".format(name, t))


def _intel_const_names(prefix, n):
    """SIMA constituent standard names under Intel. gdb 8.2 misreads Intel's
    array descriptor for the module pointer-array const_props (it returns
    bogus bounds), so const_props(i)%... cannot be evaluated by name; walk it
    by hand instead.

    Intel array descriptor (LP64): base_addr at +0, element length at +8 (the
    only two fields we rely on -- not the rank/dim header, whose offsets are
    less stable). const_props is a rank-1 contiguous array of
    constituent_prop_ptr_t, so element i sits at base + i*elem_len and its
    %prop pointer is at the component's offset within that element. %prop
    points to a constituent-properties derived type whose var_std_name is
    character(len=:),allocatable, stored inline by Intel as
    {char *data; size_t len} -- data at the component's offset, length in the
    next 8-byte word. Component offsets come from gdb's type info."""
    desc_addr = int(gdb.parse_and_eval("&{}const_props".format(prefix)))
    base = _read_u64(desc_addr)
    elem_len = _read_i64(desc_addr + 8)
    if base == 0 or elem_len <= 0:
        raise gdb.error("implausible const_props descriptor: base={:#x} "
                        "elem_len={}".format(base, elem_len))
    # peel pointer/array layers to the constituent_prop_ptr_t element type
    et = gdb.parse_and_eval("{}const_props".format(prefix)).type.strip_typedefs()
    while et.code in (gdb.TYPE_CODE_PTR, gdb.TYPE_CODE_ARRAY):
        et = et.target().strip_typedefs()
    prop_off, prop_t = _field_offset(et, "prop")
    vsn_off, _vsn_t = _field_offset(prop_t.target(), "var_std_name")
    names = []
    for i in range(n):
        prop_ptr = _read_u64(base + i * elem_len + prop_off)
        if prop_ptr == 0:
            names.append("<null prop>")
            continue
        data_ptr = _read_u64(prop_ptr + vsn_off)
        length = _read_i64(prop_ptr + vsn_off + 8)
        if data_ptr == 0 or not (0 < length <= 1024):
            names.append("<unreadable: len={}>".format(length))
            continue
        raw = bytes(gdb.selected_inferior().read_memory(data_ptr, length))
        names.append(raw.decode("latin-1").strip())
    return names


def _static_addr(name):
    """Address of a (minimal/linker) symbol via 'info address' text --
    gdb 16.2's Python API has no lookup_minimal_symbol."""
    try:
        out = gdb.execute("info address {}".format(name), to_string=True)
    except gdb.error:
        return None
    m = re.search(r"0x[0-9a-fA-F]+", out)
    return int(m.group(0), 16) if m else None


def _gnu_const_names(prefix, n):
    """SIMA constituent standard names decoded from raw memory using the
    gfortran descriptor layout. Used for optimized binaries, where Fortran
    subscript evaluation (const_props(i)%...) can abort gdb 16.2 with an
    internal error on -O2 dynamic-stride DWARF; module data itself is
    unoptimized, so the raw layout is the same as debug builds'.

    The descriptor is the module variable itself, found by its LINKER name:
    gdb resolves &const_props through the DWARF data_location to the array
    DATA, not the descriptor. const_props is a rank-1 pointer array of
    constituent_prop_ptr_t; element i sits at descriptor base + i*stride,
    its %prop pointer at the component's offset within the element. %prop
    points to a properties type whose var_std_name is
    character(len=:),allocatable: a data pointer at the component offset,
    with the length in the hidden sibling member _var_std_name_length
    gfortran appends to the derived type."""
    desc_addr = _static_addr("__cam_constituents_MOD_const_props")
    if desc_addr is None:
        raise gdb.error("no __cam_constituents_MOD_const_props symbol")
    v = gdb.parse_and_eval("&{}const_props".format(prefix))
    base = _read_u64(desc_addr + _D_BASE)
    elem_len = _read_i64(desc_addr + _D_ELEM_LEN)
    if base == 0 or elem_len <= 0:
        raise gdb.error("implausible const_props descriptor: base={:#x} "
                        "elem_len={}".format(base, elem_len))
    stride = _read_i64(desc_addr + _D_DIM0) * elem_len
    # peel pointer/array layers to the constituent_prop_ptr_t element type
    et = v.type.strip_typedefs()
    while et.code in (gdb.TYPE_CODE_PTR, gdb.TYPE_CODE_ARRAY):
        et = et.target().strip_typedefs()
    prop_off, prop_t = _field_offset(et, "prop")
    vsn_off, _vsn_t = _field_offset(prop_t.target(), "var_std_name")
    try:
        len_off, _ = _field_offset(prop_t.target(), "_var_std_name_length")
    except gdb.error:
        len_off = None
    names = []
    for i in range(n):
        prop_ptr = _read_u64(base + i * stride + prop_off)
        if prop_ptr == 0:
            names.append("<null prop>")
            continue
        data_ptr = _read_u64(prop_ptr + vsn_off)
        if len_off is not None:
            length = _read_i64(prop_ptr + len_off)
        else:
            length = _read_i64(prop_ptr + vsn_off + 8)
        if data_ptr == 0 or not (0 < length <= 1024):
            names.append("<unreadable: len={}>".format(length))
            continue
        raw = bytes(gdb.selected_inferior().read_memory(data_ptr, length))
        names.append(raw.decode("latin-1").strip())
    return names


def dump_constituents():
    try:
        if ROLE == "cam":
            # character(len=16) :: cnst_name(pcnst) in module constituents
            arr = None
            for expr in ("constituents::cnst_name",
                         "__constituents_MOD_cnst_name", "cnst_name"):
                try:
                    arr = gdb.parse_and_eval(expr)
                    break
                except gdb.error:
                    continue
            if arr is None:
                raise gdb.error("cnst_name not found by any spelling")
            total = arr.type.sizeof
            # element sizeof gives the character length regardless of how
            # gdb orders the array-of-strings dimensions
            width = gdb.parse_and_eval(expr + "(1)").type.sizeof
            n = total // width
            raw = bytes(gdb.selected_inferior().read_memory(
                int(arr.address), total))
            names = [raw[i * width:(i + 1) * width].decode("latin-1").strip()
                     for i in range(n)]
        else:
            # cam_constituents module: const_props(:) + num_constituents.
            # Linker-name spellings can hit "unknown type" (gdb 16.2), so
            # try the Fortran module:: syntax first.
            n = None
            base = None
            for prefix in ("cam_constituents::", "__cam_constituents_MOD_",
                           ""):
                try:
                    n = int(gdb.parse_and_eval(prefix + "num_constituents"))
                    base = prefix
                    break
                except gdb.error:
                    continue
            if n is None:
                raise gdb.error(
                    "num_constituents not found by any spelling")
            if compiler() == "intel":
                # gdb 8.2 can't index Intel's pointer-array descriptor; decode
                # const_props from raw memory.
                names = _intel_const_names(base, n)
            elif OPTIMIZED[0]:
                # -O2 subscript evaluation can abort gdb; decode raw instead.
                names = _gnu_const_names(base, n)
            else:
                names = []
                for i in range(1, n + 1):
                    try:
                        names.append(_string_at(
                            "{}const_props({})%prop%var_std_name".format(
                                base, i)))
                    except Exception as exc:
                        names.append("<unreadable: {}>".format(exc))
        MANIFEST["constituents"] = names
        note("captured {} constituent names".format(len(names)))
    except Exception as exc:
        note("constituent name capture FAILED: {}\n{}".format(
            exc, traceback.format_exc()))


# --------------------------------------------------------------------------
# drive loop
# --------------------------------------------------------------------------

class EntryBP(gdb.Breakpoint):
    def __init__(self, spec, scheme, target):
        super(EntryBP, self).__init__(spec)
        self.scheme = scheme    # SDF scheme name (the capture/differ tag)
        self.target = target    # symbol the breakpoint is actually planted on


class StepBP(gdb.Breakpoint):
    """Sentinel on cam_run1: fires once at the start of every timestep."""
    pass


LAST_STOP = [None]
DONE = [False]
CURRENT_STEP = [0]


def _qualified(short):
    """The gdb-demangled, module-qualified spelling (module::short) of a
    Fortran procedure, found by listing the symbol table. gfortran emits a
    DWARF name equal to the bare 'short' so 'break short' resolves directly;
    Intel mangles to module_mp_short_, which gdb DEMANGLES to module::short --
    a valid breakpoint spec, but the bare name does not resolve. So we look up
    the qualified spelling by the short name. Compiler-agnostic (gdb shows the
    same module::proc form for both); returns None if absent or not
    module-scoped. Works pre-run off the static symbol table."""
    try:
        out = gdb.execute("info functions \\b{}\\b".format(short),
                          to_string=True)
    except gdb.error:
        return None
    m = re.search(r"\b([A-Za-z_]\w*::{})\b".format(short), out)
    return m.group(1) if m else None


def _make_bp(cls, specs, *args):
    """Create a breakpoint, trying each spec. A spec may be a callable
    (evaluated only if every earlier spec failed): _qualified sweeps the
    whole symbol table per call, so the gfortran fast path must not pay for
    the Intel fallback. The Python API ignores 'set breakpoint pending off'
    and silently creates PENDING breakpoints for missing symbols (calibrated
    on gdb 16.2), so check and delete."""
    for spec in specs:
        if callable(spec):
            spec = spec()
        if not spec:
            continue
        try:
            bp = cls(spec, *args)
        except gdb.error:
            continue
        if getattr(bp, "pending", False):
            bp.delete()
            continue
        return bp
    return None


def _on_stop(ev):
    LAST_STOP[0] = ev


def _on_exit(_ev):
    DONE[0] = True


def setup():
    gdb.execute("set pagination off")
    gdb.execute("set confirm off")
    gdb.execute("set breakpoint pending off")
    gdb.execute("set width unlimited")
    # don't let gdb print frame args at every stop: slow and floods the
    # log (and pre-'next', explicit-shape dummies have garbage bounds)
    gdb.execute("set print frame-arguments none")
    gdb.events.stop.connect(_on_stop)
    gdb.events.exited.connect(_on_exit)

    plans = []
    for s in SCHEMES:
        portable = PORTABLE.get(s)
        if portable:
            # Compare at the shared portable subroutine that both wrappers
            # call. Its module name differs from the symbol, so the
            # __mod_MOD_proc linker guess can't be formed; the bare DWARF name
            # (gfortran) and the demangled module::proc (Intel) suffice.
            plans.append((s, portable,
                          (portable, lambda p=portable: _qualified(p))))
        else:
            # gfortran bare DWARF name first (the calibrated Derecho fast
            # path), then its linker name, then the module-qualified spelling
            # gdb demangles Intel symbols to (module_mp_<scheme>_run_ ->
            # module::run).
            plans.append((s, s + "_run",
                          (s + "_run", "__{0}_MOD_{0}_run".format(s),
                           lambda s=s: _qualified(s + "_run"))))

    probe_optimized([specs for _s, _t, specs in plans])
    resolved = 0
    for s, target, specs in plans:
        if OPTIMIZED[0]:
            bp = _make_entry_bp_optimized(specs, s, target)
        else:
            bp = _make_bp(EntryBP, specs, s, target)
        if bp is None:
            MANIFEST["breakpoints"][s] = "missing"
        else:
            MANIFEST["breakpoints"][s] = "{} ({})".format(
                bp.location, target) if OPTIMIZED[0] else bp.location
            resolved += 1
    note("{}/{} scheme entry points resolved".format(resolved, len(SCHEMES)))

    step_bp = _make_bp(StepBP, ("cam_run1", "__cam_comp_MOD_cam_run1",
                                lambda: _qualified("cam_run1")))
    if step_bp is None:
        note("WARNING: cam_run1 sentinel not found; hits will not be "
             "timestep-tagged and the run will not auto-terminate")
    return resolved


def main():
    if setup() == 0:
        note("no breakpoints resolved; not running")
        return
    const_done = False
    gdb.execute("run")
    while not DONE[0]:
        ev = LAST_STOP[0]
        LAST_STOP[0] = None
        if ev is None:
            note("stopped without a stop event; aborting")
            break
        entry = None
        step = None
        for b in getattr(ev, "breakpoints", None) or []:
            if isinstance(b, EntryBP):
                entry = b
            elif isinstance(b, StepBP):
                step = b
        if step is not None:
            CURRENT_STEP[0] += 1
            if KILL_AFTER and CURRENT_STEP[0] > KILL_AFTER:
                note("collected {} timesteps; terminating run".format(
                    KILL_AFTER))
                gdb.execute("kill")
                break
            note("timestep {} begins".format(CURRENT_STEP[0]))
        elif entry is None:
            note("non-breakpoint stop ({}); backtrace follows".format(
                getattr(ev, "stop_signal", "?")))
            try:
                gdb.execute("bt 25")
            except gdb.error:
                pass
            break
        elif OPTIMIZED[0]:
            # Optimized binary: the breakpoint is at the raw entry
            # instruction, the one PC where the ABI argument slots are
            # guaranteed; capture immediately (no advance -- moving even one
            # instruction may clobber the argument registers).
            frame = gdb.newest_frame()
            try:
                rec = handle_entry_optimized(entry.scheme, frame)
                ExitBP(frame, rec)
            except Exception as exc:
                note("entry capture failed for {}: {}".format(
                    entry.scheme, exc))
            if not const_done:
                dump_constituents()
                const_done = True
        else:
            # gdb's prologue skip stops inside the argument-descriptor setup,
            # where dummies raise "Location address is not set" and stepping
            # can't escape (the setup cycles). Advance to the first body
            # statement (errmsg=''), where every dummy is live, then capture.
            # Nothing but compiler-generated setup runs between entry and that
            # statement, so 'advance' lands in this same call with inputs
            # untouched (errmsg=''/errflg=0 don't modify physical args).
            frame = gdb.newest_frame()
            body = _first_body_line(entry.scheme, frame)
            if body is not None:
                try:
                    gdb.execute("advance {}:{}".format(
                        os.path.basename(body[0]), body[1]))
                    frame = gdb.newest_frame()
                    # Safety: if advance hit the frame-exit stop (target
                    # was unreachable), we're now in the caller -- fall
                    # back to capturing from wherever we are, but warn.
                    fn = frame.function()
                    if fn and entry.target not in fn.name:
                        note("{}: advance escaped frame (now in {}); "
                             "skipping capture".format(
                                 entry.scheme, fn.name))
                        try:
                            gdb.execute("continue")
                        except gdb.error:
                            break
                        continue
                except gdb.error as exc:
                    note("{}: advance to body line {} failed: {}".format(
                        entry.scheme, body[1], exc))
            else:
                # No body marker found: fall back to the historical single
                # 'next', which suffices for short argument lists.
                note("{}: no body marker found; single-next fallback".format(
                    entry.scheme))
                try:
                    gdb.execute("next")
                except gdb.error:
                    pass
                frame = gdb.newest_frame()
            try:
                rec = handle_entry(entry.scheme, frame)
                ExitBP(frame, rec)
            except Exception as exc:
                note("entry capture failed for {}: {}".format(
                    entry.scheme, exc))
            if not const_done:
                dump_constituents()
                const_done = True
        try:
            gdb.execute("continue")
        except gdb.error:
            break


try:
    main()
finally:
    with open(os.path.join(OUT, "manifest.json"), "w") as f:
        json.dump(MANIFEST, f, indent=1)
    note("manifest written: {} hits".format(len(MANIFEST["hits"])))
