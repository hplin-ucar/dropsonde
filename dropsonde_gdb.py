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

CFG = json.load(open(os.environ["DROPSONDE_CONFIG"]))
ROLE = CFG["role"]            # "cam" | "sima"
OUT = CFG["out_dir"]
SCHEMES = CFG["schemes"]      # unique scheme names, suite order
KILL_AFTER = CFG.get("kill_after_steps", 0)  # 0 = run to natural exit

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

# Upper bound on 'next' steps spent materializing dummy-argument descriptors
# at scheme entry (see the drive loop). Schemes with long argument lists
# (e.g. park_macrophysics, ~50 dummies) spread gfortran's descriptor wiring
# across many declaration lines; this caps the search so an unreadable arg
# can't stall the run. Readiness is normally reached well within this.
ARG_SETUP_MAX_STEPS = 96


def note(msg):
    MANIFEST["notes"].append(msg)
    gdb.write("dropsonde[{}]: {}\n".format(ROLE, msg))


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


def read_array(plan):
    """Read the logical array contents as contiguous subscript-major bytes."""
    extents = plan["extents"]
    strides = plan["strides"]
    elsize = plan["elsize"]
    total = 1
    for e in extents:
        total *= e
    if total == 0:
        return b""
    span = sum(s * (e - 1) for s, e in zip(strides, extents)) + elsize
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
    fn = "{:05d}_{}_h{}_{}_{}.bin".format(FILE_IDX[0], scheme, hit, arg, phase)
    with open(os.path.join(OUT, fn), "wb") as f:
        f.write(data)
    return fn


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


def _capturable_array_dtype(sym):
    """Element dtype if this dummy is a numeric array we capture, else None
    (scalars, character, and derived-type arrays such as const_props).

    Uses the static symbol type so the classification holds even before the
    descriptor is materialized -- so a derived-type arg, which never becomes
    capturable, is excluded from the readiness gate rather than stalling it.
    """
    try:
        t = sym.type.strip_typedefs()
    except gdb.error:
        return None
    if is_character(t):
        return None
    depth = 0
    while t.code == gdb.TYPE_CODE_ARRAY:
        t = t.target().strip_typedefs()
        depth += 1
    if depth == 0:
        return None
    return dtype_of(t)


def _args_ready(frame):
    """True once every numeric-array dummy resolves to a base element
    address, i.e. gfortran has finished wiring its descriptor. Scalars,
    character and derived-type args are ignored so they never gate."""
    for sym in arg_symbols(frame):
        if _capturable_array_dtype(sym) is None:
            continue
        try:
            val = sym.value(frame)
            dims, _base = array_dims(val.type)
            if DIM_ORDER == "reversed":
                dims = dims[::-1]
            elem_addr(sym.name, [d[0] for d in dims])
        except Exception:
            return False
    return True


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
    for sym in arg_symbols(frame):
        name = sym.name
        info = {"kind": "skipped"}
        rec["args"][name] = info
        try:
            val = sym.value(frame)
            t = val.type.strip_typedefs()
            if t.code == gdb.TYPE_CODE_ARRAY and not is_character(t):
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
    def __init__(self, spec, scheme):
        super(EntryBP, self).__init__(spec)
        self.scheme = scheme


class StepBP(gdb.Breakpoint):
    """Sentinel on cam_run1: fires once at the start of every timestep."""
    pass


LAST_STOP = [None]
DONE = [False]
CURRENT_STEP = [0]


def _make_bp(cls, specs, *args):
    """Create a breakpoint, trying each spec. The Python API ignores
    'set breakpoint pending off' and silently creates PENDING breakpoints
    for missing symbols (calibrated on gdb 16.2), so check and delete."""
    for spec in specs:
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

    resolved = 0
    for s in SCHEMES:
        bp = _make_bp(EntryBP, (s + "_run", "__{0}_MOD_{0}_run".format(s)),
                      s)
        if bp is None:
            MANIFEST["breakpoints"][s] = "missing"
        else:
            MANIFEST["breakpoints"][s] = bp.location
            resolved += 1
    note("{}/{} scheme entry points resolved".format(resolved, len(SCHEMES)))

    step_bp = _make_bp(StepBP, ("cam_run1", "__cam_comp_MOD_cam_run1"))
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
        else:
            # Step over gfortran's compiler-generated dummy-argument setup
            # before capturing. Bounds/descriptors for explicit- and
            # assumed-shape dummies are wired in entry code attributed to the
            # declaration lines, past gdb's post-prologue stop, so the args
            # are not readable there yet. A single 'next' covers short
            # argument lists; large schemes (e.g. park_macrophysics, ~50
            # dummies) spread this across many source lines, leaving later
            # descriptors unset ("Location address is not set"). Step until
            # every numeric-array dummy resolves -- which lands on the first
            # executable statement (errmsg=''), still scheme entry, before
            # any input is touched -- capped, and bailing if we ever leave
            # the scheme frame, so an unreadable arg can't run past entry.
            frame = gdb.newest_frame()
            steps = 0
            while steps < ARG_SETUP_MAX_STEPS:
                try:
                    gdb.execute("next")
                except gdb.error:
                    break
                steps += 1
                frame = gdb.newest_frame()
                if entry.scheme not in (frame.name() or ""):
                    note("{}: left scheme frame after {} setup steps; "
                         "capture may be partial".format(entry.scheme, steps))
                    break
                if _args_ready(frame):
                    break
            if steps > 1:
                note("{}: {} setup steps to materialize args{}".format(
                    entry.scheme, steps,
                    "" if _args_ready(frame) else " (cap hit; partial)"))
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
