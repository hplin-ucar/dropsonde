# differ.py -- offline comparison of dropsonde dumps.
#
# Reads <out>/cam/manifest.json and <out>/sima/manifest.json plus the raw
# .bin argument dumps and prints a divergence report. Bit-for-bit equality
# is the bytes fast path; only mismatching arrays are unpacked.
#
# Python 3.6 compatible; stdlib only. Standalone use:
#   python3 differ.py <out_dir>

import json
import os
import shutil
import sys
from array import array
from datetime import datetime

# CAM short name -> CCPP standard name, used to match constituent indices
# between CAM's cnst_name(:) and CAM-SIMA's registered standard names.
# Verified against CAM-SIMA src/data/registry.xml; extend as ports grow.
SPECIAL_CONST_MAP = {
    "Q": "water_vapor_mixing_ratio_wrt_moist_air_and_condensed_water",
    "CLDLIQ":
        "cloud_liquid_water_mixing_ratio_wrt_moist_air_and_condensed_water",
    "CLDICE":
        "cloud_ice_mixing_ratio_wrt_moist_air_and_condensed_water",
    "RAINQM": "rain_mixing_ratio_wrt_moist_air_and_condensed_water",
    "SNOWQM": "snow_mixing_ratio_wrt_moist_air_and_condensed_water",
    "GRAUQM": "graupel_water_mixing_ratio_wrt_moist_air_and_condensed_water",
    "NUMLIQ": "mass_number_concentration_of_cloud_liquid_wrt_moist_air_and_condensed_water",
    "NUMRAI": "mass_number_concentration_of_rain_wrt_moist_air_and_condensed_water",
    "NUMICE":
        "mass_number_concentration_of_ice_wrt_moist_air_and_condensed_water",
    "NUMSNO": "mass_number_concentration_of_snow_wrt_moist_air_and_condensed_water",
    "NUMGRA": "mass_number_concentration_of_graupel_wrt_moist_air_and_condensed_water"
}

TYPECODE = {"f8": "d", "f4": "f", "i8": "q", "i4": "i", "i2": "h", "i1": "b"}
ELSIZE = {"f8": 8, "f4": 4, "i8": 8, "i4": 4, "i2": 2, "i1": 1}

# Not compared: errmsg/errflg are uninitialized at entry on the CAM side
# (callers pass stack garbage) and any nonzero errflg aborts the run
# instantly anyway; iulog is each model's own log unit number.
SKIP_ARGS = {"errmsg", "errflg", "iulog"}

# --- portable MAM aerosol routines (modal_aero_gasaerexch/rename) -----------
# These are called from two different constituent-index conventions. CAM passes
# the mozart vmr/vmrcw arrays (gas_pcnst wide, loffset=imozart-1) with pointer
# arrays in full pcnst-space; SIMA passes one packed constituent array
# (loffset=0). Inside the routine every constituent is reached as
# array(:,:, pointer - loffset). So the constituent axis of these args cannot be
# compared byte-for-byte: it must be realigned per species by shifting each
# model's own loffset (see the _portable_* helpers). A hit is treated as
# MAM-convention when it carries a `loffset` argument on both sides.

# constituent-indexed args, tagged by which per-model array they live in:
#   'q'    = interstitial (CAM mozart vmr)
#   'qqcw' = cloud-borne  (CAM separate vmrcw array; SIMA packs it alongside q)
PORTABLE_CONST_ARG_ROLE = {
    "q": "q", "dqdt": "q", "dqdt_other": "q", "dotendrn": "q",
    "qsrflx": "q", "qsrflx_gaexch": "q", "dqdt_rnpos": "q", "dotend": "q",
    "qqcw": "qqcw", "dqqcwdt": "qqcw", "dqqcwdt_other": "qqcw",
    "dotendqqcwrn": "qqcw", "qqcwsrflx": "qqcw",
}

# Index-space / sentinel args: their values are in each model's own constituent
# space (pcnst vs packed) or are unused-slot fill, so they are expected to
# differ and carry no physics signal. Noted, never counted as divergences.
PORTABLE_CONVENTION_ARGS = {
    "loffset", "num_q",
    "lmassptr_amode", "lmassptrcw_amode", "numptr_amode", "numptrcw_amode",
    "lspecfrma_renamexf", "lspecfrmc_renamexf",
    "lspectooa_renamexf", "lspectooc_renamexf",
    "modeptr_stracoar",
}

# Mode/species physical-constant tables, shape (nspec_max, ntot_amode). Unused
# (species,mode) slots hold different padding fill in each model; the valid
# slots (a real species) should match. Compared only over slots where both
# models have a positive lmassptr_amode entry.
PORTABLE_METADATA_ARGS = {"specmw_amode", "specdens_amode"}


def load_role(outdir, role):
    d = os.path.join(outdir, role)
    path = os.path.join(d, "manifest.json")
    if not os.path.isfile(path):
        return None
    with open(path) as f:
        man = json.load(f)
    man["_dir"] = d
    return man


def blob(man, fn):
    with open(os.path.join(man["_dir"], fn), "rb") as f:
        return f.read()


def unpack(data, dt):
    a = array(TYPECODE[dt])
    a.frombytes(data)
    return a


def fmt_val(x):
    if isinstance(x, float):
        return "{:.17g}".format(x)
    return str(x)


def lin_to_sub(li, extents, los):
    subs = []
    rem = li
    for e, lo in zip(extents, los):
        subs.append(lo + rem % e)
        rem //= e
    return "(" + ",".join(str(s) for s in subs) + ")"


def diff_stats(va, vb):
    """(count, (absdiff, linear_index, a, b) of worst differing element)."""
    n = min(len(va), len(vb))
    count = 0
    worst = None
    for k in range(n):
        x, y = va[k], vb[k]
        if x == y or (x != x and y != y):
            continue
        count += 1
        d = abs(x - y)
        if d != d:
            d = float("inf")
        if worst is None or d > worst[0]:
            worst = (d, k, x, y)
    return count, worst


def array_diff_text(sdata, cdata, dt, extents, los):
    """None if identical, else a human-readable stats line."""
    if sdata == cdata:
        return None
    va = unpack(sdata, dt)
    vb = unpack(cdata, dt)
    count, worst = diff_stats(va, vb)
    if count == 0:
        return "byte-level difference only (padding?)"
    d, li, x, y = worst
    return ("{}/{} elements differ, max |diff| {:.3e} at {}: "
            "sima={} cam={}".format(count, len(va), d,
                                    lin_to_sub(li, extents, los),
                                    fmt_val(x), fmt_val(y)))


def _exit_diff_written_only(sdata, cdata, s_entry, c_entry, dt, extents,
                            los):
    """Like array_diff_text, but ignore elements left untouched (entry ==
    exit, bitwise) by the scheme in BOTH models: intent(out) args are not
    always fully defined, and untouched elements just echo each caller's
    prior buffer contents. Returns (text_or_None, n_ignored)."""
    es = ELSIZE[dt]
    n = min(len(sdata), len(cdata)) // es
    keep = []
    ignored = 0
    for k in range(n):
        a = k * es
        b = a + es
        if sdata[a:b] == cdata[a:b]:
            continue
        if s_entry[a:b] == sdata[a:b] and c_entry[a:b] == cdata[a:b]:
            ignored += 1
            continue
        keep.append(k)
    if not keep:
        return None, ignored
    va = unpack(sdata, dt)
    vb = unpack(cdata, dt)
    count = 0
    worst = None
    for k in keep:
        x, y = va[k], vb[k]
        if x == y or (x != x and y != y):
            continue
        count += 1
        d = abs(x - y)
        if d != d:
            d = float("inf")
        if worst is None or d > worst[0]:
            worst = (d, k, x, y)
    if worst is None:
        return None, ignored
    d, li, x, y = worst
    extra = (" (+{} elements never written in either model, "
             "not compared)".format(ignored)) if ignored else ""
    return ("{}/{} elements differ{}, max |diff| {:.3e} at {}: "
            "sima={} cam={}".format(count, n, extra, d,
                                    lin_to_sub(li, extents, los),
                                    fmt_val(x), fmt_val(y)), ignored)


def build_const_map(cam_names, sima_names):
    """(pairs, cam_unmatched, sima_unmatched); pairs are 0-based
    (cam_idx, sima_idx, cam_name, sima_name)."""
    sima_lookup = {n: i for i, n in enumerate(sima_names)}
    pairs = []
    cam_un = []
    used = set()
    for i, cn in enumerate(cam_names):
        j = None
        std = SPECIAL_CONST_MAP.get(cn)
        if std is not None and std in sima_lookup:
            j = sima_lookup[std]
        elif cn in sima_lookup:
            j = sima_lookup[cn]
        elif "cnst_" + cn in sima_lookup:
            # constituents auto-registered from a snapshot file keep
            # their netcdf variable name (cnst_<CAM short name>) as
            # their standard name
            j = sima_lookup["cnst_" + cn]
        else:
            cands = [k for k, n in enumerate(sima_names)
                     if cn.lower() == n or cn.lower() in n.split("_")]
            if len(cands) == 1:
                j = cands[0]
        if j is None or j in used:
            cam_un.append((i, cn))
        else:
            used.add(j)
            pairs.append((i, j, cn, sima_names[j]))
    sima_un = [(j, n) for j, n in enumerate(sima_names) if j not in used]
    return pairs, cam_un, sima_un


def _disp_scheme(scheme, portable):
    """Display label for a scheme. When the scheme is compared at a shared
    portable subroutine (dropsonde:portable SDF annotation), show
    'scheme -> portable_sub' so the report makes clear what is actually being
    compared; otherwise just the scheme name."""
    sym = (portable or {}).get(scheme)
    return "{} -> {}".format(scheme, sym) if sym else scheme


class Reporter(object):
    def __init__(self, portable=None):
        self.n_diffs = 0
        self.first_printed = False
        self.alignment_suspect = False
        self.first_pair_seen = False
        self.cur_scheme = None
        self.by_scheme = {}  # scheme -> [entry_diffs, exit_diffs]
        self.portable = portable or {}  # scheme -> portable_sub, for labels

    def diff(self, scheme, hit, phase, arg, text, extra_lines=None):
        self.n_diffs += 1
        if self.cur_scheme:
            c = self.by_scheme.setdefault(self.cur_scheme, [0, 0])
            c[0 if phase == "entry" else 1] += 1
        tag = "INPUTS DIFFER" if phase == "entry" else "OUTPUTS DIFFER"
        head = "  [{}] {} (hit {}) arg {}: {}".format(
            tag, scheme, hit, arg, text)
        if not self.first_printed:
            self.first_printed = True
            print("")
            print("FIRST DIVERGENCE " + "-" * 47)
            print(head)
            for ln in extra_lines or []:
                print("    " + ln)
            print("-" * 64)
        else:
            print(head)
            for ln in extra_lines or []:
                print("    " + ln)

    def note(self, text):
        print("  [note] " + text)


def kdesc(info):
    why = info.get("why")
    k = info.get("kind")
    return "{} ({})".format(k, why) if why else str(k)


def _species_slab(data, elsize, extents, axis, si):
    """Bytes for the si-th (0-based) species along `axis` of a Fortran-order
    (column-major) dump. Last-axis species are a contiguous block (fast path);
    an interior axis (e.g. qsrflx (ncol,nspec,nsrflx)) is a strided gather."""
    inner = 1
    for e in extents[:axis]:
        inner *= e
    block = inner * elsize                      # bytes per species step
    if axis == len(extents) - 1:
        return data[si * block:(si + 1) * block]
    outer = 1
    for e in extents[axis + 1:]:
        outer *= e
    period = block * extents[axis]              # bytes per step of the next axis
    base = si * block
    return b"".join(data[base + o * period:base + o * period + block]
                    for o in range(outer))


def _const_axis(sext, cext, sima_names):
    """Axis carrying the constituent index for a portable MAM arg, or None.
    Anchored on SIMA, whose packed array length == number of registered
    constituents; every other axis must agree between the two models."""
    if len(sext) != len(cext) or not sima_names:
        return None
    for k in range(len(sext)):
        if (sext[k] == len(sima_names) and
                sext[:k] == cext[:k] and sext[k + 1:] == cext[k + 1:]):
            return k
    return None


def _int_arg(man, args, name):
    """Decoded integer entry values of a pointer arg, or None if not captured."""
    info = args.get(name)
    if not info or not info.get("entry_file"):
        return None
    return unpack(blob(man, info["entry_file"]), info["dtype"])


def _cloudborne_pairs(srec, crec, sman, cman, loff_s, loff_c, valid,
                      sima_names):
    """Cloud-borne (qqcw) species pairs (sima_local, cam_local, name, name) in
    each model's local index space. CAM keeps cloud-borne in a separate vmrcw
    array with no registered constituent name, so alignment is driven by the
    captured mode/species pointer arrays keyed by (species,mode) slot -- the
    same logical MAM table in both models -- rather than the constituent map.
    lmassptrcw_amode/numptrcw_amode hold pcnst-space indices; local = ptr -
    loffset. Species names come from SIMA's registered cloud-borne (*_c*)
    constituents for reporting."""
    pairs = []

    def add(sptr, cptr):
        # pointers are 1-based pcnst-space indices; the routine reads
        # array(:,:, ptr - loffset), so the 0-based local index is
        # ptr - loffset - 1.
        if sptr is None or cptr is None:
            return
        sj, ci = sptr - loff_s - 1, cptr - loff_c - 1
        if sj < 0 or ci < 0:
            return
        nm = sima_names[sptr - 1] if 0 < sptr <= len(sima_names) else "?"
        pairs.append((sj, ci, nm, nm))

    # mass species: valid (species,mode) slots only (padding holds fill)
    lms = _int_arg(sman, srec["args"], "lmassptrcw_amode")
    lmc = _int_arg(cman, crec["args"], "lmassptrcw_amode")
    if lms is not None and lmc is not None:
        for k in sorted(valid):
            if k < len(lms) and k < len(lmc):
                add(lms[k], lmc[k])
    # number species: one per aerosol mode
    nms = _int_arg(sman, srec["args"], "numptrcw_amode")
    nmc = _int_arg(cman, crec["args"], "numptrcw_amode")
    if nms is not None and nmc is not None:
        for k in range(min(len(nms), len(nmc))):
            if nms[k] > 0 and nmc[k] > 0:
                add(nms[k], nmc[k])
    return pairs


def _portable_ctx(srec, crec, sman, cman, const_ctx):
    """Per-hit realignment context for a MAM-convention hit, or None. Carries
    each side's loffset, the interstitial species pairs already shifted into
    per-model local (array) index space, the valid mode/species mask used for
    the physical-constant tables, and the cloud-borne (qqcw) species pairs."""
    sa = srec["args"].get("loffset")
    ca = crec["args"].get("loffset")
    if not sa or not ca:
        return None
    loff_s = sa.get("entry_value")
    loff_c = ca.get("entry_value")
    if loff_s is None or loff_c is None:
        return None
    _cam_names, sima_names, pairs = const_ctx
    # pairs are 0-based (cam_full, sima_full); local index = full - loffset
    q_local = [(sj - loff_s, ci - loff_c, cn, sn)
               for (ci, sj, cn, sn) in pairs
               if sj - loff_s >= 0 and ci - loff_c >= 0]
    # valid (species,mode) slots: a real species has lmassptr_amode > 0 in both
    valid = set()
    la = srec["args"].get("lmassptr_amode")
    lc = crec["args"].get("lmassptr_amode")
    if (la and lc and la.get("entry_file") and lc.get("entry_file") and
            la.get("dtype") == lc.get("dtype")):
        sv = unpack(blob(sman, la["entry_file"]), la["dtype"])
        cv = unpack(blob(cman, lc["entry_file"]), lc["dtype"])
        for k in range(min(len(sv), len(cv))):
            if sv[k] > 0 and cv[k] > 0:
                valid.add(k)
    qqcw_local = _cloudborne_pairs(srec, crec, sman, cman, loff_s, loff_c,
                                   valid, sima_names)
    return {"loff_s": loff_s, "loff_c": loff_c, "q_local": q_local,
            "valid": valid, "qqcw_local": qqcw_local}


def _compare_species_axis(label, hit, phase, arg, sdata, cdata, dt, sext, cext,
                          axis, local_pairs, los, rep, space):
    """Compare a constituent-indexed array per mapped species, each side
    indexed in its own local (array) space. `space` labels which set of pairs
    is used ('interstitial' or 'cloud-borne')."""
    elsize = ELSIZE[dt]
    n_s, n_c = sext[axis], cext[axis]
    res_ext = [e for i, e in enumerate(sext) if i != axis]
    res_los = [l for i, l in enumerate(los) if i != axis]
    lines, diff_names = [], []
    for sj, ci, cn, sn in local_pairs:
        if not (0 <= sj < n_s and 0 <= ci < n_c):
            continue
        ss = _species_slab(sdata, elsize, sext, axis, sj)
        cs = _species_slab(cdata, elsize, cext, axis, ci)
        if not res_ext:                        # one element per species
            txt = None
            if ss != cs:
                txt = "sima={} cam={}".format(fmt_val(unpack(ss, dt)[0]),
                                              fmt_val(unpack(cs, dt)[0]))
        else:
            txt = array_diff_text(ss, cs, dt, res_ext, res_los)
        if txt:
            diff_names.append(cn)
            lines.append("{} (cam idx {}, sima idx {}) <-> {}: {}".format(
                cn, ci + 1, sj + 1, sn, txt))
    if lines:
        rep.diff(label, hit, phase, arg,
                 "constituent-indexed array ({}, per-model local "
                 "index), {}/{} mapped species differ: {}".format(space,
                     len(lines), len(local_pairs), ", ".join(diff_names)),
                 lines)


def _compare_metadata_valid(label, hit, phase, arg, sdata, cdata, dt, ext,
                            los, valid, rep):
    """Compare a mode/species constant table over valid slots only; unused
    padding slots hold different fill in each model and are ignored."""
    va, vb = unpack(sdata, dt), unpack(cdata, dt)
    count, worst = 0, None
    for k in range(min(len(va), len(vb))):
        if k not in valid:
            continue
        x, y = va[k], vb[k]
        if x == y or (x != x and y != y):
            continue
        count += 1
        d = abs(x - y)
        if d != d:
            d = float("inf")
        if worst is None or d > worst[0]:
            worst = (d, k, x, y)
    if count:
        d, li, x, y = worst
        rep.diff(label, hit, phase, arg,
                 "{}/{} valid slots differ (unused padding ignored), max "
                 "|diff| {:.3e} at {}: sima={} cam={}".format(
                     count, len(valid), d, lin_to_sub(li, ext, los),
                     fmt_val(x), fmt_val(y)))


def _compare_arg(label, hit, phase, arg, sinfo, cinfo, sman, cman,
                 const_ctx, rep, intent=None, port=None):
    fkey = phase + "_file"
    vkey = phase + "_value"
    cam_names, sima_names, pairs = const_ctx
    kind = sinfo.get("kind")
    if kind != cinfo.get("kind"):
        if phase == "entry":
            rep.diff(label, hit, phase, arg,
                     "argument kind mismatch: sima {} vs cam {}".format(
                         kdesc(sinfo), kdesc(cinfo)))
        return
    if kind == "error":
        if phase == "entry":
            rep.note("{} (hit {}): arg {} capture errored on both sides: "
                     "sima {} / cam {}".format(label, hit, arg,
                                               kdesc(sinfo), kdesc(cinfo)))
        return

    if kind == "array":
        sfile = sinfo.get(fkey)
        cfile = cinfo.get(fkey)
        if not sfile or not cfile:
            return  # incomplete capture (noted in manifest)
        dt = sinfo["dtype"]
        if dt != cinfo.get("dtype"):
            if phase == "entry":
                rep.diff(label, hit, phase, arg,
                         "dtype mismatch: sima {} vs cam {}".format(
                             dt, cinfo.get("dtype")))
            return
        sext = sinfo["extents"]
        cext = cinfo["extents"]
        sdata = blob(sman, sfile)
        cdata = blob(cman, cfile)

        # portable MAM args: realign the constituent axis into each model's own
        # local (array) index space before comparing (see PORTABLE_* tables).
        if port is not None:
            if arg in PORTABLE_METADATA_ARGS:
                _compare_metadata_valid(label, hit, phase, arg, sdata, cdata,
                                        dt, sext, sinfo["los"], port["valid"],
                                        rep)
                return
            role = PORTABLE_CONST_ARG_ROLE.get(arg)
            local = {"q": port["q_local"],
                     "qqcw": port["qqcw_local"]}.get(role)
            if local is not None:
                axis = _const_axis(sext, cext, const_ctx[1])
                if axis is not None:
                    space = "interstitial" if role == "q" else "cloud-borne"
                    _compare_species_axis(label, hit, phase, arg, sdata, cdata,
                                          dt, sext, cext, axis, local,
                                          sinfo["los"], rep, space)
                    return

        # constituent-indexed array: last dim is the constituent index,
        # e.g. q(ncol,pver,pcnst), fluxes cflx(ncol,pcnst), or per-species
        # config like qmincg(pcnst) / do_diffusion_const(pcnst)
        rank = len(sext)
        const_arr = (rank == len(cext) and rank in (1, 2, 3) and
                     sima_names and cam_names and
                     sext[-1] == len(sima_names) and
                     cext[-1] == len(cam_names) and
                     sext[:-1] == cext[:-1])
        if const_arr:
            chunk = ELSIZE[dt]
            for e in sext[:-1]:
                chunk *= e
            lines = []
            diff_names = []
            for ci, sj, cn, sn in pairs:
                ss = sdata[sj * chunk:(sj + 1) * chunk]
                cs = cdata[ci * chunk:(ci + 1) * chunk]
                if rank == 1:
                    # one element per species: just show the two values
                    txt = None
                    if ss != cs:
                        txt = "sima={} cam={}".format(
                            fmt_val(unpack(ss, dt)[0]),
                            fmt_val(unpack(cs, dt)[0]))
                else:
                    txt = array_diff_text(ss, cs, dt, sext[:-1],
                                          sinfo["los"][:-1])
                if txt:
                    diff_names.append(cn)
                    lines.append(
                        "{} (cam idx {}, sima idx {}) <-> {}: {}".format(
                            cn, ci + 1, sj + 1, sn, txt))
            if lines:
                rep.diff(label, hit, phase, arg,
                         "constituent-indexed array, {}/{} mapped species"
                         " differ: {}".format(len(lines), len(pairs),
                                              ", ".join(diff_names)),
                         lines)
            return

        if sext != cext:
            if phase == "entry":
                rep.diff(label, hit, phase, arg,
                         "shape mismatch: sima {} vs cam {}".format(
                             sext, cext))
            return
        if (phase == "exit" and intent == "out" and sdata != cdata and
                sinfo.get("entry_file") and cinfo.get("entry_file")):
            txt, ignored = _exit_diff_written_only(
                sdata, cdata, blob(sman, sinfo["entry_file"]),
                blob(cman, cinfo["entry_file"]), dt, sext, sinfo["los"])
            if txt is None:
                if ignored:
                    rep.note("{} (hit {}): arg {}: every element the "
                             "scheme wrote matches; {} elements were "
                             "never written in either model (exit == "
                             "entry bitwise) and still hold each "
                             "caller's unrelated pre-call memory, so "
                             "they are not compared".format(
                                 label, hit, arg, ignored))
                return
            rep.diff(label, hit, phase, arg, txt)
            return
        txt = array_diff_text(sdata, cdata, dt, sext, sinfo["los"])
        if txt:
            rep.diff(label, hit, phase, arg, txt)

    elif kind in ("scalar", "char"):
        sv = sinfo.get(vkey)
        cv = cinfo.get(vkey)
        if sv is None or cv is None:
            return
        if sv != cv:
            rep.diff(label, hit, phase, arg,
                     "sima={} cam={}".format(fmt_val(sv), fmt_val(cv)))


def compare_hit(srec, crec, sman, cman, const_ctx, rep, intents=None):
    """Compare one aligned hit pair; entry args first, then exit.

    With .meta intents available, intent(out) args are not compared at
    entry (caller-side garbage) and intent(in) args not at exit. Exit
    reports are suppressed for args whose entry already differed -- the
    'outputs differ given identical inputs' signal is gone for those.
    """
    bucket = srec.get("_bucket", "<toplevel>")
    label = _disp_scheme(srec["scheme"], rep.portable)
    if bucket != "<toplevel>":
        label += " via " + bucket
    label += " [step {}]".format(srec.get("step", "?"))
    hit = srec.get("_occ", srec["hit"])
    sch_int = (intents or {}).get(srec["scheme"], {})
    rep.cur_scheme = srec["scheme"]
    port = _portable_ctx(srec, crec, sman, cman, const_ctx)
    conv = set()
    entry_diffed = set()
    for phase in ("entry", "exit"):
        suppressed = []
        for arg, sinfo in srec["args"].items():
            if arg in SKIP_ARGS:
                continue
            it = sch_int.get(arg)
            if it == "out" and phase == "entry":
                continue
            if it == "in" and phase == "exit":
                continue
            cinfo = crec["args"].get(arg)
            if cinfo is None:
                if phase == "entry":
                    rep.note("{} (hit {}): arg {} absent on CAM side".format(
                        label, hit, arg))
                continue
            if port is not None and arg in PORTABLE_CONVENTION_ARGS:
                conv.add(arg)
                continue
            if phase == "exit" and arg in entry_diffed:
                suppressed.append(arg)
                continue
            pre = rep.n_diffs
            _compare_arg(label, hit, phase, arg, sinfo, cinfo, sman, cman,
                         const_ctx, rep, it, port)
            if phase == "entry" and rep.n_diffs > pre:
                entry_diffed.add(arg)
        if phase == "exit" and suppressed:
            rep.note("{} (hit {}): exit comparison suppressed for {} args "
                     "with differing inputs: {}".format(
                         label, hit, len(suppressed),
                         ", ".join(suppressed)))
        if phase == "entry" and not rep.first_pair_seen:
            rep.first_pair_seen = True
            if rep.n_diffs > 0:
                rep.alignment_suspect = True
    if conv:
        rep.note("{} (hit {}): index-space/convention args not compared "
                 "(expected to differ between models): {}".format(
                     label, hit, ", ".join(sorted(conv))))


def shared_caller_map(cam, sima):
    """Per scheme, caller names seen in BOTH models. A scheme called from
    inside another scheme (shared atmospheric_physics code) has the same
    caller symbol in both binaries; suite-level call sites (CCPP cap vs
    CAM driver) never match and bucket to '<toplevel>'."""
    maps = []
    for man in (cam, sima):
        d = {}
        for rec in man["hits"]:
            c = rec.get("caller", "?")
            if c != "?":
                d.setdefault(rec["scheme"], set()).add(c)
        maps.append(d)
    shared = {}
    for s in set(maps[0]) | set(maps[1]):
        shared[s] = maps[0].get(s, set()) & maps[1].get(s, set())
    return shared


def tag_hits(man, shared):
    """Group hits by (scheme, step, bucket); annotate each rec with its
    bucket and occurrence index within the group. Returns the groups."""
    groups = {}
    for rec in man["hits"]:
        s = rec["scheme"]
        c = rec.get("caller", "?")
        b = c if c in shared.get(s, set()) else "<toplevel>"
        key = (s, rec.get("step", 0), b)
        lst = groups.setdefault(key, [])
        rec["_bucket"] = b
        rec["_occ"] = len(lst)
        lst.append(rec)
    return groups


def count_matching_entry_args(srec, crec, sman, cman, sch_int=None,
                              const_ctx=None):
    """(n_bitwise_identical, n_comparable) over entry captures.
    intent(out) args are excluded when intents are known -- their entry
    values are caller-side garbage. Constituent-indexed arrays are
    matched per mapped species (the raw bytes never match across the
    models' different constituent orderings)."""
    cam_names, sima_names, pairs = const_ctx or ([], [], [])
    same = 0
    total = 0
    for arg, sinfo in srec["args"].items():
        if arg in SKIP_ARGS:
            continue
        if sch_int and sch_int.get(arg) == "out":
            continue
        cinfo = crec["args"].get(arg)
        if cinfo is None or sinfo.get("kind") != cinfo.get("kind"):
            continue
        if sinfo.get("kind") == "array":
            if not sinfo.get("entry_file") or not cinfo.get("entry_file"):
                continue
            sext = sinfo["extents"]
            cext = cinfo["extents"]
            if (len(sext) == len(cext) and len(sext) in (1, 2, 3) and
                    sima_names and cam_names and
                    sext[-1] == len(sima_names) and
                    cext[-1] == len(cam_names) and
                    sext[:-1] == cext[:-1]):
                if not pairs:
                    continue  # no species mapping: cannot judge this arg
                sdata = blob(sman, sinfo["entry_file"])
                cdata = blob(cman, cinfo["entry_file"])
                chunk = ELSIZE[sinfo["dtype"]]
                for e in sext[:-1]:
                    chunk *= e
                total += 1
                if all(sdata[sj * chunk:(sj + 1) * chunk] ==
                       cdata[ci * chunk:(ci + 1) * chunk]
                       for ci, sj, _, _ in pairs):
                    same += 1
                continue
            total += 1
            if (sext == cext and
                    blob(sman, sinfo["entry_file"]) ==
                    blob(cman, cinfo["entry_file"])):
                same += 1
        elif sinfo.get("kind") in ("scalar", "char"):
            if sinfo.get("entry_value") is None:
                continue
            total += 1
            if sinfo.get("entry_value") == cinfo.get("entry_value"):
                same += 1
    return same, total


def offset_scan(first_srec, cam_g, sima_g, sman, cman, intents=None,
                const_ctx=None):
    """When alignment looks wrong, show which cam/sima steps the first
    compared sima hit actually matches, bitwise."""
    s = first_srec["scheme"]
    b = first_srec["_bucket"]
    occ = first_srec["_occ"]
    sch_int = (intents or {}).get(s, {})
    print("offset scan: entry args of {} (sima step {}) vs every dumped "
          "step:".format(s, first_srec["step"]))
    cam_steps = sorted(set(t for (ss, t, bb) in cam_g
                           if ss == s and bb == b))
    for t in cam_steps:
        lst = cam_g.get((s, t, b), [])
        if occ < len(lst):
            same, total = count_matching_entry_args(
                first_srec, lst[occ], sman, cman, sch_int, const_ctx)
            print("  vs cam  step {}: {:2d}/{} args bitwise identical"
                  .format(t, same, total))
    sima_steps = sorted(set(t for (ss, t, bb) in sima_g
                            if ss == s and bb == b
                            and t != first_srec["step"]))
    for t in sima_steps:
        lst = sima_g.get((s, t, b), [])
        if occ < len(lst):
            same, total = count_matching_entry_args(
                first_srec, lst[occ], sman, sman, sch_int, const_ctx)
            print("  vs sima step {}: {:2d}/{} args bitwise identical "
                  "(identical => snapshot record repeated!)".format(
                      t, same, total))


class _Tee(object):
    """Write-through to several streams (console + report.txt archive)."""
    def __init__(self, *streams):
        self.streams = streams

    def write(self, s):
        for st in self.streams:
            st.write(s)

    def flush(self):
        for st in self.streams:
            st.flush()


def _archive_report(outdir, report_path):
    """Copy report.txt to a sibling reports/ directory (next to <out>, so it
    survives emptying the dump directory) under a timestamped, SDF-named
    filename. Returns the destination path, or None if it couldn't be
    written."""
    try:
        with open(os.path.join(outdir, "suite.json")) as f:
            sdf = json.load(f).get("sdf", "")
    except (IOError, ValueError):
        sdf = ""
    sdf_name = os.path.splitext(os.path.basename(sdf))[0] or "unknown"
    stamp = datetime.now().strftime("%Y-%m-%d_%H%M%S")
    reports_dir = os.path.join(os.path.dirname(outdir), "reports")
    dest = os.path.join(
        reports_dir, "dropsonde_report_{}_{}.txt".format(sdf_name, stamp))
    try:
        os.makedirs(reports_dir, exist_ok=True)
        shutil.copyfile(report_path, dest)
    except (IOError, OSError) as e:
        print("(could not archive report to {}: {})".format(dest, e))
        return None
    return dest


def report(outdir, suite_order, steps, intents=None):
    """Run the comparison; everything printed is also archived to
    <out>/report.txt so it can be re-read (or regenerated after editing
    differ.py) without re-running the models: python3 differ.py <out>. A
    timestamped, SDF-named copy is also kept in a sibling reports/ directory."""
    path = os.path.join(outdir, "report.txt")
    old = sys.stdout
    try:
        with open(path, "w") as f:
            sys.stdout = _Tee(old, f)
            ok = _report(outdir, suite_order, steps, intents)
    finally:
        sys.stdout = old
    print("(report archived to {})".format(path))
    dest = _archive_report(outdir, path)
    if dest:
        print("(report copied to {})".format(dest))
    return ok


def _report(outdir, suite_order, steps, intents=None):
    print("")
    print("=" * 64)
    print("dropsonde report  ({})".format(outdir))
    print("=" * 64)
    try:
        with open(os.path.join(outdir, "suite.json")) as f:
            suite_meta = json.load(f)
    except (IOError, ValueError):
        suite_meta = {}
    portable = suite_meta.get("portable", {}) or {}

    def disp(s):
        return _disp_scheme(s, portable)

    ver = suite_meta.get("version")
    if ver:
        print("dropsonde version: {}".format(ver))
    sdf_path = suite_meta.get("sdf")
    if sdf_path:
        print("SDF:       {}".format(sdf_path))
    cam_case = suite_meta.get("cam_case")
    if cam_case:
        print("CAM case:  {}".format(cam_case))
    sima_case = suite_meta.get("sima_case")
    if sima_case:
        print("SIMA case: {}".format(sima_case))
    print("steps:     {}".format(steps))

    cam = load_role(outdir, "cam")
    sima = load_role(outdir, "sima")
    missing = [r for r, m in (("cam", cam), ("sima", sima)) if m is None]
    if missing:
        print("ERROR: no manifest for: {} (gdb run failed? see gdb.log)"
              .format(", ".join(missing)))
        return False

    for man, role in ((cam, "cam"), (sima, "sima")):
        for note_ in man.get("notes", []):
            if "FAILED" in note_ or "WARNING" in note_:
                print("  {} gdb note: {}".format(
                    role, note_.splitlines()[0]))

    uniq = []
    for s in suite_order:
        if s not in uniq:
            uniq.append(s)

    # --- breakpoint coverage -------------------------------------------
    compared, sima_only, cam_only, neither = [], [], [], []
    for s in uniq:
        c_ok = cam["breakpoints"].get(s, "missing") != "missing"
        s_ok = sima["breakpoints"].get(s, "missing") != "missing"
        if c_ok and s_ok:
            compared.append(s)
        elif s_ok:
            sima_only.append(s)
        elif c_ok:
            cam_only.append(s)
        else:
            neither.append(s)
    print("schemes: {} compared".format(len(compared)))
    if sima_only:
        print("  SIMA-only (no CAM symbol, not compared): " +
              ", ".join(sima_only))
    if cam_only:
        print("  CAM-only (no SIMA symbol!): " + ", ".join(cam_only))
    if neither:
        print("  unresolved in both: " + ", ".join(neither))

    # --- constituent mapping -------------------------------------------
    cam_names = cam.get("constituents") or []
    sima_names = sima.get("constituents") or []
    pairs, cam_un, sima_un = build_const_map(cam_names, sima_names)
    print("")
    print("constituent mapping (CAM {} species, SIMA {}):".format(
        len(cam_names), len(sima_names)))
    for ci, sj, cn, sn in pairs:
        print("  cam q(:,:,{:2d}) {:16s} <-> sima q(:,:,{:2d}) {}".format(
            ci + 1, cn, sj + 1, sn))
    for i, cn in cam_un:
        print("  cam q(:,:,{:2d}) {:16s} <-> (unmatched, not compared)"
              .format(i + 1, cn))
    for j, sn in sima_un:
        print("  sima q(:,:,{:2d}) {} (unmatched, not compared)".format(
            j + 1, sn))
    if not cam_names or not sima_names:
        print("  (capture incomplete on {} side; constituent-indexed "
              "arrays compared element-wise only)".format(
                  "CAM" if not cam_names else "SIMA"))
    if cam_un and sima_un:
        print("")
        print("  *** WARNING: {} CAM and {} SIMA constituents are "
              "UNMATCHED on both sides.".format(len(cam_un), len(sima_un)))
        print("  *** They are not name-matched, so they are not compared "
              "via the constituent map (portable")
        print("  *** MAM cloud-borne *_c* species are an expected exception: "
              "they carry no CAM constituent")
        print("  *** name and are realigned separately from the pointer "
              "arrays). Otherwise this often means a")
        print("  *** missing SPECIAL_CONST_MAP entry (differ.py) or a "
              "duplicate standard name (see example 3 in docs/).")
        print("  *** Unmatched CAM:  {}".format(
            ", ".join(cn for _, cn in cam_un)))
        print("  *** Unmatched SIMA: {}".format(
            ", ".join(sn for _, sn in sima_un)))
    const_ctx = (cam_names, sima_names, pairs)

    # --- alignment ------------------------------------------------------
    # Hits are tagged with the timestep (cam_run1 sentinel) and bucketed
    # by caller: sima step t pairs with cam step t+1, occurrence-by-
    # occurrence within each (scheme, step, bucket).
    if not any(r.get("step", 0) >= 1 for r in sima["hits"]):
        print("")
        print("ERROR: no timestep tags on SIMA hits -- the cam_run1 "
              "sentinel did not resolve or never fired (see gdb.log).")
        return False
    if not any(r.get("step", 0) >= 1 for r in cam["hits"]):
        print("")
        print("ERROR: no timestep tags on CAM hits -- the cam_run1 "
              "sentinel did not resolve or never fired (see gdb.log).")
        return False

    shared = shared_caller_map(cam, sima)
    cam_g = tag_hits(cam, shared)
    sima_g = tag_hits(sima, shared)

    def bucket_parent(b):
        # 'bretherton_park_diff::bretherton_park_diff_run' -> the
        # compared scheme whose _run subroutine made the nested call
        name = b.split("::")[-1]
        if name.endswith("_run"):
            name = name[:-4]
        return name if name in compared else None

    print("")
    print("call alignment (hits per compared step; sima step t pairs "
          "with cam step t+1;")
    print("indented rows are calls made from inside the parent scheme):")
    width = max(len(disp(s)) for s in compared) + 2

    def vmark(n_s, n_c):
        if n_s == n_c:
            return ""
        if n_c == 0:
            return "  <-- no CAM hits: not compared"
        if n_s == 0:
            return "  <-- no SIMA hits: not compared"
        return "  <-- hit counts differ (see notes)"

    def row(label, n_s, n_c, depth):
        print("  {}{:<{w}} sima x{:<3d} cam x{:<3d}{}".format(
            "    " * depth, label, n_s, n_c, vmark(n_s, n_c),
            w=max(1, width - 4 * depth)).rstrip())

    for t in range(1, steps + 1):
        if steps > 1:
            print("  step {}:".format(t))
        info = {}  # scheme -> {"top": (n_s, n_c)|None, "via": [...]}
        kids = {}  # parent scheme -> set of nested callees
        for s in compared:
            buckets = set(b for (ss, tt, b) in sima_g
                          if ss == s and tt == t)
            buckets |= set(b for (ss, tt, b) in cam_g
                           if ss == s and tt == t + 1)
            for b in sorted(buckets):
                n_s = len(sima_g.get((s, t, b), []))
                n_c = len(cam_g.get((s, t + 1, b), []))
                e = info.setdefault(s, {"top": None, "via": []})
                if b == "<toplevel>":
                    e["top"] = (n_s, n_c)
                else:
                    p = bucket_parent(b)
                    e["via"].append((b, p, n_s, n_c))
                    if p:
                        kids.setdefault(p, set()).add(s)

        printed = set()

        def emit_nested(parent, depth):
            for c in compared:
                if c not in kids.get(parent, ()):
                    continue
                new = False
                for (b, p, n_s, n_c) in info[c]["via"]:
                    if p == parent and (c, b) not in printed:
                        printed.add((c, b))
                        row(disp(c), n_s, n_c, depth)
                        new = True
                if new:
                    emit_nested(c, depth + 1)

        for s in compared:
            e = info.get(s)
            if e is None:
                continue
            if e["top"] is not None:
                printed.add((s, "<toplevel>"))
                row(disp(s), e["top"][0], e["top"][1], 0)
                emit_nested(s, 1)
            for (b, p, n_s, n_c) in e["via"]:
                if p is None and (s, b) not in printed:
                    printed.add((s, b))
                    row("{} (via {})".format(disp(s), b), n_s, n_c, 0)
        # nested rows whose parent never printed a row of its own
        for s in compared:
            if s not in info:
                continue
            for (b, p, n_s, n_c) in info[s]["via"]:
                if (s, b) not in printed:
                    printed.add((s, b))
                    row("{} (via {})".format(disp(s), b), n_s, n_c, 0)
        never = [s for s in compared if s not in info]
        if never:
            print("  (not called in compared steps: {})".format(
                ", ".join(disp(s) for s in never)))

    # --- SDF-order pairing support ----------------------------------------
    # CAM runs the full model timestep (all parameterizations), so a
    # utility scheme like geopotential_temp is called from physics_update
    # after EVERY parameterization -- but SIMA only runs the SDF suite, so
    # it calls geopotential_temp once. To match, find the CAM call that
    # follows the SDF predecessor in CAM's execution order: that's the
    # physics_update call right after the parameterization we're debugging.
    for idx, rec in enumerate(cam["hits"]):
        rec["_exec_idx"] = idx
    compared_set = set(compared)

    def sdf_pred_of(scheme):
        """Last compared scheme before scheme in the SDF."""
        last = None
        for s in suite_order:
            if s == scheme:
                return last
            if s in compared_set:
                last = s
        return None

    def find_cam_by_sdf_order(scheme, cam_step, cam_list, used):
        """Find the cam_list entry that follows the SDF predecessor in
        CAM's execution order. Returns (cam_list_idx, hit_record,
        predecessor_name) or None."""
        pred = sdf_pred_of(scheme)
        if pred is None:
            return None
        # last execution-order hit of pred at cam_step
        pred_pos = -1
        for rec in cam["hits"]:
            if rec.get("step") == cam_step and rec["scheme"] == pred:
                pred_pos = rec["_exec_idx"]
        if pred_pos < 0:
            return None
        # next hit of our scheme after pred_pos
        for idx in range(pred_pos + 1, len(cam["hits"])):
            rec = cam["hits"][idx]
            if rec.get("step") != cam_step:
                continue
            if rec["scheme"] != scheme:
                continue
            for k, c in enumerate(cam_list):
                if c["_exec_idx"] == idx and k not in used:
                    return k, c, pred
            break
        return None

    # --- stream-order comparison ----------------------------------------
    print("")
    print("comparison (execution order):")
    rep = Reporter(portable)
    n_pairs = 0
    first_srec = None
    matched_cam = {}  # group key -> set of cam occurrences already paired
    for srec in sima["hits"]:
        t = srec.get("step", 0)
        if t < 1 or t > steps:
            continue
        key = (srec["scheme"], t + 1, srec["_bucket"])
        cam_list = cam_g.get(key, [])
        n_s = len(sima_g.get((srec["scheme"], t, srec["_bucket"]), []))
        if len(cam_list) == n_s:
            # unambiguous: pair occurrence-by-occurrence
            crec = cam_list[srec["_occ"]]
        elif not cam_list:
            continue  # cam never calls it; reported above
        else:
            # Hit counts differ: CAM runs the full model timestep, so
            # utility schemes (geopotential_temp, update_dry_static_energy)
            # are called from physics_update after every parameterization;
            # SIMA only runs the SDF suite, so it calls them once.
            sch_int = (intents or {}).get(srec["scheme"], {})
            used = matched_cam.setdefault(key, set())
            crec = None

            # Strategy 1: SDF-order pairing -- find the CAM call that
            # follows the last SDF predecessor in CAM's execution order
            result = find_cam_by_sdf_order(
                srec["scheme"], t + 1, cam_list, used)
            if result is not None:
                k, crec, pred = result
                used.add(k)
                rep.note("{} [step {}] (hit {}): paired with cam "
                         "occurrence {} of {} (SDF order: follows {} "
                         "in CAM execution)".format(
                             disp(srec["scheme"]), t, srec["_occ"],
                             k, len(cam_list), pred))
            else:
                # Strategy 2: bitwise input match (fallback)
                closest = None
                for k, cand in enumerate(cam_list):
                    if k in used:
                        continue
                    same, total = count_matching_entry_args(
                        srec, cand, sima, cam, sch_int, const_ctx)
                    if total and same == total:
                        crec = cand
                        used.add(k)
                        rep.note(
                            "{} [step {}] (hit {}): paired with cam "
                            "occurrence {} of {} (bitwise input match)"
                            .format(disp(srec["scheme"]), t, srec["_occ"],
                                    k, len(cam_list)))
                        break
                    if closest is None or same > closest[1]:
                        closest = (k, same, total)
                if crec is None:
                    why = ("all {} candidates already paired".format(
                           len(cam_list)) if closest is None else
                           "closest: occurrence {}, {}/{} args match"
                           .format(closest[0], closest[1], closest[2]))
                    rep.note(
                        "{} [step {}] (hit {}): NO cam hit with matching "
                        "inputs among {} candidates ({}); not compared"
                        .format(disp(srec["scheme"]), t, srec["_occ"],
                                len(cam_list), why))
                    continue
        if first_srec is None:
            first_srec = srec
        compare_hit(srec, crec, sima, cam, const_ctx, rep, intents)
        n_pairs += 1

    # --- summary ----------------------------------------------------------
    print("")
    if rep.alignment_suspect:
        print("*** WARNING: the very first compared scheme already has "
              "input differences.")
        print("*** Timestep alignment or initial conditions are suspect; "
              "check snapshot setup")
        print("*** before believing any divergence below the first "
              "scheme.")
        if first_srec is not None:
            offset_scan(first_srec, cam_g, sima_g, sima, cam, intents,
                        const_ctx)
        print("")
    if rep.n_diffs == 0:
        n_args = sum(len(r["args"]) for r in sima["hits"])
        print("No differences found in any subroutines!")
        print("({} hit pairs compared, {} sima arg captures, byte-for-byte"
              " identical)".format(n_pairs, n_args))
        for man, role in ((cam, "cam"), (sima, "sima")):
            for note in man.get("notes", []):
                if "FAILED" in note or "failed" in note:
                    print("  {} note: {}".format(role, note))
        return True
    print("{} differing comparisons across {} hit pairs.".format(
        rep.n_diffs, n_pairs))
    clean = [s for s in compared if s not in rep.by_scheme]
    print("per-scheme: {}/{} compared schemes bit-for-bit; differing:"
          .format(len(clean), len(compared)))
    for s in compared:
        if s in rep.by_scheme:
            ne, nx = rep.by_scheme[s]
            print("  {}: {} input / {} output diffs".format(disp(s), ne, nx))
    print("Raw dumps and entry-time addresses are in the manifests for "
          "manual gdb follow-up.")
    return False


if __name__ == "__main__":
    if len(sys.argv) != 2:
        sys.exit("usage: python3 differ.py <out_dir>")
    outdir = sys.argv[1]
    with open(os.path.join(outdir, "suite.json")) as f:
        meta = json.load(f)
    ok = report(outdir, meta["schemes"], meta["steps"],
                meta.get("intents"))
    sys.exit(0 if ok else 1)
