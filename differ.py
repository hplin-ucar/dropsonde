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
import sys
from array import array

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


def load_role(outdir, role):
    d = os.path.join(outdir, role)
    path = os.path.join(d, "manifest.json")
    if not os.path.isfile(path):
        return None
    man = json.load(open(path))
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
    extra = " (ignoring {} untouched in both models)".format(ignored) \
        if ignored else ""
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


class Reporter(object):
    def __init__(self):
        self.n_diffs = 0
        self.first_printed = False
        self.alignment_suspect = False
        self.first_pair_seen = False
        self.cur_scheme = None
        self.by_scheme = {}  # scheme -> [entry_diffs, exit_diffs]

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

    def note(self, text):
        print("  [note] " + text)


def kdesc(info):
    why = info.get("why")
    k = info.get("kind")
    return "{} ({})".format(k, why) if why else str(k)


def _compare_arg(label, hit, phase, arg, sinfo, cinfo, sman, cman,
                 const_ctx, rep, intent=None):
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
            for ci, sj, cn, sn in pairs:
                ss = sdata[sj * chunk:(sj + 1) * chunk]
                cs = cdata[ci * chunk:(ci + 1) * chunk]
                txt = array_diff_text(ss, cs, dt, sext[:-1],
                                      sinfo["los"][:-1])
                if txt:
                    lines.append("{} <-> {}: {}".format(cn, sn, txt))
            if lines:
                rep.diff(label, hit, phase, arg,
                         "constituent-indexed array, {} mapped species "
                         "differ".format(len(lines)), lines)
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
                    rep.note("{} (hit {}): arg {}: exit differs only in "
                             "{} elements untouched by the scheme in both "
                             "models (partially-defined intent(out)); "
                             "ignored".format(label, hit, arg, ignored))
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
    label = srec["scheme"]
    if bucket != "<toplevel>":
        label += " via " + bucket
    label += " [step {}]".format(srec.get("step", "?"))
    hit = srec.get("_occ", srec["hit"])
    sch_int = (intents or {}).get(srec["scheme"], {})
    rep.cur_scheme = srec["scheme"]
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
            if phase == "exit" and arg in entry_diffed:
                suppressed.append(arg)
                continue
            pre = rep.n_diffs
            _compare_arg(label, hit, phase, arg, sinfo, cinfo, sman, cman,
                         const_ctx, rep, it)
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


def report(outdir, suite_order, steps, intents=None):
    print("")
    print("=" * 64)
    print("dropsonde report  ({})".format(outdir))
    print("=" * 64)

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
    print("")
    print("alignment (sima step t <=> cam step t+1; nested calls "
          "bucketed by caller):")
    for s in compared:
        per_bucket = {}
        problems = []
        for t in range(1, steps + 1):
            buckets = set(b for (ss, tt, b) in sima_g
                          if ss == s and tt == t)
            buckets |= set(b for (ss, tt, b) in cam_g
                           if ss == s and tt == t + 1)
            for b in sorted(buckets):
                n_s = len(sima_g.get((s, t, b), []))
                n_c = len(cam_g.get((s, t + 1, b), []))
                per_bucket[b] = per_bucket.get(b, 0) + min(n_s, n_c)
                if n_s != n_c:
                    how = ("paired by bitwise input match"
                           if n_s and n_c else "extra hits not compared")
                    problems.append(
                        "    MISMATCH step {} [{}]: sima {} vs cam {} "
                        "hits ({})".format(t, b, n_s, n_c, how))
        if not per_bucket:
            print("  {}: never called within compared steps".format(s))
            continue
        desc = ", ".join("{} x{}".format(b, n)
                         for b, n in sorted(per_bucket.items()))
        print("  {}: {}".format(s, desc))
        for p in problems:
            print(p)

    # --- stream-order comparison ----------------------------------------
    print("")
    print("comparison (execution order):")
    rep = Reporter()
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
            # hit counts differ (e.g. CAM's host calls geopotential_temp
            # from many sites): pair with the cam hit whose comparable
            # inputs ALL match bitwise, if there is one.
            sch_int = (intents or {}).get(srec["scheme"], {})
            used = matched_cam.setdefault(key, set())
            crec = None
            closest = None
            for k, cand in enumerate(cam_list):
                if k in used:
                    continue
                same, total = count_matching_entry_args(
                    srec, cand, sima, cam, sch_int, const_ctx)
                if total and same == total:
                    crec = cand
                    used.add(k)
                    rep.note("{} [step {}] (hit {}): paired with cam "
                             "occurrence {} of {} (bitwise input match)"
                             .format(srec["scheme"], t, srec["_occ"],
                                     k, len(cam_list)))
                    break
                if closest is None or same > closest[1]:
                    closest = (k, same, total)
            if crec is None:
                why = ("all {} candidates already paired".format(
                       len(cam_list)) if closest is None else
                       "closest: occurrence {}, {}/{} args match".format(
                           closest[0], closest[1], closest[2]))
                rep.note("{} [step {}] (hit {}): NO cam hit with matching "
                         "inputs among {} candidates ({}); not compared"
                         .format(srec["scheme"], t, srec["_occ"],
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
            print("  {}: {} input / {} output diffs".format(s, ne, nx))
    print("Raw dumps and entry-time addresses are in the manifests for "
          "manual gdb follow-up.")
    return False


if __name__ == "__main__":
    if len(sys.argv) != 2:
        sys.exit("usage: python3 differ.py <out_dir>")
    outdir = sys.argv[1]
    meta = json.load(open(os.path.join(outdir, "suite.json")))
    ok = report(outdir, meta["schemes"], meta["steps"],
                meta.get("intents"))
    sys.exit(0 if ok else 1)
