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

# Not compared: uninitialized at entry on the CAM side (callers pass
# stack garbage), and any nonzero errflg aborts the run instantly anyway.
SKIP_ARGS = {"errmsg", "errflg"}


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

    def diff(self, scheme, hit, phase, arg, text, extra_lines=None):
        self.n_diffs += 1
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


def compare_hit(srec, crec, sman, cman, const_ctx, rep):
    """Compare one aligned hit pair; entry args first, then exit."""
    bucket = srec.get("_bucket", "<toplevel>")
    label = srec["scheme"]
    if bucket != "<toplevel>":
        label += " via " + bucket
    label += " [step {}]".format(srec.get("step", "?"))
    hit = srec.get("_occ", srec["hit"])
    cam_names, sima_names, pairs = const_ctx
    for phase, fkey, vkey in (("entry", "entry_file", "entry_value"),
                              ("exit", "exit_file", "exit_value")):
        for arg, sinfo in srec["args"].items():
            if arg in SKIP_ARGS:
                continue
            cinfo = crec["args"].get(arg)
            if cinfo is None:
                if phase == "entry":
                    rep.note("{} (hit {}): arg {} absent on CAM side".format(
                        label, hit, arg))
                continue
            kind = sinfo.get("kind")
            if kind != cinfo.get("kind"):
                if phase == "entry":
                    rep.diff(label, hit, phase, arg,
                             "argument kind mismatch: sima {} vs cam "
                             "{}".format(kdesc(sinfo), kdesc(cinfo)))
                continue
            if kind == "error":
                if phase == "entry":
                    rep.note("{} (hit {}): arg {} capture errored on both "
                             "sides: sima {} / cam {}".format(
                                 label, hit, arg, kdesc(sinfo),
                                 kdesc(cinfo)))
                continue

            if kind == "array":
                sfile = sinfo.get(fkey)
                cfile = cinfo.get(fkey)
                if not sfile or not cfile:
                    continue  # incomplete capture (noted in manifest)
                dt = sinfo["dtype"]
                if dt != cinfo.get("dtype"):
                    if phase == "entry":
                        rep.diff(label, hit, phase, arg,
                                 "dtype mismatch: sima {} vs cam {}".format(
                                     dt, cinfo.get("dtype")))
                    continue
                sext = sinfo["extents"]
                cext = cinfo["extents"]
                sdata = blob(sman, sfile)
                cdata = blob(cman, cfile)

                const_arr = (len(sext) == 3 and len(cext) == 3 and
                             sima_names and cam_names and
                             sext[2] == len(sima_names) and
                             cext[2] == len(cam_names) and
                             sext[:2] == cext[:2])
                if const_arr:
                    chunk = sext[0] * sext[1] * ELSIZE[dt]
                    lines = []
                    for ci, sj, cn, sn in pairs:
                        ss = sdata[sj * chunk:(sj + 1) * chunk]
                        cs = cdata[ci * chunk:(ci + 1) * chunk]
                        txt = array_diff_text(ss, cs, dt, sext[:2],
                                              sinfo["los"][:2])
                        if txt:
                            lines.append("{} <-> {}: {}".format(cn, sn, txt))
                    if lines:
                        rep.diff(label, hit, phase, arg,
                                 "constituent-indexed array, {} mapped "
                                 "species differ".format(len(lines)), lines)
                    continue

                if sext != cext:
                    if phase == "entry":
                        rep.diff(label, hit, phase, arg,
                                 "shape mismatch: sima {} vs cam {}".format(
                                     sext, cext))
                    continue
                txt = array_diff_text(sdata, cdata, dt, sext, sinfo["los"])
                if txt:
                    rep.diff(label, hit, phase, arg, txt)

            elif kind in ("scalar", "char"):
                sv = sinfo.get(vkey)
                cv = cinfo.get(vkey)
                if sv is None or cv is None:
                    continue
                if sv != cv:
                    rep.diff(label, hit, phase, arg,
                             "sima={} cam={}".format(fmt_val(sv),
                                                     fmt_val(cv)))
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


def report(outdir, suite_order, steps):
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
                    problems.append(
                        "    MISMATCH step {} [{}]: sima {} vs cam {} "
                        "hits (extra hits not compared)".format(
                            t, b, n_s, n_c))
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
    for srec in sima["hits"]:
        t = srec.get("step", 0)
        if t < 1 or t > steps:
            continue
        key = (srec["scheme"], t + 1, srec["_bucket"])
        cam_list = cam_g.get(key, [])
        if srec["_occ"] >= len(cam_list):
            continue  # count mismatch, already reported above
        compare_hit(srec, cam_list[srec["_occ"]], sima, cam, const_ctx, rep)
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
    print("Raw dumps and entry-time addresses are in the manifests for "
          "manual gdb follow-up.")
    return False


if __name__ == "__main__":
    if len(sys.argv) != 2:
        sys.exit("usage: python3 differ.py <out_dir>")
    outdir = sys.argv[1]
    meta = json.load(open(os.path.join(outdir, "suite.json")))
    ok = report(outdir, meta["schemes"], meta["steps"])
    sys.exit(0 if ok else 1)
