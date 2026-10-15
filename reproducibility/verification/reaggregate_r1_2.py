#!/usr/bin/env python3
"""P2AB B1 re-aggregation (manuscript R1.1). Read-only on production data; writes JSON + Markdown to --out.

Usage:  python3 reaggregate_r1_1.py --root <P2AB_constraint_ablation> --out <new dir outside PRODUCTION/FINAL/FIGURE_DATA>
Needs numpy only when the optional witness bundles are present (section 10); without them that section is SKIPPED and reported as such.
Exit codes: 0 every check passes; 1 integrity/validation/aggregate check failed; 2 usage error, missing required input, or numpy missing while bundles are present.
Provenance of each check:  RECOMPUTED = derived here from the stored rows/arrays;  INHERITED = value stored by the frozen production run or
earlier audits and only read/compared here (SAFE decisions, |Crit| counts, master hash list). The script does not re-run any search and
is not an independent verification of the search engine or of the SAFE gate.
"""
import argparse, collections as C, csv, datetime, hashlib, json, platform, sys
from pathlib import Path

LAB = {"PURE_FOUND": "FOUND", "NO_PURE_EXHAUSTED": "EXH", "NO_PURE_DEPTH_LIMITED": "DL", "CAP_EXCEEDED": "CAP"}
CELLS = ("C00", "C10", "C01P", "C11P")
B110 = ("BraTS2021_00110", "k2xy@106,111,162")
FIG = "FIGURE_DATA/shard_1/FIG/BraTS2021_00110__k2xy_106_111_162/"
REQ_BUNDLES = [FIG + n for n in ("ctx.npz", "C00-F0__G26.npz", "C00-F1__G26.npz")]
REL = ["PRODUCTION/shard_0/BLOCKS.csv", "PRODUCTION/shard_1/BLOCKS.csv", "PRODUCTION/shard_0/PATHS.csv", "PRODUCTION/shard_1/PATHS.csv", "FINAL/REUSE_MERGE_LOG.csv"]
MASTER = "FINAL/P2AB_B1_MASTER_SHA256.txt"
REQ_BLOCK_COLS = ("system", "grid", "cohort", "case_id", "block", "status", "n_states", "n_max", "L", "max_depth", "depth_limited", "lstar", "level_counts",
                  "net_raw_fg_removed", "rawfg_rem_steps", "n_repeated_flip_voxels")
REQ_PATH_COLS = ("system", "grid", "case_id", "block", "step", "type", "z", "y", "x", "safe", "loc_ok", "fo_path_ok", "fo_hist_ok", "fid2_ok")
EXP = {("F0", "G26"): [(666, 0, 0, 14), (5, 630, 45, 0), (636, 0, 4, 40), (0, 631, 49, 0)],
       ("F0", "G14"): [(661, 0, 19, 0), (0, 641, 39, 0), (588, 0, 92, 0), (0, 642, 38, 0)],
       ("F1", "G26"): [(72, 0, 606, 2), (0, 630, 50, 0), (13, 207, 365, 95), (0, 673, 7, 0)],
       ("F1", "G14"): [(15, 0, 665, 0), (0, 641, 39, 0), (3, 460, 217, 0), (0, 677, 3, 0)],
       ("F2", "G26"): [(15, 205, 460, 0), (0, 669, 11, 0), (13, 207, 365, 95), (0, 673, 7, 0)],
       ("F2", "G14"): [(3, 458, 219, 0), (0, 677, 3, 0), (3, 460, 217, 0), (0, 677, 3, 0)]}
R, I = "RECOMPUTED", "INHERITED"


class Abort(Exception):
    def __init__(self, code, msg): super().__init__(msg); self.code = code


def sha(p): return hashlib.sha256(Path(p).read_bytes()).hexdigest()
def iv(x): return int(x) if str(x).strip() not in ("", "None") else 0


def rd(p, cols):
    with open(p, newline="", encoding="utf-8") as f:
        rr = csv.DictReader(f); miss = [c for c in cols if c not in (rr.fieldnames or [])]
        if miss: raise Abort(1, f"{p.name}: required columns missing {miss}")
        return list(rr)


def run(root, out):
    # ---- 1. input presence and master-list integrity (hard failures)
    miss = [r for r in REL + [MASTER] if not (root / r).is_file()]
    have_b = all((root / r).is_file() for r in REQ_BUNDLES)
    part_b = (not have_b) and any((root / r).is_file() for r in REQ_BUNDLES)
    if part_b: raise Abort(2, "witness bundles present only partially: supply all three or none")
    if miss: raise Abort(2, f"missing required inputs: {miss}")
    master = {}
    for ln in (root / MASTER).read_text(encoding="utf-8").splitlines():
        if ln.strip(): h, p = ln.split(None, 1); master[p.strip().lstrip("./")] = h
    inputs, bad = [], []
    for r in REL + (REQ_BUNDLES if have_b else []):
        h = sha(root / r); m = master.get(r)
        inputs.append(dict(path=r, sha256=h, bytes=(root / r).stat().st_size, master_sha256=m, master_match=(m == h)));
        if m != h: bad.append(f"{r}: {'not listed in master' if m is None else 'SHA256 differs from master'}")
    if bad: raise Abort(1, "master SHA256 verification failed: " + "; ".join(bad))
    np = None
    if have_b:
        try: import numpy as np
        except ImportError: raise Abort(2, "numpy is required when witness bundles are present")

    BL = rd(root / REL[0], REQ_BLOCK_COLS) + rd(root / REL[1], REQ_BLOCK_COLS); PA = rd(root / REL[2], REQ_PATH_COLS) + rd(root / REL[3], REQ_PATH_COLS)
    RM = rd(root / REL[4], ("system", "case_id", "block", "own_label", "prior_label"))
    chk = []
    def check(name, ok, detail="", prov=R): chk.append(dict(check=name, ok=bool(ok), provenance=prov, detail=str(detail)))

    # ---- 2. path-attribute validation (hard failure on any violation)
    K = {(r["system"], r["grid"], r["case_id"], r["block"]): r for r in BL}
    P = C.defaultdict(list); pv = []
    for i, x in enumerate(PA):
        k = (x["system"], x["grid"], x["case_id"], x["block"]); P[k].append(x)
        if x["type"] not in ("ADD", "REM"): pv.append((i, "type", x["type"]))
        for c in ("safe", "loc_ok", "fo_path_ok", "fo_hist_ok", "fid2_ok"):
            if x[c] not in ("0", "1"): pv.append((i, c, x[c]))
        for c in ("step", "z", "y", "x"):
            try: int(x[c])
            except ValueError: pv.append((i, c, x[c]))
        if x.get("safe") != "1": pv.append((i, "safe!=1", x.get("safe")))
        if k not in K: pv.append((i, "no BLOCKS row", k))
    for k, r in K.items():
        if r["status"] == "PURE_FOUND":
            st = sorted(iv(x["step"]) for x in P.get(k, []))
            if st != list(range(iv(r["lstar"]))): pv.append((k, "path steps != 0..lstar-1", st))
        elif k in P: pv.append((k, "path rows for non-FOUND block", r["status"]))
    check("PATHS validation: types, 0/1 flags, integer fields, SAFE=1 on every move, steps 0..L*-1 for every FOUND row, no orphan rows", not pv, f"{len(PA)} rows; {len(pv)} violations {pv[:3]}")
    if pv: chk[-1]["ok"] = False; raise_after = Abort(1, f"path attribute validation failed: {len(pv)} violations, e.g. {pv[:3]}")
    else: raise_after = None
    check("SAFE=1 for every recorded move (value stored by the production run; gate not re-evaluated here)", all(x["safe"] == "1" for x in PA), f"{len(PA)} moves", I)

    check("rows=19040", len(BL) == 19040, len(BL)); check("no ERROR rows", all(r["status"] != "ERROR" for r in BL)); check("unique keys", len(K) == len(BL))
    check("status values known", all(r["status"] in LAB for r in BL))
    lab = {k: LAB[r["status"]] for k, r in K.items() if r["status"] in LAB}

    # ---- 3. reuse merge
    eff = dict(lab); badr = 0
    for r in RM:
        k = (r["system"], "G26", r["case_id"], r["block"])
        if r["system"] not in ("C00-F1", "C00-F2") or lab.get(k) != "CAP" or r["own_label"] != "CAP" or r["prior_label"] not in ("FOUND", "DL"): badr += 1; continue
        eff[k] = r["prior_label"]
    check("reuse log: 655 rows, all applicable (prior labels are inherited from the earlier search runs)", len(RM) == 655 and badr == 0, f"{len(RM)} rows, {badr} not applicable")

    # ---- 4. Table 1 / S2 and cohort-wise Table S3
    def cnt(sysn, g, src, coh=None):
        c = C.Counter(src[k] for k in K if k[0] == sysn and k[1] == g and (coh is None or K[k]["cohort"] == coh)); return tuple(c.get(x, 0) for x in ("FOUND", "EXH", "DL", "CAP"))
    sysname = lambda fid, c: f"{c}-{'F2' if (fid == 'F1' and c in ('C01P', 'C11P')) else fid}"
    tab, s3 = {}, {}
    for (fid, g), want in EXP.items():
        got = [cnt(sysname(fid, c), g, eff) for c in CELLS]; tab[f"{fid}|{g}"] = got
        check(f"Table {fid} {g} == S2", [list(x) for x in got] == [list(x) for x in want], got)
        for coh in ("C1", "C2"): s3[f"{fid}|{g}|{coh}"] = [cnt(sysname(fid, c), g, eff, coh) for c in CELLS]
    for fid, g in EXP:
        a, b, t = s3[f"{fid}|{g}|C1"], s3[f"{fid}|{g}|C2"], tab[f"{fid}|{g}"]
        check(f"Table S3 {fid} {g}: C1 + C2 == pooled Table", all(tuple(x + y for x, y in zip(a[i], b[i])) == t[i] for i in range(4)))
    nb = C.Counter(r["cohort"] for r in BL if r["system"] == "C00-F0" and r["grid"] == "G26")
    cc = {c: len({r["case_id"] for r in BL if r["cohort"] == c}) for c in ("C1", "C2")}
    check("cohort sizes: C1 614 blocks / 213 cases, C2 66 / 39, 222 distinct cases", nb == {"C1": 614, "C2": 66} and cc == {"C1": 213, "C2": 39} and len({r["case_id"] for r in BL}) == 222, f"{dict(nb)} {cc}")
    own = {f"{s}|{g}": cnt(s, g, lab) for s in ("C00-F1", "C00-F2") for g in ("G26", "G14")}

    # ---- 5. DL / CAP semantics
    dl = [r for r in BL if r["status"] == "NO_PURE_DEPTH_LIMITED"]; cp = [r for r in BL if r["status"] == "CAP_EXCEEDED"]
    dl_bad = sum(not (iv(r["max_depth"]) == iv(r["L"]) and iv(r["n_states"]) < iv(r["n_max"]) and r["depth_limited"] == "1") for r in dl); cp_bad = sum(iv(r["n_states"]) != iv(r["n_max"]) for r in cp)
    check("DL rows: max_depth==L, n_states<n_max, depth_limited=1", dl_bad == 0 and len(dl) == 2760, f"{len(dl)} DL rows, {dl_bad} violations")
    check("CAP rows: n_states==n_max", cp_bad == 0 and len(cp) == 919, f"{len(cp)} CAP rows, {cp_bad} violations")

    # ---- 6. witness composition
    def comp(sysn, g):
        Rr = [r for r in BL if r["system"] == sysn and r["grid"] == g and r["status"] == "PURE_FOUND"]
        return dict(n=len(Rr), lstar2=sum(r["lstar"] == "2" for r in Rr), net_rawfg_removed_gt0=sum(iv(r["net_raw_fg_removed"]) > 0 for r in Rr),
                    rawfg_rem_step_gt0=sum(iv(r["rawfg_rem_steps"]) > 0 for r in Rr), repeated_flip_gt0=sum(iv(r["n_repeated_flip_voxels"]) > 0 for r in Rr))
    W = {f"{s}|G26": comp(s, "G26") for s in ("C00-F0", "C00-F1", "C00-F2", "C01P-F0", "C01P-F2", "C10-F0")}
    for key, want in {"C00-F0|G26": dict(n=666, lstar2=576, net_rawfg_removed_gt0=662), "C00-F1|G26": dict(n=28, net_rawfg_removed_gt0=0, rawfg_rem_step_gt0=15, repeated_flip_gt0=16),
                      "C00-F2|G26": dict(n=15, net_rawfg_removed_gt0=0, rawfg_rem_step_gt0=0, repeated_flip_gt0=2), "C01P-F0|G26": dict(n=636, net_rawfg_removed_gt0=632, repeated_flip_gt0=0),
                      "C10-F0|G26": dict(n=5, repeated_flip_gt0=5)}.items():
        check(f"witness composition {key} (counters net_raw_fg_removed etc. are stored columns)", all(W[key][f] == v for f, v in want.items()), {f: W[key][f] for f in want}, I)
    # recomputed from PATHS: RAW-foreground removal steps are fid2_ok==0 rows of REM type
    def fidrem(sysn):
        ks = [k for k, r in K.items() if k[0] == sysn and k[1] == "G26" and r["status"] == "PURE_FOUND"]
        return sum(any(x["fid2_ok"] == "0" for x in P[k]) for k in ks), len(ks)
    check("C00-F0 / C01P-F0 G26: witnesses with >=1 step recorded as fidelity-violating (fid2_ok=0) [recomputed from PATHS]", True, {s: fidrem(s) for s in ("C00-F0", "C01P-F0", "C00-F1", "C00-F2")})

    # ---- 7. C10 vs C01P pairs
    path = lambda s, k: [(x["type"], x["z"], x["y"], x["x"]) for x in sorted(P[(s, "G26") + k], key=lambda x: iv(x["step"]))]
    pairs = []
    for k, r in sorted(K.items()):
        if k[0] == "C10-F0" and k[1] == "G26" and r["status"] == "PURE_FOUND":
            t = k[2:]; o = K[("C01P-F0", "G26") + t]; x = K[("C11P-F0", "G26") + t]; pc, po = path("C10-F0", t), path("C01P-F0", t)
            pairs.append(dict(case=t[0], block=t[1], c10_lstar=iv(r["lstar"]), c10_repeated_flip_voxels=iv(r["n_repeated_flip_voxels"]), c01p_status=LAB[o["status"]], c01p_lstar=iv(o["lstar"]) if o["status"] == "PURE_FOUND" else None,
                              identical_path=(pc == po and bool(pc)), c01p_first_step_loc_ok=(P[("C01P-F0", "G26") + t][0]["loc_ok"] if po else None), c11p_status=LAB[x["status"]], c11p_n_states=iv(x["n_states"])))
    check("C10 successes: 5; 4 also C01P FOUND, 1 C01P CAP; no identical witness", len(pairs) == 5 and sum(p["c01p_status"] == "FOUND" for p in pairs) == 4 and not any(p["identical_path"] for p in pairs), pairs)
    check("C10&C01P targets: C11P is DL", all(p["c11p_status"] == "DL" for p in pairs), [p["c11p_status"] for p in pairs])
    check("Table 2 values", [p["c10_lstar"] for p in pairs] == [6] * 5 and [p["c10_repeated_flip_voxels"] for p in pairs] == [2, 1, 1, 2, 1] and [p["c01p_lstar"] for p in pairs if p["c01p_status"] == "FOUND"] == [2] * 4
          and [p["c11p_n_states"] for p in pairs] == [84, 126, 55, 23, 157], pairs)

    # ---- 8. C00-F2 FOUND but C01P-F2 not
    diff = []
    for k, r in sorted(K.items()):
        if k[0] == "C00-F2" and k[1] == "G26" and r["status"] == "PURE_FOUND" and K[("C01P-F2",) + k[1:]]["status"] != "PURE_FOUND":
            o = K[("C01P-F2",) + k[1:]]; diff.append(dict(case=k[2], block=k[3], c00f2_lstar=iv(r["lstar"]), c01pf2_status=LAB[o["status"]], c01pf2_n_states=iv(o["n_states"]),
                                                          c00f2_path_flags=[(x["type"], x["z"], x["y"], x["x"], {c: x[c] for c in ("safe", "loc_ok", "fo_path_ok", "fo_hist_ok", "fid2_ok")}) for x in sorted(P[k], key=lambda x: iv(x["step"]))]))
    check("C00-F2 FOUND not in C01P-F2: 2 (01282 CAP, 01294 DL)", sorted((d["case"], d["c01pf2_status"]) for d in diff) == [("BraTS2021_01282", "CAP"), ("BraTS2021_01294", "DL")], [(d["case"], d["c01pf2_status"]) for d in diff])
    check("01294: C01P-F2 DL with 639,806 states", [(d["c01pf2_status"], d["c01pf2_n_states"]) for d in diff if d["case"] == "BraTS2021_01294"] == [("DL", 639806)])

    # ---- 9. BraTS2021_00110
    b110 = {f"{s}|{g}": dict(status=LAB[K[(s, g) + B110]["status"]], lstar=iv(K[(s, g) + B110]["lstar"]), n_states=iv(K[(s, g) + B110]["n_states"]), level_counts=K[(s, g) + B110]["level_counts"])
            for s in ("C00-F0", "C01P-F0", "C00-F1", "C00-F2", "C01P-F2", "C10-F1", "C11P-F2") for g in ("G26", "G14")}
    want = {"C00-F0|G26": ("FOUND", 2, 139), "C01P-F0|G26": ("FOUND", 2, 139), "C00-F1|G26": ("FOUND", 4, 52137), "C00-F2|G26": ("DL", 0, 24528), "C01P-F2|G26": ("DL", 0, 24528),
            "C10-F1|G26": ("EXH", 0, 1), "C11P-F2|G26": ("EXH", 0, 1), "C00-F1|G14": ("FOUND", 4, 3293)}
    got = {k: (b110[k]["status"], b110[k]["lstar"], b110[k]["n_states"]) for k in want}
    check("Table 3 values (00110); L* shown as 0 where not FOUND", all(got[k][0] == v[0] and got[k][2] == v[2] and (v[1] == 0 or got[k][1] == v[1]) for k, v in want.items()), got)
    check("00110 C00-F2 == C01P-F2 level counts (G26)", b110["C00-F2|G26"]["level_counts"] == b110["C01P-F2|G26"]["level_counts"] == "1|18|149|759|2678|6974|13949" and b110["C00-F2|G26"]["n_states"] == 24528, b110["C00-F2|G26"])
    f0path = [dict(step=iv(x["step"]) + 1, type=x["type"], voxel=(iv(x["z"]), iv(x["y"]), iv(x["x"])), safe=x["safe"], loc=x["loc_ok"], fo_path=x["fo_path_ok"], fid2=x["fid2_ok"]) for x in sorted(P[("C00-F0", "G26") + B110], key=lambda x: iv(x["step"]))]
    check("00110 C00-F0 G26 path: REM(106,110,162), ADD(106,111,162)", [(p["type"], p["voxel"]) for p in f0path] == [("REM", (106, 110, 162)), ("ADD", (106, 111, 162))], f0path)

    # ---- 10. witness bundles (optional; mask windows are not part of the public release)
    s1c = None
    if have_b:
        d = root / FIG; z = np.load(d / "C00-F0__G26.npz", allow_pickle=True); c = np.load(d / "ctx.npz", allow_pickle=True); zf = np.load(d / "C00-F1__G26.npz", allow_pickle=True)
        need = {"C00-F0__G26.npz": (z, ("coords", "wit_states", "wit_n_crit_local", "wit_b_critical", "wit_vidx")), "ctx.npz": (c, ("crop_lo", "raw")), "C00-F1__G26.npz": (zf, ("root_safe", "root_loc_ok"))}
        nk = [f"{f}:{k}" for f, (a, ks) in need.items() for k in ks if k not in a.files]
        if nk: raise Abort(1, f"required bundle arrays missing: {nk}")
        co = z["coords"] - c["crop_lo"][None, :]; raw = c["raw"][co[:, 0], co[:, 1], co[:, 2]].astype(int); st = z["wit_states"].astype(int); rawfg = (raw == 1) & (st[0] == 1)
        nrfr = [int(((st[i] == 0) & rawfg).sum()) for i in (1, 2)]
        s1c = dict(crit_after_stored=[int(v) for v in z["wit_n_crit_local"][1:]], nrfr_after_recomputed=nrfr, b_critical_states=[bool(v) for v in z["wit_b_critical"]],
                   sha256=dict(ctx=sha(d / "ctx.npz"), bundle=sha(d / "C00-F0__G26.npz")))
        sf, lf = zf["root_safe"].astype(bool), zf["root_loc_ok"].astype(bool)
        check("00110 at X_T: window 180 voxels, 41 SAFE, 4 voxels of B, none of B SAFE (SAFE flags are stored by the production run; counted here)", len(sf) == 180 and sf.sum() == 41 and lf.sum() == 4 and not (sf & lf).any(),
              dict(window=len(sf), safe=int(sf.sum()), loc=int(lf.sum()), safe_and_loc=int((sf & lf).sum())), I)
        check("00110 C00-F0 witness: NRFR after steps = [1,1] (recomputed from stored states and RAW)", nrfr == [1, 1], s1c)
        check("00110 C00-F0 witness: |Crit| after steps = [1,0] (stored by the production run)", s1c["crit_after_stored"] == [1, 0], s1c, I)

    else:
        chk.append(dict(check="section 10 (witness-bundle checks of BraTS2021_00110): SKIPPED, bundles not supplied", ok=True, provenance="SKIPPED", detail="mask windows are not distributed"))
    if raise_after: raise raise_after

    res = dict(generated_utc=datetime.datetime.now(datetime.timezone.utc).isoformat(timespec="seconds"), command=" ".join(sys.argv), python=platform.python_version(), numpy=(np.__version__ if np else None), bundles_checked=have_b, inputs=inputs,
               provenance_note="RECOMPUTED = derived by this script from stored rows/arrays; INHERITED = stored by the frozen production run or earlier audits and only read here. No search was re-run; not an independent verification.",
               row_selection=dict(table1="all 19,040 BLOCKS rows; effective label = own label, except C00-F1/C00-F2 G26 rows listed in REUSE_MERGE_LOG with own_label CAP take prior_label; F1 cells for C01P/C11P use the F2 rows (spec D4)",
                                  table_s3="same rows split by BLOCKS.cohort (C1/C2)", witnesses="BLOCKS rows PURE_FOUND (own run), per (system, grid)", pairs="C10-F0 G26 PURE_FOUND rows joined to C01P-F0, C11P-F0 on (case_id, block)",
                                  diff="C00-F2 G26 PURE_FOUND rows whose C01P-F2 G26 row is not PURE_FOUND", block110="case BraTS2021_00110, block k2xy@106,111,162"),
               table_effective=tab, table_s3_by_cohort=s3, own_run_counts_C00_F1_F2=own, witness_composition=W, c10_vs_c01p=pairs, c00f2_found_c01pf2_not=diff, block_00110=b110, c00_f0_00110_path=f0path, c00_f0_00110_extra=s1c,
               checks=chk, all_pass=all(x["ok"] for x in chk))
    out.mkdir(parents=True, exist_ok=True)
    (out / "REAGG_R1_2.json").write_text(json.dumps(res, indent=1, ensure_ascii=False, default=str), encoding="utf-8")
    md = ["# REAGG_R1_2 result", f"generated {res['generated_utc']} · command `{res['command']}` · python {res['python']} · numpy {res['numpy']}", "", res["provenance_note"], "", "## Inputs (SHA256, all equal to master list)"] + \
         [f"- `{i['path']}` {i['sha256']} ({i['bytes']} B)" for i in inputs] + ["", "## Checks"] + [f"- {'SKIP' if x['provenance'] == 'SKIPPED' else ('PASS' if x['ok'] else 'FAIL')} · [{x['provenance']}] {x['check']}" for x in chk]
    (out / "REAGG_R1_2.md").write_text("\n".join(md) + "\n", encoding="utf-8")
    print("\n".join(md)); print("ALL_PASS" if res["all_pass"] else "SOME_FAIL"); return 0 if res["all_pass"] else 1


def main():
    ap = argparse.ArgumentParser(); ap.add_argument("--root", required=True); ap.add_argument("--out", required=True); a = ap.parse_args()
    root, out = Path(a.root).resolve(), Path(a.out).resolve()
    if not root.is_dir(): print(f"ERROR: root not found: {root}", file=sys.stderr); return 2
    if any(p in out.parents or p == out for p in (root / "PRODUCTION", root / "FINAL", root / "FIGURE_DATA")): print("ERROR: --out must not be inside an input directory", file=sys.stderr); return 2
    try: return run(root, out)
    except Abort as e: print(f"ABORT({e.code}): {e}", file=sys.stderr); return e.code


if __name__ == "__main__":
    sys.exit(main())
