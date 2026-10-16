#!/usr/bin/env python3
"""P2AB constraint-ablation runner (Priority B).  Spec: P2AB_ABLATION_SPEC_R0.md.  The frozen search module p2ab_search.py (and through it the frozen P2B runner) is imported UNCHANGED;
only its module-level `allowed` is replaced by the cell/fidelity rule below, so the BFS (order, dedup, cap, PURE test) is byte-identical to the frozen one.  Python: `python -P -s -B`.
  run --out DIR [--targets CSV(cohort,case_id,block)] [--systems S,..] [--grids G26,G14] [--shard k/n] [--nmax N] [--limit N] [--figdata 0|1]  COMMON + --frozen DIR
system = CELL-FID, CELL in C00,C10,C01P,C01H,C11P,C11H (locality, flip-once PATH/HIST), FID in F0,F1,F2.  grid G26=(R=2,L=6), G14=(R=1,L=4).  exit 3 = provenance failure."""
import argparse, csv, hashlib, importlib.util, json, os, platform, shutil, sys, time
from pathlib import Path
import numpy as np

FROZEN_RUNNER = "5f32741b78e4024e0124080bd50e2f0fd874bf15efd8285611d87deffa5fdcd3"
sha = lambda p: hashlib.sha256(open(p, "rb").read()).hexdigest()
CELLS = {"C00": (0, None), "C10": (1, None), "C01P": (0, "P"), "C01H": (0, "H"), "C11P": (1, "P"), "C11H": (1, "H")}
GRIDS = {"G26": (2, 6), "G14": (1, 4)}
SYSTEMS = [f"{c}-{f}" for c in ("C00", "C10") for f in ("F0", "F1", "F2")] + [f"{c}-{f}" for c in ("C01P", "C01H", "C11P", "C11H") for f in ("F0", "F2")]
ALLSYS = [f"{c}-{f}" for c in CELLS for f in ("F0", "F1", "F2")]                      # every combination (F1 for flip-once cells only for differential checks)
COLS = ["system", "cell", "fid", "grid", "R", "L", "cohort", "case_id", "block", "n_max", "window_voxels", "status", "depth_limited", "frontier_exhausted", "n_states", "level_counts", "max_level", "max_depth", "n_children_eval",
        "n_moves_root", "n_bres", "lstar", "n_flip_add", "n_flip_rem", "net_add", "net_rem", "net_outside_B", "net_raw_fg_removed", "net_added_voxel_reverted", "net_new_bg_added", "net_removed_voxel_restored",
        "n_crit_end", "n_repeated_flip_voxels", "hist_violation_steps", "rawfg_rem_steps", "p2a_added_rem_steps", "first_c11_violation_step", "first_c11_violation_type", "gate_audited_new", "seconds", "error"]
NONDET = ("seconds", "gate_audited_new")
PCOLS = ["system", "grid", "case_id", "block", "step", "type", "z", "y", "x", "nmask", "safe", "loc_ok", "fo_path_ok", "fo_hist_ok", "fid2_ok"]
M = None; CUR = {}


def boot(a):
    """SHA-check and load the frozen search module, run its provenance gate (exit 3), install the cell rule."""
    global M
    p = Path(a.frozen).resolve() / "p2ab_search.py"
    if sha(p) != FROZEN_RUNNER: print("PROVENANCE_FAIL: p2ab_search.py", file=sys.stderr); sys.exit(3)
    s = importlib.util.spec_from_file_location("p2s", str(p)); M = importlib.util.module_from_spec(s); sys.modules["p2s"] = M; s.loader.exec_module(M)
    ra = M.boot(a); M.allowed = cell_allowed; return ra


def parse(system):
    cell, fid = system.split("-"); loc, fo = CELLS[cell]; return loc, fo, fid


def lc_G(lc, bits):
    Gb = np.broadcast_to(lc.G0, (len(bits),) + lc.G0.shape).copy(); Gb[:, 1:-1, 1:-1, 1:-1] = bits.reshape((len(bits),) + lc.wshape); return Gb


def cell_allowed(system, lc, P, rawW, bits0, tbl):
    """Replacement for the frozen `allowed`: SAFE & [fidelity F2] & [locality] & [flip-once]; also accumulates the union of admissible voxels over all expanded parents (figure data)."""
    loc, fo, fid = CUR["spec"]; Gb = lc.unpack(P); S = lc.safe(Gb, tbl); bits = Gb[:, 1:-1, 1:-1, 1:-1].reshape(len(P), -1); A = S
    FID = ~((bits == 1) & (rawW[None, :] == 1)) if fid == "F2" else None
    LOC = M.incrit(lc, Gb) if loc else None
    FO = (bits == bits0[None, :]) if fo == "P" else (np.where(bits == 0, rawW[None, :] == 0, rawW[None, :] == 1) if fo == "H" else None)
    for c in (FID, LOC, FO):
        if c is not None: A = A & c
    ac = CUR["acc"]; ac["ever"] |= A.any(0); ac["ever_safe"] |= S.any(0); ac["eval"] += int(A.sum()); return A


def root_decomp(lc, rawW, bits0, tbl):
    """Per window voxel at X_T: individual condition outcomes (independent of the system)."""
    P = lc.key0[None, :]; Gb = lc.unpack(P); S = lc.safe(Gb, tbl)[0]; L = M.incrit(lc, Gb)[0]
    return dict(safe=S, loc_ok=L, fo_path_ok=np.ones_like(S), fo_hist_ok=np.where(bits0 == 0, rawW == 0, rawW == 1), fid2_ok=~((bits0 == 1) & (rawW == 1)))


def analyze(lc, rawW, bits0, tbl, path):
    """Window-space replay of a witness: per-step rule outcomes and derived path statistics (independent re-evaluation by VERIFY/ is separate)."""
    P2B = M.P2B; vi = [int(np.ravel_multi_index((z - lc.lo[0], y - lc.lo[1], x - lc.lo[2]), lc.wshape)) for _, z, y, x in path]
    cur = lc.G0[1:-1, 1:-1, 1:-1].reshape(-1).copy(); states = [cur.copy()]; fl = []
    for (t, *_), v in zip(path, vi):
        Gb = lc_G(lc, cur[None, :]); S = bool(lc.safe(Gb, tbl)[0, v]); LOC = bool(M.incrit(lc, Gb)[0, v])
        fl.append((S, LOC, bool(cur[v] == bits0[v]), bool(rawW[v] == 0) if cur[v] == 0 else bool(rawW[v] == 1), not (cur[v] == 1 and rawW[v] == 1)))
        if (cur[v] == 1) != (t == "REM"): raise RuntimeError("path type/state mismatch")
        cur[v] ^= 1; states.append(cur.copy())
    fl = np.array(fl, bool).reshape(len(path), 5); cnt = np.bincount(vi, minlength=lc.nW)
    bad = np.flatnonzero(~(fl[:, 1] & fl[:, 2])); first = int(bad[0]) if bad.size else -1
    typ = "NONE" if first < 0 else ("BOTH" if not fl[first, 1] and not fl[first, 2] else "OUT_OF_LOCALITY" if not fl[first, 1] else "REFLIP")
    rem = np.array([t == "REM" for t, *_ in path]); vv = np.array(vi)
    st = np.stack(states); Gs = lc_G(lc, st); inc = M.incrit(lc, Gs); bc = lc.b_critical(Gs); nc = lc.new_crit(Gs)
    ncrit = sum((P2B.crit_mask(Gs, p) & lc.vblk[p]).reshape(len(st), -1).sum(1) for p in P2B.PLANES)
    d = dict(n_repeated_flip_voxels=int((cnt >= 2).sum()), hist_violation_steps=int((~fl[:, 3]).sum()), rawfg_rem_steps=int((rem & (rawW[vv] == 1)).sum()),
             p2a_added_rem_steps=int((rem & (bits0[vv] == 1) & (rawW[vv] == 0)).sum()), first_c11_violation_step=first, first_c11_violation_type=typ)
    return d, fl, dict(states=st, vidx=vv, is_rem=rem, flags=fl, incrit=inc, b_critical=bc, new_crit=nc, n_crit_local=ncrit)


def run_one(ra, system, grid, co, cid, bs, B, cache, nmax, figdir):
    P2B = M.P2B; R, L = GRIDS[grid]; loc, fo, fid = parse(system); t0 = time.time(); Xp, raw, C0 = cache["Xp"], cache["raw"], cache["C0"]
    r = dict.fromkeys(COLS, ""); r.update(system=system, cell=system.split("-")[0], fid=fid, grid=grid, R=R, L=L, cohort=co, case_id=cid, block=bs, n_max=nmax); prow = []
    try:
        CUR["spec"] = (loc, fo, fid); CUR["acc"] = dict(ever=None, ever_safe=None, eval=0)
        lc = P2B.Local(Xp, B, R); CUR["acc"]["ever"] = np.zeros(lc.nW, bool); CUR["acc"]["ever_safe"] = np.zeros(lc.nW, bool)
        res = M.bfs(Xp, raw, B, R, L, nmax, "F1" if fid == "F1" else "F0"); na = P2B.audit_gate(); lcn = res["level_counts"]
        r.update(status=res["status"], window_voxels=res["window"], depth_limited=res["depth_limited"], frontier_exhausted=res["frontier_exhausted"], n_states=res["n_states"], level_counts="|".join(map(str, lcn)), max_level=max(lcn),
                 max_depth=len(lcn) - 1, n_children_eval=CUR["acc"]["eval"], n_moves_root=res["n_moves_root"], n_bres=res["n_bres"], lstar=res["lstar"], gate_audited_new=na)
        lo, hi = lc.lo, lc.hi; rawW = raw[lo[0]:hi[0] + 1, lo[1]:hi[1] + 1, lo[2]:hi[2] + 1].reshape(-1).astype(np.uint8); bits0 = lc.G0[1:-1, 1:-1, 1:-1].reshape(-1); tbl = P2B.env["TBL"]
        fig = dict(status=np.array(res["status"]), level_counts=np.array(lcn), n_states=np.array(res["n_states"]), ever_allowed=CUR["acc"]["ever"], ever_safe=CUR["acc"]["ever_safe"], coords=lc.coords, wlo=lc.lo, whi=lc.hi,
                   **{"root_" + k: v for k, v in root_decomp(lc, rawW, bits0, tbl).items()})
        if res["status"] == "PURE_FOUND":
            acc = P2B.replay_witness(Xp, raw, B, res["path"], C0); d, fl, wit = analyze(lc, rawW, bits0, tbl, res["path"])
            r.update(n_flip_add=sum(p[0] == "ADD" for p in res["path"]), n_flip_rem=sum(p[0] == "REM" for p in res["path"]), n_crit_end=acc["n_crit_end"],
                     **{k: acc[k] for k in ("net_add", "net_rem", "net_outside_B", "net_raw_fg_removed", "net_added_voxel_reverted", "net_new_bg_added", "net_removed_voxel_restored")}, **d)
            prow = [dict(system=system, grid=grid, case_id=cid, block=bs, step=i, type=p[0], z=p[1], y=p[2], x=p[3], nmask=acc["nm"][i], safe=int(fl[i, 0]), loc_ok=int(fl[i, 1]), fo_path_ok=int(fl[i, 2]),
                         fo_hist_ok=int(fl[i, 3]), fid2_ok=int(fl[i, 4])) for i, p in enumerate(res["path"])]
            fig.update({"wit_" + k: v for k, v in wit.items()}); fig["wit_path"] = np.array([(0 if p[0] == "ADD" else 1, p[1], p[2], p[3]) for p in res["path"]])
        if figdir is not None and res["status"] != "CAP_EXCEEDED":
            fd = Path(figdir) / f"{cid}__{bs.replace(',', '_').replace('@', '_')}"; fd.mkdir(parents=True, exist_ok=True); tmp = fd / f".{system}__{grid}.tmp.npz"
            np.savez_compressed(tmp, **fig); os.replace(tmp, fd / f"{system}__{grid}.npz")
    except Exception as e: r.update(status="ERROR", error=f"{type(e).__name__}: {e}"[:300])
    r["seconds"] = round(time.time() - t0, 3); return r, prow


def save_ctx(ra, cid, bs, cache, figdir):
    """Per-target context (RAW, X_T crop, P2A provenance, critical blocks of X_T) shared by all systems; written once."""
    P2B = M.P2B; Xp, raw, C0 = cache["Xp"], cache["raw"], cache["C0"]; B = P2B.parse_bk(bs); lc = P2B.Local(Xp, B, 2); sh = raw.shape
    lo = np.maximum(lc.lo - 8, 0); hi = np.minimum(lc.hi + 8, np.array(sh) - 1); sl = tuple(slice(int(a), int(b) + 1) for a, b in zip(lo, hi)); xt = Xp[1:-1, 1:-1, 1:-1]
    ed = [(int(r["step"]), 0 if r["type"] == "ADD" else 1, int(r["z"]), int(r["y"]), int(r["x"])) for r in ra.edits.get(cid, [])]; ed = np.array(ed, int).reshape(-1, 5)
    inb = np.all((ed[:, 2:] >= lo) & (ed[:, 2:] <= hi), 1) if len(ed) else np.zeros(0, bool)
    cb = np.array(sorted(C0), dtype="<U32"); fd = Path(figdir) / f"{cid}__{bs.replace(',', '_').replace('@', '_')}"; fd.mkdir(parents=True, exist_ok=True)
    np.savez_compressed(fd / "ctx.npz", crop_lo=lo, crop_hi=hi, raw=raw[sl].astype(np.uint8), xt=xt[sl].astype(np.uint8), p2a_added=((xt == 1) & (raw == 0))[sl], p2a_removed=((xt == 0) & (raw == 1))[sl],
                        p2a_edits_in_crop=ed[inb], p2a_edits_n_total=np.array(len(ed)), xt_critical_blocks=cb, target_block=np.array(bs), shape=np.array(sh))


def rows(p): return list(csv.DictReader(open(p, encoding="utf-8-sig", newline="")))


def cmd_run(a, ra):
    P2B = M.P2B; k, n = (int(v) for v in a.shard.split("/")); systems = a.systems.split(","); grids = a.grids.split(",")
    for s in systems:
        if s not in ALLSYS: raise SystemExit(f"unknown system {s}")
    for g in grids:
        if g not in GRIDS: raise SystemExit(f"unknown grid {g}")
    tg = [(r["cohort"], r["case_id"], r["block"]) for r in rows(a.targets)] if a.targets else ra.targets(("C1", "C2")); tg = tg[:a.limit] if a.limit else tg; P2B.check_cases(tg); tg = tg[k::n]
    out = Path(a.out); out.mkdir(parents=True, exist_ok=True); final = out / f"shard_{k}"; part = out / f"shard_{k}.partial"
    if final.exists(): raise SystemExit(f"write-once: {final} exists")
    (part / "_j").mkdir(parents=True, exist_ok=True); fig = part / "FIG" if a.figdata else None; BL, PA, cache, t0, nres = [], [], {}, time.time(), 0
    for co, cid, bs in tg:
        B = P2B.parse_bk(bs)
        for g in grids:
            for s in systems:
                jf = part / "_j" / f"{s}__{g}__{cid}__{bs.replace(',', '_').replace('@', '_')}.json"
                if jf.exists(): j = json.load(open(jf)); b, p = j["row"], j["prow"]; nres += 1
                else:
                    if cache.get("cid") != cid:
                        Xp, raw, ok = ra.reconstruct(cid)
                        if not ok: raise RuntimeError(f"{cid}: reconstructed X_T SHA != Run A mask_sha")
                        cache.clear(); cache.update(cid=cid, Xp=Xp, raw=raw, C0=P2B.global_crit(Xp[1:-1, 1:-1, 1:-1]), ctx=set())
                    if fig is not None and bs not in cache["ctx"]: save_ctx(ra, cid, bs, cache, fig); cache["ctx"].add(bs)
                    b, p = run_one(ra, s, g, co, cid, bs, B, cache, a.nmax, fig)
                    with open(jf.with_suffix(".tmp"), "w") as f: json.dump(dict(row=b, prow=p), f)
                    os.replace(jf.with_suffix(".tmp"), jf)
                BL.append(b); PA += p
    for nm, rr, cols in (("BLOCKS.csv", BL, COLS), ("PATHS.csv", PA, PCOLS)):
        with open(part / nm, "w", newline="") as f: w = csv.DictWriter(f, cols); w.writeheader(); w.writerows(rr)
    json.dump(dict(runner_sha=sha(__file__), frozen_runner_sha=FROZEN_RUNNER, shard=a.shard, systems=systems, grids=grids, n_targets=len(tg), n_rows=len(BL), resumed=nres, nmax=a.nmax, python=platform.python_version(), numpy=np.__version__,
                   pythonhashseed=os.environ.get("PYTHONHASHSEED"), manifest_digest=P2B.env["FZ"].manifest_digest(), seconds=round(time.time() - t0, 1)), open(part / "shard_meta.json", "w"), indent=1)
    shutil.rmtree(part / "_j"); os.rename(part, final); print(f"shard {a.shard}: {len(tg)} targets, {len(BL)} rows, {time.time() - t0:.1f}s (resumed {nres})")


def main():
    ap = argparse.ArgumentParser(); sp = ap.add_subparsers(dest="cmd", required=True); s = sp.add_parser("run")
    for n in ("p1", "vendor", "p0", "table", "betti", "p4", "runa", "p2a", "config", "frozen"): s.add_argument(f"--{n}", required=True)
    s.add_argument("--out", required=True); s.add_argument("--targets", default=""); s.add_argument("--systems", default=",".join(SYSTEMS)); s.add_argument("--grids", default="G26,G14"); s.add_argument("--shard", default="0/1")
    s.add_argument("--nmax", type=int, default=1_000_000); s.add_argument("--limit", type=int, default=0); s.add_argument("--figdata", type=int, default=1)
    a = ap.parse_args(); ra = boot(a); cmd_run(a, ra)


if __name__ == "__main__": main()
