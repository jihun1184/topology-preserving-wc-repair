#!/usr/bin/env python3
"""P2AB completion engine FS: exact depth-limited BFS for the fidelity systems F1/F2 with a two-level lookahead (stored BFS levels 0..L-2, levels L-1 and L generated but not stored).
Definitions and the soundness argument: P2AB_COMPLETION_SPEC_R0.md.  The frozen runner p2ab_search.py (and through it the frozen P2B runner) is imported UNCHANGED.
  run --rows ROWS.csv --out DIR [--shard k/n] [--nstore N] [--max-children N]  COMMON(--p1 --vendor --p0 --table --betti --p4 --runa --p2a --config)
ROWS.csv columns: system,case_id,block,tier (system in F1,F2; tier A/B/C -> (R,L) from the frozen config).  exit 3 = provenance failure.  Python: `python -P -s -B`."""
import argparse, csv, hashlib, importlib.util, json, os, platform, shutil, sys, time
from pathlib import Path
import numpy as np

FROZEN_RUNNER = "5f32741b78e4024e0124080bd50e2f0fd874bf15efd8285611d87deffa5fdcd3"      # p2ab_search.py SHA (P2AB_SEARCH_PRE_RUN_SHA256.txt)
sha = lambda p: hashlib.sha256(open(p, "rb").read()).hexdigest()
M = None                                                                                  # frozen p2ab_search module
COLS = ["system", "cohort", "case_id", "block", "tier", "R", "L", "n_store_max", "window_voxels", "status", "n_states_stored", "level_counts_stored", "n_children_eval", "n_children_surv", "lstar", "n_flip_add", "n_flip_rem",
        "net_add", "net_rem", "net_outside_B", "net_raw_fg_removed", "net_added_voxel_reverted", "net_new_bg_added", "net_removed_voxel_restored", "n_crit_start", "n_crit_end", "n_resolved_total",
        "gate_audited_new", "seconds", "error"]
NONDET = ("seconds", "gate_audited_new")
PCOLS = ["system", "case_id", "block", "tier", "step", "type", "z", "y", "x", "nmask"]


def boot(a):
    """Load the frozen search module by path (SHA-checked), then run its provenance gate (exit 3 on failure)."""
    global M
    here = Path(a.frozen).resolve(); p = here / "p2ab_search.py"
    if sha(p) != FROZEN_RUNNER: print("PROVENANCE_FAIL: p2ab_search.py", file=sys.stderr); sys.exit(3)
    s = importlib.util.spec_from_file_location("p2s", str(p)); M = importlib.util.module_from_spec(s); sys.modules["p2s"] = M; s.loader.exec_module(M)
    return M.boot(a)


def bad_info(lc, Gb, bcrit):
    """Bad blocks of every state = {B if critical} + {critical blocks that X_T does not have}  (PURE <=> none).  Returns (nbad[m], cnt[m,nW]) with cnt = #bad blocks containing each window voxel."""
    P2B = M.P2B; m = len(Gb); cnt = np.zeros(Gb.shape, np.int8); nbad = bcrit.astype(np.int16)
    for p in P2B.PLANES:
        cm = P2B.crit_mask(Gb, p) & lc.vblk[p] & ~lc.c0[p]; nbad += cm.reshape(m, -1).sum(1).astype(np.int16)
        for sl in M.corner_slices(p, Gb.shape[-3:]): cnt[(slice(None),) + sl] += cm
    for z, y, x in lc.bidx: cnt[:, z, y, x] += bcrit
    return nbad, cnt[:, 1:-1, 1:-1, 1:-1].reshape(m, -1)


def bfs_fs(Xp, raw, B, R, L, nstore, system, max_children):
    """Same decision problem and discovery order as p2ab_search.bfs for F1/F2 at depth limit L, but levels L-1 and L are never stored (see spec section 3)."""
    P2B = M.P2B; tbl = P2B.env["TBL"]; lc = P2B.Local(Xp, B, R); nb, nW = lc.nb, lc.nW; CH = P2B.CH
    if L < 3: raise ValueError("L >= 3 required")
    if bool(lc.b_critical(lc.G0[None])[0]) is not True: raise RuntimeError("target block is not critical in the start state")
    lo, hi = lc.lo, lc.hi; rawW = raw[lo[0]:hi[0] + 1, lo[1]:hi[1] + 1, lo[2]:hi[2] + 1].reshape(-1).astype(np.uint8); bits0 = lc.G0[1:-1, 1:-1, 1:-1].reshape(-1)
    rawfg = (rawW == 1) & (bits0 == 1); f1 = system == "F1"
    vis = {lc.key0.tobytes()}; lk = [lc.key0[None, :]]; lp = [np.array([-1])]; lv = [np.array([-1])]; count = 1; cap = False; found = None
    for d in range(L - 2):                                                                    # stored levels 1..L-2: identical to the frozen BFS (dedup, cap, PURE test, first-PURE stop)
        par = lk[d]; nk, npar, nvox = [], [], []; n_lvl = 0
        for s in range(0, len(par), CH):
            P = par[s:s + CH]; pi, vi = np.nonzero(M.allowed(system, lc, P, rawW, bits0, tbl))
            if not pi.size: continue
            ck = P[pi] ^ lc.UNIT[vi]; kv = np.ascontiguousarray(ck).view(np.dtype((np.void, nb))).ravel(); first = np.sort(np.unique(kv, return_index=True)[1]); keep = []
            for j, b in zip(first.tolist(), kv[first].tolist()):
                if b in vis: continue
                vis.add(b); keep.append(j)
            if count + len(keep) > nstore: keep = keep[:nstore - count]; cap = True
            if keep:
                K = ck[keep]; Gn = lc.unpack(K); res = ~lc.b_critical(Gn); pi_, vi_ = pi[keep], vi[keep]
                if res.any():
                    nc = np.zeros(len(K), bool); nc[res] = lc.new_crit(Gn[res]); pure = res & ~nc
                    nr = ((Gn[:, 1:-1, 1:-1, 1:-1].reshape(len(K), -1) == 0) & rawfg[None, :]).sum(1); win = pure & (nr == 0) if f1 else pure
                    if win.any(): j0 = int(np.argmax(win)); found = (d + 1, n_lvl + j0, K[j0].copy()); K, pi_, vi_ = K[:j0 + 1], pi_[:j0 + 1], vi_[:j0 + 1]
                nk.append(K); npar.append(s + pi_); nvox.append(vi_); n_lvl += len(K); count += len(K)
            if found is not None or cap: break
        lk.append(np.concatenate(nk) if nk else np.zeros((0, nb), np.uint8)); lp.append(np.concatenate(npar) if npar else np.zeros(0, int)); lv.append(np.concatenate(nvox) if nvox else np.zeros(0, int))
        if found is not None or cap or not nk: break
    out = dict(n_states_stored=count, window=nW, level_counts=[len(k) for k in lk], n_children_eval=0, path=None, lstar=-1)

    def trace(level, idx):                                                                    # voxel steps from X_T to stored state (level, idx)
        st = []
        for dd in range(level, 0, -1): st.append(int(lv[dd][idx])); idx = int(lp[dd][idx])
        return st[::-1]
    steps = None
    if found is not None: out["status"] = "PURE_FOUND"; out["lstar"] = found[0]; steps = trace(found[0], found[1])
    elif cap: out["status"] = "CAP_EXCEEDED"
    elif len(lk) - 1 < L - 2 or len(lk[-1]) == 0: out["status"] = "NO_PURE_EXHAUSTED"                  # frontier empty at a stored level
    else:                                                                                     # lookahead over the children (depth L-1) and grandchildren (depth L) of the level L-2 states
        par = lk[L - 2]; b1 = b2 = None; nch = 0; nsv = 0; over = False; P0 = P2B.env["P0"]; OFFS = [tuple(o) for o in P0.OFFS]; okmap = {o: k for k, o in enumerate(OFFS)}
        Bw = np.array([int(np.ravel_multi_index((z - 1, y - 1, x - 1), lc.wshape)) for z, y, x in lc.bidx]); Wc = np.stack(np.unravel_index(np.arange(nW), lc.wshape), 1); inB = np.zeros(nW, bool); inB[Bw] = True
        delta = np.zeros((nW, len(Bw)), np.uint32)                                            # delta[v, j]: bit of voxel v in the 26-neighbourhood mask of B voxel j (0 if not a neighbour)
        for j, bw in enumerate(Bw.tolist()):
            for i, off in enumerate((Wc - Wc[bw]).tolist()):
                k = okmap.get(tuple(off))
                if k is not None: delta[i, j] = 1 << k
        for s in range(0, len(par), CH):
            P = par[s:s + CH]; pi, vi = np.nonzero(M.allowed(system, lc, P, rawW, bits0, tbl))
            if not pi.size: continue
            nch += len(pi); Gp = lc.unpack(P); pbits = Gp[:, 1:-1, 1:-1, 1:-1].reshape(len(P), -1)
            nmB = np.zeros((len(P), len(Bw)), np.uint32)                                      # neighbourhood masks of the B voxels in every parent
            for j, (z, y, x) in enumerate(lc.bidx):
                for k, (da, db, dc) in enumerate(OFFS): nmB[:, j] |= Gp[:, z + da, y + db, x + dc].astype(np.uint32) << np.uint32(k)
            # pre-filter (necessary condition, spec section 3): a child whose B is still critical (parent B critical, flip outside B) needs a flippable B voxel for the last flip; F1: NRFR(child) <= 1
            nmc = nmB[pi] ^ delta[vi]; P2B.NMFLAG[np.unique(nmc)] |= 1; sf = (tbl[nmc] & 7) == 7
            if system == "F2": sf &= ~((pbits[:, Bw][pi] == 1) & (rawW[Bw][None, :] == 1))
            kc = inB[vi] | sf.any(1) | ~lc.b_critical(Gp)[pi]                                 # the B-based test applies only to parents whose B is still critical (B may already be resolved)
            if f1:
                nrp = ((pbits == 0) & rawfg[None, :]).sum(1); dn = np.where(rawfg[vi], 2 * pbits[pi, vi].astype(np.int64) - 1, 0); kc &= (nrp[pi] + dn) <= 1
            pi, vi = pi[kc], vi[kc]
            if not pi.size: continue
            nsv += len(pi); ck = P[pi] ^ lc.UNIT[vi]; kv = np.ascontiguousarray(ck).view(np.dtype((np.void, nb))).ravel(); fi = np.sort(np.unique(kv, return_index=True)[1]); ck, pi, vi = ck[fi], pi[fi], vi[fi]
            Gc = lc.unpack(ck); bc = lc.b_critical(Gc); nbad, cnt = bad_info(lc, Gc, bc); bits = Gc[:, 1:-1, 1:-1, 1:-1].reshape(len(ck), -1)
            nr = ((bits == 0) & rawfg[None, :]).sum(1); pure = nbad == 0; win = pure & (nr == 0) if f1 else pure
            if win.any(): j0 = int(np.argmax(win)); b1 = (s + int(pi[j0]), int(vi[j0])); break  # PURE at depth L-1 (earliest possible here)
            if b2 is None:                                                                    # candidates for the last flip: must lie in every bad block; F1: must repair NRFR when NRFR(child)=1
                needh = nbad > 0; hit = (cnt == nbad[:, None]) & needh[:, None]; sv = (~needh) | hit.any(1)
                if f1: sv &= nr <= 1
                ix = np.flatnonzero(sv)
                if ix.size:
                    Ks = ck[ix]; C = M.allowed(system, lc, Ks, rawW, bits0, tbl) & np.where(needh[ix, None], hit[ix], True)
                    if f1: C &= np.where(nr[ix, None] == 1, (bits[ix] == 0) & rawfg[None, :], True)
                    qi, qv = np.nonzero(C)
                    if qi.size:
                        Gg = lc.unpack(Ks[qi] ^ lc.UNIT[qv]); rs = ~lc.b_critical(Gg); nc = np.zeros(len(qi), bool); nc[rs] = lc.new_crit(Gg[rs]); p2 = rs & ~nc
                        if f1: p2 &= ((Gg[:, 1:-1, 1:-1, 1:-1].reshape(len(qi), -1) == 0) & rawfg[None, :]).sum(1) == 0
                        if p2.any(): j = int(np.argmax(p2)); ci = int(ix[qi[j]]); b2 = (s + int(pi[ci]), int(vi[ci]), int(qv[j]))
            if nch > max_children: over = True; break
        out["n_children_eval"] = nch; out["n_children_surv"] = nsv
        if b1 is not None: out["status"] = "PURE_FOUND"; out["lstar"] = L - 1; steps = trace(L - 2, b1[0]) + [b1[1]]
        elif b2 is not None: out["status"] = "PURE_FOUND"; out["lstar"] = L; steps = trace(L - 2, b2[0]) + [b2[1], b2[2]]
        elif over: out["status"] = "CAP_BUDGET"
        else: out["status"] = "NO_PURE_COMPLETE"                                              # no PURE state within <= L flips in the window (exhausted-vs-depth-limited not distinguished here)
    if steps is not None:
        cur = lc.G0[1:-1, 1:-1, 1:-1].reshape(-1).copy(); path = []
        for v in steps: z, y, x = lc.coords[v].tolist(); path.append(("REM" if cur[v] else "ADD", z, y, x)); cur[v] ^= 1
        out["path"] = path
    return out


def one_row(ra, system, cid, bs, tier, cache, nstore, max_children):
    P2B = M.P2B; cfg = P2B.env["cfg"]; t = cfg["tiers"][tier]; t0 = time.time(); B = P2B.parse_bk(bs)
    if cache.get("cid") != cid:
        Xp, raw, ok = ra.reconstruct(cid)
        if not ok: raise RuntimeError(f"{cid}: reconstructed X_T SHA != Run A mask_sha")
        cache.clear(); cache.update(cid=cid, Xp=Xp, raw=raw, C0=P2B.global_crit(Xp[1:-1, 1:-1, 1:-1]))
    Xp, raw, C0 = cache["Xp"], cache["raw"], cache["C0"]; prow = []
    r = dict.fromkeys(COLS, ""); r.update(system=system, case_id=cid, block=bs, tier=tier, R=t["R"], L=t["L"], n_store_max=nstore, n_crit_start=len(C0))
    try:
        res = bfs_fs(Xp, raw, B, t["R"], t["L"], nstore, system, max_children); na = P2B.audit_gate()
        r.update(status=res["status"], window_voxels=res["window"], n_states_stored=res["n_states_stored"], level_counts_stored="|".join(map(str, res["level_counts"])), n_children_eval=res["n_children_eval"], n_children_surv=res.get("n_children_surv", ""),
                 lstar=res["lstar"], gate_audited_new=na)
        if res["status"] == "PURE_FOUND":
            acc = P2B.replay_witness(Xp, raw, B, res["path"], C0)
            r.update(n_flip_add=sum(p[0] == "ADD" for p in res["path"]), n_flip_rem=sum(p[0] == "REM" for p in res["path"]),
                     **{k: acc[k] for k in ("net_add", "net_rem", "net_outside_B", "net_raw_fg_removed", "net_added_voxel_reverted", "net_new_bg_added", "net_removed_voxel_restored", "n_crit_end", "n_resolved_total")})
            prow = [dict(system=system, case_id=cid, block=bs, tier=tier, step=i, type=p[0], z=p[1], y=p[2], x=p[3], nmask=acc["nm"][i]) for i, p in enumerate(res["path"])]
    except Exception as e: r.update(status="ERROR", error=f"{type(e).__name__}: {e}"[:300])
    r["seconds"] = round(time.time() - t0, 3); return r, prow


def main():
    ap = argparse.ArgumentParser(); sp = ap.add_subparsers(dest="cmd", required=True); s = sp.add_parser("run")
    for n in ("p1", "vendor", "p0", "table", "betti", "p4", "runa", "p2a", "config", "frozen"): s.add_argument(f"--{n}", required=True)
    s.add_argument("--rows", required=True); s.add_argument("--out", required=True); s.add_argument("--shard", default="0/1"); s.add_argument("--nstore", type=int, default=8_000_000)
    s.add_argument("--max-children", type=int, default=4_000_000_000); s.add_argument("--limit", type=int, default=0)
    a = ap.parse_args(); ra = boot(a); k, n = (int(v) for v in a.shard.split("/"))
    rows = [r for r in csv.DictReader(open(a.rows, encoding="utf-8"))]; rows = rows[:a.limit] if a.limit else rows; rows = rows[k::n]
    M.P2B.check_cases([(r.get("cohort", ""), r["case_id"], r["block"]) for r in rows])
    out = Path(a.out); out.mkdir(parents=True, exist_ok=True); final = out / f"shard_{k}"; part = out / f"shard_{k}.partial"
    if final.exists(): raise SystemExit(f"write-once: {final} exists")
    part.mkdir(parents=True, exist_ok=True); jd = part / "_j"; jd.mkdir(exist_ok=True); BL, PA, cache, t0 = [], [], {}, time.time(); nres = 0
    for r in rows:
        jf = jd / f"{r['system']}__{r['case_id']}__{r['block'].replace(',', '_').replace('@', '_')}__{r['tier']}.json"
        if jf.exists(): j = json.load(open(jf)); b, p = j["row"], j["prow"]; nres += 1
        else:
            b, p = one_row(ra, r["system"], r["case_id"], r["block"], r["tier"], cache, a.nstore, a.max_children); b["cohort"] = r.get("cohort", "")
            json.dump(dict(row=b, prow=p), open(jf.with_suffix(".tmp"), "w")); os.replace(jf.with_suffix(".tmp"), jf)
        BL.append(b); PA += p
    for nm, rr, cols in (("BLOCKS.csv", BL, COLS), ("PATHS.csv", PA, PCOLS)):
        with open(part / nm, "w", newline="") as f: w = csv.DictWriter(f, cols); w.writeheader(); w.writerows(rr)
    json.dump(dict(engine_sha=sha(__file__), frozen_runner_sha=FROZEN_RUNNER, shard=a.shard, n_rows=len(rows), resumed=nres, nstore=a.nstore, max_children=a.max_children, python=platform.python_version(),
                   numpy=np.__version__, pythonhashseed=os.environ.get("PYTHONHASHSEED"), manifest_digest=M.P2B.env["FZ"].manifest_digest(), seconds=round(time.time() - t0, 1)), open(part / "shard_meta.json", "w"), indent=1)
    shutil.rmtree(jd); os.rename(part, final); print(f"shard {a.shard}: {len(rows)} rows, {time.time() - t0:.1f}s (resumed {nres})")


if __name__ == "__main__": main()
