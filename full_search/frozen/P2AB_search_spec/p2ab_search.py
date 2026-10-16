#!/usr/bin/env python3
"""P2A-restricted exhaustive search (A_HIST, A_PATH, B_HIST, B_PATH) and fidelity search (F0 control, F1 endpoint, F2 pathwise) on the P2B targets.  Definitions: P2AB_SEARCH_SPEC_R0.md (SHA-checked first).
The frozen P2B runner is imported UNCHANGED (provenance gate, X_T reconstruction, window/halo geometry, gate lookup, gate audit, PURE/new-critical predicates).  Python: `python -P -s -B`.
  search --system S --out DIR --shard k/n [--limit N]     (write-once shard; DIR/shard_k)       merge --system S --out DIR --nshards n
COMMON = --p1 --vendor --p0 --table --betti --p4 --runa --p2a --config (identical to the P2B runner).  exit 3 = provenance failure."""
import argparse, csv, hashlib, importlib.util, json, os, platform, shutil, sys, time
from collections import defaultdict
from itertools import product
from pathlib import Path
import numpy as np

PIN_P2B = "ed683826f6bbf3c0352a946574b11f9eb1d4c1f75b511d6aac7518e38d19c62b"
SPEC = "P2AB_SEARCH_SPEC_R0.md"
SYSTEMS = ("A_HIST", "A_PATH", "B_HIST", "B_PATH", "F0", "F1", "F2")
P2B = None
COLS = ["system", "cohort", "case_id", "block", "tier", "R", "L", "n_max", "window_voxels", "status", "depth_limited", "frontier_exhausted", "n_states", "level_counts", "n_moves_root", "n_bres", "min_depth_bres", "n_pure_seen",
        "n_pure_unfaithful_seen", "min_nrfr_pure_seen", "lstar", "n_flip_add", "n_flip_rem", "net_add", "net_rem", "net_outside_B", "net_raw_fg_removed", "net_added_voxel_reverted", "net_new_bg_added", "net_removed_voxel_restored",
        "n_crit_start", "n_crit_end", "n_resolved_total", "gate_audited_new", "seconds", "error"]
NONDET = ("seconds", "gate_audited_new")
PCOLS = ["system", "case_id", "block", "tier", "step", "type", "z", "y", "x", "nmask"]
sha = lambda p: hashlib.sha256(open(p, "rb").read()).hexdigest()
RSHA = sha(__file__); CNT = {"new": 0, "resumed": 0}


def boot(a):
    """Provenance gate (exit 3), then import the frozen P2B runner unchanged."""
    global P2B
    rn = Path(a.p2a).parent / "P2B_local_edit_search" / "P2B_search_R1.py"; here = Path(__file__).resolve().parent; bad = []
    if sha(rn) != PIN_P2B: bad.append("P2B_search_R1.py")
    if sha(here / SPEC) != open(here / "P2AB_SEARCH_SPEC_R0.sha256").read().split()[0]: bad.append(SPEC)
    if bad: print("PROVENANCE_FAIL:", bad, file=sys.stderr); sys.exit(3)
    s = importlib.util.spec_from_file_location("p2b", str(rn)); P2B = importlib.util.module_from_spec(s); sys.modules["p2b"] = P2B; s.loader.exec_module(P2B); P2B.load_env(a); return P2B.RunA(a)


def corner_slices(plane, s):
    """Index tuples (over the last three axes) of the corner views used by P2B.corners(); OR-ing a per-block mask into them marks every voxel of every flagged block."""
    f, e = slice(None), slice
    if plane == "xy": return [(f, e(0, -1), e(0, -1)), (f, e(0, -1), e(1, None)), (f, e(1, None), e(0, -1)), (f, e(1, None), e(1, None))]
    if plane == "yz": return [(e(0, -1), e(0, -1), f), (e(0, -1), e(1, None), f), (e(1, None), e(0, -1), f), (e(1, None), e(1, None), f)]
    if plane == "zx": return [(e(0, -1), f, e(0, -1)), (e(0, -1), f, e(1, None)), (e(1, None), f, e(0, -1)), (e(1, None), f, e(1, None))]
    return [(e(i, s[0] - 1 + i), e(j, s[1] - 1 + j), e(k, s[2] - 1 + k)) for i, j, k in product((0, 1), repeat=3)]


def incrit(lc, Gb):
    """bool[m, nW]: window voxel lies in a critical block of the state (blocks fully inside window+halo and inside the domain; every block through a window voxel lies there)."""
    out = np.zeros(Gb.shape, bool)
    for p in P2B.PLANES:
        cm = P2B.crit_mask(Gb, p) & lc.vblk[p]
        for sl in corner_slices(p, Gb.shape[-3:]): out[(slice(None),) + sl] |= cm
    return out[:, 1:-1, 1:-1, 1:-1].reshape(len(Gb), -1)


def allowed(system, lc, P, rawW, bits0, tbl):
    """bool[m, nW] of permitted single flips of the states in key rows P (spec sections 2 and 3)."""
    Gb = lc.unpack(P); S = lc.safe(Gb, tbl); bits = Gb[:, 1:-1, 1:-1, 1:-1].reshape(len(P), -1)
    if system in ("F0", "F1"): return S
    if system == "F2": return S & ~((bits == 1) & (rawW[None, :] == 1))                                      # pathwise fidelity: never remove a RAW foreground voxel
    rule = (bits == bits0[None, :]) if system.endswith("PATH") else np.where(bits == 0, rawW[None, :] == 0, rawW[None, :] == 1)    # flip-once: path-local / RAW (history) rule
    A = S & incrit(lc, Gb) & rule                                                                           # locality + SAFE + flip-once
    if system.startswith("B_"): A &= (bits == 0) | ~(A & (bits == 0)).any(1)[:, None]                         # phase: REM only at the fixed point (no permitted ADD)
    return A


def bfs(Xp, raw, B, R, L, nmax, system, tbl=None):
    """Exhaustive FIFO BFS over window states (same discovery order, dedup, cap and PURE test as P2B bfs).  L=None: no depth limit (restricted systems; depth <= |W| since each voxel flips once)."""
    tbl = P2B.env["TBL"] if tbl is None else tbl; lc = P2B.Local(Xp, B, R); nb, nW = lc.nb, lc.nW; Lm = nW if L is None else L
    if bool(lc.b_critical(lc.G0[None])[0]) is not True: raise RuntimeError("target block is not critical in the start state")
    lo, hi = lc.lo, lc.hi; rawW = raw[lo[0]:hi[0] + 1, lo[1]:hi[1] + 1, lo[2]:hi[2] + 1].reshape(-1).astype(np.uint8); bits0 = lc.G0[1:-1, 1:-1, 1:-1].reshape(-1)
    rawfg = (rawW == 1) & (bits0 == 1); f1 = system == "F1"; n_root = int(allowed(system, lc, lc.key0[None], rawW, bits0, tbl).sum())
    vis = {lc.key0.tobytes()}; lk = [lc.key0[None, :]]; lp = [np.array([-1])]; lv = [np.array([-1])]; count = 1; nbres = minb = 0; minb = -1; cap = False; found = None; npure = nunf = 0; minnr = -1
    for d in range(Lm):
        par = lk[d]; nk, npar, nvox = [], [], []; n_lvl = 0
        for s in range(0, len(par), P2B.CH):
            P = par[s:s + P2B.CH]; pi, vi = np.nonzero(allowed(system, lc, P, rawW, bits0, tbl))
            if not pi.size: continue
            ck = P[pi] ^ lc.UNIT[vi]; kv = np.ascontiguousarray(ck).view(np.dtype((np.void, nb))).ravel(); first = np.sort(np.unique(kv, return_index=True)[1]); keep = []
            for j, b in zip(first.tolist(), kv[first].tolist()):
                if b in vis: continue
                vis.add(b); keep.append(j)
            if count + len(keep) > nmax: keep = keep[:nmax - count]; cap = True
            if keep:
                K = ck[keep]; Gn = lc.unpack(K); res = ~lc.b_critical(Gn); pi_, vi_ = pi[keep], vi[keep]
                if res.any():
                    nc = np.zeros(len(K), bool); nc[res] = lc.new_crit(Gn[res]); pure = res & ~nc
                    nr = ((Gn[:, 1:-1, 1:-1, 1:-1].reshape(len(K), -1) == 0) & rawfg[None, :]).sum(1)         # NRFR of every new state
                    win = pure & (nr == 0) if f1 else pure
                    if win.any(): j0 = int(np.argmax(win)); found = (d + 1, n_lvl + j0, K[j0].copy()); K, res, pi_, vi_, pure, nr = K[:j0 + 1], res[:j0 + 1], pi_[:j0 + 1], vi_[:j0 + 1], pure[:j0 + 1], nr[:j0 + 1]
                    nbres += int(res.sum()); npure += int(pure.sum()); nunf += int((pure & (nr > 0)).sum())
                    if pure.any(): m_ = int(nr[pure].min()); minnr = m_ if minnr < 0 else min(minnr, m_)
                    if res.any() and minb < 0: minb = d + 1
                nk.append(K); npar.append(s + pi_); nvox.append(vi_); n_lvl += len(K); count += len(K)
            if found is not None or cap: break
        lk.append(np.concatenate(nk) if nk else np.zeros((0, nb), np.uint8)); lp.append(np.concatenate(npar) if npar else np.zeros(0, int)); lv.append(np.concatenate(nvox) if nvox else np.zeros(0, int))
        if found is not None or cap or not nk: break
    D = len(lk) - 1; out = dict(n_states=count, window=nW, level_counts=[len(k) for k in lk], n_bres=nbres, min_depth_bres=minb, n_moves_root=n_root, n_pure_seen=npure, n_pure_unfaithful_seen=nunf, min_nrfr_pure_seen=minnr,
                                path=None, lstar=-1, depth_limited=0, frontier_exhausted=0)
    if found is not None:
        depth, idx, key = found; out["status"] = "PURE_FOUND"; out["lstar"] = depth; steps = []; i = idx
        for dd in range(depth, 0, -1): steps.append(int(lv[dd][i])); i = int(lp[dd][i])
        steps.reverse(); cur = lc.G0[1:-1, 1:-1, 1:-1].reshape(-1).copy(); path = []
        for v in steps: z, y, x = lc.coords[v].tolist(); path.append(("REM" if cur[v] else "ADD", z, y, x)); cur[v] ^= 1
        out["path"] = path
        if int.from_bytes(np.packbits(cur, bitorder="little").tobytes(), "little") != int.from_bytes(key.tobytes(), "little"): raise RuntimeError("path replay != endpoint key")
    elif cap: out["status"] = "CAP_EXCEEDED"
    elif D == Lm and len(lk[D]) > 0: out["status"] = "NO_PURE_DEPTH_LIMITED"; out["depth_limited"] = 1
    else: out["status"] = "NO_PURE_EXHAUSTED"; out["frontier_exhausted"] = 1
    return out


def p2b_status(r):
    """Map to the P2B enum for routing / comparison with P2B Run A."""
    return r["status"] if r["status"] in ("PURE_FOUND", "CAP_EXCEEDED") else ("B_RESOLVED_WITH_NEW_CRIT" if r["n_bres"] else "NO_PATH_EXHAUSTIVE")


def one_tier(ra, system, co, cid, bs, B, tier, cache):
    """One (target, tier) search + witness self-check -> (row, path rows)."""
    cfg = P2B.env["cfg"]; restricted = system[:2] in ("A_", "B_")
    if cache.get("cid") != cid:
        Xp, raw, ok = ra.reconstruct(cid)
        if not ok: raise RuntimeError(f"{cid}: reconstructed X_T SHA != Run A mask_sha")
        cache.clear(); cache.update(cid=cid, Xp=Xp, raw=raw, C0=P2B.global_crit(Xp[1:-1, 1:-1, 1:-1]))
    Xp, raw, C0 = cache["Xp"], cache["raw"], cache["C0"]; prow = []; t = cfg["tiers"][tier]; t0 = time.time()
    r = dict.fromkeys(COLS, ""); r.update(system=system, cohort=co, case_id=cid, block=bs, tier=tier, R=t["R"], L="none" if restricted else t["L"], n_max=cfg["n_max"], n_crit_start=len(C0))
    try:
        res = bfs(Xp, raw, B, t["R"], None if restricted else t["L"], cfg["n_max"], system); na = P2B.audit_gate()
        r.update(status=res["status"], window_voxels=res["window"], n_states=res["n_states"], level_counts="|".join(map(str, res["level_counts"])), gate_audited_new=na, lstar=res["lstar"],
                 **{k: res[k] for k in ("depth_limited", "frontier_exhausted", "n_moves_root", "n_bres", "min_depth_bres", "n_pure_seen", "n_pure_unfaithful_seen", "min_nrfr_pure_seen")})
        if res["status"] == "PURE_FOUND":
            acc = P2B.replay_witness(Xp, raw, B, res["path"], C0)
            r.update(n_flip_add=sum(p[0] == "ADD" for p in res["path"]), n_flip_rem=sum(p[0] == "REM" for p in res["path"]), **{k: acc[k] for k in ("net_add", "net_rem", "net_outside_B", "net_raw_fg_removed",
                     "net_added_voxel_reverted", "net_new_bg_added", "net_removed_voxel_restored", "n_crit_end", "n_resolved_total")})
            prow = [dict(system=system, case_id=cid, block=bs, tier=tier, step=i, type=p[0], z=p[1], y=p[2], x=p[3], nmask=acc["nm"][i]) for i, p in enumerate(res["path"])]
    except Exception as e: r.update(status="ERROR", error=f"{type(e).__name__}: {e}"[:300])
    r["seconds"] = round(time.time() - t0, 3); return r, prow


def run_block(ra, system, co, cid, bs, cache, jdir=None):
    """Tier ladder for one target.  Each finished (target, tier) is journaled (jdir) so an interrupted shard resumes without recomputation; a journal entry is only used if it carries this runner's SHA."""
    cfg = P2B.env["cfg"]; B = P2B.parse_bk(bs); restricted = system[:2] in ("A_", "B_"); out, prow = [], []; tier = "A"
    while tier:
        jf = None if jdir is None else Path(jdir) / f"{cid}__{bs.replace(',', '_').replace('@', '_')}__{tier}.json"; j = None
        if jf is not None and jf.exists():
            j = json.load(open(jf)); j = j if j.get("runner_sha") == RSHA else None
        if j is None:
            r, pr = one_tier(ra, system, co, cid, bs, B, tier, cache); CNT["new"] += 1
            if jf is not None: 
                with open(jf.with_suffix(".tmp"), "w") as f: json.dump(dict(runner_sha=RSHA, row=r, prow=pr), f)
                os.replace(jf.with_suffix(".tmp"), jf)
        else: r, pr = j["row"], j["prow"]; CNT["resumed"] += 1
        out.append(r); prow += pr
        tier = (None if r["status"] in ("PURE_FOUND", "ERROR") else {"A": "B", "B": "C"}.get(tier)) if restricted else P2B.next_tier(cfg, tier, p2b_status(r) if r["status"] != "ERROR" else "ERROR")
    return out, prow


def rows(p): return list(csv.DictReader(open(p, encoding="utf-8-sig", newline="")))


def cmd_search(a, ra):
    k, n = (int(v) for v in a.shard.split("/")); tg = ra.targets(("C1", "C2")); tg = tg[:a.limit] if a.limit else tg; tg = tg[k::n]; P2B.check_cases(tg)
    final = Path(a.out) / f"shard_{k}"; part = Path(a.out) / f"shard_{k}.partial"
    if final.exists(): raise SystemExit(f"write-once: {final} exists")
    (part / "_j").mkdir(parents=True, exist_ok=True); BL, PA, cache, t0 = [], [], {}, time.time()
    for co, cid, bs in tg: b, p = run_block(ra, a.system, co, cid, bs, cache, part / "_j"); BL += b; PA += p
    with open(part / "BLOCKS.csv", "w", newline="") as f: w = csv.DictWriter(f, COLS); w.writeheader(); w.writerows(BL)
    with open(part / "PATHS.csv", "w", newline="") as f: w = csv.DictWriter(f, PCOLS); w.writeheader(); w.writerows(PA)
    json.dump(dict(system=a.system, runner_sha=RSHA, spec_sha=sha(Path(__file__).with_name(SPEC)), shard=a.shard, n_targets=len(tg), n_tier_rows_computed=CNT["new"], n_tier_rows_resumed=CNT["resumed"],
                   pythonhashseed=os.environ.get("PYTHONHASHSEED"), python=platform.python_version(), numpy=np.__version__, manifest_digest=P2B.env["FZ"].manifest_digest(), seconds=round(time.time() - t0, 1), cfg=P2B.env["cfg"]),
              open(part / "shard_meta.json", "w"), indent=1)
    shutil.rmtree(part / "_j"); os.rename(part, final); print(f"{a.system} shard {a.shard}: {len(tg)} targets, {len(BL)} rows, {time.time() - t0:.1f}s (computed {CNT['new']}, resumed {CNT['resumed']})")


def cmd_merge(a):
    d = Path(a.out); BL, PA, metas = [], [], []
    for k in range(a.nshards):
        s = d / f"shard_{k}"
        if not (s / "shard_meta.json").exists(): raise SystemExit(f"missing shard {k}")
        BL += rows(s / "BLOCKS.csv"); PA += rows(s / "PATHS.csv"); metas.append(json.load(open(s / "shard_meta.json")))
    if len({(m["runner_sha"], m["spec_sha"], m["manifest_digest"], m["system"]) for m in metas}) != 1: raise SystemExit("shards disagree on runner/spec/manifest/system")
    key = lambda r: (r["case_id"], r["block"], r["tier"]); BL.sort(key=key); PA.sort(key=lambda r: (r["case_id"], r["block"], r["tier"], int(r["step"])))
    if len({key(r) for r in BL}) != len(BL): raise SystemExit("duplicate block-tier rows")
    for nm, rr, cols in (("BLOCKS.csv", BL, COLS), ("PATHS.csv", PA, PCOLS)):
        with open(d / nm, "x", newline="") as f: w = csv.DictWriter(f, cols); w.writeheader(); w.writerows(rr)
    json.dump(dict(shards=metas, n_block_tier_rows=len(BL), n_blocks=len({(r["case_id"], r["block"]) for r in BL}), system=a.system), open(d / "META.json", "x"), indent=1); print("merged", len(BL), "rows,", len(PA), "path rows")


def main():
    ap = argparse.ArgumentParser(); sp = ap.add_subparsers(dest="cmd", required=True); s = sp.add_parser("search"); m = sp.add_parser("merge")
    for n in ("p1", "vendor", "p0", "table", "betti", "p4", "runa", "p2a", "config"): s.add_argument(f"--{n}", required=True)
    s.add_argument("--system", required=True, choices=SYSTEMS); s.add_argument("--out", required=True); s.add_argument("--shard", default="0/1"); s.add_argument("--limit", type=int, default=0)
    m.add_argument("--system", required=True, choices=SYSTEMS); m.add_argument("--out", required=True); m.add_argument("--nshards", type=int, required=True)
    a = ap.parse_args()
    if a.cmd == "merge": return cmd_merge(a)
    ra = boot(a); cmd_search(a, ra)


if __name__ == "__main__": main()
