#!/usr/bin/env python3
"""Independent verifier for the constraint-ablation runner (spec P2AB_ABLATION_SPEC_R0.md section 8).  Re-implements the cell rules, PURE/NRFR tests and the BFS per voxel (frozen P1 gate on the full padded volume, window-crop critical sets from the frozen vendor,
python sets), sharing only the frozen provenance boot / X_T reconstruction with the runner.  Python: `python -P -s -B`.
  qind --root DIR --frozen DIR   hash-stratified rows re-solved by the independent BFS (status, n_states, level_counts, n_moves_root, n_bres, lstar, path)
  witness --root DIR             every PURE_FOUND path replayed on the full volume: gate SAFE, window, locality / flip-once / fidelity per system, endpoint PURE, accounting, derived path statistics
  chain --root DIR               monotonicity / inclusion invariants over all (target, grid) groups (no search)"""
import argparse, csv, glob, hashlib, json, sys, time
from collections import defaultdict
from pathlib import Path
import numpy as np

SEED = "20261009"; CELLS = {"C00": (0, None), "C10": (1, None), "C01P": (0, "P"), "C01H": (0, "H"), "C11P": (1, "P"), "C11H": (1, "H")}; GRIDS = {"G26": 6, "G14": 4}; RGRID = {"G26": 2, "G14": 1}
V = S = None
rd = lambda p: list(csv.DictReader(open(p, encoding="utf-8-sig", newline="")))


def parse(sy): c, f = sy.split("-"); return CELLS[c][0], CELLS[c][1], f


def load_root(root):
    BL, PA = [], defaultdict(list)
    for s in sorted(glob.glob(str(Path(root) / "shard_*"))):
        BL += rd(Path(s) / "BLOCKS.csv")
        for r in rd(Path(s) / "PATHS.csv"): PA[(r["system"], r["grid"], r["case_id"], r["block"])].append(r)
    return BL, {k: [(r["type"], int(r["z"]), int(r["y"]), int(r["x"]), r["nmask"]) for r in sorted(v, key=lambda r: int(r["step"]))] for k, v in PA.items()}


def moves2(ind, sy, crit):
    loc, fo, fid = parse(sy); Vx = ind.vox(crit) if loc else None; out = []
    for v in ind.W:
        c = ind.cur(v)
        if Vx is not None and v not in Vx: continue
        if fo == "P" and c != ind.x0[v]: continue
        if fo == "H" and not (ind.raw[v] == 0 if c == 0 else ind.raw[v] == 1): continue
        if fid == "F2" and c == 1 and ind.raw[v] == 1: continue
        if ind.P1.gate(ind.X, v)[2] == "SAFE": out.append((v, c))
    return out


def ibfs2(ind, sy, L, nmax):
    f1 = sy.endswith("F1"); root = frozenset(); cr = {root: ind.C0}; n_root = len(moves2(ind, sy, ind.C0)); vis = {root}; levels = [[root]]; par_of = {}; count, nbres, found = 1, 0, None
    for d in range(L):
        nxt = []
        for p in levels[d]:
            ind.apply(p); mv = moves2(ind, sy, cr[p]); ind.apply(p)
            for v, _ in mv:
                ch = p ^ {v}
                if ch in vis: continue
                vis.add(ch); nxt.append(ch); count += 1; par_of[ch] = (p, v)
                if count > nmax: raise RuntimeError("cap reached in independent BFS")
                ind.apply(ch); c = ind.crit(); nr = ind.nrfr(); ind.apply(ch); cr[ch] = c; res = ind.bs not in c; pure = res and not (c - ind.C0); nbres += res
                if pure and (nr == 0 or not f1): found = (d + 1, ch); break
            if found: break
        levels.append(nxt)
        if found or not nxt: break
    D = len(levels) - 1; o = dict(n_states=count, n_bres=nbres, n_moves_root=n_root, level_counts="|".join(str(len(x)) for x in levels), path="", lstar=-1)
    if found:
        o["status"] = "PURE_FOUND"; o["lstar"] = found[0]; st, ch = [], found[1]
        while ch != root: ch, v = par_of[ch]; st.append(v)
        st.reverse(); cur = dict(ind.x0); pth = []
        for v in st: pth.append((("REM" if cur[v] else "ADD"),) + v); cur[v] ^= 1
        o["path"] = ";".join(f"{t},{z},{y},{x}" for t, z, y, x in pth)
    else: o["status"] = "NO_PURE_DEPTH_LIMITED" if (D == L and len(levels[D]) > 0) else "NO_PURE_EXHAUSTED"
    return o


def replay2(Xp, raw, bs, sy, R, path, row):
    """Full-volume replay; returns list of violations (empty = pass)."""
    P2B = S.P2B; P1, RM = P2B.env["P1"], P2B.env["R"]; loc, fo, fid = parse(sy); B = P2B.parse_bk(bs); bad = []; X = Xp.copy(); x0 = Xp.copy(); C0 = P2B.global_crit(Xp[1:-1, 1:-1, 1:-1]); sh = np.array(raw.shape)
    cs = np.array(RM._block_coords(B)); lo = np.maximum(cs.min(0) - R, 0); hi = np.minimum(cs.max(0) + R, sh - 1); nflip = defaultdict(int); hv = rr = pr = 0; first = -1; ftype = "NONE"
    for i, (t, z, y, x, nmask) in enumerate(path):
        v = (z, y, x); cur = int(X[z + 1, y + 1, x + 1]); crit = P2B.global_crit(X[1:-1, 1:-1, 1:-1])
        if (cur == 0) != (t == "ADD"): bad.append(f"step{i}:type"); break
        if not (np.all(np.array(v) >= lo) and np.all(np.array(v) <= hi)): bad.append(f"step{i}:outside_window")
        n, _, verdict, _ = P1.gate(X, v)
        if verdict != "SAFE": bad.append(f"step{i}:{verdict}")
        if str(n) != nmask: bad.append(f"step{i}:nmask")
        Vx = {tuple(c) for b in crit for c in RM._block_coords(P2B.parse_bk(b))}; locok = v in Vx; pathok = nflip[v] == 0; histok = (raw[v] == 0) if cur == 0 else (raw[v] == 1)
        if loc and not locok: bad.append(f"step{i}:locality")
        if fo == "P" and not pathok: bad.append(f"step{i}:flip_once_path")
        if fo == "H" and not histok: bad.append(f"step{i}:flip_once_hist")
        if fid == "F2" and cur == 1 and raw[v] == 1: bad.append(f"step{i}:raw_rem")
        hv += not histok; rr += (cur == 1 and raw[v] == 1); pr += (cur == 1 and x0[z + 1, y + 1, x + 1] == 1 and raw[v] == 0)
        if first < 0 and not (locok and pathok): first = i; ftype = "BOTH" if (not locok and not pathok) else "OUT_OF_LOCALITY" if not locok else "REFLIP"
        X[z + 1, y + 1, x + 1] ^= 1; nflip[v] += 1
    Xe = X[1:-1, 1:-1, 1:-1]; CY = P2B.global_crit(Xe); X0 = x0[1:-1, 1:-1, 1:-1]
    if bs in CY or (CY - C0): bad.append("endpoint_not_pure")
    acc = defaultdict(int); bset = set(RM._block_coords(B))
    for z, y, x in np.argwhere(Xe != X0).tolist():
        acc["net_rem" if X0[z, y, x] else "net_add"] += 1; acc["net_outside_B"] += (z, y, x) not in bset
        acc["net_raw_fg_removed" if (X0[z, y, x] and raw[z, y, x]) else "net_added_voxel_reverted" if X0[z, y, x] else "net_new_bg_added" if not raw[z, y, x] else "net_removed_voxel_restored"] += 1
    if fid in ("F1", "F2") and acc["net_raw_fg_removed"]: bad.append("nrfr_nonzero")
    exp = dict(net_add=acc["net_add"], net_rem=acc["net_rem"], net_outside_B=acc["net_outside_B"], net_raw_fg_removed=acc["net_raw_fg_removed"], net_added_voxel_reverted=acc["net_added_voxel_reverted"], net_new_bg_added=acc["net_new_bg_added"],
               net_removed_voxel_restored=acc["net_removed_voxel_restored"], n_crit_end=len(CY), n_repeated_flip_voxels=sum(1 for c in nflip.values() if c >= 2), hist_violation_steps=hv, rawfg_rem_steps=rr, p2a_added_rem_steps=pr,
               first_c11_violation_step=first, first_c11_violation_type=ftype, lstar=len(path))
    for k, e in exp.items():
        if str(row[k]) != str(e): bad.append(f"row_{k}:{row[k]}!={e}")
    return bad


class Cache:
    def __init__(s, ra): s.ra, s.cid, s.d = ra, None, None
    def get(s, cid):
        if s.cid != cid:
            Xp, raw, ok = s.ra.reconstruct(cid)
            if not ok: raise RuntimeError(f"{cid}: X_T sha mismatch")
            s.cid, s.d = cid, (Xp, raw)
        return s.d


def emit(a, name, ok, detail, n=None):
    d = Path(a.out); d.mkdir(parents=True, exist_ok=True); p = d / f"{name}.json"
    if p.exists(): raise SystemExit(f"write-once: {p} exists")
    json.dump(dict(id=name, ok=bool(ok), detail=detail, n=n), open(p, "w"), indent=1, ensure_ascii=False); print(("PASS " if ok else "FAIL ") + name, detail, n if n is not None else ""); return ok


def hkey(r): return hashlib.sha256(f"{SEED}|{r['system']}|{r['grid']}|{r['case_id']}|{r['block']}".encode()).hexdigest()


def cmd_witness(a, ra):
    BL, PA = load_root(a.root); cache = Cache(ra); k, n = (int(v) for v in a.shard.split("/")); err = []; cnt = 0
    for r in sorted((r for r in BL if r["status"] == "PURE_FOUND"), key=lambda r: (r["system"], r["grid"], r["case_id"], r["block"]))[k::n]:
        Xp, raw = cache.get(r["case_id"]); bad = replay2(Xp, raw, r["block"], r["system"], int(r["R"]), PA[(r["system"], r["grid"], r["case_id"], r["block"])], r); cnt += 1
        if bad: err.append(f"{r['system']} {r['grid']} {r['case_id']} {r['block']}: {bad[:3]}")
    emit(a, f"witness_{a.shard.replace('/', 'of')}", not err, "; ".join(err[:5]) or f"{cnt} witnesses replayed on the full volume", len(err))


def cmd_qind(a, ra):
    BL, PA = load_root(a.root); cache = Cache(ra); k, n = (int(v) for v in a.shard.split("/")); strata = defaultdict(list)
    for r in BL:
        if r["status"] != "CAP_EXCEEDED" and int(r["n_states"]) <= a.max_states and int(r["n_states"]) * int(r["window_voxels"]) <= a.budget: strata[(r["system"], r["grid"], r["status"])].append(r)
    pick = sorted((x for v in strata.values() for x in sorted(v, key=hkey)[:a.per_stratum]), key=lambda r: (r["system"], r["grid"], r["case_id"], r["block"]))[k::n]; err = []
    for r in pick:
        Xp, raw = cache.get(r["case_id"]); ind = V.Ind(Xp, raw, r["block"], int(r["R"])); o = ibfs2(ind, r["system"], int(r["L"]), 10 ** 7)
        for kk in ("status", "n_states", "level_counts", "n_bres", "n_moves_root", "lstar"):
            if str(o[kk]) != str(r[kk]): err.append(f"{r['system']} {r['grid']} {r['case_id']} {r['block']} {kk}: indep {o[kk]} != runner {r[kk]}")
        if r["status"] == "PURE_FOUND" and o["path"] != ";".join(f"{t},{z},{y},{x}" for t, z, y, x, _ in PA[(r["system"], r["grid"], r["case_id"], r["block"])]): err.append(f"{r['system']} {r['case_id']}: path differs")
    emit(a, f"qind_{a.shard.replace('/', 'of')}", not err, "; ".join(err[:5]) or f"{len(pick)} hash-stratified rows (seed {SEED}, <= {a.per_stratum}/stratum) identical", len(err))


def cmd_chain(a, ra):
    BL, _ = load_root(a.root); g = defaultdict(dict)
    for r in BL: g[(r["case_id"], r["block"], r["grid"])][r["system"]] = r
    found = lambda r: r["status"] == "PURE_FOUND"; det = lambda r: r["status"] not in ("CAP_EXCEEDED", "ERROR"); err = []; n = 0
    chains = [("C11P", "C10"), ("C11P", "C01P"), ("C10", "C00"), ("C01P", "C00"), ("C11H", "C01H"), ("C01H", "C00")]
    for key, d in g.items():
        for fid in ("F0", "F1", "F2"):
            for lo, hi in chains:
                x, y = d.get(f"{lo}-{fid}"), d.get(f"{hi}-{fid}")
                if x is None or y is None or not (det(x) and det(y)): continue
                n += 1
                if found(x) and not found(y): err.append(f"{key} {fid}: {lo} found but {hi} not")
                if found(x) and found(y) and int(y["lstar"]) > int(x["lstar"]): err.append(f"{key} {fid}: depth {hi} {y['lstar']} > {lo} {x['lstar']}")
        for c in CELLS:
            for lo, hi in (("F2", "F1"), ("F1", "F0")):
                x, y = d.get(f"{c}-{lo}"), d.get(f"{c}-{hi}")
                if x is None or y is None or not (det(x) and det(y)): continue
                n += 1
                if found(x) and not found(y): err.append(f"{key} {c}: {lo} found but {hi} not")
            x, y = d.get(f"{c}-F1"), d.get(f"{c}-F2")
            if CELLS[c][1] and x is not None and y is not None and det(x) and det(y):
                n += 1
                if (found(x), x["lstar"]) != (found(y), y["lstar"]): err.append(f"{key} {c}: F1 != F2 under flip-once")
    emit(a, "chain", not err, "; ".join(err[:5]) or f"{n} pairwise invariants over {len(g)} (target, grid) groups hold", len(err))


def main():
    global V, S
    ap = argparse.ArgumentParser(); ap.add_argument("cmd", choices=["qind", "witness", "chain"])
    for n_ in ("p1", "vendor", "p0", "table", "betti", "p4", "runa", "p2a", "config", "frozen", "root", "out"): ap.add_argument(f"--{n_}", required=True)
    ap.add_argument("--shard", default="0/1"); ap.add_argument("--per-stratum", dest="per_stratum", type=int, default=10); ap.add_argument("--max-states", dest="max_states", type=int, default=20000); ap.add_argument("--budget", type=int, default=400000)
    a = ap.parse_args(); sys.path.insert(0, str(Path(a.frozen).resolve())); import p2ab_search_validate as V_, p2ab_search as S_; V, S = V_, S_
    ra = S.boot(a); globals()["cmd_" + a.cmd](a, ra)


if __name__ == "__main__": main()
