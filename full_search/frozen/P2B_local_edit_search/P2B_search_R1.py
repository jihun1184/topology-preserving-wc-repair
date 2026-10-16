#!/usr/bin/env python3
"""P2B R1 runner: bounded exhaustive local SAFE-edit search (spec P2B_SPEC_R1_1.md; R1.1 = R1 + routing/cap clarifications).  Python: `python -P -s -B`.
  q0      --out DIR                         start-state qualification on the 680 target blocks (C1 614 + C2 66); exit 0 iff 680/680 pass
  search  --cohort c1c2|c5 --out DIR --shard k/n [--limit N] [--approval FILE]   tier ladder per block; c1c2 requires the freeze-approval file
  merge   --out DIR --nshards n --cohort c1c2|c5                                 merge write-once shards
Frozen inputs are SHA-verified first (exit 3).  The module is also imported by the qualification tests (bfs / env / helpers)."""
import argparse, csv, hashlib, importlib.util, json, os, platform, sys, time
from collections import defaultdict
from pathlib import Path
import numpy as np

# ---- pins (Run A analysis inputs + imported P2A manifest module); everything else is pinned by P2A_frozen_R2.FROZEN -------------------------------
PIN_RUNA = {"P2A_RUN_CASES.csv": "bc146d6bfa7d7eda6e0129aefaaf4349985bd48376c1e6e9635d9acacf5450d6", "P2A_RUN_EDITS.csv": "ba801ac55a0f436d5cd06e541bba7f9f048fa30585c888d06cadc78927cfac11",
            "P2A_RUN_STUCK.csv": "58fa33dedf19cdb196fb83991d55e0b9699b126d7007a6b6b4cef4b40fb3fa5e", "P2A_RUN_META.json": "8097a8e2539f8276985c3259c8b3d093214f47692bbb8eea8bcbd7ec23767b1d"}
PIN_P2A_FROZEN = "944836057b2f6ec06083698192773c313bdea1401b5125b64ed440662efe4b10"
EXPECT = {"C1": (614, 213), "C2": (66, 39)}                                      # (blocks, cases) frozen in the P2B R1 spec (from the P2A final STUCK table)
CH = 4096                                                                        # parents expanded per vectorised chunk (memory only; no effect on results)
env = {}                                                                         # filled by load_env(): P0, P1, R (vendor), TBL, BETTI, cfg, ...
NMFLAG = None                                                                    # uint8[2^26]: bit0 = nmask used by a search, bit1 = audited against the frozen verdict


def sha(p):
    h = hashlib.sha256()
    with open(p, "rb") as f:
        for b in iter(lambda: f.read(1 << 22), b""): h.update(b)
    return h.hexdigest()


def _load(name, path):
    spec = importlib.util.spec_from_file_location(name, str(path)); m = importlib.util.module_from_spec(spec); sys.modules[name] = m; spec.loader.exec_module(m); return m


def rows(p): return list(csv.DictReader(open(p, encoding="utf-8-sig")))


def parse_bk(s):
    """'k2xy@88,71,122' -> BlockKey (same text format as the frozen CSVs); round-trip checked."""
    h, c = s.split("@"); z, y, x = (int(v) for v in c.split(",")); b = env["R"].BlockKey(int(h[1]), h[2:], z, y, x)
    if bk(b) != s: raise ValueError(f"block key round-trip failed: {s}")
    return b


bk = lambda b: f"k{b.k}{b.plane}@{b.z},{b.y},{b.x}"


def cohort_of(r):
    """P1 terminal block row -> C1..C5 (identical to the frozen P2A validator's definition)."""
    if r["cand_classes"] == "shell|shell" and r["category"] in ("RAW_SHELL_CLOSED", "FINAL_SHELL_CLOSED"):
        return {("2", "RAW_SHELL_CLOSED"): "C1", ("2", "FINAL_SHELL_CLOSED"): "C2", ("3", "RAW_SHELL_CLOSED"): "C3", ("3", "FINAL_SHELL_CLOSED"): "C4"}[(r["k"], r["category"])]
    return "C5"


# ======================================================================================================================== environment / provenance
def load_env(a):
    """Provenance gate (exit 3 on any mismatch), then import the frozen modules.  Nothing scientific is computed before this returns."""
    global NMFLAG
    fz = Path(a.p2a) / "P2A_frozen_R2.py"
    if sha(fz) != PIN_P2A_FROZEN: print("PROVENANCE_FAIL: P2A_frozen_R2.py sha", file=sys.stderr); sys.exit(3)
    FZ = _load("p2a_frozen_r2", fz); roots = dict(p1=a.p1, vendor=a.vendor, p0=a.p0, table=a.table, betti=a.betti, p4=a.p4)
    try: got = FZ.verify(roots)
    except FZ.ProvenanceError as e: print(f"PROVENANCE_FAIL: {e}", file=sys.stderr); sys.exit(3)
    for n, h in PIN_RUNA.items():
        if sha(Path(a.runa) / n) != h: print(f"PROVENANCE_FAIL: Run A {n}", file=sys.stderr); sys.exit(3)
    sys.path.insert(0, str(Path(a.vendor).resolve())); import r_plus_repair_production as R
    P0 = _load("p0_audit", a.p0); TBL = P0.load_table(a.table); BETTI = _load("cubical_betti", a.betti).cubical_betti
    P1 = _load("p1_frozen", Path(a.p1) / "P1_simple_constrained_rplus.py"); P1.setup(a.vendor, a.p0, a.table, None)
    cfg = json.load(open(a.config, encoding="utf-8"))
    for t in cfg["tiers"].values():
        if not (isinstance(t["R"], int) and isinstance(t["L"], int) and t["R"] >= 0 and t["L"] >= 1): raise SystemExit("bad tier in config")
    env.update(P0=P0, P1=P1, R=R, TBL=TBL, BETTI=BETTI, cfg=cfg, a=a, prov=got, FZ=FZ, roots=roots)
    NMFLAG = np.zeros(1 << 26, np.uint8)
    return env


def check_cases(tg):
    """Every case input used must be covered by the P4 manifest (P2A_frozen_R2.verify, case-coverage branch)."""
    try: env["FZ"].verify(env["roots"], case_ids=sorted({c for _, c, _ in tg}))
    except env["FZ"].ProvenanceError as e: print(f"PROVENANCE_FAIL: {e}", file=sys.stderr); sys.exit(3)


# ======================================================================================================================== inputs: X_T reconstruction and targets
class RunA:
    """Run A CASES / ordered EDITS / STUCK and the frozen P1 terminal table (read-only)."""
    def __init__(self, a):
        d = Path(a.runa); self.cases = {r["case_id"]: r for r in rows(d / "P2A_RUN_CASES.csv")}
        self.edits = defaultdict(list)
        for r in rows(d / "P2A_RUN_EDITS.csv"): self.edits[r["case_id"]].append(r)
        self.stuck = defaultdict(dict)
        for r in rows(d / "P2A_RUN_STUCK.csv"): self.stuck[r["case_id"]][r["block"]] = r
        self.p1t = defaultdict(dict)
        for r in rows(Path(a.p1) / "P1_TERMINAL_STUCK.csv"): self.p1t[r["case_id"]][r["block"]] = (cohort_of(r), r)
        self.p4 = Path(a.p4)

    def targets(self, cohorts):
        """[(cohort, case, block_str)] sorted by (case, block): P1-terminal block of that cohort that is still critical in the Run A final state."""
        t = [(co, c, b) for c, d in self.p1t.items() for b, (co, _) in d.items() if co in cohorts and b in self.stuck.get(c, {})]
        return sorted(t, key=lambda r: (r[1], r[2]))

    def reconstruct(self, cid):
        """pinned P4 RAW -> Run A ordered edit replay -> padded X_T (uint8).  Returns (Xp, raw, sha_ok)."""
        P0 = env["P0"]; raw = P0.load_zyx(self.p4 / "raw" / f"{cid}_WT_RAW.nii.gz"); Xp = np.pad(raw, 1).astype(np.uint8); er = self.edits.get(cid, [])
        if [int(r["step"]) for r in er] != list(range(len(er))): raise RuntimeError(f"{cid}: edit steps not contiguous")
        for r in er:
            z, y, x = int(r["z"]) + 1, int(r["y"]) + 1, int(r["x"]) + 1
            if (Xp[z, y, x] == 0) != (r["type"] == "ADD"): raise RuntimeError(f"{cid}: step {r['step']} {r['type']} on value {Xp[z, y, x]}")
            Xp[z, y, x] = 1 if r["type"] == "ADD" else 0
        return Xp, raw, P0.mask_sha(Xp[1:-1, 1:-1, 1:-1]) == self.cases[cid]["mask_sha"]


# ======================================================================================================================== vectorised local geometry
def corners(X, plane):
    """Block corner views in the frozen vendor order (a,b,c,d for k=2; the 8 corners for k=3) over the last three axes."""
    if plane == "xy": return (X[..., :, :-1, :-1], X[..., :, :-1, 1:], X[..., :, 1:, :-1], X[..., :, 1:, 1:])
    if plane == "yz": return (X[..., :-1, :-1, :], X[..., :-1, 1:, :], X[..., 1:, :-1, :], X[..., 1:, 1:, :])
    if plane == "zx": return (X[..., :-1, :, :-1], X[..., :-1, :, 1:], X[..., 1:, :, :-1], X[..., 1:, :, 1:])
    return tuple(X[..., i:X.shape[-3] - 1 + i, j:X.shape[-2] - 1 + j, k:X.shape[-1] - 1 + k] for i in (0, 1) for j in (0, 1) for k in (0, 1))   # "k3"


def crit_mask(X, plane):
    """Critical-block indicator of every block of `plane` fully inside X (same predicate as vendor _critical_keys_vectorized)."""
    c = corners(X, plane)
    if plane != "k3": return (c[0] == c[3]) & (c[1] == c[2]) & (c[0] != c[1])
    s = sum(v.astype(np.uint16) for v in c); out = np.zeros(s.shape, bool)
    for i, j in ((0, 7), (1, 6), (2, 5), (3, 4)): out |= (s == 2) & (c[i] == 1) & (c[j] == 1); out |= (s == 6) & (c[i] == 0) & (c[j] == 0)
    return out


PLANES = ("xy", "yz", "zx", "k3")


class Local:
    """Window + 1-voxel halo geometry for one (state X_T, block B, R).  G = padded X_T cropped to bbox(window)+-1; G-frame index g = u - lo + 1 (u = unpadded coordinate)."""
    def __init__(self, Xp, B, R):
        Rm = env["R"]; shape = tuple(s - 2 for s in Xp.shape); cs = np.array(Rm._block_coords(B))
        self.lo = np.maximum(cs.min(0) - R, 0); self.hi = np.minimum(cs.max(0) + R, np.array(shape) - 1); self.wshape = tuple((self.hi - self.lo + 1).tolist())
        self.nW = int(np.prod(self.wshape)); self.nb = (self.nW + 7) // 8
        self.G0 = np.ascontiguousarray(Xp[self.lo[0]:self.hi[0] + 3, self.lo[1]:self.hi[1] + 3, self.lo[2]:self.hi[2] + 3])             # padded index = u+1 -> slice lo..hi+2
        if self.G0.shape != tuple(w + 2 for w in self.wshape): raise RuntimeError("halo crop wrong (window not inside the domain)")
        valid = np.zeros(self.G0.shape, bool); u0 = -self.lo + 1                                                                         # valid <=> in-domain (u in [0, shape-1])
        valid[tuple(slice(max(0, u0[i]), min(self.G0.shape[i], shape[i] + u0[i])) for i in range(3))] = True
        self.vblk = {p: np.logical_and.reduce(corners(valid, p)) for p in PLANES}; self.c0 = {p: crit_mask(self.G0, p) & self.vblk[p] for p in PLANES}
        self.bidx = [tuple((np.array(c) - self.lo + 1).tolist()) for c in Rm._block_coords(B)]; self.B = B; self.k = B.k
        i = np.arange(self.nW); self.UNIT = np.zeros((self.nW, self.nb), np.uint8); self.UNIT[i, i >> 3] = (1 << (i & 7)).astype(np.uint8)
        self.coords = np.stack(np.unravel_index(i, self.wshape), 1) + self.lo                                                            # window voxel i -> unpadded (z,y,x), (z,y,x)-lexicographic
        self.key0 = np.packbits(self.G0[1:-1, 1:-1, 1:-1].reshape(-1), bitorder="little")                                                # bit i of the key = occupancy of window voxel i

    def unpack(self, keys):
        m = len(keys); bits = np.unpackbits(keys, axis=1, count=self.nW, bitorder="little")
        Gb = np.broadcast_to(self.G0, (m,) + self.G0.shape).copy(); Gb[:, 1:-1, 1:-1, 1:-1] = bits.reshape((m,) + self.wshape); return Gb

    def b_critical(self, Gb):
        """Is the target block still critical?  (bool[m])"""
        v = [Gb[:, z, y, x] for z, y, x in self.bidx]
        if self.k == 2: return (v[0] == v[3]) & (v[1] == v[2]) & (v[0] != v[1])
        s = sum(t.astype(np.uint16) for t in v); out = np.zeros(len(Gb), bool)
        for i, j in ((0, 7), (1, 6), (2, 5), (3, 4)): out |= (s == 2) & (v[i] == 1) & (v[j] == 1); out |= (s == 6) & (v[i] == 0) & (v[j] == 0)
        return out

    def new_crit(self, Gb):
        """Does the state have an in-domain critical block that X_T does not?  (blocks outside the local frame are unchanged by construction)  (bool[m])"""
        out = np.zeros(len(Gb), bool)
        for p in PLANES: out |= (crit_mask(Gb, p) & self.vblk[p] & ~self.c0[p]).reshape(len(Gb), -1).any(1)
        return out

    def safe(self, Gb, tbl):
        """SAFE indicator of every window voxel of every state: 26-neighbourhood bitmask (P0.OFFS bit order) -> frozen gate table; bits&7==7 <=> attachment contractible."""
        P0 = env["P0"]; m, a, b, c = Gb.shape; nm = np.zeros((m, a - 2, b - 2, c - 2), np.uint32)
        for k, (da, db, dc) in enumerate(P0.OFFS): nm |= Gb[:, 1 + da:a - 1 + da, 1 + db:b - 1 + db, 1 + dc:c - 1 + dc].astype(np.uint32) << np.uint32(k)
        nm = nm.reshape(m, -1); NMFLAG[nm] |= 1
        return (tbl[nm] & 7) == 7


def audit_gate(tbl=None):
    """Every neighbourhood bitmask any search has used is re-evaluated by the frozen verdict (attachment homology + table bits); raises if SAFE-ness differs from the vectorised lookup."""
    P0 = env["P0"]; tbl = env["TBL"] if tbl is None else tbl; todo = np.flatnonzero(NMFLAG == 1); n = 0
    for m in todo.tolist():
        att = P0.attachment(m); v = P0.verdict(att, env["TBL"], m)
        if v == "UNRESOLVED" or (v == "SAFE") != bool(int(tbl[m]) & 7 == 7): raise RuntimeError(f"GATE_AUDIT_FAIL nmask={m}: frozen verdict {v} vs lookup bits {int(tbl[m]) & 7}")
        NMFLAG[m] |= 2; n += 1
    return n


# ======================================================================================================================== the search
def bfs(Xp, B, R, L, nmax, tbl=None, dedup=True):
    """Exhaustive FIFO BFS over window states reachable by single-voxel SAFE flips, depth <= L, exact bitmask dedup.  Success = PURE (spec R1 section 4).
    Discovery order is the sequential-FIFO order (parent order, then voxel (z,y,x) order); the (nmax+1)-th distinct state aborts with CAP_EXCEEDED unless a PURE state was discovered earlier.
    `tbl`/`dedup` are test hooks (fault injection); production uses the defaults."""
    tbl = env["TBL"] if tbl is None else tbl; lc = Local(Xp, B, R); nb = lc.nb
    if bool(lc.b_critical(lc.G0[None])[0]) is not True: raise RuntimeError("target block is not critical in the start state")
    vis = {lc.key0.tobytes()}; lk = [lc.key0[None, :]]; lp = [np.array([-1])]; lv = [np.array([-1])]; count = 1; nbres = 0; minb = -1; cap = False; found = None
    for d in range(L):
        par = lk[d]; nk, npar, nvox = [], [], []; n_lvl = 0
        for s in range(0, len(par), CH):
            P = par[s:s + CH]; pi, vi = np.nonzero(lc.safe(lc.unpack(P), tbl))
            if not pi.size: continue
            ck = P[pi] ^ lc.UNIT[vi]; kv = np.ascontiguousarray(ck).view(np.dtype((np.void, nb))).ravel()
            first = np.sort(np.unique(kv, return_index=True)[1]) if dedup else np.arange(len(kv))                  # first occurrence in discovery order
            keep = []
            for j, b in zip(first.tolist(), kv[first].tolist()):
                if dedup:
                    if b in vis: continue
                    vis.add(b)
                keep.append(j)
            if count + len(keep) > nmax: keep = keep[:nmax - count]; cap = True                                  # states of rank > nmax are never admitted
            if keep:
                K = ck[keep]; Gn = lc.unpack(K); res = ~lc.b_critical(Gn); pi_, vi_ = pi[keep], vi[keep]
                if res.any():
                    nc = np.zeros(len(K), bool); nc[res] = lc.new_crit(Gn[res]); pure = res & ~nc
                    if pure.any():                                                                             # sequential semantics: nothing after the first PURE state is discovered
                        j0 = int(np.argmax(pure)); found = (d + 1, n_lvl + j0, K[j0].copy()); K, res, pi_, vi_ = K[:j0 + 1], res[:j0 + 1], pi_[:j0 + 1], vi_[:j0 + 1]
                    nbres += int(res.sum())
                    if res.any() and minb < 0: minb = d + 1
                nk.append(K); npar.append(s + pi_); nvox.append(vi_); n_lvl += len(K); count += len(K)
            if found is not None or cap: break
        lk.append(np.concatenate(nk) if nk else np.zeros((0, nb), np.uint8)); lp.append(np.concatenate(npar) if npar else np.zeros(0, int)); lv.append(np.concatenate(nvox) if nvox else np.zeros(0, int))
        if found is not None or cap or not nk: break
    out = dict(n_states=count, window=lc.nW, lo=lc.lo.tolist(), hi=lc.hi.tolist(), level_counts=[len(k) for k in lk], n_bres=nbres, min_depth_bres=minb, path=None, lstar=-1)
    if found is not None:
        depth, idx, key = found; out["status"] = "PURE_FOUND"; out["lstar"] = depth; out["complete"] = 0; steps = []; i = idx
        for dd in range(depth, 0, -1): steps.append(int(lv[dd][i])); i = int(lp[dd][i])
        steps.reverse(); cur = lc.G0[1:-1, 1:-1, 1:-1].reshape(-1).copy(); path = []
        for v in steps: z, y, x = lc.coords[v].tolist(); path.append(("REM" if cur[v] else "ADD", z, y, x)); cur[v] ^= 1
        out["path"] = path; out["endpoint_key"] = int.from_bytes(key.tobytes(), "little")
        if int.from_bytes(np.packbits(cur, bitorder="little").tobytes(), "little") != out["endpoint_key"]: raise RuntimeError("path replay != endpoint key")
    elif cap: out["status"] = "CAP_EXCEEDED"; out["complete"] = 0
    else: out["status"] = "B_RESOLVED_WITH_NEW_CRIT" if nbres else "NO_PATH_EXHAUSTIVE"; out["complete"] = 1
    return out


# ======================================================================================================================== per-block driver (tier ladder, witness self-check, endpoint accounting)
def global_crit(X): return {bk(k) for k in env["R"]._critical_keys_vectorized(np.ascontiguousarray(X))}


def replay_witness(Xp, raw, B, path, C0):
    """Runner-side check of a PURE witness with the frozen P1 gate on the full padded volume (the validator repeats this independently).  Returns accounting dict."""
    P1 = env["P1"]; X = Xp.copy(); seen = {X.tobytes()}; nm = []
    for t, z, y, x in path:
        if (X[z + 1, y + 1, x + 1] == 0) != (t == "ADD"): raise RuntimeError("witness step type/state mismatch")
        n, att, v, _ = P1.gate(X, (z, y, x))
        if v != "SAFE": raise RuntimeError(f"witness step ({t},{z},{y},{x}) is {v}")
        nm.append(n); X[z + 1, y + 1, x + 1] ^= 1
        if X.tobytes() in seen: raise RuntimeError("witness revisits a state")
        seen.add(X.tobytes())
    Xe = X[1:-1, 1:-1, 1:-1]; CY = global_crit(Xe)
    if bk(B) in CY or (CY - C0): raise RuntimeError("endpoint is not PURE")
    X0 = Xp[1:-1, 1:-1, 1:-1]; bset = set(env["R"]._block_coords(B)); acc = defaultdict(int)
    for z, y, x in np.argwhere(Xe != X0).tolist():
        acc["net_rem" if X0[z, y, x] else "net_add"] += 1; acc["net_outside_B"] += (z, y, x) not in bset
        acc["net_raw_fg_removed" if (X0[z, y, x] and raw[z, y, x]) else "net_added_voxel_reverted" if X0[z, y, x] else "net_new_bg_added" if not raw[z, y, x] else "net_removed_voxel_restored"] += 1
    acc.update(n_crit_end=len(CY), n_resolved_total=len(C0 - CY), nm=nm); return acc


COLS = ["cohort", "case_id", "block", "tier", "R", "L", "n_max", "window_voxels", "status", "n_states", "complete", "level_counts", "n_bres_states", "min_depth_bres", "lstar",
        "n_flip_add", "n_flip_rem", "net_add", "net_rem", "net_outside_B", "net_raw_fg_removed", "net_added_voxel_reverted", "net_new_bg_added", "net_removed_voxel_restored",
        "n_crit_start", "n_crit_end", "n_resolved_total", "n_resolved_other", "pure_reachable_by_B", "gate_audited_new", "seconds", "error"]
NONDET = ("seconds", "gate_audited_new")                                          # diagnostic columns excluded from the determinism comparison
PCOLS = ["case_id", "block", "tier", "step", "type", "z", "y", "x", "nmask"]


def next_tier(cfg, tier, status):
    """Tier routing (config-driven, spec R1.1 section 5.1): the next tier if `status` of `tier` is listed in routing[tier].when, else None (stop)."""
    rt = cfg["routing"].get(tier); return rt["next"] if rt and status in rt["when"] else None


def run_block(ra, co, cid, bs, cache):
    """Tier ladder for one block.  Always starts from the case's own X_T; no state is shared between blocks.  Returns (block rows, path rows)."""
    cfg = env["cfg"]; B = parse_bk(bs)
    if cache.get("cid") != cid:
        Xp, raw, ok = ra.reconstruct(cid)
        if not ok: raise RuntimeError(f"{cid}: reconstructed X_T SHA != Run A mask_sha")
        cache.clear(); cache.update(cid=cid, Xp=Xp, raw=raw, C0=global_crit(Xp[1:-1, 1:-1, 1:-1]))
    Xp, raw, C0 = cache["Xp"], cache["raw"], cache["C0"]; out, prow = [], []; tier = "A"
    while tier:
        t = cfg["tiers"][tier]; t0 = time.time(); r = dict.fromkeys(COLS, ""); r.update(cohort=co, case_id=cid, block=bs, tier=tier, R=t["R"], L=t["L"], n_max=cfg["n_max"], n_crit_start=len(C0))
        try:
            res = bfs(Xp, B, t["R"], t["L"], cfg["n_max"]); na = audit_gate()
            r.update(status=res["status"], window_voxels=res["window"], n_states=res["n_states"], complete=res["complete"], level_counts="|".join(map(str, res["level_counts"])),
                     n_bres_states=res["n_bres"], min_depth_bres=res["min_depth_bres"], lstar=res["lstar"], gate_audited_new=na)
            if res["status"] == "PURE_FOUND":
                acc = replay_witness(Xp, raw, B, res["path"], C0)
                r.update(n_flip_add=sum(p[0] == "ADD" for p in res["path"]), n_flip_rem=sum(p[0] == "REM" for p in res["path"]), n_resolved_other=acc["n_resolved_total"] - 1,
                         **{k: acc[k] for k in ("net_add", "net_rem", "net_outside_B", "net_raw_fg_removed", "net_added_voxel_reverted", "net_new_bg_added", "net_removed_voxel_restored", "n_crit_end", "n_resolved_total")})
                prow += [dict(case_id=cid, block=bs, tier=tier, step=i, type=p[0], z=p[1], y=p[2], x=p[3], nmask=acc["nm"][i]) for i, p in enumerate(res["path"])]
        except Exception as e:                                                                        # any exception becomes an ERROR row, never a silent skip
            r.update(status="ERROR", error=f"{type(e).__name__}: {e}"[:300])
        r["seconds"] = round(time.time() - t0, 3); out.append(r)
        tier = next_tier(cfg, tier, r["status"])
    pb = int(any(x["status"] == "PURE_FOUND" and x["tier"] in ("A", "B") for x in out))
    for r in out: r["pure_reachable_by_B"] = pb
    return out, prow


# ======================================================================================================================== Q0
def q0(ra, tg):
    """Start-state qualification.  Per block: X_T SHA, target critical (vendor recount), identity (cohort, Run A STUCK row GATE_ALL_UNSAFE, P1 terminal row, case critical set == STUCK rows),
    every voxel flip of the block UNSAFE under the frozen P1.gate.  Returns (per-block rows, {case: sha_ok})."""
    P1 = env["P1"]; res = []; cache = {}; ncase = {}
    for co, cid, bs in tg:
        r = dict(cohort=co, case_id=cid, block=bs, sha_ok=0, critical=0, identity=0, all_unsafe=0, detail="")
        try:
            if cache.get("cid") != cid:
                Xp, raw, ok = ra.reconstruct(cid); cache.update(cid=cid, Xp=Xp, ok=ok, C0=global_crit(Xp[1:-1, 1:-1, 1:-1])); ncase[cid] = int(ok)
            Xp = cache["Xp"]; B = parse_bk(bs); st = ra.stuck[cid][bs]; c1, p1r = ra.p1t[cid][bs]
            r["sha_ok"] = int(cache["ok"]); r["critical"] = int(bs in cache["C0"])
            r["identity"] = int(c1 == co and p1r["k"] == "2" and p1r["cand_classes"] == "shell|shell" and st["cause"] == "GATE_ALL_UNSAFE" and st["k"] == "2" and B.k == 2 and cache["C0"] == set(ra.stuck[cid]))
            vs = [P1.gate(Xp, c)[2] for c in env["R"]._block_coords(B)]; r["all_unsafe"] = int(len(vs) == 4 and all(v == "UNSAFE" for v in vs)); r["detail"] = ",".join(vs)
        except Exception as e: r["detail"] = f"{type(e).__name__}: {e}"[:200]
        res.append(r)
    return res, ncase


# ======================================================================================================================== CLI
def args_common(ap):
    for n in ("p1", "vendor", "p0", "table", "betti", "p4", "runa", "p2a", "config"): ap.add_argument(f"--{n}", required=True)


def targets_for(ra, cohort, limit):
    t = ra.targets({"c1c2": ("C1", "C2"), "c5": ("C5",)}[cohort]); return t[:limit] if limit else t


def main():
    ap = argparse.ArgumentParser(); sp = ap.add_subparsers(dest="cmd", required=True)
    q = sp.add_parser("q0"); args_common(q); q.add_argument("--out", required=True)
    s = sp.add_parser("search"); args_common(s); s.add_argument("--cohort", required=True, choices=("c1c2", "c5")); s.add_argument("--out", required=True); s.add_argument("--shard", default="0/1")
    s.add_argument("--limit", type=int, default=0); s.add_argument("--approval", default="")
    m = sp.add_parser("merge"); m.add_argument("--out", required=True); m.add_argument("--nshards", type=int, required=True); m.add_argument("--cohort", required=True, choices=("c1c2", "c5"))
    a = ap.parse_args()
    if a.cmd == "merge": return merge(a)
    load_env(a); ra = RunA(a); me = sha(__file__)
    if a.cmd == "q0":
        tg = targets_for(ra, "c1c2", 0); check_cases(tg); res, ncase = q0(ra, tg); out = Path(a.out); out.mkdir(parents=True, exist_ok=False)
        with open(out / "P2B_Q0.csv", "w", newline="") as f: w = csv.DictWriter(f, list(res[0])); w.writeheader(); w.writerows(res)
        cnt = {c: sum(r["cohort"] == c for r in res) for c in ("C1", "C2")}; cs = {c: len({r["case_id"] for r in res if r["cohort"] == c}) for c in ("C1", "C2")}
        fails = {k: sum(not r[k] for r in res) for k in ("sha_ok", "critical", "identity", "all_unsafe")}
        ok = len(res) == 680 and not any(fails.values()) and all((cnt[c], cs[c]) == EXPECT[c] for c in EXPECT)
        summ = dict(total=len(res), per_cohort_blocks=cnt, per_cohort_cases=cs, n_cases=len(ncase), n_cases_sha_ok=sum(ncase.values()), fails=fails, pass_all=bool(ok), runner_sha=me,
                    python=platform.python_version(), numpy=np.__version__)
        json.dump(summ, open(out / "P2B_Q0_SUMMARY.json", "w"), indent=1); print(json.dumps(summ)); sys.exit(0 if ok else 1)
    k, n = (int(v) for v in a.shard.split("/")); tg = targets_for(ra, a.cohort, a.limit)[k::n]; out = Path(a.out) / f"shard_{k}"; check_cases(tg)
    if a.cohort == "c1c2":                                                                                   # production gate: only after freeze approval
        ap_ = json.load(open(a.approval, encoding="utf-8")) if a.approval else None
        want = dict(runner=me, config=sha(a.config), spec=sha(Path(__file__).with_name("P2B_SPEC_R1_1.md")), validator=sha(Path(__file__).with_name("P2B_validate_R1.py")))
        if not ap_ or any(ap_.get(x) != y for x, y in want.items()): print("APPROVAL_FAIL: freeze approval file missing or SHA mismatch", file=sys.stderr); sys.exit(3)
    out.mkdir(parents=True, exist_ok=False)                                                                  # write-once: an existing shard is never overwritten
    allrows, allpaths, cache = [], [], {}; t0 = time.time()
    for co, cid, bs in tg:
        br, pr = run_block(ra, co, cid, bs, cache); allrows += br; allpaths += pr
    with open(out / "P2B_RUN_BLOCKS.csv", "w", newline="") as f: w = csv.DictWriter(f, COLS); w.writeheader(); w.writerows(allrows)
    with open(out / "P2B_RUN_PATHS.csv", "w", newline="") as f: w = csv.DictWriter(f, PCOLS); w.writeheader(); w.writerows(allpaths)
    json.dump(dict(runner_sha=me, config_sha=sha(a.config), cohort=a.cohort, shard=a.shard, n_targets=len(tg), pythonhashseed=os.environ.get("PYTHONHASHSEED"), python=platform.python_version(),
                   numpy=np.__version__, manifest_digest=env["FZ"].manifest_digest(), seconds=round(time.time() - t0, 1), cfg=env["cfg"]), open(out / "shard_meta.json", "w"), indent=1)
    print(f"shard {a.shard}: {len(tg)} targets, {len(allrows)} block-tier rows, {time.time() - t0:.1f}s")


def merge(a):
    d = Path(a.out); BL, PA, metas = [], [], []
    for k in range(a.nshards):
        s = d / f"shard_{k}"
        if not (s / "shard_meta.json").exists(): raise SystemExit(f"missing shard {k}")
        BL += rows(s / "P2B_RUN_BLOCKS.csv"); PA += rows(s / "P2B_RUN_PATHS.csv"); metas.append(json.load(open(s / "shard_meta.json")))
    if len({(m["runner_sha"], m["config_sha"], m["manifest_digest"]) for m in metas}) != 1: raise SystemExit("shards disagree on runner/config/manifest")
    key = lambda r: (r["case_id"], r["block"], r["tier"]); BL.sort(key=key); PA.sort(key=lambda r: (r["case_id"], r["block"], r["tier"], int(r["step"])))
    if len({key(r) for r in BL}) != len(BL): raise SystemExit("duplicate block-tier rows")
    for nm, rr, cols in (("P2B_RUN_BLOCKS.csv", BL, COLS), ("P2B_RUN_PATHS.csv", PA, PCOLS)):
        with open(d / nm, "x", newline="") as f: w = csv.DictWriter(f, cols); w.writeheader(); w.writerows(rr)
    json.dump(dict(shards=metas, n_block_tier_rows=len(BL), n_blocks=len({(r["case_id"], r["block"]) for r in BL}), cohort=a.cohort), open(d / "P2B_RUN_META.json", "x"), indent=1)
    print("merged", len(BL), "rows,", len(PA), "path rows")


if __name__ == "__main__": main()
