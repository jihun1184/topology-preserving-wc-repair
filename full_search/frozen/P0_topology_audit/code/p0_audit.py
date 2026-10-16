#!/usr/bin/env python3
"""P0 topology-change audit of the stored ICEIC R+ edit trajectories (TP-WC project).

Model (definition).  Foreground X ⊂ Z^3 is read as the union K(X) of closed unit cubes (the (26,6) model of ICEIC; in the
Mazo et al. 2012 language: the support of the regular binary image zeta(X) on F^3 ordered by ⊇).
Toggling one voxel p changes K by the closed cube C=C(p):  K^+ = K ∪ C,  attachment A = K ∩ C  (K = complex WITHOUT p).
A depends only on the 26 neighbours of p.  Gate:   SAFE  <=>  A is contractible  (h0=1, h1=0, h2=0)
                                                   UNSAFE <=> otherwise  (class: empty / multi / loop / shell).
Independent confirmations (all verified exhaustively on 2^26 neighbourhoods by exhaustive_gate.cpp, see table file):
  Bertrand-Malandain simple point (T26=1 and T6bar=1)   and   greedy collapse of A to a point.
Per-edit Betti deltas are exact and incremental:  d_b0 = 1 - c_fg   (c_fg = # distinct 26-components touched, forward union-find),
  d_b2 = c_bg - 1  (c_bg = # distinct 6-components of the background touching p, reverse union-find), d_chi = 1 - chi(A),
  d_b1 = d_b0 + d_b2 - d_chi.

Usage (python -I):
  p0_audit.py selftest  --table T.bin --betti-module cubical_betti.py
  p0_audit.py replay    --p4 DIR --out DIR --table T.bin [--cases ID,ID | --all] [--workers N]
  p0_audit.py scratch   --p4 DIR --betti-module cubical_betti.py --cases ID,... [--random-small N]
"""
from __future__ import annotations
import argparse, csv, hashlib, importlib.util, json, sys, time
from functools import lru_cache
from multiprocessing import Pool
from pathlib import Path
import numpy as np

# ----------------------------------------------------------------------------------------------- gate (definition)
OFFS = [(a, b, c) for a in (-1, 0, 1) for b in (-1, 0, 1) for c in (-1, 0, 1) if (a, b, c) != (0, 0, 0)]  # bit k <-> OFFS[k]
CODE = [(i // 9, i // 3 % 3, i % 3) for i in range(27)]      # face of the cube C: per axis 0=low, 1=high, 2=full interval; id 26 = C itself
DIM = [c.count(2) for c in CODE]
_ID = lambda c: 9 * c[0] + 3 * c[1] + c[2]
_sub = lambda g, f: all(CODE[f][a] == 2 or CODE[g][a] == CODE[f][a] for a in range(3))      # face g ⊆ face f
# neighbour k shares with C the face f(k) (offset -1 -> low, +1 -> high, 0 -> full); it contributes the closure of f(k) to A
CL = [sum(1 << g for g in range(26) if _sub(g, _ID(tuple(2 if d == 0 else d > 0 for d in o)))) for o in OFFS]
POW2 = (1 << np.arange(26)).astype(np.int64)


@lru_cache(maxsize=None)
def attachment(nmask: int):
    """26-bit neighbour occupancy -> (V,E,S,h0,h1,h2) of A = K ∩ C.  h1 from chi = h0-h1+h2 (A ⊂ S^2 so h_k=0 for k>2; h2=1 iff A=∂C)."""
    A = 0
    for k in range(26):
        if nmask >> k & 1: A |= CL[k]
    V = E = S = 0; par = {}
    def find(x):                                   # union-find over the (<=8) vertices of A, path halving
        while par[x] != x: par[x] = par[par[x]]; x = par[x]
        return x
    for g in range(26):
        if A >> g & 1:
            if DIM[g] == 0: V += 1; par[g] = g
            elif DIM[g] == 1: E += 1
            else: S += 1
    comps = V
    for g in range(26):
        if A >> g & 1 and DIM[g] == 1:
            ax = CODE[g].index(2); lo = list(CODE[g]); hi = list(CODE[g]); lo[ax], hi[ax] = 0, 1
            u, w = _ID(lo), _ID(hi)
            if u not in par or w not in par: raise RuntimeError("attachment not closed under faces")     # error handling: closure invariant
            ru, rw = find(u), find(w)
            if ru != rw: par[ru] = rw; comps -= 1
    h2 = int(S == 6); h1 = comps + h2 - (V - E + S)
    if h1 < 0: raise RuntimeError(f"negative h1 for mask {nmask}")
    return V, E, S, comps, h1, h2


def classify(att):
    V, E, S, h0, h1, h2 = att
    if (h0, h1, h2) == (1, 0, 0): return "contractible"
    if V == 0: return "empty"          # A = ∅ : voxel touches nothing -> new component
    if h2: return "shell"              # A = ∂C : fills a single-voxel cavity
    if h0 >= 2: return "multi"         # A disconnected
    return "loop"                      # A connected, h1 >= 1


def nmask_of(vol, z, y, x):
    """26-neighbour occupancy bitmask of voxel (z,y,x) of an array that has >=1 voxel of zero padding around the data."""
    nb = (vol[z - 1:z + 2, y - 1:y + 2, x - 1:x + 2] > 0).ravel()
    return int(np.dot(np.delete(nb, 13).astype(np.int64), POW2))


def load_table(path):
    t = np.fromfile(path, dtype=np.uint8)
    if t.size != 1 << 26: raise ValueError(f"gate table has {t.size} entries, expected 2^26")
    return t


def verdict(att, table, nmask):
    """SAFE/UNSAFE/UNRESOLVED.  UNRESOLVED iff the independent confirmations (table bits: A-homology, BM, collapse) disagree."""
    a = int(att[3:] == (1, 0, 0)); bits = int(table[nmask]) & 7; ref = (bits & 1, bits >> 1 & 1, bits >> 2 & 1)
    if len(set(ref)) != 1 or ref[0] != a: return "UNRESOLVED"
    return "SAFE" if a else "UNSAFE"


# ----------------------------------------------------------------------------------------------- helpers
def load_zyx(p):
    import nibabel as nib
    return np.ascontiguousarray(np.transpose(np.asanyarray(nib.load(str(p)).dataobj), (2, 1, 0))).astype(np.uint8)


MASK_PREFIX = b"ICEIC2027_P4_MASK_R0\0"            # same hash convention as the frozen p4_contract.mask_sha256 (read, not imported)
def mask_sha(m):
    h = hashlib.sha256(); h.update(MASK_PREFIX); h.update(np.asarray(m.shape, dtype="<u8").tobytes()); h.update(np.ascontiguousarray(m, dtype=np.uint8).tobytes())
    return h.hexdigest()


def file_sha(p):
    h = hashlib.sha256()
    with open(p, "rb") as f:
        for ch in iter(lambda: f.read(1 << 20), b""): h.update(ch)
    return h.hexdigest()


def read_edits(p4, cid):
    rows = [r for r in csv.DictReader(open(Path(p4) / "edits" / f"{cid}_WT_EDITS.csv", encoding="utf-8-sig")) if r["condition"] == "WC"]
    rows.sort(key=lambda r: int(r["selection_index"]))
    if [int(r["selection_index"]) for r in rows] != list(range(len(rows))): raise ValueError(f"{cid}: edit indices not contiguous")
    return [(int(r["z"]), int(r["y"]), int(r["x"])) for r in rows]


SIX = [(1, 0, 0), (-1, 0, 0), (0, 1, 0), (0, -1, 0), (0, 0, 1), (0, 0, -1)]
def _uf(n):
    par = list(range(n))
    def find(a):
        while par[a] != a: par[a] = par[par[a]]; a = par[a]
        return a
    def union(a, b):
        a, b = find(a), find(b)
        if a != b: par[a] = b
    return find, union


def replay(raw, edits, table):
    """Replay stored add-only edits on RAW.  Returns per-edit dicts and the final padded-free state.  Raises on any contract violation."""
    from scipy import ndimage as ndi
    Xp = np.pad(raw, 1); n = len(edits); pe = [(z + 1, y + 1, x + 1) for z, y, x in edits]
    if len(set(pe)) != n: raise ValueError("duplicate edit coordinate")
    for z, y, x in pe:
        if not (0 < z < Xp.shape[0] - 1 and 0 < y < Xp.shape[1] - 1 and 0 < x < Xp.shape[2] - 1): raise ValueError("edit outside volume")
        if Xp[z, y, x]: raise ValueError("edit is not 0->1 on the replayed state (stored trajectory inconsistent)")
    # reverse pass: background 6-components of the FINAL state; un-adding each edit merges background components
    Xf = Xp.copy()
    for e in pe: Xf[e] = 1
    bg, nb = ndi.label(Xf == 0, structure=ndi.generate_binary_structure(3, 1))
    find, union = _uf(nb + 1 + n); cb = [0] * n; nxt = nb + 1
    for j in range(n - 1, -1, -1):
        z, y, x = pe[j]; ids = {find(int(bg[z + a, y + b, x + c])) for a, b, c in SIX if bg[z + a, y + b, x + c] > 0}
        cb[j] = len(ids); bg[z, y, x] = nxt
        for i in ids: union(nxt, i)
        nxt += 1
    # forward pass: foreground 26-components; gate evaluated on the pre-edit state
    fg, nf = ndi.label(Xp > 0, structure=np.ones((3, 3, 3), dtype=int))
    find, union = _uf(nf + 1 + n); out = []; nxt = nf + 1
    for j, (z, y, x) in enumerate(pe):
        t0 = time.perf_counter(); nm = nmask_of(fg, z, y, x); att = attachment(nm); v = verdict(att, table, nm)
        ids = {find(int(fg[z + a, y + b, x + c])) for a, b, c in OFFS if fg[z + a, y + b, x + c] > 0}
        V, E, S, h0, h1, h2 = att; dchi = 1 - (V - E + S); db0 = 1 - len(ids); db2 = cb[j] - 1; db1 = db0 + db2 - dchi
        out.append(dict(edit=j, z=z - 1, y=y - 1, x=x - 1, nmask=nm, V=V, E=E, F=S, h0=h0, h1=h1, h2=h2, cls=classify(att), verdict=v,
                        dchi=dchi, db0=db0, db1=db1, db2=db2, changed=int((db0, db1, db2) != (0, 0, 0)), gate_us=round((time.perf_counter() - t0) * 1e6)))
        fg[z, y, x] = nxt
        for i in ids: union(nxt, i)
        nxt += 1
    return out, Xf[1:-1, 1:-1, 1:-1]


_TABLE = None
def _init_worker(tablepath):
    global _TABLE; _TABLE = load_table(tablepath)       # once per worker process


def case_job(args):
    cid, p4, a, betti_row = args
    table = _TABLE; t0 = time.perf_counter(); P = Path(p4); res = dict(case_id=cid)
    try:
        raw = load_zyx(P / "raw" / f"{cid}_WT_RAW.nii.gz"); wc = load_zyx(P / "wc" / f"{cid}_WT_WC.nii.gz"); edits = read_edits(P, cid)
        rows, final = replay(raw, edits, table)
        res.update(n_edits=len(edits), prov_raw_nifti=file_sha(P / "raw" / f"{cid}_WT_RAW.nii.gz") == a["raw_nifti_sha256"],
                   prov_wc_nifti=file_sha(P / "wc" / f"{cid}_WT_WC.nii.gz") == a["wc_nifti_sha256"],
                   prov_edit_log=file_sha(P / "edits" / f"{cid}_WT_EDITS.csv") == a["edit_log_sha256"],
                   prov_raw_mask=mask_sha(raw) == a["raw_mask_sha256"], prov_wc_mask=mask_sha(final) == a["wc_mask_sha256"],
                   replay_equals_stored_wc=bool(np.array_equal(final, wc)))
        if betti_row:  # endpoint check against the independently computed stored RAW/WC Betti vectors
            res["endpoint_betti_match"] = [sum(r[k] for r in rows) for k in ("db0", "db1", "db2")] == betti_row
        v = [r["verdict"] for r in rows]; ch = [r["changed"] for r in rows]
        res.update(n_safe=v.count("SAFE"), n_unsafe=v.count("UNSAFE"), n_unresolved=v.count("UNRESOLVED"), n_changed=sum(ch),
                   first_unsafe=next((i for i, x in enumerate(v) if x != "SAFE"), -1), first_change=next((i for i, x in enumerate(ch) if x), -1),
                   false_safe=sum(1 for r in rows if r["verdict"] == "SAFE" and r["changed"]),
                   betti_invisible_unsafe=sum(1 for r in rows if r["verdict"] == "UNSAFE" and not r["changed"]),
                   seconds=round(time.perf_counter() - t0, 3), error="")
        for r in rows: r["case_id"] = cid
        return res, rows
    except Exception as e:                       # report, never silently skip
        res["error"] = f"{type(e).__name__}: {e}"; return res, []


def stored_refs(p4, betti_csv):
    audit = {r["case_id"]: r for r in csv.DictReader(open(Path(p4) / "case_materialization_audit.csv", encoding="utf-8-sig"))}
    betti = {}
    for r in csv.DictReader(open(betti_csv, encoding="utf-8-sig")):
        if r["arm"] == "WC": betti[r["case_id"]] = [int(r["delta_beta0"]), int(r["delta_beta1"]), int(r["delta_beta2"])]
    return dict(audit=audit, betti=betti)


def write_csv(path, rows, cols):
    with open(path, "w", newline="", encoding="utf-8") as f:
        w = csv.DictWriter(f, fieldnames=cols, extrasaction="ignore"); w.writeheader(); w.writerows(rows)


EDIT_COLS = ["case_id", "edit", "z", "y", "x", "nmask", "V", "E", "F", "h0", "h1", "h2", "cls", "verdict", "dchi", "db0", "db1", "db2", "changed", "gate_us"]
CASE_COLS = ["case_id", "n_edits", "n_safe", "n_unsafe", "n_unresolved", "n_changed", "first_unsafe", "first_change", "false_safe", "betti_invisible_unsafe",
             "prov_raw_nifti", "prov_wc_nifti", "prov_edit_log", "prov_raw_mask", "prov_wc_mask", "replay_equals_stored_wc", "endpoint_betti_match", "seconds", "error"]


def cmd_replay(a):
    p4 = Path(a.p4).resolve(); stored = stored_refs(p4, a.betti_csv); ids = sorted(stored["audit"]) if a.all else a.cases.split(",")
    for c in ids:
        if c not in stored["audit"]: sys.exit(f"unknown case {c}")
    jobs = [(c, str(p4), stored["audit"][c], stored["betti"].get(c)) for c in ids]
    t0 = time.time(); res = []; rows = []
    with Pool(a.workers, initializer=_init_worker, initargs=(a.table,)) as pool:
        for k, (r, e) in enumerate(pool.imap(case_job, jobs, chunksize=1)):
            res.append(r); rows += e
            if (k + 1) % 50 == 0: print(f"  {k + 1}/{len(jobs)} cases, {time.time() - t0:.0f}s", flush=True)
    out = Path(a.out); out.mkdir(parents=True, exist_ok=True)
    write_csv(out / "p0_edits.csv", rows, EDIT_COLS); write_csv(out / "p0_cases.csv", res, CASE_COLS)
    print(f"wrote {len(rows)} edit rows, {len(res)} case rows, errors={sum(1 for r in res if r['error'])}, {time.time() - t0:.0f}s")


# ----------------------------------------------------------------------------------------------- self-test (required synthetic tests)
def betti_fn(path):
    spec = importlib.util.spec_from_file_location("cubical_betti", path); m = importlib.util.module_from_spec(spec); sys.modules["cubical_betti"] = m; spec.loader.exec_module(m)
    return lambda vol: m.cubical_betti(vol).beta


def cmd_selftest(a):
    beta = betti_fn(a.betti_module); table = load_table(a.table); rng = np.random.default_rng(20261007); fails = 0

    def check(name, ok):
        nonlocal fails; fails += (not ok); print(f"  [{'PASS' if ok else 'FAIL'}] {name}")

    def toggle(vol, p):
        """gate verdict/class and from-scratch Betti(K ∪ C) - Betti(K), K = vol with voxel p cleared (so it covers addition and removal)."""
        base = vol.copy(); base[p] = 0; nm = nmask_of(base, *p); att = attachment(nm)
        full = base.copy(); full[p] = 1
        return verdict(att, table, nm), classify(att), tuple(int(i) for i in np.subtract(beta(full), beta(base))), att

    def vol(*cells, shape=(9, 9, 9)):
        v = np.zeros(shape, np.uint8)
        for c in cells: v[c] = 1
        return v
    bar = [(4, 4, 3), (4, 4, 4), (4, 4, 5)]
    v, c, d, _ = toggle(vol(*bar), (4, 4, 5)); check(f"safe addition (extend a bar by one voxel): {v}/{c}, d_beta={d}", v == "SAFE" and d == (0, 0, 0))
    v, c, d, _ = toggle(vol(*bar), (4, 4, 5)); check(f"safe removal (remove end of a bar): {v}/{c}, d_beta={d}", v == "SAFE" and d == (0, 0, 0))
    ring = [(4, 3, 4), (4, 5, 4), (4, 4, 3), (4, 4, 5)]               # 4 face neighbours of the centre: closed (edge-touching) ring around a vertical tunnel
    v, c, d, att = toggle(vol(*ring, (4, 4, 4)), (4, 4, 4)); check(f"loop closing (plug a tunnel): {v}/{c}, A(V,E,F,h0,h1,h2)={att}, d_beta={d}", v == "UNSAFE" and c == "loop" and att == (8, 12, 4, 1, 1, 0) and d == (0, -1, 0))
    shell = [(z, y, x) for z in (3, 4, 5) for y in (3, 4, 5) for x in (3, 4, 5)]
    v, c, d, att = toggle(vol(*shell), (4, 4, 4)); check(f"cavity filling: {v}/{c}, A={att}, d_beta={d}", v == "UNSAFE" and c == "shell" and att == (8, 12, 6, 1, 0, 1) and d == (0, 0, -1))
    v, c, d, att = toggle(vol((4, 4, 2), (4, 4, 6), (4, 4, 4)), (4, 4, 4)); check(f"isolated addition (new component): {v}/{c}, d_beta={d}", v == "UNSAFE" and c == "empty" and d == (1, 0, 0))
    v, c, d, att = toggle(vol((4, 4, 3), (4, 4, 5), (4, 4, 4)), (4, 4, 4)); check(f"bridge between two components: {v}/{c}, d_beta={d}", v == "UNSAFE" and c == "multi" and d == (-1, 0, 0))
    # 6 Python gate == C++ table (A bit) on random patterns of varying density
    bad = 0
    for dens in (0.05, 0.2, 0.5, 0.8, 0.95):
        for _ in range(40000):
            m = int(np.dot((rng.random(26) < dens).astype(np.int64), POW2)); bad += ((table[m] & 1) != int(attachment(m)[3:] == (1, 0, 0)))
    check(f"python gate == exhaustive C++ table on 200k random patterns (mismatches={bad})", bad == 0)
    # 7 theorem check on small random volumes: gate SAFE -> Betti unchanged;  record Betti-invisible UNSAFE examples
    stats = {"SAFE/same": 0, "SAFE/changed": 0, "UNSAFE/same": 0, "UNSAFE/changed": 0}; ex = None
    for _ in range(int(a.random_small)):
        v0 = (rng.random((8, 8, 8)) < rng.uniform(0.2, 0.7)).astype(np.uint8); v0 = np.pad(v0, 1)
        p = tuple(int(i) for i in rng.integers(2, 8, 3))
        if v0[p]: continue
        nm = nmask_of(v0, *p); vd = verdict(attachment(nm), table, nm); b0 = beta(v0); w = v0.copy(); w[p] = 1; d = tuple(np.subtract(beta(w), b0))
        stats[f"{vd}/{'same' if d == (0, 0, 0) else 'changed'}"] += 1
        if vd == "UNSAFE" and d == (0, 0, 0) and ex is None: ex = (v0.copy(), p, attachment(nm))
    print("  random small volumes (gate / from-scratch Betti):", stats)
    check("no SAFE edit changes the Betti vector", stats["SAFE/changed"] == 0)
    # 8 targeted search.  Betti unchanged => d_chi=0 => chi(A)=1.  So a Betti-invisible UNSAFE edit needs a NON-contractible A with chi(A)=1,
    #   e.g. (h0,h1,h2)=(2,1,0).  Sample such neighbourhoods and embed them in random 9^3 volumes; look for from-scratch d_beta == (0,0,0).
    if ex is None:
        tries = hits = 0
        while tries < a.targeted and ex is None:
            m = int(np.dot((rng.random(26) < rng.uniform(0.1, 0.4)).astype(np.int64), POW2)); att = attachment(m)
            if att[3:] not in ((2, 1, 0), (3, 2, 0)): continue
            tries += 1; v0 = (rng.random((9, 9, 9)) < rng.uniform(0.1, 0.5)).astype(np.uint8); v0[3:6, 3:6, 3:6] = 0
            for k, o in enumerate(OFFS):
                if m >> k & 1: v0[4 + o[0], 4 + o[1], 4 + o[2]] = 1
            w = v0.copy(); w[4, 4, 4] = 1; d = tuple(int(i) for i in np.subtract(beta(w), beta(v0)))
            if d == (0, 0, 0): ex = (v0, (4, 4, 4), att)
        print(f"  targeted search (chi(A)=1, non-contractible A): {tries} embeddings tried")
    print("  Betti-equal but UNSAFE example found:", bool(ex), "" if not ex else f"voxel {ex[1]}, attachment (V,E,F,h0,h1,h2)={ex[2]}")
    if ex: np.save(Path(a.out) / "betti_invisible_example.npy", ex[0]); json.dump(dict(voxel=ex[1], attachment=ex[2], note="from-scratch d_beta=(0,0,0) although attachment is not contractible"), open(Path(a.out) / "betti_invisible_example.json", "w"))
    sys.exit(1 if fails else 0)


def cmd_scratch(a):
    """Validate the incremental per-step Betti deltas against the frozen from-scratch implementation, step by step."""
    beta = betti_fn(a.betti_module); p4 = Path(a.p4).resolve(); stored = stored_refs(p4, a.betti_csv); table = load_table(a.table)
    ids = a.cases.split(",") if a.cases else []
    if a.random_small:
        rng = np.random.default_rng(7); pool = [c for c in sorted(stored["audit"]) if 0 < int(stored["audit"][c]["wc_edit_count"]) <= a.max_edits]
        ids += list(rng.choice(pool, a.random_small, replace=False))
    tot = bad = 0
    for cid in ids:
        raw = load_zyx(p4 / "raw" / f"{cid}_WT_RAW.nii.gz"); edits = read_edits(p4, cid); rows, _ = replay(raw, edits, table); X = raw.copy(); prev = beta(X)
        for r, (z, y, x) in zip(rows, edits):
            X[z, y, x] = 1; cur = beta(X); d = tuple(int(i) for i in np.subtract(cur, prev)); ok = d == (r["db0"], r["db1"], r["db2"]); tot += 1; bad += (not ok)
            if not ok: print("MISMATCH", cid, r["edit"], d, (r["db0"], r["db1"], r["db2"]))
            prev = cur
        print(f"  {cid}: {len(edits)} steps checked, RAW beta unknown-> final {prev}", flush=True)
    print(f"from-scratch step check: {tot} steps, mismatches={bad}"); sys.exit(1 if bad else 0)


if __name__ == "__main__":
    ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter); sub = ap.add_subparsers(dest="cmd", required=True)
    for name in ("selftest", "replay", "scratch"):
        s = sub.add_parser(name); s.add_argument("--table", required=True); s.add_argument("--out", default=".")
        s.add_argument("--betti-module"); s.add_argument("--p4"); s.add_argument("--betti-csv"); s.add_argument("--cases", default="")
        s.add_argument("--all", action="store_true"); s.add_argument("--workers", type=int, default=2)
        s.add_argument("--random-small", type=int, default=300); s.add_argument("--targeted", type=int, default=3000); s.add_argument("--max-edits", type=int, default=10)
    a = ap.parse_args()
    {"selftest": cmd_selftest, "replay": cmd_replay, "scratch": cmd_scratch}[a.cmd](a)
