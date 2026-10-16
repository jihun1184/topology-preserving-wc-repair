#!/usr/bin/env python3
"""P1: original-grid, add-only, simple-constrained R+ (TP-WC project).

Question.  Can R+ remove every critical block using only locally simple voxel additions?

Algorithm (R+ with ONE added gate; everything else is the frozen production R+ code, imported, not re-implemented)
  state X_t (binary, (z,y,x)), FIFO queue of critical blocks in the frozen total order, exactly as r_plus_repair_production.repair_r_plus.
  pop critical block B;  Z(B) = zero voxels of B (lexicographic order)
  S_t(B) = { p in Z(B) : A_{X_t}(p) = K(X_t) ∩ C_p is contractible }            (P0 gate, closed-cube (26,6) model)
  S_t(B) != {}  ->  add p* = lexicographically smallest element of S_t(B)           (R+ tie rule restricted to SAFE voxels)
  S_t(B) == {}  ->  DEFER B, store fingerprint of its local state R(B) = U_{q in Z(B)} N26[q];
                    B is re-queued only when an accepted edit changes a voxel of R(B) (fingerprint changes), never otherwise.
  queue empty   ->  every critical block left is deferred with an unchanged fingerprint: TERMINAL_STUCK.
  Any zero voxel of B resolves B (frozen R+ asserts this), so the gate restricts WHICH voxel resolves B, not whether one exists.

Betti numbers are NOT used by the algorithm; they are only a secondary check (independent P0 replay, per-step delta must be 0).

Usage (python -I):
  P1_simple_constrained_rplus.py selftest     --out DIR [--cases N]
  P1_simple_constrained_rplus.py illustrative --out DIR
  P1_simple_constrained_rplus.py all          --out DIR [--workers 2] [--no-control]
  P1_simple_constrained_rplus.py variants     --out DIR --base P1_CASES_876.csv   (rule-sensitivity diagnostic on the TERMINAL_STUCK cases; not part of P1)
common: --vendor DIR(r_plus_repair_production.py)  --p0 p0_audit.py  --table gate_table_2p26.bin  --p4 DIR  --p0e p0e_candidates.csv  --p0edits p0_edits.csv
"""
from __future__ import annotations
import argparse, csv, hashlib, importlib.util, itertools, json, sys, time
from collections import Counter, defaultdict, deque
from multiprocessing import Pool
from pathlib import Path
import numpy as np

R = P0 = TABLE = None            # frozen modules / gate table, set by setup()
P0E = {}                         # (case, edit) -> row of P0 p0e_candidates.csv
ILLUSTRATIVE = ("BraTS2021_00451", "BraTS2021_01152", "BraTS2021_01563")


class InvariantError(RuntimeError):
    """Provenance / invariant / implementation-consistency failure.  Always reported as status ERROR, never as STUCK."""


def _load(name, path):
    spec = importlib.util.spec_from_file_location(name, str(path)); m = importlib.util.module_from_spec(spec); sys.modules[name] = m; spec.loader.exec_module(m); return m


def setup(vendor, p0path, tablepath, p0e=None):
    global R, P0, TABLE, P0E
    sys.path.insert(0, str(Path(vendor).resolve()))                  # production R+ imports dwc_checker(_fast) by bare name
    import r_plus_repair_production as r_mod
    R = r_mod; P0 = _load("p0_audit", p0path); TABLE = P0.load_table(tablepath)
    if p0e:
        P0E = {(r["case_id"], int(r["edit"])): r for r in csv.DictReader(open(p0e, encoding="utf-8-sig"))}


# ------------------------------------------------------------------------------------------------ small helpers
OFF27 = [(a, b, c) for a in (-1, 0, 1) for b in (-1, 0, 1) for c in (-1, 0, 1)]
bkey = lambda b: f"k{b.k}{b.plane}@{b.z},{b.y},{b.x}"
vstr = lambda c: f"{c[0]},{c[1]},{c[2]}"


def gate(Xp, p):
    """P0 gate for 0-voxel p (unpadded coords) on padded state Xp -> (nmask, attachment, verdict, class).  UNRESOLVED (table bits disagree) is an error."""
    nm = P0.nmask_of(Xp, p[0] + 1, p[1] + 1, p[2] + 1); att = P0.attachment(nm); v = P0.verdict(att, TABLE, nm)
    if v == "UNRESOLVED": raise InvariantError(f"gate implementations disagree at {p} nmask={nm}")
    return nm, att, v, P0.classify(att)


def region(zeros, shape):
    """R(B) = union of the closed 26-neighbourhoods of the zero voxels of B, restricted to the domain (outside is constant 0)."""
    s = {(z + a, y + b, x + c) for z, y, x in zeros for a, b, c in OFF27}
    return sorted(q for q in s if 0 <= q[0] < shape[0] and 0 <= q[1] < shape[1] and 0 <= q[2] < shape[2])


def fingerprint(x, reg):
    arr = np.asarray(reg, dtype=np.int16); h = hashlib.blake2b(digest_size=8)
    h.update(arr.tobytes()); h.update(x[arr[:, 0], arr[:, 1], arr[:, 2]].tobytes()); return h.hexdigest()


# ------------------------------------------------------------------------------------------------ the algorithm
def run_p1(raw, *, gate_on=True, watch=frozenset(), trace=None, cid="", select="lexmin", lifo=False):
    """Run (gated) R+.  gate_on=False reproduces the frozen production R+ exactly (control).  Returns a dict; raises InvariantError on any inconsistency.
    select/lifo are ONLY for the rule-sensitivity diagnostic (variants command); the P1 algorithm is select='lexmin', lifo=False (FIFO queue, lexicographically smallest SAFE voxel).
    select in {lexmin, lexmax, minnew, rand<seed>}: which SAFE voxel of the popped block is added."""
    rng = np.random.default_rng(int(select[4:])) if select.startswith("rand") else None
    Xp = np.pad(np.asarray(raw, np.uint8), 1); x = Xp[1:-1, 1:-1, 1:-1]; shape = x.shape            # x is a VIEW of Xp: one mutable state, two index frames
    init = R._critical_keys_vectorized(x); initial_set = set(init)
    q = deque(init); queued = set(init); ncrit = len(init)
    deferred = {}                                    # BlockKey -> (fingerprint, region)
    watchers = defaultdict(set)                      # in-domain voxel -> deferred blocks whose region contains it
    ever_def, ever_react, released, def_res_other = set(), set(), set(), set(); dcount = Counter(); born = {k: -1 for k in init}
    st = Counter(); edits, recs, H = [], [], {}
    maxcorr = int(np.count_nonzero(x == 0)); max_pops = len(init) + 600 * maxcorr + 10      # runaway guard only (a pop that re-queues is bounded by #watchers)

    def hh(k): return H.setdefault(k, dict(epis=0, own=0, other=0, own_after_defer=False, defers=0, pops=[]))
    for k in init:
        if k in watch: hh(k)["epis"] += 1

    while q:
        st["pops"] += 1
        if st["pops"] > max_pops: raise InvariantError("QUEUE_POP_BOUND_EXCEEDED")
        b = q.pop() if lifo else q.popleft(); queued.remove(b)
        if not R._is_critical(x, b): st["stale"] += 1; continue
        zeros = sorted(c for c in R._block_coords(b) if x[c] == 0)
        if not zeros: raise InvariantError(f"critical block {bkey(b)} without zero voxel")
        gz = [(c,) + gate(Xp, c) for c in zeros]                                   # (voxel, nmask, att, verdict, class)
        safe = [g for g in gz if g[3] == "SAFE"] if gate_on else gz
        reg = region(zeros, shape) if (b in watch or not safe) else None
        fp = fingerprint(x, reg) if reg is not None else None
        if b in watch: hh(b)["pops"].append((len(edits), len(safe), fp))
        zdesc = ";".join(f"{vstr(g[0])}:{g[3]}/{g[4]}" for g in gz)
        if not safe:                                                               # ---- DEFER
            deferred[b] = (fp, reg)
            for v_ in reg: watchers[v_].add(b)
            st["defer_ev"] += 1; ever_def.add(b); dcount[b] += 1
            if b in watch: hh(b)["defers"] += 1
            if trace is not None:
                trace.append(dict(case_id=cid, seq=len(trace), event="DEFER", block=bkey(b), n_crit=ncrit, rplus_voxel=vstr(zeros[0]), rplus_gate=f"{gz[0][3]}/{gz[0][4]}",
                                  zero_voxels=zdesc, safe_set="", p1_voxel="", deferred_size=len(deferred), cum_edits=len(edits), stuck_now=""))
            continue
        if select == "lexmin" or len(safe) == 1: g = safe[0]
        elif select == "lexmax": g = safe[-1]
        elif rng is not None: g = safe[int(rng.integers(len(safe)))]
        elif select == "minnew":                                                   # one-step lookahead: fewest newly created critical blocks, ties -> lexmin
            def _nn(c):
                inc_ = R._incident_block_keys(shape, c[0]); bf = {k: R._is_critical(x, k) for k in inc_}; x[c[0]] = 1
                n_ = sum(R._is_critical(x, k) and not bf[k] for k in inc_); x[c[0]] = 0; return n_
            g = min(safe, key=lambda c: (_nn(c), c[0]))
        else: raise ValueError(f"unknown select {select}")
        p, nm, att, v, cls = g
        if x[p] != 0 or p not in zeros or (gate_on and v != "SAFE"): raise InvariantError(f"bad selection {p}")
        inc = R._incident_block_keys(shape, p)
        before = {k: R._is_critical(x, k) for k in inc}
        x[p] = 1
        after = {k: R._is_critical(x, k) for k in inc}
        newc = [k for k in inc if after[k] and not before[k]]; resolved = [k for k in inc if before[k] and not after[k]]
        if b not in resolved: raise InvariantError(f"{bkey(b)} not resolved by its own edit {p}")
        n_before = ncrit; ncrit += len(newc) - len(resolved)
        react = []
        for k in list(watchers.get(p, ())):                                        # deferred blocks whose local state just changed
            fp0, reg0 = deferred.pop(k)
            for v_ in reg0: watchers[v_].discard(k)
            if k in after: crit = after[k]
            else: crit = R._is_critical(x, k)
            if crit:
                if fingerprint(x, reg0) == fp0: raise InvariantError("reactivated block with unchanged fingerprint")
                react.append(k); ever_react.add(k); st["react_ev"] += 1
        for k in resolved:
            if k != b and k in ever_def: def_res_other.add(k)
        if b in ever_def: released.add(b)
        for k in newc:
            born[k] = len(edits)                                                    # index of the accepted edit that (re)created this critical block
            if k in watch: hh(k)["epis"] += 1
        for k in resolved:
            if k in watch:
                h = hh(k)
                if k == b: h["own"] += 1; h["own_after_defer"] |= (h["defers"] > 0)
                else: h["other"] += 1
        enq = sorted({k for k in inc if after[k] and k not in queued} | set(react), key=R._sort_key)
        for k in enq:
            if k not in queued: q.append(k); queued.add(k)
        st["new_crit"] += len(newc); st["edits_new_crit"] += bool(newc)
        edits.append(p)
        recs.append(dict(case_id=cid, step=len(edits) - 1, z=p[0], y=p[1], x=p[2], block=bkey(b), rplus_voxel=vstr(zeros[0]), p1_is_rplus=int(p == zeros[0]),
                         n_zero=len(zeros), n_safe=len(safe), nmask=nm, V=att[0], E=att[1], F=att[2], h0=att[3], h1=att[4], h2=att[5], verdict=v, cls=cls,
                         gate_cert=f"{v}/{cls}/tbl{int(TABLE[nm]) & 7}", n_crit_before=n_before, n_crit_after=ncrit, n_new_crit=len(newc), n_resolved=len(resolved),
                         n_react=len(react), deferred_size=len(deferred)))
        if trace is not None:
            trace.append(dict(case_id=cid, seq=len(trace), event="ACCEPT", block=bkey(b), n_crit=n_before, rplus_voxel=vstr(zeros[0]), rplus_gate=f"{gz[0][3]}/{gz[0][4]}",
                              zero_voxels=zdesc, safe_set=";".join(vstr(s[0]) for s in safe), p1_voxel=vstr(p), V=att[0], E=att[1], F=att[2], h0=att[3], h1=att[4], h2=att[5],
                              gate=f"{v}/{cls}", resolved_blocks="|".join(bkey(k) for k in resolved), new_crit_blocks="|".join(bkey(k) for k in newc),
                              reactivated="|".join(bkey(k) for k in react), deferred_size=len(deferred), cum_edits=len(edits), stuck_now=""))

    # ---- post-conditions (invariants 5, 9): bookkeeping == independent recount
    final_crit = R._critical_keys_vectorized(x)
    if set(final_crit) != set(deferred): raise InvariantError(f"final critical set != deferred set ({len(final_crit)} vs {len(deferred)})")
    if len(final_crit) != ncrit: raise InvariantError(f"incremental critical count {ncrit} != recount {len(final_crit)}")
    stuck = []
    for b in sorted(deferred, key=R._sort_key):
        fp0, reg0 = deferred[b]
        if fingerprint(x, reg0) != fp0: raise InvariantError("stuck block fingerprint changed without reactivation")
        zs = sorted(c for c in R._block_coords(b) if x[c] == 0)
        gz = [(c,) + gate(Xp, c) for c in zs]
        if gate_on and any(g[3] == "SAFE" for g in gz): raise InvariantError(f"stuck block {bkey(b)} has a SAFE candidate")
        stuck.append((b, gz))
        if trace is not None:
            trace.append(dict(case_id=cid, seq=len(trace), event="TERMINAL_STUCK", block=bkey(b), n_crit=ncrit, zero_voxels=";".join(f"{vstr(g[0])}:{g[3]}/{g[4]}" for g in gz),
                              deferred_size=len(deferred), cum_edits=len(edits), stuck_now=1))
    n2f = sum(b.k == 2 for b in final_crit)
    return dict(edits=edits, recs=recs, final=x.copy(), Xp=Xp, init=init, initial_set=initial_set, n2i=sum(b.k == 2 for b in init), n3i=sum(b.k == 3 for b in init),
                n2f=n2f, n3f=len(final_crit) - n2f, stuck=stuck, stuck_keys=set(deferred), H=H, st=st, ever_def=ever_def, dcount=dcount, born=born, ever_react=ever_react, released=released,
                def_res_other=def_res_other)


# ------------------------------------------------------------------------------------------------ original R+ trajectory (replay, P0-E reconstruction)
def orig_states(cid, raw, edits0):
    """Replay the stored R+ trajectory; for each edit return the candidate blocks (BlockKey) reconstructed exactly as in P0-E, with alternatives' verdicts."""
    Xp = np.pad(np.asarray(raw, np.uint8), 1); x = Xp[1:-1, 1:-1, 1:-1]; out = []
    for j, p in enumerate(edits0):
        if x[p] != 0: raise InvariantError(f"{cid}: stored edit {j} not 0->1")
        nm, att, v, cls = gate(Xp, p); cands = []
        for k in R._incident_block_keys(x.shape, p):
            if not R._is_critical(x, k): continue
            zs = sorted(c for c in R._block_coords(k) if x[c] == 0)
            if zs[0] != p: continue
            alts = [gate(Xp, c)[2] for c in zs[1:]]
            cands.append((k, zs, alts))
        any_safe = int(any("SAFE" in a for _, _, a in cands)); all_safe = int(bool(cands) and all("SAFE" in a for _, _, a in cands))
        row = P0E.get((cid, j))
        if row is not None and (int(row["n_candidate_blocks"]) != len(cands) or int(row["any_safe_alt"]) != any_safe or int(row["all_safe_alt"]) != all_safe
                                or int(row["n_cand_with_safe_alt"]) != sum("SAFE" in a for _, _, a in cands) or row["verdict"] != v
                                or sorted(row["cand_kinds"].split("|")) != sorted(("cube" if k.k == 3 else k.plane) for k, _, _ in cands)):
            raise InvariantError(f"{cid} edit {j}: reconstruction disagrees with P0 p0e_candidates.csv")
        fps = {}
        if v == "UNSAFE":                                           # original local state of each candidate block (pre-edit) and whether S(B) was empty there
            for k, zs, al in cands: fps[k] = (fingerprint(x, region(zs, x.shape)), "SAFE" not in al)
        out.append(dict(edit=j, p=p, verdict=v, cls=cls, cands=cands, any_safe=any_safe, all_safe=all_safe, fps=fps)); x[p] = 1
    return out, x.copy()


def fate(k, h, stuck_keys, orig):
    """Block-level P1 fate of an original candidate block: (A|B|C, detail, met_same_local_state).
    orig = (fingerprint of R(B) at the original pre-edit state, S(B) was empty there).  The gate depends only on R(B), so an identical fingerprint with S(B)=empty must give S=empty."""
    orig_fp, orig_empty = orig
    same = bool(h) and any(fp == orig_fp for _, _, fp in h["pops"])
    if h and orig_empty and any(fp == orig_fp and ns > 0 for _, ns, fp in h["pops"]):
        raise InvariantError(f"{bkey(k)}: SAFE candidate at a local state identical to the original no-SAFE state")
    if k in stuck_keys: return "C", "TERMINAL_STUCK", same
    if h and h["own"] > 0: return "B", "RELEASED_AFTER_DEFER" if h["own_after_defer"] else "RELEASED_DIRECT", same
    if h and h["other"] > 0: return "A", "RESOLVED_BY_OTHER_EDIT", same
    if h and h["epis"] > 0: raise InvariantError(f"{bkey(k)} critical in P1 but neither resolved nor stuck")
    return "A", "NEVER_CRITICAL", same


# ------------------------------------------------------------------------------------------------ stuck-block diagnostics
def unlock_r1(Xp, b, max_zeros=12):
    """DIAGNOSTIC (necessary condition only).  Exhaustive DFS over SAFE additions of zero voxels inside the radius-1 window of B (the only voxels that can change the
    gate of a zero voxel of B).  'UNLOCKABLE': some sequence of such additions makes a zero voxel of B SAFE.  This does NOT check that the case is then solved."""
    X = Xp.copy(); sh = X.shape; bc = list(R._block_coords(b)); zs = [c for c in bc if X[c[0] + 1, c[1] + 1, c[2] + 1] == 0]
    lo = np.min(bc, 0) - 1; hi = np.max(bc, 0) + 1
    W = [(z, y, x) for z in range(lo[0], hi[0] + 1) for y in range(lo[1], hi[1] + 1) for x in range(lo[2], hi[2] + 1)
         if 0 <= z < sh[0] - 2 and 0 <= y < sh[1] - 2 and 0 <= x < sh[2] - 2 and X[z + 1, y + 1, x + 1] == 0 and (z, y, x) not in zs]
    if len(W) > max_zeros: return "WINDOW_TOO_LARGE", len(W), 0
    seen, stack = set(), [()]
    while stack:
        added = stack.pop(); key = frozenset(added)
        if key in seen: continue
        seen.add(key)
        for c in added: X[c[0] + 1, c[1] + 1, c[2] + 1] = 1
        try:
            if any(gate(X, c)[2] == "SAFE" for c in zs): return "UNLOCKABLE", len(W), len(seen)
            for c in W:
                if c not in key and gate(X, c)[2] == "SAFE": stack.append(added + (c,))
        finally:
            for c in added: X[c[0] + 1, c[1] + 1, c[2] + 1] = 0
    return "NOT_UNLOCKABLE_R1", len(W), len(seen)


def stuck_rows(cid, res, orig_added, raw):
    """Per terminal-stuck block: local pattern, per-candidate gate class (final and RAW state), exact Delta-Betti of the (forbidden) addition (P0 incremental formulas,
    single voxel), category, removal-gate diagnostic and radius-1 unlock diagnostic."""
    if not res["stuck"]: return []
    from scipy import ndimage as ndi
    Xp, x = res["Xp"], res["final"]; shape = x.shape; Xraw = np.pad(np.asarray(raw, np.uint8), 1)
    fg, _ = ndi.label(Xp > 0, structure=np.ones((3, 3, 3), int)); six = [(1, 0, 0), (-1, 0, 0), (0, 1, 0), (0, -1, 0), (0, 0, 1), (0, 0, -1)]
    nz = np.argwhere(Xp > 0); lo, hi = nz.min(0), nz.max(0); rows = []
    S6 = ndi.generate_binary_structure(3, 1)
    for b, gz in res["stuck"]:
        cand = []
        for (c, nm, att, v, cls) in gz:
            pz, py, px = c[0] + 1, c[1] + 1, c[2] + 1
            ids = {int(fg[pz + a, py + bb, px + cc]) for a, bb, cc in P0.OFFS if fg[pz + a, py + bb, px + cc] > 0}
            l = np.minimum(lo, [pz, py, px]) - 1; h = np.maximum(hi, [pz, py, px]) + 2; l = np.maximum(l, 0); h = np.minimum(h, Xp.shape)
            sub = Xp[l[0]:h[0], l[1]:h[1], l[2]:h[2]].copy(); sub[pz - l[0], py - l[1], px - l[2]] = 1
            lab, _ = ndi.label(sub == 0, structure=S6)
            cb = {int(lab[pz - l[0] + a, py - l[1] + bb, px - l[2] + cc]) for a, bb, cc in six if lab[pz - l[0] + a, py - l[1] + bb, px - l[2] + cc] > 0}
            db0 = 1 - len(ids); db2 = len(cb) - 1; dchi = 1 - (att[0] - att[1] + att[2]); db1 = db0 + db2 - dchi
            cand.append(dict(v=c, nm=nm, att=att, cls=cls, d=(db0, db1, db2), orig=int(c in orig_added), raw_cls=gate(Xraw, c)[3]))
        ones = [c for c in R._block_coords(b) if x[c] == 1]
        rem = [(c,) + gate(Xp, c)[2:] for c in ones]                              # removal gate: o in X is removed; A = K(X \ o) ∩ C_o depends on the 26 neighbours only
        raw_all_shell = all(c["raw_cls"] == "shell" for c in cand); fin_all_shell = all(c["cls"] == "shell" for c in cand); born = res["born"].get(b, -1)
        cat = "RAW_SHELL_CLOSED" if (born == -1 and b in res["initial_set"] and raw_all_shell) else "FINAL_SHELL_CLOSED" if fin_all_shell else "OPEN"
        ul = unlock_r1(Xp, b) if cat == "OPEN" else ("PERMANENT_SHELL", 0, 0)
        pat = "".join(str(int(x[c])) for c in R._block_coords(b))
        rows.append(dict(case_id=cid, block=bkey(b), k=b.k, plane=b.plane or "cube", z=b.z, y=b.y, x=b.x, initial_critical=int(b in res["initial_set"]), born_step=born, block_pattern=pat,
                         n_zero=len(gz), category=cat, cand_classes="|".join(sorted(c["cls"] for c in cand)), raw_classes="|".join(sorted(c["raw_cls"] for c in cand)),
                         cand_detail=";".join(f"{vstr(c['v'])}:{c['cls']}:(V{c['att'][0]},E{c['att'][1]},F{c['att'][2]};h{c['att'][3]}{c['att'][4]}{c['att'][5]}):dB{c['d']}:nm{c['nm']:#x}:origR+edit{c['orig']}" for c in cand),
                         cand_dbetti="|".join(str(c["d"]) for c in cand), orig_rplus_added_any=int(any(c["orig"] for c in cand)),
                         n_ones=len(ones), rem_safe_resolving=sum(1 for r_ in rem if r_[1] == "SAFE"), rem_classes="|".join(sorted(r_[2] for r_ in rem)),
                         unlock_r1=ul[0], unlock_window_zeros=ul[1], unlock_states=ul[2], n_defer_events=res["dcount"][b]))
    return rows


# ------------------------------------------------------------------------------------------------ per-case job
_G = {}
def _init_worker(vendor, p0path, tablepath, p0e, opts):
    setup(vendor, p0path, tablepath, p0e); _G.update(opts)


def case_job(args):
    cid, p4, audit = args; t0 = time.perf_counter(); P = Path(p4)
    out = dict(case=dict(case_id=cid, status="ERROR", error=""), edits=[], stuck=[], fates=[])
    try:
        raw = P0.load_zyx(P / "raw" / f"{cid}_WT_RAW.nii.gz"); edits0 = P0.read_edits(P, cid)
        prov = dict(prov_raw_mask=P0.mask_sha(raw) == audit["raw_mask_sha256"], prov_edit_log=P0.file_sha(P / "edits" / f"{cid}_WT_EDITS.csv") == audit["edit_log_sha256"])
        if not all(prov.values()): raise InvariantError(f"provenance failed {prov}")
        if len(edits0) != int(audit["wc_edit_count"]): raise InvariantError("stored edit count != audit")
        ctrl_ok = ""; rp = {}
        if _G.get("control", True):                                                # gate OFF must reproduce the stored production-R+ trajectory exactly
            c = run_p1(raw, gate_on=False, cid=cid); ctrl_ok = int(c["edits"] == edits0)
            if not ctrl_ok: raise InvariantError("gate-OFF run != stored R+ trajectory")
            rp = dict(rplus_new_crit=c["st"]["new_crit"], rplus_edits_new_crit=c["st"]["edits_new_crit"], rplus_pops=c["st"]["pops"])
        states, wc = orig_states(cid, raw, edits0)
        watch = {k for s in states if s["verdict"] == "UNSAFE" for k, _, _ in s["cands"]}
        r = run_p1(raw, gate_on=True, watch=watch, cid=cid)
        # ---- independent verification (P0 replay: forward/backward union-find, gate re-evaluated on the replayed pre-edit state)
        rows, fin = P0.replay(raw, r["edits"], TABLE) if r["edits"] else ([], np.asarray(raw, np.uint8))
        if not np.array_equal(fin, r["final"]): raise InvariantError("replayed trajectory != final mask")                   # invariant 8
        if len(set(r["edits"])) != len(r["edits"]) or np.any((raw == 1) & (r["final"] == 0)): raise InvariantError("duplicate edit / monotonicity")   # 7
        bad = [w for w in rows if w["verdict"] != "SAFE"]
        if bad: raise InvariantError(f"non-SAFE accepted edit in independent replay: {bad[0]['edit']}")
        betti = tuple(sum(w[k] for w in rows) for k in ("db0", "db1", "db2")); chg = sum(w["changed"] for w in rows)
        if chg: raise InvariantError("a SAFE accepted edit changed the Betti vector (theorem violated)")
        cen = R.census_dwc(r["final"])                                              # padded census (secondary): must equal the in-domain recount for these masks
        out["edits"] = r["recs"]
        srows = stuck_rows(cid, r, set(edits0), raw); scat = {w["block"]: w["category"] for w in srows}
        Xraw = np.pad(np.asarray(raw, np.uint8), 1)                                  # RAW-level permanent obstruction check (monotonicity lemma): initial critical block whose every zero is shell in RAW
        rawclosed = [b for b in r["init"] if all(gate(Xraw, c)[3] == "shell" for c in R._block_coords(b) if raw[c] == 0)]
        if any(b not in r["stuck_keys"] for b in rawclosed): raise InvariantError("RAW-shell-closed initial block was resolved (contradicts monotonicity lemma)")
        # ---- stress cohort / fate of every original UNSAFE state
        orig_added = set(r["edits"]); st = r["st"]
        for s in states:
            if s["verdict"] != "UNSAFE": continue
            fs = []
            for k, zs, alts in s["cands"]:
                f = fate(k, r["H"].get(k), r["stuck_keys"], s["fps"][k]); fs.append((k, f))
            cnt = Counter(f[0] for _, f in fs); n = len(fs)
            cls_any = "C" if cnt["C"] else "B" if cnt["B"] else "A"
            cls_all = "C" if cnt["C"] == n else "A" if cnt["A"] == n else "B" if cnt["B"] == n else "MIXED"
            out["fates"].append(dict(case_id=cid, orig_edit=s["edit"], voxel=vstr(s["p"]), unsafe_class=s["cls"], stress=int(s["any_safe"] == 0), n_cand=n, n_A=cnt["A"], n_B=cnt["B"], n_C=cnt["C"],
                                     fate_any=cls_any, fate_all=cls_all, cand_fates="|".join(f"{bkey(k)}:{f[0]}:{f[1]}{('/' + scat[bkey(k)]) if f[0] == 'C' else ''}:same{int(f[2])}" for k, f in fs),
                                     met_same_local_state=int(any(f[2] for _, f in fs)), voxel_added_in_p1=int(s["p"] in orig_added)))
        out["stuck"] = srows
        nst = len(r["stuck"]); status = "SOLVED" if nst == 0 else "TERMINAL_STUCK"
        if (cen.n2, cen.n3) != (r["n2f"], r["n3f"]): raise InvariantError(f"padded census ({cen.n2},{cen.n3}) != in-domain recount ({r['n2f']},{r['n3f']})")
        out["case"] = dict(case_id=cid, status=status, error="", n_rplus_edits=len(edits0), n_p1_edits=len(r["edits"]),
                           edit_inflation=(round(len(r["edits"]) / len(edits0), 4) if edits0 else ""), n_unsafe_rplus=sum(s["verdict"] == "UNSAFE" for s in states),
                           init_crit_n2=r["n2i"], init_crit_n3=r["n3i"], init_crit=r["n2i"] + r["n3i"], final_crit=nst, final_n2=r["n2f"], final_n3=r["n3f"],
                           n_deferred_events=st["defer_ev"], n_deferred_blocks=len(r["ever_def"]), n_reactivation_events=st["react_ev"], n_reactivated_blocks=len(r["ever_react"]),
                           n_released_blocks=len(r["released"]), n_raw_shell_closed_initial=len(rawclosed), n_deferred_resolved_by_other=len(r["def_res_other"]), n_terminal_stuck=nst,
                           n_new_crit=st["new_crit"], n_edits_new_crit=st["edits_new_crit"], queue_pops=st["pops"], stale_pops=st["stale"], ctrl_equals_stored_rplus=ctrl_ok,
                           replay_ok=1, betti_delta_sum=str(betti), n_betti_changing_edits=chg, **rp, **prov, seconds=round(time.perf_counter() - t0, 2))
        out["_final_sha"] = P0.mask_sha(r["final"]); out["case"]["p1_mask_sha256"] = out["_final_sha"]
    except Exception as e:                                                          # ERROR is never merged with STUCK
        out["case"]["error"] = f"{type(e).__name__}: {e}"; out["case"].setdefault("seconds", round(time.perf_counter() - t0, 2))
    return out


# ------------------------------------------------------------------------------------------------ writers / summary
def write_csv(path, rows, cols=None):
    cols = cols or (list(rows[0]) if rows else [])
    with open(path, "w", newline="", encoding="utf-8") as f:
        w = csv.DictWriter(f, fieldnames=cols, extrasaction="ignore"); w.writeheader(); w.writerows(rows)


def audit_rows(p4):
    return {r["case_id"]: r for r in csv.DictReader(open(Path(p4) / "case_materialization_audit.csv", encoding="utf-8-sig"))}


def cmd_all(a):
    out = Path(a.out); out.mkdir(parents=True, exist_ok=True); audit = audit_rows(a.p4); ids = sorted(audit) if not a.cases else a.cases.split(",")
    jobs = [(c, a.p4, audit[c]) for c in ids]; t0 = time.time(); res = []
    with Pool(a.workers, initializer=_init_worker, initargs=(a.vendor, a.p0, a.table, a.p0e, dict(control=not a.no_control))) as pool:
        for k, r in enumerate(pool.imap(case_job, jobs, chunksize=1)):
            res.append(r)
            if (k + 1) % 25 == 0: print(f"  {k + 1}/{len(jobs)} cases, {time.time() - t0:.0f}s", flush=True)
    cases = [r["case"] for r in res]; edits = [e for r in res for e in r["edits"]]; stuck = [s for r in res for s in r["stuck"]]; fates = [f for r in res for f in r["fates"]]
    CASE_COLS = ["case_id", "status", "n_rplus_edits", "n_p1_edits", "edit_inflation", "n_unsafe_rplus", "init_crit_n2", "init_crit_n3", "init_crit", "final_crit", "final_n2", "final_n3",
                 "n_deferred_events", "n_deferred_blocks", "n_reactivation_events", "n_reactivated_blocks", "n_released_blocks", "n_raw_shell_closed_initial", "n_deferred_resolved_by_other", "n_terminal_stuck",
                 "n_new_crit", "n_edits_new_crit", "rplus_new_crit", "rplus_edits_new_crit", "queue_pops", "stale_pops", "rplus_pops", "ctrl_equals_stored_rplus", "replay_ok",
                 "betti_delta_sum", "n_betti_changing_edits", "prov_raw_mask", "prov_edit_log", "p1_mask_sha256", "seconds", "error"]
    write_csv(out / "P1_CASES_876.csv", cases, CASE_COLS)
    write_csv(out / "P1_EDITS_876.csv", edits, ["case_id", "step", "z", "y", "x", "block", "rplus_voxel", "p1_is_rplus", "n_zero", "n_safe", "nmask", "V", "E", "F", "h0", "h1", "h2", "verdict", "cls",
                                                  "gate_cert", "n_crit_before", "n_crit_after", "n_new_crit", "n_resolved", "n_react", "deferred_size"])
    FATE_COLS = ["case_id", "orig_edit", "voxel", "unsafe_class", "stress", "n_cand", "n_A", "n_B", "n_C", "fate_any", "fate_all", "cand_fates", "met_same_local_state", "voxel_added_in_p1"]
    write_csv(out / "P1_UNSAFE_FATE_ALL.csv", fates, FATE_COLS); write_csv(out / "P1_STRESS_2413.csv", [f for f in fates if f["stress"]], FATE_COLS)
    write_csv(out / "P1_TERMINAL_STUCK.csv", stuck, ["case_id", "block", "k", "plane", "z", "y", "x", "initial_critical", "born_step", "block_pattern", "n_zero", "category", "cand_classes", "raw_classes",
                                                    "cand_detail", "cand_dbetti", "orig_rplus_added_any", "n_ones", "rem_safe_resolving", "rem_classes", "unlock_r1", "unlock_window_zeros",
                                                    "unlock_states", "n_defer_events"])
    print(f"wrote {len(cases)} cases, {len(edits)} edits, {len(fates)} unsafe-state fates, {len(stuck)} stuck blocks, {time.time() - t0:.0f}s")
    write_summary(out, cases, edits, fates, stuck)


def write_summary(out, cases, edits, fates, stuck):
    ok = [c for c in cases if c["status"] != "ERROR"]; S = Counter(c["status"] for c in cases); wi = [c for c in ok if c["n_rplus_edits"]]
    st = [f for f in fates if f["stress"]]; fc = lambda key, rows: dict(Counter(f[key] for f in rows))
    infl_solved = [c for c in wi if c["status"] == "SOLVED"]
    sm = dict(cases=len(cases), status=dict(S), n_cases_with_rplus_edits=len(wi), n_p1_edits=len(edits), n_rplus_edits=sum(c["n_rplus_edits"] for c in ok),
              accepted_edit_verdicts=dict(Counter(e["verdict"] for e in edits)), accepted_edit_classes=dict(Counter(e["cls"] for e in edits)),
              n_betti_changing_accepted_edits=sum(c["n_betti_changing_edits"] for c in ok), ctrl_equals_stored_rplus=dict(Counter(str(c["ctrl_equals_stored_rplus"]) for c in ok)),
              total_init_crit=sum(c["init_crit"] for c in ok), total_final_crit=sum(c["final_crit"] for c in ok),
              terminal_stuck_cases=S["TERMINAL_STUCK"], terminal_stuck_blocks=len(stuck), stuck_by_k=dict(Counter(s["k"] for s in stuck)),
              stuck_initial_critical=sum(s["initial_critical"] for s in stuck), stuck_category=dict(Counter(s["category"] for s in stuck)),
              stuck_unlock_r1=dict(Counter(s["unlock_r1"] for s in stuck)), stuck_removal_safe_resolving={f"{k[0]}|rem_safe>0={k[1]}": v for k, v in Counter((s["category"], int(s["rem_safe_resolving"]) > 0) for s in stuck).items()},
              cases_with_raw_shell_closed=sum(1 for c in ok if int(c["n_raw_shell_closed_initial"]) > 0), raw_shell_closed_blocks=sum(int(c["n_raw_shell_closed_initial"]) for c in ok),
              edit_inflation_global_all=(sum(c["n_p1_edits"] for c in wi) / max(1, sum(c["n_rplus_edits"] for c in wi))),
              edit_inflation_global_solved=(sum(c["n_p1_edits"] for c in infl_solved) / max(1, sum(c["n_rplus_edits"] for c in infl_solved))),
              deferred=dict(events=sum(c["n_deferred_events"] for c in ok), blocks=sum(c["n_deferred_blocks"] for c in ok), reactivation_events=sum(c["n_reactivation_events"] for c in ok),
                            reactivated_blocks=sum(c["n_reactivated_blocks"] for c in ok), released_blocks=sum(c["n_released_blocks"] for c in ok),
                            resolved_by_other=sum(c["n_deferred_resolved_by_other"] for c in ok), cases_with_deferral=sum(c["n_deferred_blocks"] > 0 for c in ok)),
              new_crit=dict(p1_total=sum(c["n_new_crit"] for c in ok), p1_edits_creating=sum(c["n_edits_new_crit"] for c in ok), rplus_total=sum(c.get("rplus_new_crit", 0) for c in ok),
                            rplus_edits_creating=sum(c.get("rplus_edits_new_crit", 0) for c in ok)),
              p1_edit_choice=dict(equal_to_rplus_voxel=sum(e["p1_is_rplus"] for e in edits), differs=sum(1 - e["p1_is_rplus"] for e in edits)),
              unsafe_states=len(fates), stress_states=len(st), stress_fate_any=fc("fate_any", st), stress_fate_all=fc("fate_all", st),
              stress_block_level=dict(A=sum(f["n_A"] for f in st), B=sum(f["n_B"] for f in st), C=sum(f["n_C"] for f in st)), stress_n_cand_hist=dict(Counter(f["n_cand"] for f in st)),
              stress_voxel_added_in_p1=sum(f["voxel_added_in_p1"] for f in st), stress_met_same_local_state=sum(f["met_same_local_state"] for f in st),
              unsafe_fate_any=fc("fate_any", fates), nonstress_fate_any=fc("fate_any", [f for f in fates if not f["stress"]]),
              errors=[dict(case_id=c["case_id"], error=c["error"]) for c in cases if c["status"] == "ERROR"])
    json.dump(sm, open(Path(out) / "P1_SUMMARY.json", "w"), indent=1, ensure_ascii=False); print(json.dumps(sm, indent=1, ensure_ascii=False)[:6000])


def cmd_illustrative(a):
    out = Path(a.out); out.mkdir(parents=True, exist_ok=True); audit = audit_rows(a.p4); P = Path(a.p4); allrows = []
    for cid in ILLUSTRATIVE:
        raw = P0.load_zyx(P / "raw" / f"{cid}_WT_RAW.nii.gz"); edits0 = P0.read_edits(P, cid); tr = []
        states, _ = orig_states(cid, raw, edits0)
        r = run_p1(raw, gate_on=True, trace=tr, cid=cid, watch={k for s in states if s["verdict"] == "UNSAFE" for k, _, _ in s["cands"]})
        rows, fin = P0.replay(raw, r["edits"], TABLE)
        un = [s for s in states if s["verdict"] == "UNSAFE"]
        print(f"{cid}: R+ edits={len(edits0)} (UNSAFE {len(un)}) | P1 edits={len(r['edits'])} deferred_events={r['st']['defer_ev']} stuck={len(r['stuck'])} "
              f"init_crit={len(r['init'])} replay_all_SAFE={all(w['verdict']=='SAFE' for w in rows)} betti_delta={tuple(sum(w[k] for w in rows) for k in ('db0','db1','db2'))}")
        for s in un:
            fs = [fate(k, r["H"].get(k), r["stuck_keys"], s["fps"][k]) for k, _, _ in s["cands"]]
            print(f"   original UNSAFE edit {s['edit']} {vstr(s['p'])} {s['cls']} any_safe_alt={s['any_safe']} -> voxel added in P1: {s['p'] in set(r['edits'])}; candidate-block fates {[(f[0], f[1]) for f in fs]}")
        for t in tr: t["case_id"] = cid
        allrows += tr
    cols = ["case_id", "seq", "event", "block", "n_crit", "rplus_voxel", "rplus_gate", "zero_voxels", "safe_set", "p1_voxel", "V", "E", "F", "h0", "h1", "h2", "gate", "resolved_blocks",
            "new_crit_blocks", "reactivated", "deferred_size", "cum_edits", "stuck_now"]
    write_csv(out / "P1_ILLUSTRATIVE.csv", allrows, cols); print("wrote", len(allrows), "trace rows")


# ------------------------------------------------------------------------------------------------ rule-sensitivity diagnostic (NOT part of P1)
VARIANTS = [("lexmax", False), ("minnew", False), ("lexmin", True)] + [(f"rand{i}", False) for i in range(1, 6)]
def variant_job(args):
    cid, p4 = args; raw = P0.load_zyx(Path(p4) / "raw" / f"{cid}_WT_RAW.nii.gz"); out = dict(case_id=cid)
    for sel, lifo in VARIANTS:
        try:
            r = run_p1(raw, select=sel, lifo=lifo, cid=cid); out[f"{sel}{'+lifo' if lifo else ''}"] = len(r["stuck"])
        except Exception as e: out[f"{sel}{'+lifo' if lifo else ''}"] = f"ERROR {type(e).__name__}: {e}"
    return out


def cmd_variants(a):
    base = {r["case_id"]: r for r in csv.DictReader(open(a.base, encoding="utf-8-sig"))}; ids = sorted(c for c, r in base.items() if r["status"] == "TERMINAL_STUCK")
    t0 = time.time(); rows = []
    with Pool(a.workers, initializer=_init_worker, initargs=(a.vendor, a.p0, a.table, a.p0e, {})) as pool:
        for k, r in enumerate(pool.imap(variant_job, [(c, a.p4) for c in ids], chunksize=1)):
            r["baseline"] = int(base[r["case_id"]]["n_terminal_stuck"]); r["n_raw_shell_closed_initial"] = int(base[r["case_id"]]["n_raw_shell_closed_initial"]); rows.append(r)
            if (k + 1) % 50 == 0: print(f"  {k + 1}/{len(ids)} {time.time() - t0:.0f}s", flush=True)
    cols = ["case_id", "baseline", "n_raw_shell_closed_initial"] + [f"{s}{'+lifo' if l else ''}" for s, l in VARIANTS]
    write_csv(Path(a.out) / "P1_VARIANTS_DIAGNOSTIC.csv", rows, cols)
    num = lambda r: [r[c] for c in cols[3:] if isinstance(r[c], int)]
    best = [min([r["baseline"]] + num(r)) for r in rows]
    sm = dict(stuck_cases=len(rows), errors=sum(isinstance(r[c], str) for r in rows for c in cols[3:]), solved_by_some_variant=sum(b == 0 for b in best),
              total_stuck_blocks_baseline=sum(r["baseline"] for r in rows), total_stuck_blocks_best_per_case=sum(best),
              per_variant_solved={c: sum(r[c] == 0 for r in rows) for c in cols[3:]}, per_variant_total_stuck={c: sum(r[c] for r in rows if isinstance(r[c], int)) for c in cols[3:]},
              raw_shell_closed_cases_solved_by_variant=sum(1 for r, b in zip(rows, best) if r["n_raw_shell_closed_initial"] > 0 and b == 0),
              best_le_raw_shell_closed=sum(1 for r, b in zip(rows, best) if b >= r["n_raw_shell_closed_initial"]), n_cases_best_lt_baseline=sum(b < r["baseline"] for r, b in zip(rows, best)))
    json.dump(sm, open(Path(a.out) / "P1_VARIANTS_SUMMARY.json", "w"), indent=1); print(json.dumps(sm, indent=1))


# ------------------------------------------------------------------------------------------------ self-test
def cmd_selftest(a):
    cb = _load("cubical_betti", a.betti).cubical_betti; beta = lambda X: cb(X) if X.any() else type("Z", (), {"beta": (0, 0, 0)})(); rng = np.random.default_rng(20261007); fails = []; log = []
    def say(s): print(s, flush=True); log.append(s)
    def check(name, ok):
        say(f"  [{'PASS' if ok else 'FAIL'}] {name}")
        if not ok: fails.append(name)
    say("# P1 self-test\n")
    # 1. random small volumes: gate-OFF == frozen R+, P1 invariants, per-step from-scratch Betti, determinism, behaviours exercised
    from scipy import ndimage as ndi
    beh = Counter(); n = a.vols; stepchk = bad_step = bad_ctrl = bad_det = bad_solved = 0
    for i in range(n):
        shp = tuple(int(s) for s in rng.integers(5, 11, 3)); d = rng.uniform(0.2, 0.7)
        v = (rng.random(shp) < d).astype(np.uint8)
        if i % 2: v = (ndi.uniform_filter(v.astype(float), 3, mode="constant") > 0.45).astype(np.uint8)       # smoothed blobs: fewer, more structured defects
        r_prod, aud = R.repair_r_plus(v, return_audit=True); c = run_p1(v, gate_on=False); bad_ctrl += (c["edits"] != list(aud.edits)) or not np.array_equal(c["final"], r_prod)
        r1 = run_p1(v); r2 = run_p1(v); bad_det += (r1["edits"] != r2["edits"] or not np.array_equal(r1["final"], r2["final"]))
        X = np.pad(v, 1); prev = beta(X).beta
        for (z, y, x_) in r1["edits"]:
            X[z + 1, y + 1, x_ + 1] = 1; cur = beta(X).beta; stepchk += 1; bad_step += (cur != prev); prev = cur
        stuck = len(r1["stuck"]); beh["stuck"] += stuck > 0; beh["solved"] += stuck == 0; beh["deferred"] += r1["st"]["defer_ev"] > 0; beh["released"] += len(r1["released"]) > 0
        beh["reactivated"] += len(r1["ever_react"]) > 0; beh["def_resolved_by_other"] += len(r1["def_res_other"]) > 0; beh["p1_differs_from_rplus"] += any(not e["p1_is_rplus"] for e in r1["recs"])
        beh["rplus_changes_betti"] += beta(np.pad(r_prod, 1)).beta != beta(np.pad(v, 1)).beta
        beh["new_crit_by_p1"] += r1["st"]["new_crit"] > 0
        if stuck == 0 and beta(np.pad(r1["final"], 1)).beta != beta(np.pad(v, 1)).beta: bad_solved += 1
    say(f"1. {n} random volumes (5..10^3, density 0.2-0.7, half smoothed)")
    check(f"gate-OFF run == frozen production R+ (edits and mask) on all {n} volumes (mismatches={bad_ctrl})", bad_ctrl == 0)
    check(f"P1 is deterministic (two runs identical; mismatches={bad_det})", bad_det == 0)
    check(f"every accepted P1 edit leaves the from-scratch Betti vector unchanged ({stepchk} steps, violations={bad_step})", bad_step == 0)
    check(f"SOLVED volumes end with the raw Betti vector (violations={bad_solved})", bad_solved == 0)
    say(f"   behaviours exercised (volumes): {dict(beh)}")
    for k in ("stuck", "solved", "deferred", "released", "p1_differs_from_rplus", "rplus_changes_betti"):
        check(f"behaviour '{k}' exercised at least once", beh[k] > 0)
    # 2. hand-made (P0 selftest fixtures): a plugged tunnel and a filled cavity must be UNSAFE; a bar extension must be SAFE
    ring = np.zeros((9, 9, 9), np.uint8)
    for c in [(4, 3, 4), (4, 5, 4), (4, 4, 3), (4, 4, 5)]: ring[c] = 1
    nm, att, vv, cl = gate(np.pad(ring, 1), (4, 4, 4)); check(f"tunnel-plug voxel: {vv}/{cl} A={att}", vv == "UNSAFE" and cl == "loop" and att == (8, 12, 4, 1, 1, 0))
    shell = np.zeros((9, 9, 9), np.uint8); shell[3:6, 3:6, 3:6] = 1; shell[4, 4, 4] = 0
    nm, att, vv, cl = gate(np.pad(shell, 1), (4, 4, 4)); check(f"cavity-fill voxel: {vv}/{cl} A={att}", vv == "UNSAFE" and cl == "shell" and att == (8, 12, 6, 1, 0, 1))
    bar = np.zeros((9, 9, 9), np.uint8); bar[4, 4, 3:5] = 1
    nm, att, vv, cl = gate(np.pad(bar, 1), (4, 4, 5)); check(f"bar-extension voxel: {vv}/{cl}", vv == "SAFE")
    # 3. real cases: gate-OFF reproduces the stored trajectory, independent P0 replay accepts the P1 trajectory
    audit = audit_rows(a.p4); ids = sorted(audit); pick = list(ILLUSTRATIVE) + [str(c) for c in rng.choice(ids, a.n_real, replace=False)]
    for cid in pick:
        res = case_job((cid, a.p4, audit[cid])); c = res["case"]; res2 = case_job((cid, a.p4, audit[cid]))
        check(f"{cid}: second run identical (mask sha, edit rows, stuck rows)", res2["case"].get("p1_mask_sha256") == c.get("p1_mask_sha256") and res2["edits"] == res["edits"] and res2["stuck"] == res["stuck"])
        check(f"{cid}: status={c['status']} ctrl==stored R+ ({c['ctrl_equals_stored_rplus']}) P1 edits={c.get('n_p1_edits')} vs R+ {c.get('n_rplus_edits')} stuck={c.get('n_terminal_stuck')}"
              f" betti_delta={c.get('betti_delta_sum')} err='{c['error']}'", c["status"] != "ERROR" and c["ctrl_equals_stored_rplus"] == 1 and c["n_betti_changing_edits"] == 0)
    say(f"\nRESULT: {'ALL PASS' if not fails else 'FAILURES: ' + '; '.join(fails)}")
    Path(a.out).mkdir(parents=True, exist_ok=True); open(Path(a.out) / "P1_SELFTEST.md", "w", encoding="utf-8").write("\n".join(log) + "\n")
    sys.exit(1 if fails else 0)


if __name__ == "__main__":
    ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter); ap.add_argument("cmd", choices=["selftest", "illustrative", "all", "variants"])
    ap.add_argument("--out", default="."); ap.add_argument("--vendor", required=True); ap.add_argument("--p0", required=True); ap.add_argument("--table", required=True)
    ap.add_argument("--p4", required=True); ap.add_argument("--p0e", required=True); ap.add_argument("--base", default=""); ap.add_argument("--betti"); ap.add_argument("--workers", type=int, default=2)
    ap.add_argument("--no-control", action="store_true"); ap.add_argument("--cases", default=""); ap.add_argument("--n-real", type=int, default=6); ap.add_argument("--vols", type=int, default=400); a = ap.parse_args()
    setup(a.vendor, a.p0, a.table, a.p0e); _G.update(control=not a.no_control)
    {"selftest": cmd_selftest, "illustrative": cmd_illustrative, "all": cmd_all, "variants": cmd_variants}[a.cmd](a)
