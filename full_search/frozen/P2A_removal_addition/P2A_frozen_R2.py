#!/usr/bin/env python3
"""P2A R2 frozen-input manifest + verifier (spec R2 §3, T0/T9a).  Imported by the runner, the validator and the test harness; read-only on every input.
verify() raises ProvenanceError on ANY mismatch / missing file; callers must abort BEFORE producing any scientific output (T9a)."""
import hashlib, json
from pathlib import Path

class ProvenanceError(RuntimeError): pass

# (label, root, relative path, sha256).  roots: p1 = P1 folder, vendor = production R+ / dwc vendor folder, p0 = p0_audit.py, table = gate table, betti = cubical_betti.py, p4 = extracted P4 materialization.
FROZEN = [
 ("P1_simple_constrained_rplus.py", "p1", "P1_simple_constrained_rplus.py", "a27d1b6b00e83437a547db2f1f50455b62581cbd3f857a66820988f4623686b6"),
 ("P1_SUMMARY.json", "p1", "P1_SUMMARY.json", "4bd772be032cedef5171051b78c762b97a030b6b1ba7245b71e9240bbf9691c2"),
 ("P1_TERMINAL_STUCK.csv", "p1", "P1_TERMINAL_STUCK.csv", "348df726f2c783ba0e041c53f26cb10d3757db78f45c539e544755cf38e86d81"),
 ("P1_CASES_876.csv", "p1", "P1_CASES_876.csv", "9fe9e58e7f355b1d31c3fa334b953b66bdb5dff4f1dbeed225aec5409483f730"),
 ("P1_EDITS_876.csv", "p1", "P1_EDITS_876.csv", "ea7a512e1714aa7d266e30e4637d425fc97b8232a239c1eb8658d254e044354d"),
 ("vendor/r_plus_repair_production.py", "vendor", "r_plus_repair_production.py", "878b483f300415f6a3bde2169ce6630de890842c8a79326d99a0f50a2ecb255b"),
 ("vendor/dwc_checker.py", "vendor", "dwc_checker.py", "d0ed22c221fc85f1713a47c4a3f98934755e14584616b8fdd3494b205f69a51f"),
 ("vendor/dwc_checker_fast.py", "vendor", "dwc_checker_fast.py", "292af628823183d03e5348aca9b587eead0ff5a89f4a462052309ee0558cab61"),
 ("p0_audit.py", "p0", "", "8d549b628b6e5e9af688d8192e88644a15356e3c42cdeba761d7f6e011c8eb58"),
 ("gate_table_2p26.bin", "table", "", "59407fef2acb5b06059a13f19092047d172bbddd49bc85f811d822293b4c4830"),     # regenerated-and-pinned P2A artifact (not asserted to be P0's original full SHA)
 ("cubical_betti.py", "betti", "", "d87571d09f828288a9660c65b6b841d4f7939f70d9c9b1b082bfc99bde5a5b64"),
 ("P4 MANIFEST_SHA256.txt", "p4", "MANIFEST_SHA256.txt", "86d1d96441f49d5a8df6ad708973eb167be2e883868524d9172696e4b355a99e"),            # content of the pinned zip 784fc692…; its 3,512 entries are verified below
 ("P4 case_materialization_audit.csv", "p4", "case_materialization_audit.csv", "ed33ffc1ea4dec0a89ec3760a7aa7584eb61d18fbe323818f75ce71227109621"),
]
P4_MANIFEST_ENTRIES = 3512


def sha(p):
    h = hashlib.sha256()
    with open(p, "rb") as f:
        for b in iter(lambda: f.read(1 << 22), b""): h.update(b)
    return h.hexdigest()


def manifest_digest():
    """Digest of the pinned table above (recorded in run metadata so two runs can be shown to share the same frozen contract)."""
    return hashlib.sha256(json.dumps(FROZEN, sort_keys=True).encode()).hexdigest()


def verify(roots, case_ids=None, check_p4_files=True):
    """roots: dict p1/vendor/p0/table/betti/p4 -> path.  Returns {label: sha}.  Raises ProvenanceError (message lists every failure)."""
    bad, got = [], {}
    for lab, root, rel, exp in FROZEN:
        p = Path(roots[root]) / rel if rel else Path(roots[root])
        try: g = sha(p)
        except Exception as e: bad.append(f"{lab}: unreadable ({type(e).__name__}: {e})"); continue
        got[lab] = g
        if g != exp: bad.append(f"{lab}: sha {g} != pinned {exp}")
    if check_p4_files and not bad:                                  # every file listed in the P4 manifest, byte level (3,512 entries)
        p4 = Path(roots["p4"]); lines = [l.rstrip("\n") for l in open(p4 / "MANIFEST_SHA256.txt", encoding="utf-8") if l.strip()]; listed = set()
        if len(lines) != P4_MANIFEST_ENTRIES: bad.append(f"P4 manifest has {len(lines)} entries, expected {P4_MANIFEST_ENTRIES}")
        for l in lines:
            h, rel = l.split("  ", 1); rel = rel[2:] if rel.startswith("./") else rel; listed.add(rel)
            try:
                if sha(p4 / rel) != h: bad.append(f"P4 file {rel}: sha mismatch")
            except Exception as e: bad.append(f"P4 file {rel}: unreadable ({type(e).__name__})")
        for c in (case_ids or []):                                  # every case input must be covered by the manifest
            for rel in (f"raw/{c}_WT_RAW.nii.gz", f"edits/{c}_WT_EDITS.csv"):
                if rel not in listed: bad.append(f"case file not covered by P4 manifest: {rel}")
    if bad: raise ProvenanceError("; ".join(bad[:12]) + (f" (+{len(bad) - 12} more)" if len(bad) > 12 else ""))
    return got
