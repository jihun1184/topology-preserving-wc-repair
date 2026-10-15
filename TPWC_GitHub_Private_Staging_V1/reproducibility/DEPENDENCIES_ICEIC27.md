# Dependencies and reproduction boundary of this minimal package

This document describes core_code as shipped, not the larger OR4.

Available: verification/reaggregate_r1_2.py, the five non-image production/reuse CSV files under inputs/, their original master hash list, the S4 attachment cell ledger and unchanged S4 verifier/source identifiers.

Unavailable: results/; verification/run_reaggregate.sh; verification/regenerate_split.py; split_contract/; safe_gate/; search_modules/; note1_checks/; figure_sources/; witness or replay NPZ. These may exist in the earlier, larger candidate but are not included here. There is no claim of split regeneration, SAFE-table generation, figure reproduction or a complete search rerun from this archive.

## Withheld upstream inputs and modules (lineage records, not release permissions)

|Name|Recorded role|Distribution here|
|---|---|---|
|r_plus_repair_production.py|ICEIC-related R+ production repair / start-state lineage|Not included|
|cubical_betti.py|Attachment homology audit in the search stack|Not included|
|p0_audit.py and P1_simple_constrained_rplus.py|Earlier gate loader/stage dependencies|Not included|
|RAW masks and recorded R+ edits|Reconstruction of X_T from BraTS|Not included|
|Full search specs/configuration and SAFE lookup table|Provenance-controlled full search stack|Not included|
|Original witness/context NPZ|Independent stored-state replay/hash checks|Not included|

The earlier OR4 dependency record identifies r_plus_repair_production.py by SHA-256 878b483f300415f6a3bde2169ce6630de890842c8a79326d99a0f50a2ecb255b and cubical_betti.py by d87571d09f828288a9660c65b6b841d4f7939f70d9c9b1b082bfc99bde5a5b64. These establish recorded lineage, not authorship, ownership or permission to distribute. Their current source files were not newly audited here.

Recounting stored results is executable from this package. Re-running the full search is not. Public access, reviewer-only access or embargoed access to withheld code is undecided; no such access service is asserted.
