# Recorded-result and S4 verification subpackage

This subdirectory verifies supplied non-image records; it does not execute a
search. Version1.0.0's separate `../full_search/` adapter and root README describe
the new core-search capability and externally supplied input requirements.

Run with normal Python (not `python -O`), from this directory:

```text
python verify_package.py
python verification/reaggregate_r1_2.py --root inputs --out verification_output
python s4_records/verify_s4_records.py
```

Default checks use the standard library. They recount 19,040 stored outcomes
and recompute boundary ranks/Betti/SAFE for the three supplied S4 cell ledgers.
They do not rerun those 19,040 searches or extract the ledgers from source NPZ.
The original absent-NPZ checks remain SKIPPED. The optional historical project
mode and `close_evidence_FROZEN.py` need the historical project and are not the
portable core-search command.

Original scientific code/record bytes are preserved. FILE_SHA256.json and
SHA256SUMS.txt exclude themselves and update only README/rights-document hashes;
the other18 digests are unchanged. The release-level manifest verifies all files.
The historical DEPENDENCIES_ICEIC27.md is retained as provenance for this original
subpackage, not as the version1.0.0 root capability statement.

Code/docs use MIT; author-generated non-image records use CC BY4.0 under exact
root LICENSE_SCOPE.json. Provider masks/NPZ are absent and not relicensed. See
root DATA_ACCESS.md, DATA_LICENSE.md and THIRD_PARTY_NOTICES.md.
