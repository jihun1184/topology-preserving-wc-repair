# Table S4 local attachment records

Run: python s4_records/verify_s4_records.py from core_code, or python verify_s4_records.py in this directory.

Included: ATTACHMENT_RECALCULATION.json/.csv; ATTACHMENT_CELL_LEDGER.json; ORIGINAL_WITNESS_PATH_ROWS.json; PRESERVED_28_RESULT_ROWS.json; verify_s4_records.py; close_evidence_FROZEN.py; RECORD_SHA256.json.

The verifier checks the source-code digest, boundary composition zero, ranks, all three frozen Betti vectors and stored SAFE flags. SAFE is recomputed as beta == [1,0,0]; any inconsistent or non-boolean stored flag raises a fatal SAFE mismatch error. It hashes close_evidence_FROZEN.py as an identified source file and never executes that historical manuscript-building script.

The frozen probe source and FROZEN_PROBE_REPLAY.log are not included in this minimal candidate. Their earlier identifiers/results remain in the original R4J/OR4 records. The two original NPZ files and replay NPZ are not provided; --project-root requires external original inputs and NumPy.

State digests refer to C-order uint8 bytes in stored coordinates order. Ledger-only validation does not freshly establish those digests against NPZ arrays. It verifies the local cellular boundary calculation.

This local candidate does not replace OR3, alter S4 or grant a public-release licence.
