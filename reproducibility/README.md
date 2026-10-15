# Minimal R4J reproducibility candidate

Status: INTERNAL_VALIDATED_FREEZE_V1. Internal result re-aggregation is approved; public release is HOLD. Not uploaded, not licensed and not a complete search rerun package.

Run from the extracted core_code directory with normal Python 3 (do not use python -O):

    python verify_package.py
    python verification/reaggregate_r1_2.py --root inputs --out verification_output
    python s4_records/verify_s4_records.py

All three scripts exist in this package. The two scientific verification commands in their default modes require only the Python standard library. verification_output is a generated directory, not a shipped input or missing script.

inputs/PRODUCTION/shard_0 and shard_1 each contain BLOCKS.csv and PATHS.csv. inputs/FINAL contains REUSE_MERGE_LOG.csv and P2AB_B1_MASTER_SHA256.txt. They are unchanged non-image records from the existing OR4. There is no results/ directory in this minimal package.

FILE_SHA256.json and SHA256SUMS.txt identify the shipped files. verify_package.py checks FILE_SHA256.json with standard Python, without requiring sha256sum on Windows. The two root hash lists exclude themselves to avoid self-referential hashes; the delivery record gives archive hashes.

The re-aggregation reads the 19,040 recorded result rows and checks tables against the frozen stored outcomes. It does not independently reproduce the searches or re-evaluate the recorded SAFE flags. It explicitly reports the unavailable witness-NPZ checks as SKIPPED.

The S4 verifier computes boundary ranks and Betti vectors from the preserved cell ledger and checks each stored SAFE flag against the recomputed condition beta == [1,0,0]. This verifies the three fixed-voxel attachment results, not a new witness search. Source NPZ, original RAW masks and starting-state generation are absent. The verifier's optional --project-root mode needs the original source project and NumPy and is not a self-contained public-package command.

s4_records/close_evidence_FROZEN.py is provided solely as the byte-identical source identified by the S4 digest. The verifier hashes it; it does not execute it. Its historical build command and project-specific paths require the original project and are not instructions to run this package.

Search modules, SAFE exhaustive-check sources/scripts, split regeneration, Note1 finite-check scripts and figure-rendering sources are excluded from this minimal candidate. They remain preserved in the earlier candidate/OR4. No run_reaggregate.sh, regenerate_split.py, note1_checks/ or figure_sources/ execution is advertised here.

Read DEPENDENCIES_ICEIC27.md for withheld dependencies and RIGHTS_AND_RELEASE.md for permission and licence decisions. Neither file grants publication rights. No GitHub/Zenodo upload or DOI reservation has occurred. Code Availability in the manuscript remains unchanged; its final wording awaits the actual release decision.
