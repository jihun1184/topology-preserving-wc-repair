# Minimal R4J reproducibility candidate

Status: the selected minimal recount and SAFE-consistent S4 materials are now in a public GitHub repository. The maintainer reports that the author approval required for this public scope is complete. No specific code/data license has been confirmed for this documentation update. This is not a complete search rerun package.

The internal validation status describes the preserved source package. Git line-ending conversion prevents a generic checkout from inheriting its byte-level PASS status. Read the repository-root README transport-integrity note before running the checks; v1.0.0 requires a verified byte-preserving delivery artifact.

Run from the extracted core_code directory with normal Python 3 (do not use python -O):

    python verify_package.py
    python verification/reaggregate_r1_2.py --root inputs --out verification_output
    python s4_records/verify_s4_records.py

All three scripts exist in this package. The two scientific verification commands in their default modes require only the Python standard library. verification_output is a generated directory, not a shipped input or missing script.

inputs/PRODUCTION/shard_0 and shard_1 each contain BLOCKS.csv and PATHS.csv. inputs/FINAL contains REUSE_MERGE_LOG.csv and P2AB_B1_MASTER_SHA256.txt. They are unchanged non-image records from the existing OR4. There is no results/ directory in this minimal package.

FILE_SHA256.json and SHA256SUMS.txt identify the shipped files. verify_package.py checks FILE_SHA256.json with standard Python, without requiring sha256sum on Windows. The two package hash lists exclude themselves to avoid self-referential hashes. This public documentation update changes only the README and rights-document entries; the other 18 entries retain their preserved digests. No new public release archive hash is asserted.

The re-aggregation reads the 19,040 recorded result rows and checks tables against the frozen stored outcomes. It does not independently reproduce the searches or re-evaluate the recorded SAFE flags. It explicitly reports the unavailable witness-NPZ checks as SKIPPED.

The S4 verifier computes boundary ranks and Betti vectors from the preserved cell ledger and checks each stored SAFE flag against the recomputed condition beta == [1,0,0]. This verifies the three fixed-voxel attachment results, not a new witness search. Source NPZ, original RAW masks and starting-state generation are absent. The verifier's optional --project-root mode needs the original source project and NumPy and is not a self-contained public-package command.

s4_records/close_evidence_FROZEN.py is provided solely as the byte-identical source identified by the S4 digest. The verifier hashes it; it does not execute it. Its historical build command and project-specific paths require the original project and are not instructions to run this package.

Search modules, SAFE exhaustive-check sources/scripts, split regeneration, Note1 finite-check scripts and figure-rendering sources are excluded from this minimal candidate. They remain preserved in the earlier candidate/OR4. No run_reaggregate.sh, regenerate_split.py, note1_checks/ or figure_sources/ execution is advertised here.

Read DEPENDENCIES_ICEIC27.md for unavailable dependencies and RIGHTS_AND_RELEASE.md for the current public scope and license candidates. GitHub publication has occurred. No particular software/data license, Zenodo deposit, DOI or v1.0.0 release is established by this update. The manuscript remains unchanged; any future availability wording must distinguish public recorded-result verification from unavailable full search, RAW-to-X_T generation and voxel replay.
