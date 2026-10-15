# Topology-Preserving Well-Composed Repair — Reproducibility Materials

> **Private staging repository. Not approved for public release.**
>
> This repository is a curated copy of the internally verified *minimal recount + SAFE-consistent S4* candidate supplied for the JMIV manuscript. It is not the full search implementation, and inclusion in this staging copy does not grant a software or data license.

## Scope

| Reproduced calculation | Supported by this repository? |
| --- | --- |
| Verify supplied files against recorded SHA-256 digests | Yes |
| Recount tabulated outcomes from 19,040 **recorded** rows | Yes |
| Recompute boundary ranks, Betti vectors, and SAFE flags for the three supplied S4 ledgers | Yes |
| Rerun the full bounded search or reconstruct RAW-to-`X_T` starting states | **No** |
| Rerun voxel/NPZ witness replay | **No** — source NPZ arrays are not included |

No BraTS MRI, original segmentation masks, or witness/replay NPZ files are included.

## Quick start

Install Python 3.9 or later. No external package is required for the following default checks. Run with **normal Python**, not `python -O`:

```bash
cd reproducibility
python verify_package.py
python verification/reaggregate_r1_2.py --root inputs --out verification_output
python s4_records/verify_s4_records.py
```

These commands verify supplied records; they **do not regenerate the search or infer SAFE flags from original voxel data**. For details and dependencies, read [`reproducibility/README.md`](reproducibility/README.md) and [`reproducibility/DEPENDENCIES_ICEIC27.md`](reproducibility/DEPENDENCIES_ICEIC27.md).

## Provenance and integrity

The 22 files in `reproducibility/` are copied without scientific changes from `R4K_INTERNAL_REVIEW_CORRECTED_V1.zip` → `r4k_review/ESM_4.zip` → `or4_candidate/minimal_recount/core_code/`.

Original hash manifests have not been regenerated or edited. `verify_package.py` checks the original `FILE_SHA256.json`.

## Access and licensing

This is a **private review-stage package**. Copyright ownership, author/institutional permissions, and separate code/data licenses remain undecided. Do not switch the repository to Public, assign a license, or attach the full ESM_4.zip until those decisions are authorized.

See [`reproducibility/RIGHTS_AND_RELEASE.md`](reproducibility/RIGHTS_AND_RELEASE.md).
