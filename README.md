# Topology-Preserving Well-Composed Repair — Reproducibility Materials

> **Public research repository.**
>
> This repository is a curated copy of the internally verified *minimal recount + SAFE-consistent S4* candidate supplied for the JMIV manuscript. The repository is publicly accessible, and the maintainer reports that the author approval needed for this public scope is complete. It is not the full search implementation. No specific code or data license has been confirmed in the materials reviewed for this documentation update.

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

The commands below describe the minimal package checks. Before interpreting their results, read the transport-integrity note below.

Install Python 3.9 or later. No external package is required for the following default checks. Run with **normal Python**, not `python -O`:

```bash
cd reproducibility
python verify_package.py
python verification/reaggregate_r1_2.py --root inputs --out verification_output
python s4_records/verify_s4_records.py
```

These commands verify supplied records; they **do not regenerate the search or infer SAFE flags from original voxel data**. For details and dependencies, read [`reproducibility/README.md`](reproducibility/README.md) and [`reproducibility/DEPENDENCIES_ICEIC27.md`](reproducibility/DEPENDENCIES_ICEIC27.md).

## Provenance and integrity

The initial 22-file `reproducibility/` payload was selected from `R4K_INTERNAL_REVIEW_CORRECTED_V1.zip` → `r4k_review/ESM_4.zip` → `or4_candidate/minimal_recount/core_code/`.

Scientific code and recorded calculation inputs/results are unchanged by this documentation update. The two package hash lists are updated only for the revised `README.md` and `RIGHTS_AND_RELEASE.md`; their other 18 recorded digests are preserved. Earlier frozen packages remain unchanged.

### Transport-integrity limitation before v1.0.0

At the initial public commit `6099c69623567a22c2e14a85683082e1e2361a35`, 17 of the 20 Git-blob digests differ from the preserved package manifest because the original payload contains mixed LF/CRLF line endings. A Windows checkout with `core.autocrlf=true` instead mismatches three entries, including `close_evidence_FROZEN.py`; the package and S4 checks fail in that checkout. A default clone or GitHub-generated source archive must therefore not be described as the byte-identical verified payload.

This documentation update does not rewrite code/data or replace their scientific digests. Before publishing v1.0.0, resolve byte-preserving delivery of the approved payload and verify the exact download artifact. See the [release checklist](PRIVATE_RELEASE_CHECKLIST.md).

## Access and licensing

The repository is **Public**. Author approval for the currently published scope has been reported by the maintainer. This does not identify a particular license, copyright holder, contributor agreement or institutional determination. No code/data license is applied by this update.

Separate code and generated-data license choices remain to be confirmed. Additional material, including the full ESM_4.zip, is outside the current scope unless separately approved. A v1.0.0 release will be considered only after final review; no DOI or release date is asserted.

See [`reproducibility/RIGHTS_AND_RELEASE.md`](reproducibility/RIGHTS_AND_RELEASE.md).
