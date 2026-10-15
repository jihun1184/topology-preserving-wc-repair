# Public repository and v1.0.0 final review checklist

The filename is retained for link continuity from private staging. This is now a public-repository release checklist.

## Recorded status

- [x] The repository is Public.
- [x] The maintainer reports that author approval required for the current published scope is complete.
- [x] Public documentation distinguishes recorded-result recount and S4 ledger/SAFE verification from unavailable full search, RAW-to-X_T generation and voxel replay.
- This checklist does not certify an unprovided institutional decision, contributor agreement or specific license.
- No GitHub Release is created by this documentation update.

## Final changes to confirm

- [ ] Approve an explicit code license and its exact file coverage (MIT candidate; Apache-2.0 alternative), or explicitly decide and describe any remaining unlicensed scope.
- [ ] Separately approve generated-data and documentation terms (CC BY 4.0 candidate for eligible generated records), rights-holder names and attribution.
- [ ] Confirm the v1.0.0 scope: the current 22-file minimal payload plus public documentation; no full ESM_4.zip, upstream search stack, RAW masks or NPZ inputs without separate approval.
- [ ] Resolve byte-preserving delivery. At initial commit 6099c69623567a22c2e14a85683082e1e2361a35, raw Git blobs mismatch 17/20 preserved package digests. A Windows core.autocrlf=true checkout mismatches three entries. The default package and S4 checks fail there. Do not relabel this as PASS or silently replace scientific digests.
- [ ] Review the final documentation diff, exact target commit, release title/notes and asset inventory; do not create v1.0.0 before final review.

## Validation procedure for the exact release candidate

1. Compare the final repository tree with the initial public commit: code, calculation inputs/results and scientific hash references must be unchanged unless a separate transport correction is explicitly reviewed. For this documentation update, only the README/rights entries in the two package manifests change.
2. Prepare the precise delivery artifact using the approved preserved file bytes. Validate it after extraction in a clean directory. Independently test any GitHub-generated source ZIP/tarball; a working-tree test does not certify an archive.
3. With ordinary Python (not python -O), run from reproducibility/:

       python verify_package.py
       python verification/reaggregate_r1_2.py --root inputs --out verification_output
       python s4_records/verify_s4_records.py

   Require the 20 manifest entries, 19,040-row recount and all three S4 boundary/Betti/SAFE states to pass. NPZ replay remains unavailable/SKIPPED. These are existing verification checks, not a new search or experiment.
4. In a disposable copy only, invert the Y0 SAFE field in s4_records/ATTACHMENT_RECALCULATION.json and run the S4 verifier directly. Require a nonzero exit and a SAFE mismatch for Y0. Do not edit or ship the negative-control copy.
5. Record per-file SHA-256 and the exact asset SHA-256 outside the hashed artifact; retain earlier frozen packages and historical scientific digests. Check for unintended image/NPZ inputs, internal review documents and transient outputs.
6. After author final review, create the tag v1.0.0 at the reviewed commit, attach only approved verified assets, and publish the release. Then verify the actual download hashes and record the real release URL/date. A GitHub release does not itself establish a Zenodo DOI. [GitHub release procedure](https://docs.github.com/en/repositories/releasing-projects-on-github/managing-releases-in-a-repository)

No manuscript, proof, Figure or Table change is part of this checklist.
