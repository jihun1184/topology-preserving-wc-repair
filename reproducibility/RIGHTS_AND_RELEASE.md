# Public scope, rights and license decisions

## Current status

This repository is Public. The maintainer has reported that the author approval needed to publish the current minimal recount and SAFE-consistent S4 scope is complete. The earlier internal-review statements that publication had not been authorized are superseded for this scope.

This report does not identify the approving parties, a particular institutional decision, an executed contributor agreement, copyright ownership allocations or a chosen code/data license. Those details are not inferred from repository visibility, manuscript authorship, CRediT, funding or file hashes.

## Scope of the public materials

| Unit | Current provision | License decision |
| --- | --- | --- |
| Minimal recount, package-integrity and S4 verification scripts; frozen code-identity source | Included in the 22-file reproducibility payload; no new scientific code changes | Confirm exact files, rights holders and code license |
| Non-image production/reuse CSV and S4 ledger, Betti, path and hash records | Included; recorded values unchanged | Confirm separate generated-data terms and attribution |
| Public documentation and metadata | Included; documentation updated for public status | Confirm whether code or separate documentation terms apply |
| Full ESM_4.zip, upstream/full-search stack, RAW masks, starting-state inputs and witness/context/replay NPZ | Not included | Any future provision requires a separate scope and rights decision |

The full search, RAW-to-X_T reconstruction and voxel replay are unavailable here. Running the ledger-based S4 verifier does not establish any of those capabilities. No access-on-request service is promised.

## License candidates — proposals only

- **Code: MIT** is the first candidate for concise permissive reuse with retention of copyright and license notices. [MIT text](https://choosealicense.com/licenses/mit/)
- **Code alternative: Apache-2.0** if the rights holders prefer an explicit contributor patent grant and its associated conditions. [Authoritative text](https://www.apache.org/licenses/LICENSE-2.0)
- **Generated non-image records: CC BY 4.0** is a candidate if the relevant rights can be licensed and source-data obligations permit it. Attribution must identify the approved creators/source and changes. It does not grant rights in absent BraTS images or NPZ arrays. [Official terms](https://creativecommons.org/licenses/by/4.0/)

No candidate is adopted by this document. No LICENSE file, SPDX declaration or license badge is added. Author approval of public access is recorded separately from approval of a reuse license.

## Release and availability

v1.0.0 is pending final review. Public access has occurred; a GitHub Release, Zenodo deposit, DOI and release date are not established by this update. The repository URL is https://github.com/jihun1184/topology-preserving-wc-repair.

Before release, confirm license choices and their file coverage, review the precise tagged commit and assets, and resolve the line-ending transport issue described in the root README. Do not represent an automatically generated GitHub source archive as the original verified payload without testing its exact bytes. See [the final review checklist](../PRIVATE_RELEASE_CHECKLIST.md).
