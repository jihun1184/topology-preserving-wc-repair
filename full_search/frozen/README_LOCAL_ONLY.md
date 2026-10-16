# Local reviewer code kit

Preserved existing sources and non-image provenance records. This is an author-review candidate, not a public release or a license grant. No image masks/NPZ, credentials, prior manuscript, or new scientific algorithm is included.

The complete reproduction executed in this audit is documented in ../rerun_core.py, ../RERUN_CORE_RESULTS.json, ../ASSET_SHA256_INVENTORY.json and ../../../FULL_SEARCH_REPRODUCIBILITY_CLOSURE.md. Sources in this directory are byte copies. They expect the same relative project structure and the arguments recorded in RERUN_CORE_RESULTS.json. `vendor` and `cubical_betti.py` must be supplied through the matching arguments, not renamed inside historical modules.

Python 3.12.14, NumPy 2.3.5, nibabel 5.4.2 were used. Install those libraries in a separate environment. Original historical search used another recorded environment; agreement was tested here. The raw masks and 3512-file materialization manifest payload must be obtained through author-approved access and the BraTS provider terms; copying this kit alone cannot reconstruct those inputs. Full P4 archive SHA and individual00110 file SHA appear in the inventory. For verified frozen inputs, run original p2ab_ablate.py through its `run` CLI or the audit adapter; `--nmax 2000000`, G26 R2 L6, G14 R1 L4.

The public repository remains a smaller recount/S4 package. This local kit has not been uploaded or submitted and does not close external reviewer access by itself.
