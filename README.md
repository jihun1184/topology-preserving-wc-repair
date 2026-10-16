# Topology-preserving well-composed repair: version 1.0.0

Research software for bounded repair reachability under SAFE single-voxel moves.
This version adds preserved search sources and an executable single-case
reproduction adapter to the earlier recorded-results and S4 ledger package.

| Operation | Supported scope |
|---|---|
| Integrity checking and aggregation | The supplied 19,040 recorded result rows; not 19,040 new search runs |
| S4 calculations | Boundary ranks, Betti vectors and SAFE flags from three non-image cell ledgers |
| Actual core search | BraTS2021_00110: C00-F1/G26, C00-F2/G26, C01P-F2/G26 and C00-F1/G14, using your authorized input |
| Input reconstruction | Replay frozen ordered Run A edits on that RAW WT mask to reconstruct XT |
| Full-cohort stages/search | Original sources/records preserved; the adapter does not rerun 876-case stages or all 19,040 searches |

No BraTS images, segmentation masks, RAW masks, voxel-context NPZ or private
cache files are distributed. Expected records and actually executed searches
are distinguished.

## Setup

Use CPython 3.12. Validated: Python 3.12.14 on Windows 11; other platforms have
not been tested. Create a private environment:

```text
python -m venv /private/jmiv-venv
/private/jmiv-venv/bin/python -m pip install -r requirements-lock.txt
```

On Windows use `C:/private/jmiv-venv/Scripts/python.exe` instead of `bin/python`.
The lock specifies NumPy 2.3.5, nibabel 5.4.2 and validated transitive packages.
Dependencies are not bundled. Run normal Python, not `python -O`.

In **every** command below, replace `python` with that environment's full Python
path (`/private/jmiv-venv/bin/python` or
`C:/private/jmiv-venv/Scripts/python.exe`). Creating/installing a venv does not
activate it; using the explicit interpreter avoids running a system Python.

## Verify non-image package without BraTS data

From the extracted release root:

```text
python verify_release.py
cd reproducibility
python verify_package.py
python verification/reaggregate_r1_2.py --root inputs --out verification_output
python s4_records/verify_s4_records.py
```

The three existing checks do not run a search or derive cell ledgers from masks.
Original scientific code/data bytes are preserved. Revised rights documentation
is reflected in the package's two hash lists.

## Execute core search with your own authorized input

Obtain BraTS 2021 Task 1 training segmentation `BraTS2021_00110` through the
provider's access process. See [DATA_ACCESS.md](DATA_ACCESS.md). The adapter
uses labels1,2,4 for whole tumour and transposes NIfTI XYZ axes to ZYX. It rejects
a mask whose shape and canonical array digest differ from the frozen contract.
Equivalent recompressed NIfTI containers are accepted.

From the extracted release root, with output/cache directories **outside** it:

```text
python -B full_search/reproduce_core.py --segmentation /private/BraTS2021_00110_seg.nii.gz --out /private/jmiv-results --work-dir /private/jmiv-cache
```

Alternatively provide an authorized binary WT mask:

```text
python -B full_search/reproduce_core.py --raw-wt /private/BraTS2021_00110_WT_RAW.nii.gz --out /private/jmiv-results --work-dir /private/jmiv-cache
```

Use a new `--out` each time; completed results are not overwritten. The WT route
was tested with the frozen input in a newly installed isolated environment.
The multiclass option implements the specified label union; an independently
downloaded original multiclass segmentation was not obtained in this validation.
Neither option downloads, uploads or redistributes the provider data.

Expected output (states include root and newly discovered depth-L states):

```text
C00-F1 G26 PURE_FOUND 52137 1|41|881|13361|37853
C00-F2 G26 NO_PURE_DEPTH_LIMITED 24528 1|18|149|759|2678|6974|13949
C01P-F2 G26 NO_PURE_DEPTH_LIMITED 24528 1|18|149|759|2678|6974|13949
C00-F1 G14 PURE_FOUND 3293 1|13|112|775|2392
PASS: all four original-engine results and witnesses match the frozen expected records
```

The preserved engine uses Nmax=2,000,000, G26=(R2,L6), G14=(R1,L4). Both F1 paths
have length4. F2 is completed bounded absence through depth6, not unbounded
exhaustion; the cap is not reached. Deterministic columns and witness/path-rule
flags must match `full_search/expected_results.json`. Timings and cumulative gate
audit counts are excluded from exact equality. Failure exits nonzero.

## Scientific and provenance scope

`full_search/frozen/` contains byte-identical historical sources/records checked
against `FROZEN_SOURCE_SHA256.json`. BFS, geometry, SAFE, admissible rules,
witness replay and goals are unchanged. `reproduce_core.py` is a new input/release
adapter, not a new search algorithm. It regenerates the pinned SAFE table using
the original classifier and requires SHA-256
`59407fef2acb5b06059a13f19092047d172bbddd49bc85f811d822293b4c4830`.

The adapter uses a **selected-case provenance profile**. It validates frozen
files, the supplied RAW and reconstructed XT; it does not call the historical
cohort bootstrap requiring 3,512 P4 files. Do not describe it as verifying or
rerunning that whole archive. The original bootstrap is preserved in the sources.

XT is reconstructed with original `RunA.reconstruct`, a privately materialized WT
NIfTI and frozen ordered edits, followed by an exact array-digest check. Stage1/2
are not rerun to generate those records anew. The private mask cache is outside
the source tree and absent from release assets.

## Licensing and version identification

Code/software documentation: [MIT](LICENSE). Author-generated non-image records:
[CC BY 4.0](DATA_LICENSE.md); exact coverage in `LICENSE_SCOPE.json`. External
data/dependencies retain their terms. See [THIRD_PARTY_NOTICES.md](THIRD_PARTY_NOTICES.md).

Use the versioned asset and SHA256SUMS/RELEASE_PROVENANCE.json to identify exact
bytes, tag and commit. Before final author approval this tree/local tag is only
a candidate; a prospective Release URL does not establish publication.
