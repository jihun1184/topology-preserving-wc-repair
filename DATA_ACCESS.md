# External input access

Start at the original [BraTS 2021 provider page](https://www.med.upenn.edu/cbica/brats2021/),
Registration/Data Request and applicable usage/citation conditions. Obtain Task1
training segmentations, not Task2 DICOM, hidden challenge validation labels or
another year's relabeled data. The release does not register users or accept
provider terms on their behalf; access remains subject to the provider's process.

Select `BraTS2021_00110/BraTS2021_00110_seg.nii.gz`. Historical file SHA-256:
`9dde5c40e32df29acbc2cfd9fef367cc7dea71ad34054cb156e524167ee0def0`.
The required canonical **WT array hash**, unaffected by gzip recompression, is
`8571e2efb9581ab0de6b407db9009d627ccbe927f845c768197ccd262c4bb2cd`;
ZYX shape (155,240,240), foreground90740. WT is the union of labels1,2,4. All
these conditions are checked before search. Do not disable the digest check.

An authorized binary WT NIfTI is also accepted. Historical WT file SHA:
`9ae059f0171623b0848f2caa5d3ea727078e151627bf460d1825c99f693e36e9`.
Scientific identity depends on its canonical array digest, not the compressed
container. Input/output/cache paths remain private; output/cache must be outside
the source tree. One WT NIfTI is generated in the local cache for original replay.

No MRI, segmentation, WT or context/replay arrays are included. Non-image case
IDs, edit coordinates, outcome rows, cell ledgers and hashes are generated records.
The repository's license does not grant rights in external data or replace the
provider agreement. Cite Baid et al.(2021), Menze et al.(2015), Bakas et al.(2017)
as required by the provider and already cited in the manuscript.

Validation used the project's existing WT input, supplied through
a command-line path to isolated software. Independent provider registration or
new download of the multiclass segmentation was not performed. An authorized
researcher with the matching array can execute the four core searches without
the historical project or full 876-case materialization archive.

Input possession and a matching hash do not independently certify the holder's
data-use rights; lawful provider access remains each user's responsibility.
