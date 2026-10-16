# Provenance and third-party dependencies

The `full_search/frozen/vendor` name comes from the historical project layout.
Its files implement the authors' R+ and critical-block procedures; calling this
folder vendor does not establish a third-party license or assign new ownership.
All preserved code bytes and records are inventoried. Publication of this
candidate requires the authors' approval of the stated rights and license scope.

The mathematical background is acknowledged in the manuscript and preserved
source/specification comments. Implementing or citing a published mathematical
method is distinguished from reproducing another party's source code.

Python, NumPy, nibabel, packaging and typing-extensions are external dependencies,
not copied into the release. Their original licenses apply. BraTS provider data
are external inputs and are not licensed by the repository's MIT/CC BY notices.

Historical README/rights notes retained inside frozen source bundles describe
their earlier state. The current root README, LICENSE, DATA_LICENSE.md and
LICENSE_SCOPE.json define this version's scope. Prior Git commits are preserved;
this release does not rewrite their publication or licensing history.

Frozen metadata can retain original machine paths (for example /sessions/...)
as historical provenance. They are not installation instructions or credential
material. The portable adapter takes user-provided input/output paths and does
not require those machines or directories. Frozen bytes are preserved for hashes.
