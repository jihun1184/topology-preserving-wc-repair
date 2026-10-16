"""Betti numbers of the T1 closed-voxel cubical complex.

Run the independent fixture and boundary-matrix checks with
    python -m unittest -v test_cubical_betti

Only NumPy is required. This implementation crops to the occupied bounding
box, counts the complex's cells by vectorized Boolean unions, and counts
foreground/background components by row-run union-find. It never unpickles a
checkpoint or reads a dataset. See T2_BETTI_IMPLEMENTATION_20261001.md for
the mathematical derivation and production limits.
"""

from __future__ import annotations

from dataclasses import dataclass
from itertools import combinations, product

import numpy as np


@dataclass(frozen=True)
class CubicalBetti:
    beta0: int
    beta1: int
    beta2: int
    chi: int
    # (vertices, edges, faces, cubes) of the closed cubical complex.
    cell_counts: tuple[int, int, int, int]
    # Occupied bounding-box size, before the one-voxel background collar.
    bbox_shape: tuple[int, int, int]

    @property
    def beta(self) -> tuple[int, int, int]:
        return (self.beta0, self.beta1, self.beta2)


def _binary_volume(mask: np.ndarray) -> np.ndarray:
    arr = np.asarray(mask)
    if arr.ndim != 3 or any(size == 0 for size in arr.shape):
        raise ValueError("mask must be a nonempty 3-D array")
    if arr.dtype.kind not in "biuf" or not np.all(np.isfinite(arr)):
        raise ValueError("mask must contain finite binary numeric values")
    if not np.all((arr == 0) | (arr == 1)):
        raise ValueError("mask must contain only 0 or 1")
    return arr.astype(np.bool_, copy=False)


def _crop_occupied(mask: np.ndarray) -> np.ndarray:
    """Remove only empty outer slabs; translation preserves homology."""
    # Axis profiles avoid storing one coordinate triplet per foreground voxel.
    # np.nonzero would allocate hundreds of MB for a nearly full BraTS volume.
    bounds = []
    for axis in range(3):
        profile = np.any(mask, axis=tuple(other for other in range(3) if other != axis))
        occupied = np.flatnonzero(profile)
        bounds.append(slice(int(occupied[0]), int(occupied[-1]) + 1))
    return mask[tuple(bounds)]


def _cell_count(mask: np.ndarray, moving_axes: tuple[int, ...]) -> int:
    """Count unique cells of one orientation, with all cube faces included.

    A k-cell varies in k coordinate directions. For each occupied 3-cube,
    its 2^(3-k) choices of fixed-coordinate endpoints are ORed into a
    Boolean lattice. OR removes duplicate shared cells exactly once.
    """
    fixed_axes = tuple(axis for axis in range(3) if axis not in moving_axes)
    out_shape = tuple(n + (axis in fixed_axes) for axis, n in enumerate(mask.shape))
    cells = np.zeros(out_shape, dtype=np.bool_)
    for endpoints in product((0, 1), repeat=len(fixed_axes)):
        offset = dict(zip(fixed_axes, endpoints))
        target = tuple(slice(offset.get(axis, 0), offset.get(axis, 0) + n)
                       for axis, n in enumerate(mask.shape))
        cells[target] |= mask
    return int(np.count_nonzero(cells))


def _cell_counts(mask: np.ndarray) -> tuple[int, int, int, int]:
    counts = []
    for dim in range(4):
        counts.append(sum(_cell_count(mask, axes)
                          for axes in combinations(range(3), dim)))
    return tuple(counts)  # type: ignore[return-value]


def _row_runs(row: np.ndarray) -> list[tuple[int, int]]:
    """Half-open x intervals of ones in a fixed (z,y) row."""
    if not np.any(row):
        return []
    if np.all(row):
        return [(0, row.size)]
    padded = np.empty(row.size + 2, dtype=np.bool_)
    padded[0] = padded[-1] = False
    padded[1:-1] = row
    changes = np.flatnonzero(padded[1:] != padded[:-1])
    return [(int(a), int(b)) for a, b in zip(changes[0::2], changes[1::2])]


class _UnionFind:
    def __init__(self) -> None:
        self.parent: list[int] = []
        self.size: list[int] = []
        self.components = 0

    def add(self) -> int:
        index = len(self.parent)
        self.parent.append(index)
        self.size.append(1)
        self.components += 1
        return index

    def root(self, index: int) -> int:
        parent = self.parent
        while index != parent[index]:
            parent[index] = parent[parent[index]]
            index = parent[index]
        return index

    def union(self, a: int, b: int) -> None:
        a, b = self.root(a), self.root(b)
        if a == b:
            return
        if self.size[a] < self.size[b]:
            a, b = b, a
        self.parent[b] = a
        self.size[a] += self.size[b]
        self.components -= 1


def _join_rows(
    current: list[tuple[int, int, int]],
    earlier: list[tuple[int, int, int]],
    x_reach: int,
    groups: _UnionFind,
) -> None:
    """Union runs when their x intervals are within x_reach voxels."""
    j = 0
    for begin, stop, label in current:
        while j < len(earlier) and earlier[j][1] + x_reach <= begin:
            j += 1
        k = j
        while k < len(earlier) and earlier[k][0] < stop + x_reach:
            groups.union(label, earlier[k][2])
            k += 1


def _component_count(mask: np.ndarray, adjacency: int) -> int:
    if adjacency not in (6, 26):
        raise ValueError("adjacency must be 6 or 26")
    depth, height, _ = mask.shape
    groups = _UnionFind()
    previous_plane: list[list[tuple[int, int, int]]] = [[] for _ in range(height)]
    x_reach = 0 if adjacency == 6 else 1

    for z in range(depth):
        plane: list[list[tuple[int, int, int]]] = []
        for y in range(height):
            runs = [(a, b, groups.add()) for a, b in _row_runs(mask[z, y])]
            if y:
                _join_rows(runs, plane[y - 1], x_reach, groups)
            if z:
                neighbors = (y,) if adjacency == 6 else range(max(0, y - 1), min(height, y + 2))
                for old_y in neighbors:
                    _join_rows(runs, previous_plane[old_y], x_reach, groups)
            plane.append(runs)
        previous_plane = plane
    return groups.components


def cubical_betti(mask: np.ndarray) -> CubicalBetti:
    """Return beta_0,beta_1,beta_2 for the union of closed foreground cubes.

    beta_0 is the number of 26-connected foreground components. By
    Alexander duality, beta_2 counts bounded complementary components;
    these correspond to 6-connected background components after adding a
    one-voxel exterior collar, minus the exterior component. Finally the
    Euler-Poincare identity gives beta_1 = beta_0 + beta_2 - chi.
    """
    binary = _binary_volume(mask)
    if not np.any(binary):
        return CubicalBetti(0, 0, 0, 0, (0, 0, 0, 0), (0, 0, 0))
    cropped = _crop_occupied(binary)
    cells = _cell_counts(cropped)
    chi = cells[0] - cells[1] + cells[2] - cells[3]
    beta0 = _component_count(cropped, 26)
    # The added empty collar includes all exterior paths. It is needed even
    # when the caller's original array happens to be all foreground.
    background = ~np.pad(cropped, 1, mode="constant", constant_values=False)
    beta2 = _component_count(background, 6) - 1
    beta1 = beta0 + beta2 - chi
    if beta2 < 0 or beta1 < 0:
        raise RuntimeError("inconsistent cubical invariants; inspect connectivity and cell counts")
    return CubicalBetti(beta0, beta1, beta2, chi, cells, tuple(cropped.shape))
