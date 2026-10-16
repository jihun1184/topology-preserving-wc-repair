from __future__ import annotations

from dataclasses import dataclass
from typing import FrozenSet, Tuple

import numpy as np


K2Key = Tuple[str, int, int, int]   # (shared-edge axis, z, y, x)
K3Key = Tuple[int, int, int]        # canonical 2x2x2 block anchor (z, y, x)


@dataclass(frozen=True)
class DWCCensus:
    n2: int
    n3: int
    m2: int
    m3: int
    cc2: FrozenSet[K2Key]
    cc3: FrozenSet[K3Key]

    @property
    def defect_burden(self) -> int:
        return self.n2 + self.n3

    @property
    def normalized_defect_burden(self) -> float:
        denom = self.m2 + self.m3
        return 0.0 if denom == 0 else self.defect_burden / denom

    @property
    def is_dwc(self) -> bool:
        return self.defect_burden == 0


def _binary_array(a: np.ndarray, shape: tuple[int, ...]) -> np.ndarray:
    a = np.asarray(a)
    if a.shape != shape:
        raise ValueError(f"expected shape {shape}, got {a.shape}")
    if not np.isin(a, (0, 1, False, True)).all():
        raise ValueError("input must be binary")
    return a.astype(np.uint8, copy=False)


def is_critical_k2(block: np.ndarray) -> bool:
    """Boutry k=2 critical configuration on one canonical 2x2 block."""
    b = _binary_array(block, (2, 2))
    return bool(
        b[0, 0] == b[1, 1]
        and b[0, 1] == b[1, 0]
        and b[0, 0] != b[0, 1]
    )


def is_critical_k3(block: np.ndarray) -> bool:
    """Boutry k=3 critical configuration on one canonical 2x2x2 block."""
    b = _binary_array(block, (2, 2, 2))
    s = int(b.sum())
    if s not in (2, 6):
        return False
    minority = 1 if s == 2 else 0
    pts = np.argwhere(b == minority)
    return bool(len(pts) == 2 and np.all(np.abs(pts[0] - pts[1]) == 1))


def _mixed(block: np.ndarray) -> bool:
    return int(block.min()) != int(block.max())


def census_dwc(mask: np.ndarray) -> DWCCensus:
    """
    Count unique spatial critical blocks for a finite binary 3D mask.

    Contract:
      * array order is (z, y, x);
      * outside the finite volume is background (0);
      * one-voxel zero padding realizes that extension;
      * each geometric k=2 or k=3 block is keyed and counted once.
    """
    x = np.asarray(mask)
    if x.ndim != 3:
        raise ValueError(f"expected 3D mask, got ndim={x.ndim}")
    if not np.isin(x, (0, 1, False, True)).all():
        raise ValueError("mask must be binary")
    x = x.astype(np.uint8, copy=False)

    nz, ny, nx = x.shape
    p = np.pad(x, 1, mode="constant", constant_values=0)

    cc2: set[K2Key] = set()
    cc3: set[K3Key] = set()
    m2 = 0
    m3 = 0

    # Shared edge parallel to z: 2x2 pattern in (y, x), fixed z.
    for zp in range(1, nz + 1):
        for yp in range(0, ny + 1):
            for xp in range(0, nx + 1):
                b = p[zp, yp:yp + 2, xp:xp + 2]
                m2 += int(_mixed(b))
                if is_critical_k2(b):
                    cc2.add(("z", zp - 1, yp - 1, xp - 1))

    # Shared edge parallel to y: 2x2 pattern in (z, x), fixed y.
    for yp in range(1, ny + 1):
        for zp in range(0, nz + 1):
            for xp in range(0, nx + 1):
                b = p[zp:zp + 2, yp, xp:xp + 2]
                m2 += int(_mixed(b))
                if is_critical_k2(b):
                    cc2.add(("y", zp - 1, yp - 1, xp - 1))

    # Shared edge parallel to x: 2x2 pattern in (z, y), fixed x.
    for xp in range(1, nx + 1):
        for zp in range(0, nz + 1):
            for yp in range(0, ny + 1):
                b = p[zp:zp + 2, yp:yp + 2, xp]
                m2 += int(_mixed(b))
                if is_critical_k2(b):
                    cc2.add(("x", zp - 1, yp - 1, xp - 1))

    # k=3 blocks.
    for zp in range(0, nz + 1):
        for yp in range(0, ny + 1):
            for xp in range(0, nx + 1):
                b = p[zp:zp + 2, yp:yp + 2, xp:xp + 2]
                m3 += int(_mixed(b))
                if is_critical_k3(b):
                    cc3.add((zp - 1, yp - 1, xp - 1))

    return DWCCensus(
        n2=len(cc2),
        n3=len(cc3),
        m2=m2,
        m3=m3,
        cc2=frozenset(cc2),
        cc3=frozenset(cc3),
    )
