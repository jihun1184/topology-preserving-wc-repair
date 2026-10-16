from __future__ import annotations

from dataclasses import dataclass
import numpy as np


@dataclass(frozen=True)
class DWCFastCensus:
    n2: int
    n3: int
    m2: int
    m3: int

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


def _validate(mask: np.ndarray) -> np.ndarray:
    x = np.asarray(mask)
    if x.ndim != 3:
        raise ValueError(f"expected 3D mask, got ndim={x.ndim}")
    if not np.isin(x, (0, 1, False, True)).all():
        raise ValueError("mask must be binary")
    return x.astype(np.uint8, copy=False)


def _count_k2(a, b, c, d) -> tuple[int, int]:
    # Binary mixed-label exposure: not all 0 and not all 1.
    s = (a.astype(np.uint8, copy=False) + b + c + d)
    m = int(np.count_nonzero((s > 0) & (s < 4)))
    n = int(np.count_nonzero((a == d) & (b == c) & (a != b)))
    return n, m


def census_dwc_fast(mask: np.ndarray) -> DWCFastCensus:
    """Vectorized count-equivalent implementation of frozen P2-1 census_dwc.

    Same finite-volume convention: one-voxel zero padding; outside=background.
    Returns counts only (N2,N3,M2,M3), not the explicit critical-block key sets.
    """
    x = _validate(mask)
    nz, ny, nx = x.shape
    p = np.pad(x, 1, mode="constant", constant_values=0)

    n2 = 0
    m2 = 0

    # Shared edge parallel to z; shape (nz, ny+1, nx+1).
    a = p[1:nz+1, 0:ny+1, 0:nx+1]
    b = p[1:nz+1, 0:ny+1, 1:nx+2]
    c = p[1:nz+1, 1:ny+2, 0:nx+1]
    d = p[1:nz+1, 1:ny+2, 1:nx+2]
    n, m = _count_k2(a,b,c,d); n2 += n; m2 += m

    # Shared edge parallel to y; shape (nz+1, ny, nx+1).
    a = p[0:nz+1, 1:ny+1, 0:nx+1]
    b = p[0:nz+1, 1:ny+1, 1:nx+2]
    c = p[1:nz+2, 1:ny+1, 0:nx+1]
    d = p[1:nz+2, 1:ny+1, 1:nx+2]
    n, m = _count_k2(a,b,c,d); n2 += n; m2 += m

    # Shared edge parallel to x; shape (nz+1, ny+1, nx).
    a = p[0:nz+1, 0:ny+1, 1:nx+1]
    b = p[0:nz+1, 1:ny+2, 1:nx+1]
    c = p[1:nz+2, 0:ny+1, 1:nx+1]
    d = p[1:nz+2, 1:ny+2, 1:nx+1]
    n, m = _count_k2(a,b,c,d); n2 += n; m2 += m

    # k=3 blocks; shape (nz+1, ny+1, nx+1).
    corners = (
        p[0:nz+1,0:ny+1,0:nx+1], p[0:nz+1,0:ny+1,1:nx+2],
        p[0:nz+1,1:ny+2,0:nx+1], p[0:nz+1,1:ny+2,1:nx+2],
        p[1:nz+2,0:ny+1,0:nx+1], p[1:nz+2,0:ny+1,1:nx+2],
        p[1:nz+2,1:ny+2,0:nx+1], p[1:nz+2,1:ny+2,1:nx+2],
    )
    s = np.zeros(corners[0].shape, dtype=np.uint8)
    for q in corners:
        s += q
    m3 = int(np.count_nonzero((s > 0) & (s < 8)))
    crit = np.zeros(s.shape, dtype=bool)
    for i,j in ((0,7),(1,6),(2,5),(3,4)):
        crit |= ((s == 2) & (corners[i] == 1) & (corners[j] == 1))
        crit |= ((s == 6) & (corners[i] == 0) & (corners[j] == 0))
    n3 = int(np.count_nonzero(crit))

    return DWCFastCensus(n2=n2,n3=n3,m2=m2,m3=m3)
