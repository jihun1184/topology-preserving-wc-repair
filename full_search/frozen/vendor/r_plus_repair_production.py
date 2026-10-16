from __future__ import annotations

from collections import deque
from dataclasses import dataclass
import hashlib
from typing import Iterable, Literal, NamedTuple

import numpy as np

from dwc_checker import is_critical_k2, is_critical_k3
from dwc_checker_fast import census_dwc_fast as census_dwc


Plane = Literal["xy", "yz", "zx"]
_PLANE_RANK = {"xy": 0, "yz": 1, "zx": 2}


class BlockKey(NamedTuple):
    k: int
    plane: str
    z: int
    y: int
    x: int


@dataclass(frozen=True)
class RepairAudit:
    edits: tuple[tuple[int, int, int], ...]
    initial_n2: int
    initial_n3: int
    queue_pops: int
    stale_pops: int
    propagated_new_blocks: tuple[BlockKey, ...]
    output_sha256: str

    @property
    def edit_count(self) -> int:
        return len(self.edits)


class RepairFailureError(RuntimeError):
    pass


def _fail(reason: str) -> "None":
    raise RepairFailureError(f"REPAIR_FAILURE_EXCEPTION:{reason}")


def _validate_binary_3d(mask: np.ndarray) -> np.ndarray:
    x = np.asarray(mask)
    if x.ndim != 3:
        raise ValueError(f"expected 3D binary array in (z,y,x) order, got ndim={x.ndim}")
    if not np.isin(x, (0, 1, False, True)).all():
        raise ValueError("input must be binary")
    return x.astype(np.uint8, copy=False)


def _sort_key(b: BlockKey) -> tuple[int, int, int, int, int]:
    # Frozen total order across dimensions/orientations:
    # (z, y, x), then k=2 before k=3, then xy < yz < zx.
    return (
        b.z,
        b.y,
        b.x,
        0 if b.k == 2 else 1,
        _PLANE_RANK.get(b.plane, 0),
    )


def _block_coords(b: BlockKey) -> tuple[tuple[int, int, int], ...]:
    z, y, x = b.z, b.y, b.x
    if b.k == 3:
        return tuple(
            (z + dz, y + dy, x + dx)
            for dz in (0, 1)
            for dy in (0, 1)
            for dx in (0, 1)
        )
    if b.k != 2:
        raise ValueError(f"unsupported block dimension k={b.k}")
    if b.plane == "xy":
        return ((z, y, x), (z, y, x + 1), (z, y + 1, x), (z, y + 1, x + 1))
    if b.plane == "yz":
        return ((z, y, x), (z, y + 1, x), (z + 1, y, x), (z + 1, y + 1, x))
    if b.plane == "zx":
        return ((z, y, x), (z, y, x + 1), (z + 1, y, x), (z + 1, y, x + 1))
    raise ValueError(f"unsupported k=2 plane {b.plane!r}")


def _block_array(x: np.ndarray, b: BlockKey) -> np.ndarray:
    z, y, x0 = b.z, b.y, b.x
    if b.k == 3:
        return x[z:z + 2, y:y + 2, x0:x0 + 2]
    if b.plane == "xy":
        return x[z, y:y + 2, x0:x0 + 2]
    if b.plane == "yz":
        return x[z:z + 2, y:y + 2, x0]
    if b.plane == "zx":
        return x[z:z + 2, y, x0:x0 + 2]
    raise ValueError(f"unsupported block {b}")


def _is_critical(x: np.ndarray, b: BlockKey) -> bool:
    a = _block_array(x, b)
    return is_critical_k2(a) if b.k == 2 else is_critical_k3(a)


def _all_block_keys(shape: tuple[int, int, int]) -> Iterable[BlockKey]:
    nz, ny, nx = shape
    for z in range(nz):
        for y in range(max(0, ny - 1)):
            for x in range(max(0, nx - 1)):
                yield BlockKey(2, "xy", z, y, x)
    for z in range(max(0, nz - 1)):
        for y in range(max(0, ny - 1)):
            for x in range(nx):
                yield BlockKey(2, "yz", z, y, x)
    for z in range(max(0, nz - 1)):
        for y in range(ny):
            for x in range(max(0, nx - 1)):
                yield BlockKey(2, "zx", z, y, x)
    for z in range(max(0, nz - 1)):
        for y in range(max(0, ny - 1)):
            for x in range(max(0, nx - 1)):
                yield BlockKey(3, "", z, y, x)


def _critical_keys_vectorized(x: np.ndarray) -> list[BlockKey]:
    """Detect all current in-domain critical blocks, then return frozen total order."""
    nz, ny, nx = x.shape
    out: list[BlockKey] = []

    if ny >= 2 and nx >= 2:
        a = x[:, :-1, :-1]
        b = x[:, :-1, 1:]
        c = x[:, 1:, :-1]
        d = x[:, 1:, 1:]
        mask = (a == d) & (b == c) & (a != b)
        out.extend(BlockKey(2, "xy", int(z), int(y), int(xx)) for z, y, xx in np.argwhere(mask))

    if nz >= 2 and ny >= 2:
        a = x[:-1, :-1, :]
        b = x[:-1, 1:, :]
        c = x[1:, :-1, :]
        d = x[1:, 1:, :]
        mask = (a == d) & (b == c) & (a != b)
        out.extend(BlockKey(2, "yz", int(z), int(y), int(xx)) for z, y, xx in np.argwhere(mask))

    if nz >= 2 and nx >= 2:
        a = x[:-1, :, :-1]
        b = x[:-1, :, 1:]
        c = x[1:, :, :-1]
        d = x[1:, :, 1:]
        mask = (a == d) & (b == c) & (a != b)
        out.extend(BlockKey(2, "zx", int(z), int(y), int(xx)) for z, y, xx in np.argwhere(mask))

    if nz >= 2 and ny >= 2 and nx >= 2:
        corners = [
            x[:-1, :-1, :-1], x[:-1, :-1, 1:],
            x[:-1, 1:, :-1],  x[:-1, 1:, 1:],
            x[1:, :-1, :-1],  x[1:, :-1, 1:],
            x[1:, 1:, :-1],   x[1:, 1:, 1:],
        ]
        s = sum(c.astype(np.uint16, copy=False) for c in corners)
        opposite_pairs = ((0, 7), (1, 6), (2, 5), (3, 4))
        crit = np.zeros_like(s, dtype=bool)
        for i, j in opposite_pairs:
            crit |= ((s == 2) & (corners[i] == 1) & (corners[j] == 1))
            crit |= ((s == 6) & (corners[i] == 0) & (corners[j] == 0))
        out.extend(BlockKey(3, "", int(z), int(y), int(xx)) for z, y, xx in np.argwhere(crit))

    out.sort(key=_sort_key)
    return out


def _incident_block_keys(shape: tuple[int, int, int], p: tuple[int, int, int]) -> list[BlockKey]:
    nz, ny, nx = shape
    pz, py, px = p
    keys: set[BlockKey] = set()

    for y0 in (py - 1, py):
        for x0 in (px - 1, px):
            if 0 <= pz < nz and 0 <= y0 < ny - 1 and 0 <= x0 < nx - 1:
                keys.add(BlockKey(2, "xy", pz, y0, x0))

    for z0 in (pz - 1, pz):
        for y0 in (py - 1, py):
            if 0 <= z0 < nz - 1 and 0 <= y0 < ny - 1 and 0 <= px < nx:
                keys.add(BlockKey(2, "yz", z0, y0, px))

    for z0 in (pz - 1, pz):
        for x0 in (px - 1, px):
            if 0 <= z0 < nz - 1 and 0 <= py < ny and 0 <= x0 < nx - 1:
                keys.add(BlockKey(2, "zx", z0, py, x0))

    for z0 in (pz - 1, pz):
        for y0 in (py - 1, py):
            for x0 in (px - 1, px):
                if 0 <= z0 < nz - 1 and 0 <= y0 < ny - 1 and 0 <= x0 < nx - 1:
                    keys.add(BlockKey(3, "", z0, y0, x0))

    return sorted(keys, key=_sort_key)


def _solve_binary_cc(x: np.ndarray, b: BlockKey) -> tuple[int, int, int]:
    """
    Binary specialization of Boutry et al. (2015) SOLVECC.

    For a binary critical block, Algorithm 1 raises one maximal point from the
    lower-valued side to the upper value. In {0,1}, every admissible correction
    is therefore 0 -> 1. The frozen tie rule selects the lexicographically
    smallest admissible zero voxel.
    """
    if not _is_critical(x, b):
        _fail("SOLVECC_CALLED_ON_NONCRITICAL_BLOCK")
    candidates = sorted(c for c in _block_coords(b) if int(x[c]) == 0)
    if not candidates:
        _fail("NO_BINARY_INCREASING_CANDIDATE")
    p = candidates[0]
    x[p] = 1
    if _is_critical(x, b):
        _fail("LOCAL_CORRECTION_DID_NOT_RESOLVE_BLOCK")
    return p


def array_sha256(x: np.ndarray) -> str:
    a = np.ascontiguousarray(x.astype(np.uint8, copy=False))
    h = hashlib.sha256()
    h.update(b"ICEIC2027_RPLUS_ARRAY_R0\0")
    h.update(np.asarray(a.shape, dtype=np.int64).tobytes(order="C"))
    h.update(a.tobytes(order="C"))
    return h.hexdigest()


def repair_r_plus(mask: np.ndarray, *, return_audit: bool = False) -> np.ndarray | tuple[np.ndarray, RepairAudit]:
    """Deterministic same-domain binary increasing DWC repair."""
    original = _validate_binary_3d(mask)
    x = original.copy()

    initial_keys = _critical_keys_vectorized(x)
    initial_set = set(initial_keys)
    initial_n2 = sum(b.k == 2 for b in initial_keys)
    initial_n3 = sum(b.k == 3 for b in initial_keys)

    q = deque(initial_keys)
    queued = set(initial_keys)
    edits: list[tuple[int, int, int]] = []
    propagated_new: set[BlockKey] = set()
    stale_pops = 0
    queue_pops = 0

    max_corrections = int(np.count_nonzero(x == 0))
    # Each correction can consider at most 12 k=2 + 8 k=3 incident blocks.
    max_queue_pops = len(initial_keys) + 20 * max_corrections

    while q:
        queue_pops += 1
        if queue_pops > max_queue_pops:
            _fail("QUEUE_POP_BOUND_EXCEEDED")

        b = q.popleft()
        queued.remove(b)
        if not _is_critical(x, b):
            stale_pops += 1
            continue

        p = _solve_binary_cc(x, b)
        edits.append(p)
        if len(edits) > max_corrections:
            _fail("MONOTONE_CORRECTION_BOUND_EXCEEDED")

        for nb in _incident_block_keys(x.shape, p):
            if _is_critical(x, nb) and nb not in queued:
                q.append(nb)
                queued.add(nb)
                if nb not in initial_set:
                    propagated_new.add(nb)

    if np.any((original == 1) & (x == 0)):
        _fail("STRICT_MONOTONICITY_VIOLATION")

    final = census_dwc(x)
    if final.n2 != 0 or final.n3 != 0:
        _fail(f"POSTCONDITION_DWC_FAILED:N2={final.n2}:N3={final.n3}")

    audit = RepairAudit(
        edits=tuple(edits),
        initial_n2=initial_n2,
        initial_n3=initial_n3,
        queue_pops=queue_pops,
        stale_pops=stale_pops,
        propagated_new_blocks=tuple(sorted(propagated_new, key=_sort_key)),
        output_sha256=array_sha256(x),
    )
    return (x, audit) if return_audit else x
