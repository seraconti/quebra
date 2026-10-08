"""One permutation set, shared by every check that needs one.

Two reasons this is a shared object rather than a private detail of each check:

- **Cost.** The bench runs ~324k replicates and each check wants ~999 permutations. Drawing
  the index matrix once per replicate instead of once per check per replicate is the
  difference between a 40-minute run and a 3-hour one.
- **Paired comparisons.** The bench compares checks against each other cell by cell
  (criterion 4) and Arm B against Arm C at rho=0. When the checks see the SAME
  permutations, those differences are paired, and a paired difference has a smaller
  Monte Carlo error than the difference of two independent estimates. Sharing is not just
  cheaper, it makes the comparison sharper.

**Permutation is WITHIN segment, never across.** Each segment's `tau` is a function of
that segment's own gaps, so reordering gaps inside a segment leaves `tau` and the multiset
`{x_ij}` fixed and the conditional null is exact under iid gaps. Moving a gap between
segments changes both `tau_j` values and destroys that exactness. Hence the block
structure.
"""

from __future__ import annotations

from collections.abc import Sequence
from dataclasses import dataclass

import numpy as np

DEFAULT_N_PERM = 999


@dataclass(frozen=True)
class PermutationSet:
    """A block-structured index matrix over the concatenated gap vector.

    `indices` has shape `(n_perm, sum(sizes))`; row `b` is a permutation that maps each
    segment's slice onto itself. `sizes` records the block layout so a check cannot
    accidentally slice segments out in a different order than they were concatenated in.
    """

    sizes: tuple[int, ...]
    indices: np.ndarray

    @property
    def n_perm(self) -> int:
        return int(self.indices.shape[0])

    @property
    def total(self) -> int:
        return int(sum(self.sizes))

    def apply(self, values: np.ndarray) -> np.ndarray:
        """`(n_perm, total)` matrix of `values` reordered by each permutation."""
        values = np.asarray(values, dtype=float)
        if len(values) != self.total:
            raise ValueError(
                f"permutation set is over {self.total} values; got {len(values)}"
            )
        return values[self.indices]

    def blocks(self) -> list[tuple[int, int]]:
        """`[start, stop)` column ranges, one per segment, in concatenation order."""
        out: list[tuple[int, int]] = []
        off = 0
        for size in self.sizes:
            out.append((off, off + size))
            off += size
        return out


def block_permutations(
    sizes: list[int] | tuple[int, ...],
    n_perm: int = DEFAULT_N_PERM,
    rng: np.random.Generator | None = None,
) -> PermutationSet:
    """Draw `n_perm` within-segment permutations.

    `rng` is REQUIRED. Defaulting it to `np.random.default_rng()` would make every
    permutation p-value a fresh random variable: three consecutive calls on identical input
    returned p = 0.3860 / 0.4040 / 0.3790. That is fatal downstream, because
    `Job.build_identity` folds only the job-file hash and the dataset hashes - two runs
    producing opposite verdicts would share an identity, and the reuse rule would serve
    whichever ran first. Callers pass a generator built from an integer `seed` that travels
    as a step kwarg, so the seed reaches the provenance label.
    """
    if n_perm < 1:
        raise ValueError(f"n_perm must be >= 1; got {n_perm}")
    if any(int(size) < 0 for size in sizes):
        raise ValueError(f"segment sizes must be non-negative; got {list(sizes)}")
    if rng is None:
        raise ValueError(
            "block_permutations requires an explicit rng. None would draw OS "
            "entropy, which made permutation p-values irreproducible across runs while "
            "the run identity stayed unchanged. Build one from an integer seed: "
            "np.random.default_rng(seed)."
        )
    generator = rng
    columns: list[np.ndarray] = []
    offset = 0
    for size in sizes:
        size = int(size)
        block = np.tile(np.arange(size), (n_perm, 1))
        if size > 1:
            # permuted() shuffles each ROW independently, which is what a per-replicate
            # permutation means. shuffle(axis=1) would apply one common permutation.
            block = generator.permuted(block, axis=1)
        columns.append(block + offset)
        offset += size
    indices = np.hstack(columns) if columns else np.zeros((n_perm, 0), dtype=np.intp)
    return PermutationSet(sizes=tuple(int(s) for s in sizes), indices=indices)


def check_permuted(
    perm: PermutationSet | None, permuted: np.ndarray | None, caller: str
) -> None:
    """Guard a caller-supplied `permuted` matrix against its `PermutationSet`.

    `statistic_batch` takes an optional pre-gathered `(B, total)` matrix so the battery can
    build it once. Unguarded, that is the one route by which the shared-permutation
    optimisation could silently corrupt every number: a matrix from a DIFFERENT record was
    measured returning `p = 1.0` where the correct answer is 0.705, and one block too narrow
    truncates silently while the denominator uses the nominal width. Both now raise.

    What this does NOT catch, stated because the guard should not be trusted further than
    it goes: a matrix of the right shape built from a different record of the same segment
    layout. Detecting that needs a content hash on every call, which costs more than the
    gather it was meant to save. `battery.run_battery` is the only caller that passes
    `permuted`, and it builds it from the same `perm` three lines earlier.
    """
    if permuted is None:
        return
    if perm is None:
        raise ValueError(
            f"{caller} got a `permuted` matrix without the `perm` it was built from; "
            "there is then nothing to validate it against."
        )
    expected = (perm.n_perm, sum(perm.sizes))
    if permuted.shape != expected:
        raise ValueError(
            f"{caller} got a permuted matrix of shape {permuted.shape}, but its "
            f"permutation set implies {expected}. A mismatched matrix silently produces a "
            "finite, plausible, wrong statistic."
        )


def permutation_p_value(observed: float, null: np.ndarray) -> float:
    """Upper-tail permutation p-value, tie-corrected.

    `p = (1 + #{null >= observed}) / (1 + B)`. The `+1`s and the `>=` are load-bearing,
    not cosmetic: this project's duration vectors are heavily tied (at the 3 us threshold
    only 43 distinct values across 655 windows; at 5-7 us nearly every window is one read
    long), and a tied permuted statistic is evidence FOR the null, not against it. The
    naive `#{null > obs}/B` would count ties as refutations and report an
    anti-conservative p on exactly the cells where the data is thinnest.

    Including the observed value in its own reference set is also what makes the test
    exactly valid at finite B rather than valid only as B grows.
    """
    null = np.asarray(null, dtype=float)
    if null.size == 0:
        raise ValueError("permutation p-value needs at least one null draw")
    n_bad = int((~np.isfinite(null)).sum())
    if n_bad:
        raise ValueError(
            f"{n_bad} of {null.size} permuted statistics are not finite. A permutation "
            "of the observed gaps cannot make a well-posed statistic undefined, so this "
            "is a bug in the statistic, not a property of the data."
        )
    if not np.isfinite(observed):
        raise ValueError(f"observed statistic must be finite; got {observed}")
    return float((1 + int((null >= observed).sum())) / (1 + null.size))


def two_sided_p_value(observed: float, null: np.ndarray) -> float:
    """Tie-corrected permutation p-value for a signed statistic, via `|.|`.

    C1's Lewis-Robinson statistic is signed - negative for a decreasing trend, positive
    for an increasing one - and the check is two-sided, so the absolute value is the
    statistic whose upper tail matters. Taking `|.|` of BOTH the observed value and the
    null draws keeps the permutation argument intact, since `|.|` is applied identically
    to every member of the reference set.
    """
    return permutation_p_value(abs(observed), np.abs(np.asarray(null, dtype=float)))


def resolve_perm(
    sizes: Sequence[int],
    perm: PermutationSet | None,
    n_perm: int,
    rng: np.random.Generator | None,
) -> PermutationSet:
    """Build the permutation set if the caller did not supply one, or check the one it did.

    SPEC 0008 R8.5a. This exact block was copied into six call sites - C1, C2, CvM, C5, C6 and
    `battery` - in two spellings that differed only in whether `sizes` was already computed.
    Six copies of a guard is the shape `AGENTS.md` section 4 warns about: a fix applied to one
    comparison and not its twin is not a fix.

    It is deliberately NOT a precondition. The `if` branch CONSTRUCTS the missing argument and
    only the `elif` validates, so a decorator that runs before the body could express half of
    it at most. That asymmetry is why R8.5b evaluates a contract library against a different
    invariant instead of this one.

    Raising on a mismatch rather than rebuilding is the load-bearing half: a caller that passed
    a set blocked for different segment sizes is testing a different record, and silently
    rebuilding would turn that into a plausible, wrong p-value.
    """
    wanted = tuple(int(size) for size in sizes)
    if perm is None:
        # `block_permutations` is what refuses `rng=None`; not duplicated here, so there is
        # one message for that failure rather than two that can drift apart.
        return block_permutations(wanted, n_perm, rng)
    if perm.sizes != wanted:
        raise ValueError(
            f"permutation set is blocked as {perm.sizes} but the segments are {wanted}"
        )
    return perm
