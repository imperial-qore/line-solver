"""
Distinct permutations of a multiset, in the row order the state builders expect.

Equivalent to ``unique(perms(v),'rows')`` without materialising the n!
intermediate: a vector with repeated entries costs its multinomial count
n!/prod(c_i!) rather than n!.

ROW ORDER IS PART OF THE CONTRACT. These rows go straight into the local state
space, whose first row (after the reversal the builders apply) is the default
initial state, so reordering them moves which state a chain starts in. This
module lists them in ASCENDING LEXICOGRAPHIC order of the sorted multiset, for
every input including the all-distinct case.

That differs from the MATLAB twin ``matlab/util/multiset_perms.m``, which for an
all-distinct vector defers to MATLAB's ``perms`` and so lists that one case in
REVERSE lexicographic order. The divergence is pre-existing and deliberate here:
the native state spaces were enumerated and their goldens taken against the
ascending listing, so normalising the two would re-baseline them. Do not
"fix" one to match the other without re-generating both sides.

THE CAP. A station holding many jobs of many classes has a factorial-sized
buffer ordering: 5 classes with 3 jobs each is 14!/(2!*3!^4) = 33.6M rows, which
does not fit and does not finish. Above ``MAX_MULTISET_PERMS`` rows this returns
the single sorted representative instead of enumerating, which the callers treat
as one equivalent ordering; it perturbs only the initial transient.
"""

from math import factorial

__all__ = ["multiset_perms", "MAX_MULTISET_PERMS"]

# Row budget above which a single representative ordering is returned instead of
# the full enumeration. Keep in step with the native SSA seeding, which assumes
# a class-ascending representative row.
MAX_MULTISET_PERMS = 5000


def multiset_perms(vec):
    """
    All distinct permutations of ``vec``, ascending lexicographic.

    Args:
        vec: iterable of comparable elements, repeats allowed.

    Returns:
        list of tuples, one per distinct permutation. Empty list for empty
        input. A single sorted tuple when the enumeration would exceed
        ``MAX_MULTISET_PERMS`` rows.

    Examples:
        >>> multiset_perms([1, 1, 2])
        [(1, 1, 2), (1, 2, 1), (2, 1, 1)]
        >>> multiset_perms([])
        []
    """
    items = sorted(vec)
    n = len(items)
    if n == 0:
        return []

    # Multiplicity of each distinct value, and hence the exact row count. The
    # product of binomials keeps this an integer without going through n!.
    counts = []
    run = 1
    for i in range(1, n):
        if items[i] == items[i - 1]:
            run += 1
        else:
            counts.append(run)
            run = 1
    counts.append(run)

    if len(counts) == 1:
        # one distinct value: every arrangement is the same one
        return [tuple(items)]

    rows = factorial(n)
    for c in counts:
        rows //= factorial(c)
    if rows > MAX_MULTISET_PERMS:
        return [tuple(items)]

    # Knuth TAOCP 7.2.1.2 Algorithm L: starting from the smallest arrangement,
    # step to the next one in lexicographic order until the sequence is
    # exhausted. It visits each distinct permutation exactly once, so repeated
    # values never produce a duplicate row and no de-duplication pass is needed.
    out = [tuple(items)]
    while True:
        # rightmost position whose successor is strictly larger
        j = n - 2
        while j >= 0 and items[j] >= items[j + 1]:
            j -= 1
        if j < 0:
            break
        # rightmost value strictly larger than items[j], then reverse the tail
        k = n - 1
        while items[j] >= items[k]:
            k -= 1
        items[j], items[k] = items[k], items[j]
        items[j + 1:] = reversed(items[j + 1:])
        out.append(tuple(items))

    return out
