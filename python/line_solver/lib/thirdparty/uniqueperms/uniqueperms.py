"""
Generate all unique permutations of a vector with possibly duplicate elements.

Ported from MATLAB to Python. Original MATLAB implementation by:
    Author: John D'Errico
    e-mail: woodchips@rochester.rr.com
    Release: 1.0
    Release date: 2/25/08

The result is equivalent to unique(perms(vec), 'rows') in MATLAB,
but generated directly via a recursive multiset permutation algorithm
that avoids the O(n!) overhead of generating all permutations first.

Includes overflow protection: when the number of unique permutations
exceeds MAX_UNIQUE_PERMS, a single sorted representative permutation
is returned to prevent combinatorial explosion.
"""

from math import factorial

# Maximum number of unique permutations to generate before falling back to
# a single representative state. Prevents combinatorial explosion for large queues.
MAX_UNIQUE_PERMS = 5000


def uniqueperms(vec):
    """
    Generate all unique permutations of a vector with possibly replicate elements.

    Efficient multiset permutation algorithm matching MATLAB's uniqueperms.
    Avoids the O(n!) overhead of set(permutations(vec)) by directly generating
    only unique permutations.

    For single-value vectors (e.g., [1,1,1,...]), returns immediately with one row.
    For large multisets where the number of unique permutations exceeds
    MAX_UNIQUE_PERMS, returns a single representative (sorted) permutation
    to prevent combinatorial explosion.

    Args:
        vec: List or 1D array of elements (possibly with repeats)

    Returns:
        List of tuples, each a unique permutation

    Example:
        >>> uniqueperms([1, 1, 1, 2, 2])
        [(1, 1, 1, 2, 2), (1, 1, 2, 1, 2), (1, 1, 2, 2, 1), ...]
    """
    vec = list(vec)

    if len(vec) == 0:
        return []

    # Count multiplicities
    counts = {}
    for v in vec:
        counts[v] = counts.get(v, 0) + 1

    unique_vals = sorted(counts.keys())

    # Single unique value: only one permutation
    if len(unique_vals) == 1:
        return [tuple(vec)]

    # All elements unique: just generate all permutations
    n = len(vec)
    if n == len(unique_vals):
        # Estimate is n! which we check against the cap
        try:
            num_perms = factorial(n)
        except (OverflowError, ValueError):
            num_perms = MAX_UNIQUE_PERMS + 1

        if num_perms > MAX_UNIQUE_PERMS:
            return [tuple(sorted(vec))]
        return uniqueperms_recursive(sorted(vec))

    # 2 or more unique elements, at least one replicate
    # Estimate number of unique permutations: n! / (a1! * a2! * ... * ak!)
    try:
        num_perms = factorial(n)
        for c in counts.values():
            num_perms //= factorial(c)
    except (OverflowError, ValueError):
        num_perms = MAX_UNIQUE_PERMS + 1

    if num_perms > MAX_UNIQUE_PERMS:
        # Too many permutations -- return single sorted representative
        return [tuple(sorted(vec))]

    # Generate unique permutations recursively (matching MATLAB uniqueperms algorithm)
    return uniqueperms_recursive(sorted(vec))


def uniqueperms_recursive(vec):
    """
    Recursive multiset permutation generator.

    Generates all unique permutations of the given vector by choosing each
    distinct element as the first element, then recursively permuting
    the remainder.

    This directly mirrors the MATLAB recursive algorithm:
        for each unique element u in vec:
            remove first occurrence of u from vec -> rest
            prepend u to each uniqueperms(rest)

    Args:
        vec: List of elements (should be sorted for canonical ordering)

    Returns:
        List of tuples, each a unique permutation
    """
    if len(vec) <= 1:
        return [tuple(vec)]

    result = []
    seen = set()
    for i, val in enumerate(vec):
        if val in seen:
            continue
        seen.add(val)
        rest = vec[:i] + vec[i + 1:]
        for perm in uniqueperms_recursive(rest):
            result.append((val,) + perm)

    return result
