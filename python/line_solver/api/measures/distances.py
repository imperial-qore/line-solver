"""
Statistical Distance and Divergence Measures.

Native Python implementations of distance and divergence measures
for comparing probability distributions.

References:
    Cha, S.-H., "Comprehensive Survey on Distance/Similarity Measures between
    Probability Density Functions." International Journal of Mathematical
    Models and Methods in Applied Sciences, 2007.
"""

import numpy as np
from typing import Union

# Type alias for array-like inputs
ArrayLike = Union[np.ndarray, list]


def _to_array(x: ArrayLike) -> np.ndarray:
    """Convert input to numpy array and flatten."""
    return np.asarray(x, dtype=np.float64).ravel()


# =============================================================================
# Information-theoretic measures (Shannon's entropy family)
# =============================================================================

def ms_entropy(x: ArrayLike) -> float:
    """
    Compute entropy H(x) of a discrete variable.

    Args:
        x: Array of discrete values.

    Returns:
        Shannon entropy in bits.
    """
    x = _to_array(x)
    if len(x) == 0:
        return 0.0

    # Count occurrences of each unique value
    _, counts = np.unique(x.astype(int), return_counts=True)
    n = len(x)
    probs = counts / n

    # Calculate entropy
    entropy = -np.sum(probs * np.log2(probs + 1e-15))
    return max(0.0, float(entropy))


def ms_jointentropy(x: ArrayLike, y: ArrayLike) -> float:
    """
    Compute joint entropy H(x,y) of two discrete variables.

    Args:
        x: First array of discrete values.
        y: Second array of discrete values.

    Returns:
        Joint Shannon entropy in bits.
    """
    x = _to_array(x)
    y = _to_array(y)

    if len(x) == 0 or len(y) == 0:
        return 0.0

    n = len(x)
    # Create joint pairs
    pairs = list(zip(x.astype(int), y.astype(int)))

    # Count occurrences of each unique pair
    from collections import Counter
    counts = Counter(pairs)
    probs = np.array(list(counts.values())) / n

    # Calculate joint entropy
    entropy = -np.sum(probs * np.log2(probs + 1e-15))
    return max(0.0, float(entropy))


def ms_condentropy(x: ArrayLike, y: ArrayLike) -> float:
    """
    Compute conditional entropy H(x|y).

    Args:
        x: Target variable.
        y: Conditioning variable.

    Returns:
        Conditional entropy H(x|y) in bits.
    """
    return ms_jointentropy(x, y) - ms_entropy(y)


def ms_mutinfo(x: ArrayLike, y: ArrayLike) -> float:
    """
    Compute mutual information I(x;y).

    Args:
        x: First variable.
        y: Second variable.

    Returns:
        Mutual information in bits.
    """
    return ms_entropy(x) + ms_entropy(y) - ms_jointentropy(x, y)


def ms_nmi(x: ArrayLike, y: ArrayLike) -> float:
    """
    Compute normalized mutual information.

    Args:
        x: First variable.
        y: Second variable.

    Returns:
        NMI in [0, 1].
    """
    hx = ms_entropy(x)
    hy = ms_entropy(y)

    if hx == 0 or hy == 0:
        return 0.0

    mi = ms_mutinfo(x, y)
    return float(2 * mi / (hx + hy))


def ms_nvi(x: ArrayLike, y: ArrayLike) -> float:
    """
    Compute normalized variation of information.

    Args:
        x: First variable.
        y: Second variable.

    Returns:
        NVI distance.
    """
    hxy = ms_jointentropy(x, y)
    if hxy == 0:
        return 0.0

    mi = ms_mutinfo(x, y)
    return float(1.0 - mi / hxy)


def ms_kullbackleibler(P: ArrayLike, Q: ArrayLike) -> float:
    """
    Kullback-Leibler divergence D_KL(P||Q).

    Measures information lost when Q approximates P.

    Args:
        P: True probability distribution.
        Q: Approximating distribution.

    Returns:
        KL divergence (non-negative, unbounded).
    """
    P = _to_array(P)
    Q = _to_array(Q)

    if len(P) != len(Q):
        raise ValueError("Distributions must have the same size")

    kl = 0.0
    for i in range(len(P)):
        if P[i] > 0:
            if Q[i] == 0:
                return np.inf
            kl += P[i] * np.log(P[i] / Q[i])

    return float(kl)


def ms_relatentropy(P: ArrayLike, Q: ArrayLike) -> float:
    """
    Relative entropy (same as Kullback-Leibler divergence).

    Args:
        P: True probability distribution.
        Q: Approximating distribution.

    Returns:
        Relative entropy.
    """
    return ms_kullbackleibler(P, Q)


def ms_jensenshannon(P: ArrayLike, Q: ArrayLike) -> float:
    """
    Jensen-Shannon divergence (symmetric KL divergence).

    Args:
        P: First probability distribution.
        Q: Second probability distribution.

    Returns:
        JS divergence in [0, log(2)].
    """
    P = _to_array(P)
    Q = _to_array(Q)

    if len(P) != len(Q):
        raise ValueError("Distributions must have the same size")

    sum1 = 0.0
    sum2 = 0.0

    for i in range(len(P)):
        m = (P[i] + Q[i]) / 2.0
        if P[i] > 0 and m > 0:
            sum1 += P[i] * np.log(P[i] / m)
        if Q[i] > 0 and m > 0:
            sum2 += Q[i] * np.log(Q[i] / m)

    return float(0.5 * (sum1 + sum2))


def ms_jeffreys(P: ArrayLike, Q: ArrayLike) -> float:
    """
    Jeffreys divergence (symmetric KL).

    J(P,Q) = KL(P||Q) + KL(Q||P)

    Args:
        P: First probability distribution.
        Q: Second probability distribution.

    Returns:
        Jeffreys divergence.
    """
    return ms_kullbackleibler(P, Q) + ms_kullbackleibler(Q, P)


def ms_kdivergence(P: ArrayLike, Q: ArrayLike) -> float:
    """
    K-divergence.

    Args:
        P: First probability distribution.
        Q: Second probability distribution.

    Returns:
        K-divergence.
    """
    P = _to_array(P)
    Q = _to_array(Q)

    if len(P) != len(Q):
        raise ValueError("Distributions must have the same size")

    result = 0.0
    for i in range(len(P)):
        if P[i] > 0:
            m = (P[i] + Q[i]) / 2.0
            if m > 0:
                result += P[i] * np.log(P[i] / m)

    return float(result)


def ms_topsoe(P: ArrayLike, Q: ArrayLike) -> float:
    """
    Topsoe divergence.

    Args:
        P: First probability distribution.
        Q: Second probability distribution.

    Returns:
        Topsoe divergence.
    """
    return ms_kdivergence(P, Q) + ms_kdivergence(Q, P)


def ms_jensendifference(P: ArrayLike, Q: ArrayLike) -> float:
    """
    Jensen difference divergence.

    Args:
        P: First probability distribution.
        Q: Second probability distribution.

    Returns:
        Jensen difference.
    """
    P = _to_array(P)
    Q = _to_array(Q)

    if len(P) != len(Q):
        raise ValueError("Distributions must have the same size")

    result = 0.0
    for i in range(len(P)):
        m = (P[i] + Q[i]) / 2.0
        hp = -P[i] * np.log(P[i] + 1e-15) if P[i] > 0 else 0.0
        hq = -Q[i] * np.log(Q[i] + 1e-15) if Q[i] > 0 else 0.0
        hm = -m * np.log(m + 1e-15) if m > 0 else 0.0
        result += hm - (hp + hq) / 2.0

    return float(result)


def ms_taneja(P: ArrayLike, Q: ArrayLike) -> float:
    """
    Taneja divergence.

    Args:
        P: First probability distribution.
        Q: Second probability distribution.

    Returns:
        Taneja divergence.
    """
    P = _to_array(P)
    Q = _to_array(Q)

    if len(P) != len(Q):
        raise ValueError("Distributions must have the same size")

    result = 0.0
    for i in range(len(P)):
        if P[i] > 0 and Q[i] > 0:
            m = (P[i] + Q[i]) / 2.0
            result += m * np.log(m / np.sqrt(P[i] * Q[i]))

    return float(result)


# =============================================================================
# Minkowski family (Lp norms)
# =============================================================================

def ms_minkowski(P: ArrayLike, Q: ArrayLike, p: float = 2.0) -> float:
    """
    Minkowski distance (Lp norm).

    Args:
        P: First distribution.
        Q: Second distribution.
        p: Order parameter (p >= 1).

    Returns:
        Minkowski distance.
    """
    P = _to_array(P)
    Q = _to_array(Q)

    if len(P) != len(Q):
        raise ValueError("Distributions must have the same size")
    if p < 1:
        raise ValueError("Order parameter p must be >= 1")

    return float(np.sum(np.abs(P - Q) ** p) ** (1.0 / p))


def ms_euclidean(P: ArrayLike, Q: ArrayLike) -> float:
    """
    Euclidean distance (L2 norm).

    Args:
        P: First distribution.
        Q: Second distribution.

    Returns:
        Euclidean distance.
    """
    P = _to_array(P)
    Q = _to_array(Q)

    if len(P) != len(Q):
        raise ValueError("Distributions must have the same size")

    return float(np.sqrt(np.sum((P - Q) ** 2)))


def ms_squaredeuclidean(P: ArrayLike, Q: ArrayLike) -> float:
    """
    Squared Euclidean distance.

    Args:
        P: First distribution.
        Q: Second distribution.

    Returns:
        Squared Euclidean distance.
    """
    P = _to_array(P)
    Q = _to_array(Q)

    if len(P) != len(Q):
        raise ValueError("Distributions must have the same size")

    return float(np.sum((P - Q) ** 2))


def ms_cityblock(P: ArrayLike, Q: ArrayLike) -> float:
    """
    City block (Manhattan) distance (L1 norm).

    Args:
        P: First distribution.
        Q: Second distribution.

    Returns:
        Manhattan distance.
    """
    P = _to_array(P)
    Q = _to_array(Q)

    if len(P) != len(Q):
        raise ValueError("Distributions must have the same size")

    return float(np.sum(np.abs(P - Q)))


def ms_chebyshev(P: ArrayLike, Q: ArrayLike) -> float:
    """
    Chebyshev distance (L-infinity norm).

    Args:
        P: First distribution.
        Q: Second distribution.

    Returns:
        Chebyshev distance.
    """
    P = _to_array(P)
    Q = _to_array(Q)

    if len(P) != len(Q):
        raise ValueError("Distributions must have the same size")

    return float(np.max(np.abs(P - Q)))


def ms_sorensen(P: ArrayLike, Q: ArrayLike) -> float:
    """
    Sorensen distance.

    Args:
        P: First distribution.
        Q: Second distribution.

    Returns:
        Sorensen distance.
    """
    P = _to_array(P)
    Q = _to_array(Q)

    if len(P) != len(Q):
        raise ValueError("Distributions must have the same size")

    num = np.sum(np.abs(P - Q))
    den = np.sum(P + Q)

    return float(num / den) if den > 0 else 0.0


def ms_gower(P: ArrayLike, Q: ArrayLike) -> float:
    """
    Gower distance.

    Args:
        P: First distribution.
        Q: Second distribution.

    Returns:
        Gower distance.
    """
    P = _to_array(P)
    Q = _to_array(Q)

    if len(P) != len(Q):
        raise ValueError("Distributions must have the same size")

    return float(np.sum(np.abs(P - Q)) / len(P))


def ms_soergel(P: ArrayLike, Q: ArrayLike) -> float:
    """
    Soergel distance.

    Args:
        P: First distribution.
        Q: Second distribution.

    Returns:
        Soergel distance.
    """
    P = _to_array(P)
    Q = _to_array(Q)

    if len(P) != len(Q):
        raise ValueError("Distributions must have the same size")

    num = np.sum(np.abs(P - Q))
    den = np.sum(np.maximum(P, Q))

    return float(num / den) if den > 0 else 0.0


def ms_lorentzian(P: ArrayLike, Q: ArrayLike) -> float:
    """
    Lorentzian distance.

    Args:
        P: First distribution.
        Q: Second distribution.

    Returns:
        Lorentzian distance.
    """
    P = _to_array(P)
    Q = _to_array(Q)

    if len(P) != len(Q):
        raise ValueError("Distributions must have the same size")

    return float(np.sum(np.log1p(np.abs(P - Q))))


def ms_canberra(P: ArrayLike, Q: ArrayLike) -> float:
    """
    Canberra distance.

    Args:
        P: First distribution.
        Q: Second distribution.

    Returns:
        Canberra distance.
    """
    P = _to_array(P)
    Q = _to_array(Q)

    if len(P) != len(Q):
        raise ValueError("Distributions must have the same size")

    result = 0.0
    for i in range(len(P)):
        den = abs(P[i]) + abs(Q[i])
        if den > 0:
            result += abs(P[i] - Q[i]) / den

    return float(result)


def ms_avgl1linfty(P: ArrayLike, Q: ArrayLike) -> float:
    """
    Average of L1 and L-infinity distances.

    Args:
        P: First distribution.
        Q: Second distribution.

    Returns:
        Average L1 and L-infinity distance.
    """
    l1 = ms_cityblock(P, Q)
    linf = ms_chebyshev(P, Q)
    return float((l1 + linf) / 2.0)


# =============================================================================
# Inner product family
# =============================================================================

def ms_product(P: ArrayLike, Q: ArrayLike) -> float:
    """
    Inner product (dot product).

    Args:
        P: First distribution.
        Q: Second distribution.

    Returns:
        Inner product.
    """
    P = _to_array(P)
    Q = _to_array(Q)

    if len(P) != len(Q):
        raise ValueError("Distributions must have the same size")

    return float(np.dot(P, Q))


def ms_cosine(P: ArrayLike, Q: ArrayLike) -> float:
    """
    Cosine distance (1 - cosine similarity).

    Args:
        P: First distribution.
        Q: Second distribution.

    Returns:
        Cosine distance.
    """
    P = _to_array(P)
    Q = _to_array(Q)

    if len(P) != len(Q):
        raise ValueError("Distributions must have the same size")

    dot = np.dot(P, Q)
    norm_p = np.sqrt(np.sum(P ** 2))
    norm_q = np.sqrt(np.sum(Q ** 2))

    if norm_p == 0 or norm_q == 0:
        return 1.0

    return float(1.0 - dot / (norm_p * norm_q))


def ms_jaccard(P: ArrayLike, Q: ArrayLike) -> float:
    """
    Jaccard distance.

    Args:
        P: First distribution.
        Q: Second distribution.

    Returns:
        Jaccard distance.
    """
    P = _to_array(P)
    Q = _to_array(Q)

    if len(P) != len(Q):
        raise ValueError("Distributions must have the same size")

    dot = np.dot(P, Q)
    sum_p2 = np.sum(P ** 2)
    sum_q2 = np.sum(Q ** 2)

    den = sum_p2 + sum_q2 - dot
    return float(1.0 - dot / den) if den > 0 else 0.0


def ms_dice(P: ArrayLike, Q: ArrayLike) -> float:
    """
    Dice distance.

    Args:
        P: First distribution.
        Q: Second distribution.

    Returns:
        Dice distance.
    """
    P = _to_array(P)
    Q = _to_array(Q)

    if len(P) != len(Q):
        raise ValueError("Distributions must have the same size")

    dot = np.dot(P, Q)
    sum_p2 = np.sum(P ** 2)
    sum_q2 = np.sum(Q ** 2)

    den = sum_p2 + sum_q2
    return float(1.0 - 2 * dot / den) if den > 0 else 0.0


def ms_kumarhassebrook(P: ArrayLike, Q: ArrayLike) -> float:
    """
    Kumar-Hassebrook distance.

    Args:
        P: First distribution.
        Q: Second distribution.

    Returns:
        Kumar-Hassebrook distance.
    """
    P = _to_array(P)
    Q = _to_array(Q)

    if len(P) != len(Q):
        raise ValueError("Distributions must have the same size")

    dot = np.dot(P, Q)
    sum_p2 = np.sum(P ** 2)
    sum_q2 = np.sum(Q ** 2)

    den = sum_p2 + sum_q2 - dot
    return float(dot / den) if den > 0 else 0.0


def ms_harmonicmean(P: ArrayLike, Q: ArrayLike) -> float:
    """
    Harmonic mean similarity.

    Args:
        P: First distribution.
        Q: Second distribution.

    Returns:
        Harmonic mean similarity.
    """
    P = _to_array(P)
    Q = _to_array(Q)

    if len(P) != len(Q):
        raise ValueError("Distributions must have the same size")

    result = 0.0
    for i in range(len(P)):
        den = P[i] + Q[i]
        if den > 0:
            result += P[i] * Q[i] / den

    return float(2 * result)


# =============================================================================
# Fidelity and squared chord family
# =============================================================================

def ms_fidelity(P: ArrayLike, Q: ArrayLike) -> float:
    """
    Fidelity (Bhattacharyya coefficient).

    Args:
        P: First distribution.
        Q: Second distribution.

    Returns:
        Fidelity coefficient in [0, 1].
    """
    P = _to_array(P)
    Q = _to_array(Q)

    if len(P) != len(Q):
        raise ValueError("Distributions must have the same size")

    return float(np.sum(np.sqrt(P * Q)))


def ms_bhattacharyya(P: ArrayLike, Q: ArrayLike) -> float:
    """
    Bhattacharyya distance.

    Args:
        P: First distribution.
        Q: Second distribution.

    Returns:
        Bhattacharyya distance.
    """
    bc = ms_fidelity(P, Q)
    return float(-np.log(bc + 1e-15))


def ms_hellinger(P: ArrayLike, Q: ArrayLike) -> float:
    """
    Hellinger distance.

    Args:
        P: First distribution.
        Q: Second distribution.

    Returns:
        Hellinger distance in [0, sqrt(2)].
    """
    P = _to_array(P)
    Q = _to_array(Q)

    if len(P) != len(Q):
        raise ValueError("Distributions must have the same size")

    diff = np.sqrt(P) - np.sqrt(Q)
    return float(np.sqrt(2.0 * np.sum(diff ** 2)))


def ms_matusita(P: ArrayLike, Q: ArrayLike) -> float:
    """
    Matusita distance.

    Args:
        P: First distribution.
        Q: Second distribution.

    Returns:
        Matusita distance.
    """
    P = _to_array(P)
    Q = _to_array(Q)

    if len(P) != len(Q):
        raise ValueError("Distributions must have the same size")

    diff = np.sqrt(P) - np.sqrt(Q)
    return float(np.sqrt(np.sum(diff ** 2)))


def ms_squaredchord(P: ArrayLike, Q: ArrayLike) -> float:
    """
    Squared chord distance.

    Args:
        P: First distribution.
        Q: Second distribution.

    Returns:
        Squared chord distance.
    """
    P = _to_array(P)
    Q = _to_array(Q)

    if len(P) != len(Q):
        raise ValueError("Distributions must have the same size")

    diff = np.sqrt(P) - np.sqrt(Q)
    return float(np.sum(diff ** 2))


# =============================================================================
# Chi-squared family
# =============================================================================

def ms_pearsonchisquared(P: ArrayLike, Q: ArrayLike) -> float:
    """
    Pearson chi-squared distance.

    Args:
        P: First distribution.
        Q: Second distribution.

    Returns:
        Pearson chi-squared statistic.
    """
    P = _to_array(P)
    Q = _to_array(Q)

    if len(P) != len(Q):
        raise ValueError("Distributions must have the same size")

    result = 0.0
    for i in range(len(P)):
        if Q[i] > 0:
            result += (P[i] - Q[i]) ** 2 / Q[i]

    return float(result)


def ms_neymanchisquared(P: ArrayLike, Q: ArrayLike) -> float:
    """
    Neyman chi-squared distance.

    Args:
        P: First distribution.
        Q: Second distribution.

    Returns:
        Neyman chi-squared statistic.
    """
    P = _to_array(P)
    Q = _to_array(Q)

    if len(P) != len(Q):
        raise ValueError("Distributions must have the same size")

    result = 0.0
    for i in range(len(P)):
        if P[i] > 0:
            result += (P[i] - Q[i]) ** 2 / P[i]

    return float(result)


def ms_squaredchisquared(P: ArrayLike, Q: ArrayLike) -> float:
    """
    Squared chi-squared distance.

    Args:
        P: First distribution.
        Q: Second distribution.

    Returns:
        Squared chi-squared.
    """
    P = _to_array(P)
    Q = _to_array(Q)

    if len(P) != len(Q):
        raise ValueError("Distributions must have the same size")

    result = 0.0
    for i in range(len(P)):
        den = P[i] + Q[i]
        if den > 0:
            result += (P[i] - Q[i]) ** 2 / den

    return float(result)


def ms_probsymmchisquared(P: ArrayLike, Q: ArrayLike) -> float:
    """
    Probabilistic symmetric chi-squared distance.

    Args:
        P: First distribution.
        Q: Second distribution.

    Returns:
        Probabilistic symmetric chi-squared.
    """
    return 2.0 * ms_squaredchisquared(P, Q)


def ms_divergence(P: ArrayLike, Q: ArrayLike) -> float:
    """
    Divergence (chi-squared type).

    Args:
        P: First distribution.
        Q: Second distribution.

    Returns:
        Divergence measure.
    """
    P = _to_array(P)
    Q = _to_array(Q)

    if len(P) != len(Q):
        raise ValueError("Distributions must have the same size")

    result = 0.0
    for i in range(len(P)):
        den = (P[i] + Q[i]) ** 2
        if den > 0:
            result += (P[i] - Q[i]) ** 2 / den

    return float(2 * result)


def ms_additivesymmetricchisquared(P: ArrayLike, Q: ArrayLike) -> float:
    """
    Additive symmetric chi-squared distance.

    Args:
        P: First distribution.
        Q: Second distribution.

    Returns:
        Additive symmetric chi-squared.
    """
    P = _to_array(P)
    Q = _to_array(Q)

    if len(P) != len(Q):
        raise ValueError("Distributions must have the same size")

    result = 0.0
    for i in range(len(P)):
        prod = P[i] * Q[i]
        if prod > 0:
            result += (P[i] - Q[i]) ** 2 * (P[i] + Q[i]) / prod

    return float(result)


def ms_clark(P: ArrayLike, Q: ArrayLike) -> float:
    """
    Clark distance.

    Args:
        P: First distribution.
        Q: Second distribution.

    Returns:
        Clark distance.
    """
    P = _to_array(P)
    Q = _to_array(Q)

    if len(P) != len(Q):
        raise ValueError("Distributions must have the same size")

    result = 0.0
    for i in range(len(P)):
        den = abs(P[i]) + abs(Q[i])
        if den > 0:
            result += ((P[i] - Q[i]) / den) ** 2

    return float(np.sqrt(result))


# =============================================================================
# Intersection family
# =============================================================================

def ms_intersection(P: ArrayLike, Q: ArrayLike) -> float:
    """
    Intersection distance.

    Args:
        P: First distribution.
        Q: Second distribution.

    Returns:
        Intersection distance.
    """
    P = _to_array(P)
    Q = _to_array(Q)

    if len(P) != len(Q):
        raise ValueError("Distributions must have the same size")

    return float(np.sum(np.minimum(P, Q)))


def ms_czekanowski(P: ArrayLike, Q: ArrayLike) -> float:
    """
    Czekanowski distance.

    Args:
        P: First distribution.
        Q: Second distribution.

    Returns:
        Czekanowski distance.
    """
    P = _to_array(P)
    Q = _to_array(Q)

    if len(P) != len(Q):
        raise ValueError("Distributions must have the same size")

    intersection = np.sum(np.minimum(P, Q))
    return float(1.0 - 2 * intersection / np.sum(P + Q))


def ms_motyka(P: ArrayLike, Q: ArrayLike) -> float:
    """
    Motyka distance.

    Args:
        P: First distribution.
        Q: Second distribution.

    Returns:
        Motyka distance.
    """
    P = _to_array(P)
    Q = _to_array(Q)

    if len(P) != len(Q):
        raise ValueError("Distributions must have the same size")

    intersection = np.sum(np.minimum(P, Q))
    total = np.sum(P + Q)

    return float(1.0 - intersection / total) if total > 0 else 0.0


def ms_kulczynskis(P: ArrayLike, Q: ArrayLike) -> float:
    """
    Kulczynski similarity.

    Args:
        P: First distribution.
        Q: Second distribution.

    Returns:
        Kulczynski similarity.
    """
    P = _to_array(P)
    Q = _to_array(Q)

    if len(P) != len(Q):
        raise ValueError("Distributions must have the same size")

    intersection = np.sum(np.minimum(P, Q))
    diff = np.sum(np.abs(P - Q))

    return float(intersection / diff) if diff > 0 else np.inf


def ms_kulczynskid(P: ArrayLike, Q: ArrayLike) -> float:
    """
    Kulczynski distance.

    Args:
        P: First distribution.
        Q: Second distribution.

    Returns:
        Kulczynski distance.
    """
    s = ms_kulczynskis(P, Q)
    return float(1.0 / s) if s > 0 and not np.isinf(s) else 0.0


def ms_ruzicka(P: ArrayLike, Q: ArrayLike) -> float:
    """
    Ruzicka similarity.

    Args:
        P: First distribution.
        Q: Second distribution.

    Returns:
        Ruzicka similarity.
    """
    P = _to_array(P)
    Q = _to_array(Q)

    if len(P) != len(Q):
        raise ValueError("Distributions must have the same size")

    intersection = np.sum(np.minimum(P, Q))
    union = np.sum(np.maximum(P, Q))

    return float(intersection / union) if union > 0 else 0.0


def ms_tanimoto(P: ArrayLike, Q: ArrayLike) -> float:
    """
    Tanimoto distance.

    Args:
        P: First distribution.
        Q: Second distribution.

    Returns:
        Tanimoto distance.
    """
    P = _to_array(P)
    Q = _to_array(Q)

    if len(P) != len(Q):
        raise ValueError("Distributions must have the same size")

    diff = np.sum(np.abs(P - Q))
    max_sum = np.sum(np.maximum(P, Q))

    return float(diff / max_sum) if max_sum > 0 else 0.0


def ms_wavehegdes(P: ArrayLike, Q: ArrayLike) -> float:
    """
    Wave-Hedges distance.

    Args:
        P: First distribution.
        Q: Second distribution.

    Returns:
        Wave-Hedges distance.
    """
    P = _to_array(P)
    Q = _to_array(Q)

    if len(P) != len(Q):
        raise ValueError("Distributions must have the same size")

    result = 0.0
    for i in range(len(P)):
        max_val = max(abs(P[i]), abs(Q[i]))
        if max_val > 0:
            result += abs(P[i] - Q[i]) / max_val

    return float(result)


def ms_kumarjohnson(P: ArrayLike, Q: ArrayLike) -> float:
    """
    Kumar-Johnson distance.

    Args:
        P: First distribution.
        Q: Second distribution.

    Returns:
        Kumar-Johnson distance.
    """
    P = _to_array(P)
    Q = _to_array(Q)

    if len(P) != len(Q):
        raise ValueError("Distributions must have the same size")

    result = 0.0
    for i in range(len(P)):
        prod = P[i] * Q[i]
        if prod > 0:
            result += (P[i] ** 2 - Q[i] ** 2) ** 2 / (2 * prod ** 1.5)

    return float(result)


# =============================================================================
# Empirical distribution tests
# =============================================================================

def ms_kolmogorov_smirnov(X: ArrayLike, Y: ArrayLike) -> float:
    """
    Kolmogorov-Smirnov distance between two empirical distributions.

    Args:
        X: First sample.
        Y: Second sample.

    Returns:
        KS distance (max absolute CDF difference).
    """
    X = _to_array(X)
    Y = _to_array(Y)

    # Remove NaN values
    X = X[~np.isnan(X)]
    Y = Y[~np.isnan(Y)]

    if len(X) == 0 or len(Y) == 0:
        return 0.0

    # Sort both arrays
    X_sorted = np.sort(X)
    Y_sorted = np.sort(Y)

    # Combine and sort all values
    combined = np.sort(np.concatenate([X, Y]))

    # Compute empirical CDFs
    nx, ny = len(X), len(Y)
    max_dist = 0.0

    for val in combined:
        ecdf_x = np.sum(X_sorted <= val) / nx
        ecdf_y = np.sum(Y_sorted <= val) / ny
        max_dist = max(max_dist, abs(ecdf_x - ecdf_y))

    return float(max_dist)


def ms_kuiper(X: ArrayLike, Y: ArrayLike) -> float:
    """
    Kuiper distance between two empirical distributions.

    Args:
        X: First sample.
        Y: Second sample.

    Returns:
        Kuiper distance (sum of max positive and negative CDF differences).
    """
    X = _to_array(X)
    Y = _to_array(Y)

    # Remove NaN values
    X = X[~np.isnan(X)]
    Y = Y[~np.isnan(Y)]

    if len(X) == 0 or len(Y) == 0:
        return 0.0

    # Sort both arrays
    X_sorted = np.sort(X)
    Y_sorted = np.sort(Y)

    # Combine and sort all values
    combined = np.sort(np.concatenate([X, Y]))

    # Compute empirical CDFs
    nx, ny = len(X), len(Y)
    max_pos = 0.0
    max_neg = 0.0

    for val in combined:
        ecdf_x = np.sum(X_sorted <= val) / nx
        ecdf_y = np.sum(Y_sorted <= val) / ny
        diff = ecdf_x - ecdf_y
        max_pos = max(max_pos, diff)
        max_neg = max(max_neg, -diff)

    return float(max_pos + max_neg)


def ms_cramer_von_mises(X: ArrayLike, Y: ArrayLike) -> float:
    """
    Cramer-von Mises distance between two empirical distributions.

    Args:
        X: First sample.
        Y: Second sample.

    Returns:
        Cramer-von Mises statistic.
    """
    X = _to_array(X)
    Y = _to_array(Y)

    # Remove NaN values
    X = X[~np.isnan(X)]
    Y = Y[~np.isnan(Y)]

    if len(X) == 0 or len(Y) == 0:
        return 0.0

    nx, ny = len(X), len(Y)
    n = nx + ny

    # Compute ranks
    combined = np.concatenate([X, Y])
    ranks = np.argsort(np.argsort(combined)) + 1

    # Separate ranks
    ranks_x = ranks[:nx]
    ranks_y = ranks[nx:]

    # Compute statistic
    u = nx * np.sum((ranks_x - np.arange(1, nx + 1)) ** 2)
    u += ny * np.sum((ranks_y - np.arange(1, ny + 1)) ** 2)

    return float(u / (n * nx * ny))


def ms_anderson_darling(X: ArrayLike, Y: ArrayLike) -> float:
    """
    Anderson-Darling distance between two empirical distributions.

    Args:
        X: First sample.
        Y: Second sample.

    Returns:
        Anderson-Darling statistic.
    """
    X = _to_array(X)
    Y = _to_array(Y)

    # Remove NaN values
    X = X[~np.isnan(X)]
    Y = Y[~np.isnan(Y)]

    if len(X) == 0 or len(Y) == 0:
        return 0.0

    nx, ny = len(X), len(Y)
    n = nx + ny

    # Combine and sort
    combined = np.sort(np.concatenate([X, Y]))

    # Compute weighted CDF differences
    result = 0.0
    for i, val in enumerate(combined[:-1]):
        ecdf_x = np.sum(X <= val) / nx
        ecdf_y = np.sum(Y <= val) / ny
        ecdf_comb = (i + 1) / n

        if ecdf_comb > 0 and ecdf_comb < 1:
            weight = 1.0 / (ecdf_comb * (1 - ecdf_comb))
            result += weight * (ecdf_x - ecdf_y) ** 2

    return float(result * nx * ny / n)


def ms_wasserstein(X: ArrayLike, Y: ArrayLike) -> float:
    """
    Wasserstein distance (Earth Mover's Distance) between two empirical distributions.

    Args:
        X: First sample.
        Y: Second sample.

    Returns:
        1-Wasserstein distance.
    """
    X = _to_array(X)
    Y = _to_array(Y)

    # Remove NaN values
    X = X[~np.isnan(X)]
    Y = Y[~np.isnan(Y)]

    if len(X) == 0 or len(Y) == 0:
        return 0.0

    nx, ny = len(X), len(Y)

    # Create combined array with weights
    combined = []
    for x in X:
        combined.append((x, 1.0 / nx))
    for y in Y:
        combined.append((y, -1.0 / ny))

    # Sort by value
    combined.sort(key=lambda p: p[0])

    # Compute Wasserstein distance
    distance = 0.0
    cumulative_weight = 0.0

    for i in range(len(combined) - 1):
        cumulative_weight += combined[i][1]
        width = combined[i + 1][0] - combined[i][0]
        distance += abs(cumulative_weight) * width

    return float(distance)


__all__ = [
    # Entropy measures
    'ms_entropy',
    'ms_jointentropy',
    'ms_condentropy',
    'ms_mutinfo',
    'ms_nmi',
    'ms_nvi',
    # Information-theoretic divergences
    'ms_kullbackleibler',
    'ms_relatentropy',
    'ms_jensenshannon',
    'ms_jeffreys',
    'ms_kdivergence',
    'ms_topsoe',
    'ms_jensendifference',
    'ms_taneja',
    # Minkowski family
    'ms_minkowski',
    'ms_euclidean',
    'ms_squaredeuclidean',
    'ms_cityblock',
    'ms_chebyshev',
    'ms_sorensen',
    'ms_gower',
    'ms_soergel',
    'ms_lorentzian',
    'ms_canberra',
    'ms_avgl1linfty',
    # Inner product family
    'ms_product',
    'ms_cosine',
    'ms_jaccard',
    'ms_dice',
    'ms_kumarhassebrook',
    'ms_harmonicmean',
    # Fidelity family
    'ms_fidelity',
    'ms_bhattacharyya',
    'ms_hellinger',
    'ms_matusita',
    'ms_squaredchord',
    # Chi-squared family
    'ms_pearsonchisquared',
    'ms_neymanchisquared',
    'ms_squaredchisquared',
    'ms_probsymmchisquared',
    'ms_divergence',
    'ms_additivesymmetricchisquared',
    'ms_clark',
    # Intersection family
    'ms_intersection',
    'ms_czekanowski',
    'ms_motyka',
    'ms_kulczynskis',
    'ms_kulczynskid',
    'ms_ruzicka',
    'ms_tanimoto',
    'ms_wavehegdes',
    'ms_kumarjohnson',
    # Empirical tests
    'ms_kolmogorov_smirnov',
    'ms_kuiper',
    'ms_cramer_von_mises',
    'ms_anderson_darling',
    'ms_wasserstein',
]
