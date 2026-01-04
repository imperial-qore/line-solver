"""
Metrics and similarity/distance measures API wrappers for LINE solver.

This module provides Python wrappers for the MS (Metrics and Similarity)
functions available in the Java/Kotlin backend.
"""

import jpype
from .. import jlineMatrixFromArray, jlineMatrixToArray


def ms_additivesymmetricchisquared(p, q):
    """
    Computes the additive symmetric chi-squared divergence between distributions.

    Args:
        p: First probability distribution.
        q: Second probability distribution.

    Returns:
        float: The additive symmetric chi-squared divergence.
    """
    result = jpype.JPackage('jline').api.ms.Ms_additivesymmetricchisquaredKt.ms_additivesymmetricchisquared(
        jlineMatrixFromArray(p), jlineMatrixFromArray(q)
    )
    return float(result)


def ms_adtest(data1, data2):
    """
    Performs Anderson-Darling test for comparing two samples.

    Args:
        data1: First data sample.
        data2: Second data sample.

    Returns:
        dict: Test statistic and p-value.
    """
    result = jpype.JPackage('jline').api.ms.Ms_adtestKt.ms_adtest(
        jlineMatrixFromArray(data1), jlineMatrixFromArray(data2)
    )
    return {
        'statistic': float(result.getStatistic()),
        'pvalue': float(result.getPvalue())
    }


def ms_avgl1linfty(p, q):
    """
    Computes the average of L1 and L-infinity norms.

    Args:
        p: First vector.
        q: Second vector.

    Returns:
        float: The average of L1 and L-infinity distances.
    """
    result = jpype.JPackage('jline').api.ms.Ms_avgl1linftyKt.ms_avgl1linfty(
        jlineMatrixFromArray(p), jlineMatrixFromArray(q)
    )
    return float(result)


def ms_bhattacharyya(p, q):
    """
    Computes the Bhattacharyya distance between distributions.

    Args:
        p: First probability distribution.
        q: Second probability distribution.

    Returns:
        float: The Bhattacharyya distance.
    """
    result = jpype.JPackage('jline').api.ms.Ms_bhattacharyyaKt.ms_bhattacharyya(
        jlineMatrixFromArray(p), jlineMatrixFromArray(q)
    )
    return float(result)


def ms_canberra(p, q):
    """
    Computes the Canberra distance between vectors.

    Args:
        p: First vector.
        q: Second vector.

    Returns:
        float: The Canberra distance.
    """
    result = jpype.JPackage('jline').api.ms.Ms_canberraKt.ms_canberra(
        jlineMatrixFromArray(p), jlineMatrixFromArray(q)
    )
    return float(result)


def ms_chebyshev(p, q):
    """
    Computes the Chebyshev (L-infinity) distance between vectors.

    Args:
        p: First vector.
        q: Second vector.

    Returns:
        float: The Chebyshev distance.
    """
    result = jpype.JPackage('jline').api.ms.Ms_chebyshevKt.ms_chebyshev(
        jlineMatrixFromArray(p), jlineMatrixFromArray(q)
    )
    return float(result)


def ms_chisquared(p, q):
    """
    Computes the chi-squared statistic between distributions.

    Args:
        p: First probability distribution.
        q: Second probability distribution.

    Returns:
        float: The chi-squared statistic.
    """
    result = jpype.JPackage('jline').api.ms.Ms_chisquaredKt.ms_chisquared(
        jlineMatrixFromArray(p), jlineMatrixFromArray(q)
    )
    return float(result)


def ms_cityblock(p, q):
    """
    Computes the city block (Manhattan/L1) distance between vectors.

    Args:
        p: First vector.
        q: Second vector.

    Returns:
        float: The city block distance.
    """
    result = jpype.JPackage('jline').api.ms.Ms_cityblockKt.ms_cityblock(
        jlineMatrixFromArray(p), jlineMatrixFromArray(q)
    )
    return float(result)


def ms_clark(p, q):
    """
    Computes the Clark distance between distributions.

    Args:
        p: First probability distribution.
        q: Second probability distribution.

    Returns:
        float: The Clark distance.
    """
    result = jpype.JPackage('jline').api.ms.Ms_clarkKt.ms_clark(
        jlineMatrixFromArray(p), jlineMatrixFromArray(q)
    )
    return float(result)


def ms_condentropy(p, q):
    """
    Computes the conditional entropy H(P|Q).

    Args:
        p: First probability distribution.
        q: Second probability distribution.

    Returns:
        float: The conditional entropy.
    """
    result = jpype.JPackage('jline').api.ms.Ms_condentropyKt.ms_condentropy(
        jlineMatrixFromArray(p), jlineMatrixFromArray(q)
    )
    return float(result)


def ms_cosine(p, q):
    """
    Computes the cosine similarity between vectors.

    Args:
        p: First vector.
        q: Second vector.

    Returns:
        float: The cosine similarity (1 - cosine distance).
    """
    result = jpype.JPackage('jline').api.ms.Ms_cosineKt.ms_cosine(
        jlineMatrixFromArray(p), jlineMatrixFromArray(q)
    )
    return float(result)


def ms_cramer_von_mises(data1, data2):
    """
    Performs Cramer-von Mises test for comparing two samples.

    Args:
        data1: First data sample.
        data2: Second data sample.

    Returns:
        dict: Test statistic and p-value.
    """
    result = jpype.JPackage('jline').api.ms.Ms_cramer_von_misesKt.ms_cramer_von_mises(
        jlineMatrixFromArray(data1), jlineMatrixFromArray(data2)
    )
    return {
        'statistic': float(result.getStatistic()),
        'pvalue': float(result.getPvalue())
    }


def ms_czekanowski(p, q):
    """
    Computes the Czekanowski (Dice/Sorensen) coefficient.

    Args:
        p: First vector.
        q: Second vector.

    Returns:
        float: The Czekanowski coefficient.
    """
    result = jpype.JPackage('jline').api.ms.Ms_czekanowskiKt.ms_czekanowski(
        jlineMatrixFromArray(p), jlineMatrixFromArray(q)
    )
    return float(result)


def ms_dice(p, q):
    """
    Computes the Dice coefficient between sets.

    Args:
        p: First set (binary vector).
        q: Second set (binary vector).

    Returns:
        float: The Dice coefficient.
    """
    result = jpype.JPackage('jline').api.ms.Ms_diceKt.ms_dice(
        jlineMatrixFromArray(p), jlineMatrixFromArray(q)
    )
    return float(result)


def ms_divergence(p, q):
    """
    Computes the divergence measure between distributions.

    Args:
        p: First probability distribution.
        q: Second probability distribution.

    Returns:
        float: The divergence measure.
    """
    result = jpype.JPackage('jline').api.ms.Ms_divergenceKt.ms_divergence(
        jlineMatrixFromArray(p), jlineMatrixFromArray(q)
    )
    return float(result)


def ms_entropy(p):
    """
    Computes the Shannon entropy of a distribution.

    Args:
        p: Probability distribution.

    Returns:
        float: The Shannon entropy.
    """
    result = jpype.JPackage('jline').api.ms.Ms_entropyKt.ms_entropy(
        jlineMatrixFromArray(p)
    )
    return float(result)


def ms_euclidean(p, q):
    """
    Computes the Euclidean (L2) distance between vectors.

    Args:
        p: First vector.
        q: Second vector.

    Returns:
        float: The Euclidean distance.
    """
    result = jpype.JPackage('jline').api.ms.Ms_euclideanKt.ms_euclidean(
        jlineMatrixFromArray(p), jlineMatrixFromArray(q)
    )
    return float(result)


def ms_fidelity(p, q):
    """
    Computes the fidelity (Bhattacharyya coefficient) between distributions.

    Args:
        p: First probability distribution.
        q: Second probability distribution.

    Returns:
        float: The fidelity measure.
    """
    result = jpype.JPackage('jline').api.ms.Ms_fidelityKt.ms_fidelity(
        jlineMatrixFromArray(p), jlineMatrixFromArray(q)
    )
    return float(result)


def ms_gower(p, q):
    """
    Computes the Gower distance for mixed-type data.

    Args:
        p: First vector.
        q: Second vector.

    Returns:
        float: The Gower distance.
    """
    result = jpype.JPackage('jline').api.ms.Ms_gowerKt.ms_gower(
        jlineMatrixFromArray(p), jlineMatrixFromArray(q)
    )
    return float(result)


def ms_harmonicmean(p, q):
    """
    Computes the harmonic mean similarity.

    Args:
        p: First vector.
        q: Second vector.

    Returns:
        float: The harmonic mean similarity.
    """
    result = jpype.JPackage('jline').api.ms.Ms_harmonicmeanKt.ms_harmonicmean(
        jlineMatrixFromArray(p), jlineMatrixFromArray(q)
    )
    return float(result)


def ms_hellinger(p, q):
    """
    Computes the Hellinger distance between distributions.

    Args:
        p: First probability distribution.
        q: Second probability distribution.

    Returns:
        float: The Hellinger distance.
    """
    result = jpype.JPackage('jline').api.ms.Ms_hellingerKt.ms_hellinger(
        jlineMatrixFromArray(p), jlineMatrixFromArray(q)
    )
    return float(result)


def ms_intersection(p, q):
    """
    Computes the intersection (histogram intersection) between distributions.

    Args:
        p: First probability distribution.
        q: Second probability distribution.

    Returns:
        float: The intersection value.
    """
    result = jpype.JPackage('jline').api.ms.Ms_intersectionKt.ms_intersection(
        jlineMatrixFromArray(p), jlineMatrixFromArray(q)
    )
    return float(result)


def ms_jaccard(p, q):
    """
    Computes the Jaccard similarity coefficient.

    Args:
        p: First set (binary vector).
        q: Second set (binary vector).

    Returns:
        float: The Jaccard similarity coefficient.
    """
    result = jpype.JPackage('jline').api.ms.Ms_jaccardKt.ms_jaccard(
        jlineMatrixFromArray(p), jlineMatrixFromArray(q)
    )
    return float(result)


def ms_jeffreys(p, q):
    """
    Computes the Jeffreys divergence (symmetric KL divergence).

    Args:
        p: First probability distribution.
        q: Second probability distribution.

    Returns:
        float: The Jeffreys divergence.
    """
    result = jpype.JPackage('jline').api.ms.Ms_jeffreysKt.ms_jeffreys(
        jlineMatrixFromArray(p), jlineMatrixFromArray(q)
    )
    return float(result)


def ms_jensendifference(p, q):
    """
    Computes the Jensen difference divergence.

    Args:
        p: First probability distribution.
        q: Second probability distribution.

    Returns:
        float: The Jensen difference divergence.
    """
    result = jpype.JPackage('jline').api.ms.Ms_jensendifferenceKt.ms_jensendifference(
        jlineMatrixFromArray(p), jlineMatrixFromArray(q)
    )
    return float(result)


def ms_jensenshannon(p, q):
    """
    Computes the Jensen-Shannon divergence.

    Args:
        p: First probability distribution.
        q: Second probability distribution.

    Returns:
        float: The Jensen-Shannon divergence.
    """
    result = jpype.JPackage('jline').api.ms.Ms_jensenshannonKt.ms_jensenshannon(
        jlineMatrixFromArray(p), jlineMatrixFromArray(q)
    )
    return float(result)


def ms_jointentropy(p, q):
    """
    Computes the joint entropy H(P,Q).

    Args:
        p: First probability distribution.
        q: Second probability distribution.

    Returns:
        float: The joint entropy.
    """
    result = jpype.JPackage('jline').api.ms.Ms_jointentropyKt.ms_jointentropy(
        jlineMatrixFromArray(p), jlineMatrixFromArray(q)
    )
    return float(result)


def ms_kdivergence(p, q):
    """
    Computes the K-divergence between distributions.

    Args:
        p: First probability distribution.
        q: Second probability distribution.

    Returns:
        float: The K-divergence.
    """
    result = jpype.JPackage('jline').api.ms.Ms_kdivergenceKt.ms_kdivergence(
        jlineMatrixFromArray(p), jlineMatrixFromArray(q)
    )
    return float(result)


def ms_kolmogorov_smirnov(data1, data2):
    """
    Performs Kolmogorov-Smirnov test for comparing two samples.

    Args:
        data1: First data sample.
        data2: Second data sample.

    Returns:
        dict: Test statistic and p-value.
    """
    result = jpype.JPackage('jline').api.ms.Ms_kolmogorov_smirnovKt.ms_kolmogorov_smirnov(
        jlineMatrixFromArray(data1), jlineMatrixFromArray(data2)
    )
    return {
        'statistic': float(result.getStatistic()),
        'pvalue': float(result.getPvalue())
    }


def ms_kuiper(data1, data2):
    """
    Performs Kuiper test for comparing two samples.

    Args:
        data1: First data sample.
        data2: Second data sample.

    Returns:
        dict: Test statistic and p-value.
    """
    result = jpype.JPackage('jline').api.ms.Ms_kuiperKt.ms_kuiper(
        jlineMatrixFromArray(data1), jlineMatrixFromArray(data2)
    )
    return {
        'statistic': float(result.getStatistic()),
        'pvalue': float(result.getPvalue())
    }


def ms_kulczynskid(p, q):
    """
    Computes the Kulczynski D similarity.

    Args:
        p: First vector.
        q: Second vector.

    Returns:
        float: The Kulczynski D similarity.
    """
    result = jpype.JPackage('jline').api.ms.Ms_kulczynskidKt.ms_kulczynskid(
        jlineMatrixFromArray(p), jlineMatrixFromArray(q)
    )
    return float(result)


def ms_kulczynskis(p, q):
    """
    Computes the Kulczynski S similarity.

    Args:
        p: First vector.
        q: Second vector.

    Returns:
        float: The Kulczynski S similarity.
    """
    result = jpype.JPackage('jline').api.ms.Ms_kulczynskisKt.ms_kulczynskis(
        jlineMatrixFromArray(p), jlineMatrixFromArray(q)
    )
    return float(result)


def ms_kullbackleibler(p, q):
    """
    Computes the Kullback-Leibler divergence KL(P||Q).

    Args:
        p: First probability distribution.
        q: Second probability distribution.

    Returns:
        float: The KL divergence.
    """
    result = jpype.JPackage('jline').api.ms.Ms_kullbackleiblerKt.ms_kullbackleibler(
        jlineMatrixFromArray(p), jlineMatrixFromArray(q)
    )
    return float(result)


def ms_kumarhassebrook(p, q):
    """
    Computes the Kumar-Hassebrook similarity.

    Args:
        p: First vector.
        q: Second vector.

    Returns:
        float: The Kumar-Hassebrook similarity.
    """
    result = jpype.JPackage('jline').api.ms.Ms_kumarhassebrookKt.ms_kumarhassebrook(
        jlineMatrixFromArray(p), jlineMatrixFromArray(q)
    )
    return float(result)


def ms_kumarjohnson(p, q):
    """
    Computes the Kumar-Johnson divergence.

    Args:
        p: First probability distribution.
        q: Second probability distribution.

    Returns:
        float: The Kumar-Johnson divergence.
    """
    result = jpype.JPackage('jline').api.ms.Ms_kumarjohnsonKt.ms_kumarjohnson(
        jlineMatrixFromArray(p), jlineMatrixFromArray(q)
    )
    return float(result)


def ms_lorentzian(p, q):
    """
    Computes the Lorentzian distance.

    Args:
        p: First vector.
        q: Second vector.

    Returns:
        float: The Lorentzian distance.
    """
    result = jpype.JPackage('jline').api.ms.Ms_lorentzianKt.ms_lorentzian(
        jlineMatrixFromArray(p), jlineMatrixFromArray(q)
    )
    return float(result)


def ms_matusita(p, q):
    """
    Computes the Matusita distance.

    Args:
        p: First probability distribution.
        q: Second probability distribution.

    Returns:
        float: The Matusita distance.
    """
    result = jpype.JPackage('jline').api.ms.Ms_matusitaKt.ms_matusita(
        jlineMatrixFromArray(p), jlineMatrixFromArray(q)
    )
    return float(result)


def ms_minkowski(p, q, order=2):
    """
    Computes the Minkowski distance of given order.

    Args:
        p: First vector.
        q: Second vector.
        order: Order of the Minkowski distance (default 2).

    Returns:
        float: The Minkowski distance.
    """
    result = jpype.JPackage('jline').api.ms.Ms_minkowskiKt.ms_minkowski(
        jlineMatrixFromArray(p), jlineMatrixFromArray(q), float(order)
    )
    return float(result)


def ms_motyka(p, q):
    """
    Computes the Motyka similarity.

    Args:
        p: First vector.
        q: Second vector.

    Returns:
        float: The Motyka similarity.
    """
    result = jpype.JPackage('jline').api.ms.Ms_motykaKt.ms_motyka(
        jlineMatrixFromArray(p), jlineMatrixFromArray(q)
    )
    return float(result)


def ms_mutinfo(p, q):
    """
    Computes the mutual information I(P;Q).

    Args:
        p: First probability distribution.
        q: Second probability distribution.

    Returns:
        float: The mutual information.
    """
    result = jpype.JPackage('jline').api.ms.Ms_mutinfoKt.ms_mutinfo(
        jlineMatrixFromArray(p), jlineMatrixFromArray(q)
    )
    return float(result)


def ms_neymanchisquared(p, q):
    """
    Computes the Neyman chi-squared statistic.

    Args:
        p: First probability distribution.
        q: Second probability distribution.

    Returns:
        float: The Neyman chi-squared statistic.
    """
    result = jpype.JPackage('jline').api.ms.Ms_neymanchisquaredKt.ms_neymanchisquared(
        jlineMatrixFromArray(p), jlineMatrixFromArray(q)
    )
    return float(result)


def ms_nmi(p, q):
    """
    Computes the normalized mutual information.

    Args:
        p: First probability distribution.
        q: Second probability distribution.

    Returns:
        float: The normalized mutual information.
    """
    result = jpype.JPackage('jline').api.ms.Ms_nmiKt.ms_nmi(
        jlineMatrixFromArray(p), jlineMatrixFromArray(q)
    )
    return float(result)


def ms_nvi(p, q):
    """
    Computes the normalized variation of information.

    Args:
        p: First probability distribution.
        q: Second probability distribution.

    Returns:
        float: The normalized variation of information.
    """
    result = jpype.JPackage('jline').api.ms.Ms_nviKt.ms_nvi(
        jlineMatrixFromArray(p), jlineMatrixFromArray(q)
    )
    return float(result)


def ms_pearsonchisquared(p, q):
    """
    Computes the Pearson chi-squared statistic.

    Args:
        p: First probability distribution.
        q: Second probability distribution.

    Returns:
        float: The Pearson chi-squared statistic.
    """
    result = jpype.JPackage('jline').api.ms.Ms_pearsonchisquaredKt.ms_pearsonchisquared(
        jlineMatrixFromArray(p), jlineMatrixFromArray(q)
    )
    return float(result)


def ms_probsymmchisquared(p, q):
    """
    Computes the probability symmetric chi-squared divergence.

    Args:
        p: First probability distribution.
        q: Second probability distribution.

    Returns:
        float: The probability symmetric chi-squared divergence.
    """
    result = jpype.JPackage('jline').api.ms.Ms_probsymmchisquaredKt.ms_probsymmchisquared(
        jlineMatrixFromArray(p), jlineMatrixFromArray(q)
    )
    return float(result)


def ms_product(p, q):
    """
    Computes the inner product similarity.

    Args:
        p: First vector.
        q: Second vector.

    Returns:
        float: The inner product.
    """
    result = jpype.JPackage('jline').api.ms.Ms_productKt.ms_product(
        jlineMatrixFromArray(p), jlineMatrixFromArray(q)
    )
    return float(result)


def ms_relatentropy(p, q):
    """
    Computes the relative entropy (KL divergence).

    Args:
        p: First probability distribution.
        q: Second probability distribution.

    Returns:
        float: The relative entropy.
    """
    result = jpype.JPackage('jline').api.ms.Ms_relatentropyKt.ms_relatentropy(
        jlineMatrixFromArray(p), jlineMatrixFromArray(q)
    )
    return float(result)


def ms_ruzicka(p, q):
    """
    Computes the Ruzicka similarity.

    Args:
        p: First vector.
        q: Second vector.

    Returns:
        float: The Ruzicka similarity.
    """
    result = jpype.JPackage('jline').api.ms.Ms_ruzickaKt.ms_ruzicka(
        jlineMatrixFromArray(p), jlineMatrixFromArray(q)
    )
    return float(result)


def ms_soergel(p, q):
    """
    Computes the Soergel distance.

    Args:
        p: First vector.
        q: Second vector.

    Returns:
        float: The Soergel distance.
    """
    result = jpype.JPackage('jline').api.ms.Ms_soergelKt.ms_soergel(
        jlineMatrixFromArray(p), jlineMatrixFromArray(q)
    )
    return float(result)


def ms_sorensen(p, q):
    """
    Computes the Sorensen (Dice) coefficient.

    Args:
        p: First vector.
        q: Second vector.

    Returns:
        float: The Sorensen coefficient.
    """
    result = jpype.JPackage('jline').api.ms.Ms_sorensenKt.ms_sorensen(
        jlineMatrixFromArray(p), jlineMatrixFromArray(q)
    )
    return float(result)


def ms_squaredchord(p, q):
    """
    Computes the squared chord distance.

    Args:
        p: First probability distribution.
        q: Second probability distribution.

    Returns:
        float: The squared chord distance.
    """
    result = jpype.JPackage('jline').api.ms.Ms_squaredchordKt.ms_squaredchord(
        jlineMatrixFromArray(p), jlineMatrixFromArray(q)
    )
    return float(result)


def ms_squaredeuclidean(p, q):
    """
    Computes the squared Euclidean distance.

    Args:
        p: First vector.
        q: Second vector.

    Returns:
        float: The squared Euclidean distance.
    """
    result = jpype.JPackage('jline').api.ms.Ms_squaredeuclideanKt.ms_squaredeuclidean(
        jlineMatrixFromArray(p), jlineMatrixFromArray(q)
    )
    return float(result)


def ms_taneja(p, q):
    """
    Computes the Taneja divergence.

    Args:
        p: First probability distribution.
        q: Second probability distribution.

    Returns:
        float: The Taneja divergence.
    """
    result = jpype.JPackage('jline').api.ms.Ms_tanejaKt.ms_taneja(
        jlineMatrixFromArray(p), jlineMatrixFromArray(q)
    )
    return float(result)


def ms_tanimoto(p, q):
    """
    Computes the Tanimoto (extended Jaccard) coefficient.

    Args:
        p: First vector.
        q: Second vector.

    Returns:
        float: The Tanimoto coefficient.
    """
    result = jpype.JPackage('jline').api.ms.Ms_tanimotoKt.ms_tanimoto(
        jlineMatrixFromArray(p), jlineMatrixFromArray(q)
    )
    return float(result)


def ms_topsoe(p, q):
    """
    Computes the Topsoe divergence.

    Args:
        p: First probability distribution.
        q: Second probability distribution.

    Returns:
        float: The Topsoe divergence.
    """
    result = jpype.JPackage('jline').api.ms.Ms_topsoeKt.ms_topsoe(
        jlineMatrixFromArray(p), jlineMatrixFromArray(q)
    )
    return float(result)


def ms_wasserstein(p, q):
    """
    Computes the Wasserstein (earth mover's) distance.

    Args:
        p: First probability distribution.
        q: Second probability distribution.

    Returns:
        float: The Wasserstein distance.
    """
    result = jpype.JPackage('jline').api.ms.Ms_wassersteinKt.ms_wasserstein(
        jlineMatrixFromArray(p), jlineMatrixFromArray(q)
    )
    return float(result)


def ms_wavehegdes(p, q):
    """
    Computes the Wave Hedges distance.

    Args:
        p: First vector.
        q: Second vector.

    Returns:
        float: The Wave Hedges distance.
    """
    result = jpype.JPackage('jline').api.ms.Ms_wavehegdesKt.ms_wavehegdes(
        jlineMatrixFromArray(p), jlineMatrixFromArray(q)
    )
    return float(result)