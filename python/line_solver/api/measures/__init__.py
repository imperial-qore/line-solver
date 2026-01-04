"""
Statistical Distance and Divergence Measures.

Native Python implementations for comparing probability distributions
using various distance metrics and divergence measures.

Includes:
    - Information-theoretic: Entropy, KL divergence, Jensen-Shannon
    - Minkowski family: Euclidean, Manhattan, Chebyshev
    - Inner product: Cosine, Jaccard, Dice
    - Fidelity: Bhattacharyya, Hellinger
    - Chi-squared: Pearson, Neyman
    - Empirical tests: Kolmogorov-Smirnov, Wasserstein
"""

from .distances import (
    # Entropy measures
    ms_entropy,
    ms_jointentropy,
    ms_condentropy,
    ms_mutinfo,
    ms_nmi,
    ms_nvi,
    # Information-theoretic divergences
    ms_kullbackleibler,
    ms_relatentropy,
    ms_jensenshannon,
    ms_jeffreys,
    ms_kdivergence,
    ms_topsoe,
    ms_jensendifference,
    ms_taneja,
    # Minkowski family
    ms_minkowski,
    ms_euclidean,
    ms_squaredeuclidean,
    ms_cityblock,
    ms_chebyshev,
    ms_sorensen,
    ms_gower,
    ms_soergel,
    ms_lorentzian,
    ms_canberra,
    ms_avgl1linfty,
    # Inner product family
    ms_product,
    ms_cosine,
    ms_jaccard,
    ms_dice,
    ms_kumarhassebrook,
    ms_harmonicmean,
    # Fidelity family
    ms_fidelity,
    ms_bhattacharyya,
    ms_hellinger,
    ms_matusita,
    ms_squaredchord,
    # Chi-squared family
    ms_pearsonchisquared,
    ms_neymanchisquared,
    ms_squaredchisquared,
    ms_probsymmchisquared,
    ms_divergence,
    ms_additivesymmetricchisquared,
    ms_clark,
    # Intersection family
    ms_intersection,
    ms_czekanowski,
    ms_motyka,
    ms_kulczynskis,
    ms_kulczynskid,
    ms_ruzicka,
    ms_tanimoto,
    ms_wavehegdes,
    ms_kumarjohnson,
    # Empirical tests
    ms_kolmogorov_smirnov,
    ms_kuiper,
    ms_cramer_von_mises,
    ms_anderson_darling,
    ms_wasserstein,
)

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
