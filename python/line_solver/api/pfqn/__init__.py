"""
Product-form queueing network (PFQN) algorithms.

Native Python implementations of analytical algorithms for product-form
queueing networks, including Mean Value Analysis (MVA), normalizing constant
methods, and various approximation techniques.

Key algorithms:
    pfqn_mva: Standard Mean Value Analysis
    pfqn_ca: Convolution Algorithm
    pfqn_nc: Normalizing Constant methods
    pfqn_bs: Balanced System analysis
    pfqn_aql: Approximate queue lengths
"""

from .mva import (
    pfqn_mva,
    pfqn_mva_single_class,
    pfqn_bs,
    pfqn_aql,
    pfqn_sqni,
    pfqn_qd,
    pfqn_qdlin,
    pfqn_qli,
    pfqn_fli,
    pfqn_bsfcfs,
    pfqn_joint,
)

from .nc import (
    pfqn_ca,
    pfqn_nc,
    pfqn_panacea,
    pfqn_propfair,
    pfqn_ls,
)

from .linearizer import (
    pfqn_linearizer,
    pfqn_gflinearizer,
    pfqn_egflinearizer,
    SchedStrategy,
)

from .mvald import (
    pfqn_mvald,
    pfqn_mvams,
)

from .mixed import (
    pfqn_mvamx,
)

from .bounds import (
    pfqn_xzabalow,
    pfqn_xzabaup,
    pfqn_qzgblow,
    pfqn_qzgbup,
    pfqn_xzgsblow,
    pfqn_xzgsbup,
)

from .asymptotic import (
    pfqn_le,
    pfqn_cub,
    pfqn_mci,
    pfqn_grnmol,
    pfqn_le_fpi,
    pfqn_le_fpiZ,
    pfqn_le_hessian,
    pfqn_le_hessianZ,
)

from .ncld import (
    pfqn_ncld,
    pfqn_gld,
    pfqn_gldsingle,
    pfqn_mushift,
    pfqn_comomrm_ld,
    pfqn_fnc,
    PfqnNcResult,
    PfqnComomrmLdResult,
    PfqnFncResult,
)


from .replicas import (
    pfqn_unique,
    pfqn_expand,
    pfqn_combine_mi,
    PfqnUniqueResult,
)

from .utils import (
    pfqn_lldfun,
    pfqn_mu_ms,
    pfqn_nc_sanitize,
    pfqn_cdfun,
    pfqn_ljdfun,
    factln,
    factln_vec,
    softmin,
    oner,
    multichoose,
    matchrow,
)

from .comom import (
    pfqn_comom,
    pfqn_comomrm,
    pfqn_comomrm_orig,
    pfqn_comomrm_ms,
    pfqn_procomom2,
    ComomResult,
)

from .quadrature import (
    pfqn_mmint2,
    pfqn_mmint2_gausslegendre,
    pfqn_mmint2_gausslaguerre,
    pfqn_mmsample2,
    logsumexp,
)

from .schmidt import (
    pfqn_schmidt,
    pfqn_schmidt_ext,
    SchmidtResult,
    pprod,
    hashpop,
)

from .recal import (
    pfqn_recal,
)

from .mvaldmx import (
    pfqn_mvaldmx,
    pfqn_mvaldmx_ec,
    pfqn_mvaldms,
)

from .linearizerms import (
    pfqn_linearizerms,
    pfqn_linearizermx,
    pfqn_conwayms,
)

from .ljd import (
    ljd_linearize,
    ljd_delinearize,
    ljcd_interpolate,
    infradius_h,
    infradius_hnorm,
)

from .kt import (
    pfqn_kt,
)

from .ab_amva import (
    pfqn_ab_amva,
    pfqn_ab_core,
    AbAmvaResult,
)

from .rd import (
    pfqn_rd,
    RdOptions,
    RdResult,
)

from .laplace import (
    pfqn_nrl,
    pfqn_nrp,
    pfqn_lap,
    laplaceapprox,
    num_hess,
)

from .stdf import (
    pfqn_stdf,
    pfqn_stdf_heur,
)

__all__ = [
    # MVA algorithms
    'pfqn_mva',
    'pfqn_mva_single_class',
    'pfqn_bs',
    'pfqn_aql',
    'pfqn_sqni',
    'pfqn_qd',
    'pfqn_qdlin',
    'pfqn_qli',
    'pfqn_fli',
    'pfqn_bsfcfs',
    'pfqn_joint',
    # Normalizing constant algorithms
    'pfqn_ca',
    'pfqn_nc',
    'pfqn_panacea',
    'pfqn_propfair',
    'pfqn_ls',
    # Linearizer algorithms
    'pfqn_linearizer',
    'pfqn_gflinearizer',
    'pfqn_egflinearizer',
    'SchedStrategy',
    # Load-dependent MVA
    'pfqn_mvald',
    'pfqn_mvams',
    # Mixed MVA
    'pfqn_mvamx',
    # Bounds
    'pfqn_xzabalow',
    'pfqn_xzabaup',
    'pfqn_qzgblow',
    'pfqn_qzgbup',
    'pfqn_xzgsblow',
    'pfqn_xzgsbup',
    # Asymptotic methods
    'pfqn_le',
    'pfqn_cub',
    'pfqn_mci',
    'pfqn_grnmol',
    'pfqn_le_fpi',
    'pfqn_le_fpiZ',
    'pfqn_le_hessian',
    'pfqn_le_hessianZ',
    # Load-dependent NC algorithms
    'pfqn_ncld',
    'pfqn_gld',
    'pfqn_gldsingle',
    'pfqn_mushift',
    'pfqn_comomrm_ld',
    'pfqn_fnc',
    'PfqnNcResult',
    'PfqnComomrmLdResult',
    'PfqnFncResult',
    # Replica consolidation
    'pfqn_unique',
    'pfqn_expand',
    'pfqn_combine_mi',
    'PfqnUniqueResult',
    # Utility functions
    'pfqn_lldfun',
    'pfqn_mu_ms',
    'pfqn_nc_sanitize',
    'pfqn_cdfun',
    'pfqn_ljdfun',
    'factln',
    'factln_vec',
    'softmin',
    'oner',
    'multichoose',
    'matchrow',
    # COMOM methods
    'pfqn_comom',
    'pfqn_comomrm',
    'pfqn_comomrm_orig',
    'pfqn_comomrm_ms',
    'pfqn_procomom2',
    'ComomResult',
    # Quadrature methods
    'pfqn_mmint2',
    'pfqn_mmint2_gausslegendre',
    'pfqn_mmint2_gausslaguerre',
    'pfqn_mmsample2',
    'logsumexp',
    # Schmidt's exact MVA
    'pfqn_schmidt',
    'pfqn_schmidt_ext',
    'SchmidtResult',
    'pprod',
    'hashpop',
    # RECAL method
    'pfqn_recal',
    # Load-dependent mixed MVA
    'pfqn_mvaldmx',
    'pfqn_mvaldmx_ec',
    'pfqn_mvaldms',
    # Multi-server and mixed linearizers
    'pfqn_linearizerms',
    'pfqn_linearizermx',
    'pfqn_conwayms',
    # LJD indexing
    'ljd_linearize',
    'ljd_delinearize',
    'ljcd_interpolate',
    'infradius_h',
    'infradius_hnorm',
    # Knessl-Tier expansion
    'pfqn_kt',
    # Akyildiz-Bolch AMVA
    'pfqn_ab_amva',
    'pfqn_ab_core',
    'AbAmvaResult',
    # Reduced Decomposition
    'pfqn_rd',
    'RdOptions',
    'RdResult',
    # Laplace approximation methods
    'pfqn_nrl',
    'pfqn_nrp',
    'pfqn_lap',
    'laplaceapprox',
    'num_hess',
    # Sojourn time distribution
    'pfqn_stdf',
    'pfqn_stdf_heur',
]
