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

from .sens import (
    pfqn_sens,
    PfqnSens,
)

from .sens_mva import (
    pfqn_sens_mva,
    PfqnSensMva,
)

from .sens_mvaldmx import (
    pfqn_sens_ldmx_ec,
    pfqn_sens_mvaldmx,
    PfqnSensMvaldmx,
)

from .sens_mom import (
    pfqn_sens_mom,
    PfqnSensMom,
)

from .sens_respt import (
    pfqn_sens_respt,
    PfqnSensRespt,
)

from .sens_linearizer import (
    pfqn_sens_linearizer,
    PfqnSensLinearizer,
)

from .momlin import (
    pfqn_momlin,
    PfqnMomlin,
)

from .qlen_moments import pfqn_qlen_joint_moments
from .nc import (
    pfqn_ca,
    pfqn_is,
    pfqn_nc,
    pfqn_nc_resolved_method,
    pfqn_panacea,
    pfqn_propfair,
    pfqn_ls,
    pfqn_clw,
    pfqn_clw_lld,
)

from .linearizer import (
    pfqn_linearizer,
    pfqn_gflinearizer,
    pfqn_egflinearizer,
    SchedStrategy,
)

from .cftp import (
    pfqn_cftp,
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
    pfqn_mwrbb,
)

from .bound_hierarchies import (
    pfqn_pbh,
    pfqn_pbk,
    pfqn_bjbk,
    pfqn_cbh,
    pfqn_mcub,
    pfqn_ssd,
    pfqn_sib,
    pfqn_ldbcmp,
)

from .marie import (
    pfqn_marie,
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
    pfqn_ld_is,
    pfqn_gld,
    pfqn_gldsingle,
    pfqn_mushift,
    pfqn_comomrm_ld,
    pfqn_fnc,
    PfqnNcResult,
    PfqnComomrmLdResult,
    PfqnFncResult,
)

from .ncldmx import (
    pfqn_ncldmx,
    PfqnNcldmxResult,
)


from .conv import (
    pfqn_conv,
    solver_nc_conv,
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
    pfqn_procomom,
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

from .dac import (
    pfqn_dac,
)

from .recal import (
    pfqn_recal,
)

from .mvac import (
    pfqn_mvac,
)

from .mvacld import (
    pfqn_mvacld,
)

from .mvaldmx import (
    pfqn_mvaldmx,
    pfqn_ldmx_ec,
    pfqn_mvaldms,
)

from .linearizerms import (
    pfqn_linearizerms,
    pfqn_conwayms,
)
# pfqn_linearizermx lives in its own module (the MATLAB-faithful version with the
# QN0 warm-start argument); linearizerms.py no longer defines a second copy.
from .linearizermx import pfqn_linearizermx

from .ljd import (
    ljd_linearize,
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
    'pfqn_is',
    'pfqn_nc',
    'pfqn_qlen_joint_moments',
    'pfqn_nc_resolved_method',
    'pfqn_panacea',
    'pfqn_propfair',
    'pfqn_ls',
    'pfqn_clw',
    'pfqn_clw_lld',
    # Linearizer algorithms
    'pfqn_linearizer',
    'pfqn_gflinearizer',
    'pfqn_egflinearizer',
    'pfqn_momlin',
    'PfqnMomlin',
    'SchedStrategy',
    # Sensitivity analysis
    'pfqn_sens',
    'PfqnSens',
    # Exact queue-length moment recursions
    'pfqn_sens_mva',
    'PfqnSensMva',
    'pfqn_sens_ldmx_ec',
    'pfqn_sens_mvaldmx',
    'PfqnSensMvaldmx',
    # Strelen higher-moment analysis (up to order three)
    'pfqn_sens_mom',
    'PfqnSensMom',
    'pfqn_sens_respt',
    'PfqnSensRespt',
    'pfqn_sens_linearizer',
    'PfqnSensLinearizer',
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
    'pfqn_mwrbb',
    # Hierarchical / multiserver / LD bound methods
    'pfqn_pbh',
    'pfqn_pbk',
    'pfqn_bjbk',
    'pfqn_cbh',
    'pfqn_mcub',
    'pfqn_ssd',
    'pfqn_sib',
    'pfqn_ldbcmp',
    'pfqn_marie',
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
    'pfqn_ld_is',
    'pfqn_ncldmx',
    'PfqnNcldmxResult',
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
    'pfqn_procomom',
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
    'pfqn_dac',
    'pfqn_recal',
    # MVAC (mean value analysis by chain)
    'pfqn_mvac',
    'pfqn_mvacld',
    # Load-dependent mixed MVA
    'pfqn_mvaldmx',
    'pfqn_ldmx_ec',
    'pfqn_mvaldms',
    # Multi-server and mixed linearizers
    'pfqn_linearizerms',
    'pfqn_linearizermx',
    'pfqn_conwayms',
    # LJD indexing
    'ljd_linearize',
    'infradius_h',
    'infradius_hnorm',
    # Knessl-Tier expansion
    'pfqn_kt',
    # Akyildiz-Bolch AMVA
    'pfqn_ab_amva',
    'pfqn_ab_core',
    'AbAmvaResult',
    # Reduction Heuristic
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
    # Coupling from the past (exact sampler)
    'pfqn_cftp',
]
