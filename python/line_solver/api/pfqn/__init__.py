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
    pfqn_mva_ilock,
    pfqn_mva_single_class,
    pfqn_bs,
    pfqn_aql,
    pfqn_sqni,
    pfqn_qli,
    pfqn_fli,
    pfqn_joint,
    pfqn_jointmarg,
)

from .cntol import (
    pfqn_cntol,
    is_cntol,
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

from .busyp import (
    pfqn_busyp,
)

from .busyp_multiclass import (
    pfqn_busyp_multiclass,
)

from .busyp_clw import (
    pfqn_busyp_clw,
)

from .respt_ps import (
    pfqn_respt_ps_moments,
    PfqnResptPsMoments,
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
    pfqn_perm,
)

from .qsa import (
    pfqn_qsa,
)

from .lcp import (
    pfqn_lcp,
)

from .chow import (
    pfqn_chow,
)

from .pam import (
    pfqn_pam,
)

from .clust import (
    pfqn_clust,
)

from .dmlin import (
    pfqn_dmlin,
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
    pfqn_mvams_ilock,
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
    pfqn_harel_bounds,
    pfqn_harel_lb,
    pfqn_harel_ub,
)

from .bound_hierarchies import (
    pfqn_pbh,
    pfqn_looping,
    pfqn_pbk,
    pfqn_bjbk,
    pfqn_cbh,
    pfqn_mcub,
    pfqn_ssd,
    pfqn_sib,
    pfqn_ldbcmp,
    pfqn_scb,
    pfqn_scbgap,
    pfqn_usumbound,
    pfqn_minclasses,
)

from .explicit import (
    pfqn_explicit,
)
from .explicit_ld import (
    pfqn_explicit_ld,
)
from .rgf import (
    pfqn_rgf,
    pfqn_rgfmc,
)

from .gerasimov import (
    pfqn_gerasimov,
)

from .nonintegral import (
    pfqn_dnc,
    pfqn_nintmva,
)

from .interval import (
    pfqn_mva_interval,
    PfqnMvaIntervalResult,
)

from .tay import (
    pfqn_tay,
)

from .scat import (
    pfqn_scat,
)

from .robustness import (
    pfqn_hst,
)

from .sjn import (
    pfqn_mvasjn,
    pfqn_amvasjn,
    SjnOptions,
)
from .marie import (
    pfqn_marie,
)

from .asymptotic import (
    pfqn_le,
    pfqn_lekt,
    pfqn_lekt_route,
    pfqn_ble,
    pfqn_aghq,
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
    pfqn_panaceald,
    pfqn_ld_is,
    pfqn_gld,
    pfqn_gldsingle,
    pfqn_lldsingle,
    pfqn_lld,
    pfqn_xia,
    pfqn_mushift,
    pfqn_comomrm_ld,
    pfqn_fnc,
    PfqnNcResult,
    PfqnComomrmLdResult,
    PfqnFncResult,
)

from .oi import (
    pfqn_ncoi,
    pfqn_ncjd,
    pfqn_oi_fnc,
    pfqn_oi_insvc,
    pfqn_mvaoi,
    pfqn_mvajd,
    pfqn_mvaoi_marg,
)

from .clwjd import (
    pfqn_clwoi,
    pfqn_clwjd,
)

from .pas import (
    pfqn_pas_is,
    pfqn_pas_nc,
    pas_placement,
    pas_swap2order,
)

from .ncldmx import (
    pfqn_ncldmx,
    PfqnNcldmxResult,
)


from .conv import (
    pfqn_conv,
    solver_nc_conv,
)

from .manjunath import pfqn_manjunath

from .replicas import (
    pfqn_unique,
    pfqn_expand,
    pfqn_combine_mi,
    PfqnUniqueResult,
)

from .sdr import (
    pfqn_sdrcoeff,
    pfqn_sdrprob,
    pfqn_sdr,
    pfqn_sdrvisits,
    pfqn_sdrmva,
)

from .qdamva import pfqn_qdamva
from .qdlin import pfqn_qdlin
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
    multichoosecon,
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
    pfqn_bkt,
)

from .mcmc import (
    pfqn_mcmc,
    PfqnMcmcResult,
)
from .bk import (
    pfqn_bk,
    pfqn_bkue,
    pfqn_bklc,
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

from .nre import (
    pfqn_nre,
    pfqn_nre_full,
    PfqnNreResult,
)

from .laplace import (
    pfqn_nrl,
    pfqn_nrp,
    pfqn_lap,
    laplaceapprox,
    num_hess,
)

from .cyclet import (
    pfqn_cyclet_ofree,
)

from .stdf import (
    pfqn_stdf,
    pfqn_stdf_heur,
)

__all__ = [
    # MVA algorithms
    'pfqn_mva',
    'pfqn_mva_ilock',
    'pfqn_mva_single_class',
    'pfqn_bs',
    'pfqn_cntol',
    'is_cntol',
    'pfqn_aql',
    'pfqn_sqni',
    'pfqn_qdlin',
    'pfqn_qli',
    'pfqn_fli',
    'pfqn_joint',
    'pfqn_jointmarg',
    # Normalizing constant algorithms
    'pfqn_ca',
    'pfqn_manjunath',
    'pfqn_is',
    'pfqn_nc',
    'pfqn_qlen_joint_moments',
    'pfqn_nc_resolved_method',
    'pfqn_panacea',
    'pfqn_panaceald',
    'pfqn_propfair',
    'pfqn_ls',
    'pfqn_mcmc',
    'PfqnMcmcResult',
    'pfqn_clw',
    'pfqn_clw_lld',
    'pfqn_perm',
    # Linearizer algorithms
    'pfqn_qsa',
    'pfqn_lcp',
    'pfqn_chow',
    'pfqn_looping',
    'pfqn_pam',
    'pfqn_clust',
    'pfqn_dmlin',
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
    'pfqn_busyp',
    'pfqn_busyp_multiclass',
    'pfqn_busyp_clw',
    'PfqnSensRespt',
    # Processor-sharing sojourn-time moments (Mitra-Morrison)
    'pfqn_respt_ps_moments',
    'PfqnResptPsMoments',
    'pfqn_sens_linearizer',
    'PfqnSensLinearizer',
    # Load-dependent MVA
    'pfqn_mvald',
    'pfqn_mvams',
    'pfqn_mvams_ilock',
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
    'pfqn_harel_bounds',
    'pfqn_harel_lb',
    'pfqn_harel_ub',
    # Hierarchical / multiserver / LD bound methods
    'pfqn_pbh',
    'pfqn_looping',
    'pfqn_pbk',
    'pfqn_bjbk',
    'pfqn_cbh',
    'pfqn_mcub',
    'pfqn_ssd',
    'pfqn_sib',
    'pfqn_ldbcmp',
    'pfqn_scb',
    'pfqn_scbgap',
    'pfqn_usumbound',
    'pfqn_minclasses',
    'pfqn_explicit',
    'pfqn_explicit_ld',
    'pfqn_rgf',
    'pfqn_rgfmc',
    'pfqn_gerasimov',
    'pfqn_dnc',
    'pfqn_nintmva',
    'pfqn_mva_interval',
    'PfqnMvaIntervalResult',
    'pfqn_tay',
    'pfqn_scat',
    'pfqn_hst',
    'pfqn_marie',
    # Asymptotic methods
    'pfqn_le',
    'pfqn_lekt',
    'pfqn_lekt_route',
    'pfqn_ble',
    'pfqn_aghq',
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
    'pfqn_lldsingle',
    'pfqn_lld',
    'pfqn_xia',
    'pfqn_mushift',
    'pfqn_comomrm_ld',
    'pfqn_fnc',
    'pfqn_ncoi',
    'pfqn_ncjd',
    'pfqn_clwoi',
    'pfqn_clwjd',
    'pfqn_oi_fnc',
    'pfqn_oi_insvc',
    'pfqn_pas_is',
    'pfqn_pas_nc',
    'pas_placement',
    'pas_swap2order',
    'pfqn_mvaoi',
    'pfqn_mvajd',
    'pfqn_mvaoi_marg',
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
    'pfqn_qdamva',
    'pfqn_mu_ms',
    'pfqn_nc_sanitize',
    'pfqn_cdfun',
    'pfqn_sdrcoeff',
    'pfqn_sdrprob',
    'pfqn_sdr',
    'pfqn_sdrvisits',
    'pfqn_sdrmva',
    'factln',
    'factln_vec',
    'softmin',
    'oner',
    'multichoose',
    'multichoosecon',
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
    # Knessl-Tier expansion, and its Stirling-remainder correction
    'pfqn_kt',
    'pfqn_bkt',
    # Birman-Kogan saddle point, uniform expansion and load concealment
    'pfqn_bk',
    'pfqn_bkue',
    'pfqn_bklc',
    # Akyildiz-Bolch AMVA
    'pfqn_ab_amva',
    'pfqn_ab_core',
    'AbAmvaResult',
    # Reduction Heuristic
    'pfqn_rd',
    'RdOptions',
    'RdResult',
    # Laplace approximation methods
    'pfqn_nre',
    'pfqn_nre_full',
    'PfqnNreResult',
    'pfqn_nrl',
    'pfqn_nrp',
    'pfqn_lap',
    'laplaceapprox',
    'num_hess',
    # Sojourn time distribution
    'pfqn_cyclet_ofree',
    'pfqn_stdf',
    'pfqn_stdf_heur',
    # Coupling from the past (exact sampler)
    'pfqn_cftp',
    # Shortest-job-next stations (Kant 1992)
    'pfqn_mvasjn',
    'pfqn_amvasjn',
    'SjnOptions',
]
