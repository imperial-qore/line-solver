
"""
LINE Solver API Functions

This module provides direct access to low-level LINE algorithms and analytical
methods. These functions offer fine-grained control over specific computations
and are primarily intended for advanced users and researchers.

The API is organized into several domains:

- **cache**: Cache modeling algorithms (miss ratios, replacement policies)
- **ctmc/dtmc**: Continuous and discrete-time Markov chain analysis  
- **pfqn**: Product-form queueing network algorithms (MVA, NC, etc.)
- **mam**: Matrix-analytic methods for MAP/PH distributions
- **npfqn**: Non-product-form queueing network approximations
- **qsys**: Single queueing system analysis (M/M/1, M/G/1, etc.)
- **sn**: Stochastic network utilities and analysis
- **polling**: Polling system analysis
- **lossn/lsn**: Loss network and layered stochastic network analysis
- **trace**: Trace analysis utilities

Most users should use the high-level solver classes (SolverMVA, SolverJMT, etc.)
rather than calling these API functions directly.
"""

from .cache import (
    cache_mva, cache_prob_asy, cache_gamma_lp, cache_spm, cache_xi_fp,
    cache_miss_spm, cache_prob_erec, cache_prob_fpi, cache_prob_spm,
    cache_erec, cache_t_hlru, cache_t_lrum, cache_t_lrum_map, cache_ttl_hlru,
    cache_ttl_lrua, cache_ttl_lrum, cache_ttl_lrum_map, cache_ttl_tree, cache_xi_bvh,
    cache_miss, cache_mva_miss, cache_miss_asy, cache_erec_aux, cache_par,
    cache_t_hlru_aux, cache_t_lrum_aux,
    # Deprecated aliases
    cache_rayint, cache_miss_rayint, cache_prob_rayint,
    # New functions
    cache_gamma, cache_miss_fpi, cache_rrm_meanfield_ode,
    # Importance sampling functions
    cache_is, cache_prob_is, cache_miss_is
)

from .ctmc import (
    ctmc_uniformization, ctmc_timereverse, ctmc_makeinfgen, ctmc_solve,
    ctmc_transient, ctmc_simulate, ctmc_rand, ctmc_ssg, ctmc_stochcomp,
    ctmc_ssg_reachability, ctmc_randomization
)

from .dtmc import (
    dtmc_solve, dtmc_stochcomp, dtmc_timereverse, dtmc_makestochastic, dtmc_rand,
    dtmc_simulate, dtmc_isfeasible
)

from .mc import (
    ctmc_makeinfgen, ctmc_solve, ctmc_transient, ctmc_simulate, ctmc_randomization,
    ctmc_uniformization, ctmc_stochcomp, ctmc_timereverse, ctmc_rand,
    dtmc_solve as mc_dtmc_solve, dtmc_makestochastic, dtmc_isfeasible, dtmc_simulate as mc_dtmc_simulate,
    dtmc_rand as mc_dtmc_rand, dtmc_stochcomp as mc_dtmc_stochcomp, dtmc_timereverse as mc_dtmc_timereverse,
    # New functions
    ctmc_courtois, ctmc_kms, ctmc_multi, ctmc_pseudostochcomp, ctmc_relsolve,
    ctmc_solve_reducible, ctmc_stmonotone, ctmc_takahashi, ctmc_testpf_kolmogorov,
    dtmc_solve_reducible, dtmc_uniformization
)

from .pfqn import (
    pfqn_ca, pfqn_panacea, pfqn_bs, pfqn_mva, pfqn_aql,
    pfqn_mvald, pfqn_mvaldms, pfqn_mvaldmx, pfqn_mvams, pfqn_mvamx,
    pfqn_nc, pfqn_gld, pfqn_gldsingle, pfqn_comomrm,
    pfqn_linearizer, pfqn_linearizerms, pfqn_linearizerpp, pfqn_linearizermx,
    pfqn_kt, pfqn_recal,
    pfqn_cub, pfqn_mmint2, pfqn_ls, pfqn_rd,
    pfqn_fnc, pfqn_propfair, pfqn_xia,
    pfqn_xzabalow, pfqn_xzabaup, pfqn_xzgsblow, pfqn_xzgsbup,
    pfqn_conwayms, pfqn_egflinearizer, pfqn_gflinearizer, pfqn_gld_complex,
    pfqn_gldsingle_complex, pfqn_le_hessian, pfqn_le_hessianZ, pfqn_lldfun,
    pfqn_mci, pfqn_mmint2_gausslegendre, pfqn_mmsample2, pfqn_mushift,
    pfqn_cdfun, pfqn_nca, pfqn_ncld, pfqn_pff_delay, pfqn_sqni,
    pfqn_qzgblow, pfqn_qzgbup, pfqn_nc_sanitize, pfqn_comomrm_ld, pfqn_mvaldmx_ec,
    pfqn_nrl, pfqn_nrp, pfqn_stdf, pfqn_stdf_heur, pfqn_conwayms_core,
    pfqn_conwayms_estimate, pfqn_conwayms_forwardmva, pfqn_mu_ms_gnaux,
    # New functions
    pfqn_ab, pfqn_le, pfqn_le_fpi, pfqn_le_fpiz, pfqn_le_hessianz,
    pfqn_mom, pfqn_mu_ms, pfqn_procomom2, pfqn_schmidt,
    # LCFS queueing network functions
    pfqn_lcfsqn_ca, pfqn_lcfsqn_mva, pfqn_lcfsqn_nc,
    # Replicated stations support
    pfqn_unique, pfqn_expand, pfqn_combine_mi
)

from .mam import (
    map_pie, map_mean, map_var, map_scv, map_skew, map_moment, map_lambda,
    map_acf, map_acfc, map_idc, map_gamma, map_gamma2, map_cdf, map_piq,
    map_embedded, map_count_mean, map_count_var, map_varcount,
    map2_fit, aph_fit, aph2_fit, aph2_fitall, aph2_adjust, mmpp2_fit, mmpp2_fit1,
    mmap_mixture_fit, mmap_mixture_fit_mmap, mamap2m_fit_gamma_fb_mmap, mamap2m_fit_gamma_fb,
    map_exponential, map_erlang, map_hyperexp, map_scale, map_normalize,
    map_timereverse, map_mark, map_infgen,
    map_super, map_sum, map_sumind, map_checkfeasible, map_isfeasible,
    map_feastol, map_largemap, aph2_assemble, ph_reindex, map_rand, map_randn,
    mmap_lambda, mmap_count_mean, mmap_count_var, mmap_count_idc, mmap_idc,
    mmap_sigma2, mmap_exponential, mmap_mixture, mmap_super, mmap_super_safe,
    mmap_compress, mmap_normalize, mmap_scale, mmap_timereverse, mmap_hide,
    mmap_shorten, mmap_maps, mmap_pc, mmap_forward_moment, mmap_backward_moment,
    mmap_cross_moment, mmap_sample, mmap_rand,
    map_sample, aph_rand, randp,
    mmap_count_lambda, mmap_isfeasible, mmap_mark,
    aph_bernstein, map_jointpdf_derivative, map_ccdf_derivative,
    qbd_R, qbd_R_logred, qbd_rg,
    map_pdf, map_prob, map_joint, map_mixture, map_max, map_renewal,
    map_stochcomp,
    qbd_mapmap1, qbd_raprap1, qbd_bmapbmap1, qbd_setupdelayoff,
    map_kurt, mmap_sigma2_cell, amap2_adjust_gamma, amap2_fitall_gamma,
    mmpp2_fit_mu00, mmpp2_fit_mu11, mmpp2_fit_q01, mmpp2_fit_q10,
    assess_compression_quality, compress_adaptive, compress_autocorrelation,
    compress_spectral, compress_with_quality_control,
    # New MAM functions added to fill API gaps
    amap2_fit_gamma, amap2_fit_gamma_map, amap2_fit_gamma_trace,
    aph2_fit_map, aph2_fit_trace, aph_simplify,
    qbd_r, qbd_r_logred,
    # Additional missing MAP/MMAP/MMPP functions from manual
    map_block, map_feasblock, map_kpc, map_pntiter, map_pntquad,
    mmap_embedded, mmap_pie, mmap_issym, mmap_modulate, mmap_mixture_order2, mmap_sum,
    mmpp2_fitc, mmpp2_fitc_approx, mmpp_rand,
    # MAPQN functions
    mapqn_bnd_lr, mapqn_bnd_lr_mva, mapqn_bnd_lr_pf, mapqn_bnd_qr, mapqn_bnd_qr_delay,
    mapqn_bnd_qr_ld, mapqn_lpmodel, mapqn_parameters, mapqn_parameters_factory,
    mapqn_qr_bounds_bas, mapqn_qr_bounds_rsrd,
    # M3PP functions
    m3pp22_fitc_approx_cov, m3pp22_fitc_approx_cov_multiclass, m3pp22_interleave_fitc,
    m3pp2m_fitc, m3pp2m_fitc_approx, m3pp2m_fitc_approx_ag, m3pp2m_fitc_approx_ag_multiclass,
    m3pp2m_interleave, m3pp_interleave_fitc, m3pp_interleave_fitc_theoretical,
    m3pp_interleave_fitc_trace, m3pp_rand, m3pp_superpos_fitc, m3pp_superpos_fitc_theoretical,
    # MAMAP functions
    mamap22_fit_gamma_fs_trace, mamap22_fit_multiclass, mamap2m_coefficients,
    mamap2m_fit, mamap2m_fit_fb_multiclass, mamap2m_fit_mmap, mamap2m_fit_trace,
    # MAPH functions
    maph2m_fit, maph2m_fit_mmap, maph2m_fit_multiclass, maph2m_fit_trace,
    maph2m_fitc_approx, maph2m_fitc_theoretical
)

from .npfqn import (
    npfqn_nonexp_approx, npfqn_traffic_merge, npfqn_traffic_merge_cs, npfqn_traffic_split_cs
)

from .qsys import (
    qsys_mm1, qsys_mmk, qsys_gm1, qsys_mg1, qsys_gig1_approx_lin,
    qsys_gig1_approx_kk, qsys_gig1_approx_whitt, qsys_gig1_approx_allencunneen,
    qsys_gig1_approx_heyman, qsys_gig1_approx_kobayashi, qsys_gig1_approx_marchal,
    qsys_gig1_ubnd_kingman, qsys_gigk_approx, qsys_gigk_approx_kingman,
    # New QSYS functions added to fill API gaps
    qsys_gg1, qsys_gig1_approx_gelenbe, qsys_gig1_approx_kimura, qsys_gig1_approx_klb,
    qsys_gig1_approx_myskja, qsys_gig1_approx_myskja2, qsys_gig1_lbnd,
    qsys_gigk_approx_cosmetatos, qsys_gigk_approx_whitt, qsys_mg1k_loss,
    qsys_mmkk, qsys_mmm, qsys_mminf, qsys_mginf, qsys_mm1k_loss, qsys_mg1k_loss_mgs,
    # M/G/1 with non-preemptive priorities
    qsys_mg1_prio,
    # MAP/PH queue analyzers (BUTools integration)
    qsys_mapg1, qsys_mapg1_cv, qsys_mapmap1, qsys_mapmc, qsys_mapm1,
    qsys_mapph1, qsys_phph1, qsys_phm1, qsys_phmc,
    # MAP/D/c queue analyzers
    qsys_mapdc, qsys_mapd1, qsys_phdc,
    # Batch arrivals
    qsys_mxm1,
    # Size-based scheduling
    qsys_mg1_srpt, qsys_mg1_psjf, qsys_mg1_fb, qsys_mg1_lrpt
)

from .lossn import lossn_erlangfp

from .sn import (
    sn_deaggregate_chain_results, sn_get_arvr_from_tput, sn_get_demands_chain,
    sn_get_node_arvr_from_tput, sn_get_node_tput_from_tput, sn_get_product_form_chain_params,
    sn_get_product_form_params, sn_get_residt_from_respt, sn_get_state_aggr,
    sn_is_state_valid, sn_refresh_visits, sn_has_class_switching, sn_has_fork_join,
    sn_has_load_dependence, sn_has_multi_server, sn_has_priorities, sn_has_product_form,
    sn_has_closed_classes, sn_has_open_classes, sn_has_mixed_classes, sn_has_multi_chain,
    sn_is_closed_model, sn_is_open_model, sn_is_mixed_model,
    sn_has_product_form_not_het_fcfs, sn_print_routing_matrix, sn_has_multi_class,
    # New SN functions added to fill API gaps
    sn_has_dps, sn_has_dps_prio, sn_has_fcfs, sn_has_fractional_populations,
    sn_has_gps, sn_has_gps_prio, sn_has_hol, sn_has_homogeneous_scheduling,
    sn_has_inf, sn_has_lcfs, sn_has_lcfspr, sn_has_polling, sn_has_ps,
    sn_has_rr, sn_has_siro, sn_has_sjf, sn_has_slc, sn_has_snc,
    sn_has_srpt, sn_has_state_dependence, sn_print, sn_summary, sn_validate,
    sn_get_arv_r_from_tput, sn_get_node_arv_r_from_tput, sn_has_lcfs_pr,
    sn_has_multi_class_fcfs, sn_has_multi_class_heter_exp_fcfs,
    sn_has_multi_class_heter_fcfs, sn_is_population_model, sn_rtnodes_to_rtorig
)

from .polling import (
    polling_qsys_1limited, polling_qsys_exhaustive, polling_qsys_gated
)


from .lsn import lsn_max_multiplicity

from .trace import (
    trace_mean, trace_var, trace_skew, mtrace_mean,
    mtrace_backward_moment, mtrace_bootstrap, mtrace_count, mtrace_cov,
    mtrace_cross_moment, mtrace_forward_moment, mtrace_iat2counts,
    mtrace_joint, mtrace_moment, mtrace_moment_simple, mtrace_pc,
    mtrace_sigma, mtrace_sigma2, mtrace_split, mtrace_summary
)

from .wf import (
    wf_analyzer, wf_auto_integration, wf_branch_detector,
    wf_loop_detector, wf_parallel_detector, wf_pattern_updater,
    wf_sequence_detector
)

from .ms import (
    ms_additivesymmetricchisquared, ms_adtest, ms_avgl1linfty, ms_bhattacharyya,
    ms_canberra, ms_chebyshev, ms_chisquared, ms_cityblock, ms_clark,
    ms_condentropy, ms_cosine, ms_cramer_von_mises, ms_czekanowski, ms_dice,
    ms_divergence, ms_entropy, ms_euclidean, ms_fidelity, ms_gower,
    ms_harmonicmean, ms_hellinger, ms_intersection, ms_jaccard, ms_jeffreys,
    ms_jensendifference, ms_jensenshannon, ms_jointentropy, ms_kdivergence,
    ms_kolmogorov_smirnov, ms_kuiper, ms_kulczynskid, ms_kulczynskis,
    ms_kullbackleibler, ms_kumarhassebrook, ms_kumarjohnson, ms_lorentzian,
    ms_matusita, ms_minkowski, ms_motyka, ms_mutinfo, ms_neymanchisquared,
    ms_nmi, ms_nvi, ms_pearsonchisquared, ms_probsymmchisquared, ms_product,
    ms_relatentropy, ms_ruzicka, ms_soergel, ms_sorensen, ms_squaredchord,
    ms_squaredeuclidean, ms_taneja, ms_tanimoto, ms_topsoe, ms_wasserstein,
    ms_wavehegdes
)

__all__ = [
    'cache_mva', 'cache_prob_asy', 'cache_gamma_lp', 'cache_spm', 'cache_xi_fp',
    'cache_miss_spm', 'cache_prob_erec', 'cache_prob_fpi', 'cache_prob_spm',
    'cache_rayint', 'cache_miss_rayint', 'cache_prob_rayint',  # deprecated
    'cache_erec', 'cache_t_hlru', 'cache_t_lrum', 'cache_t_lrum_map', 'cache_ttl_hlru',
    'cache_ttl_lrua', 'cache_ttl_lrum', 'cache_ttl_lrum_map', 'cache_ttl_tree', 'cache_xi_bvh',
    'cache_miss', 'cache_mva_miss', 'cache_miss_asy', 'cache_erec_aux', 'cache_par',
    'cache_t_hlru_aux', 'cache_t_lrum_aux',
    # New cache functions
    'cache_gamma', 'cache_miss_fpi', 'cache_rrm_meanfield_ode',
    # Importance sampling functions
    'cache_is', 'cache_prob_is', 'cache_miss_is',
    'ctmc_uniformization', 'ctmc_timereverse', 'ctmc_makeinfgen', 'ctmc_solve',
    'ctmc_transient', 'ctmc_simulate', 'ctmc_rand', 'ctmc_ssg', 'ctmc_stochcomp',
    'ctmc_ssg_reachability', 'ctmc_randomization',
    'dtmc_solve', 'dtmc_stochcomp', 'dtmc_timereverse', 'dtmc_makestochastic', 'dtmc_rand',
    'dtmc_simulate', 'dtmc_isfeasible',
    'ctmc_makeinfgen', 'ctmc_solve', 'ctmc_transient', 'ctmc_simulate', 'ctmc_randomization',
    'ctmc_uniformization', 'ctmc_stochcomp', 'ctmc_timereverse', 'ctmc_rand',
    'mc_dtmc_solve', 'dtmc_makestochastic', 'dtmc_isfeasible', 'mc_dtmc_simulate',
    'mc_dtmc_rand', 'mc_dtmc_stochcomp', 'mc_dtmc_timereverse',
    # New MC functions
    'ctmc_courtois', 'ctmc_kms', 'ctmc_multi', 'ctmc_pseudostochcomp', 'ctmc_relsolve',
    'ctmc_solve_reducible', 'ctmc_stmonotone', 'ctmc_takahashi', 'ctmc_testpf_kolmogorov',
    'dtmc_solve_reducible', 'dtmc_uniformization',
    'pfqn_ca', 'pfqn_panacea', 'pfqn_bs', 'pfqn_mva', 'pfqn_aql',
    'pfqn_mvald', 'pfqn_mvaldms', 'pfqn_mvaldmx', 'pfqn_mvams', 'pfqn_mvamx',
    'pfqn_nc', 'pfqn_gld', 'pfqn_gldsingle', 'pfqn_comomrm',
    'pfqn_linearizer', 'pfqn_linearizerms', 'pfqn_linearizerpp', 'pfqn_linearizermx',
    'pfqn_kt', 'pfqn_recal',
    'pfqn_cub', 'pfqn_mmint2', 'pfqn_ls', 'pfqn_rd',
    'pfqn_fnc', 'pfqn_propfair', 'pfqn_xia',
    'pfqn_xzabalow', 'pfqn_xzabaup', 'pfqn_xzgsblow', 'pfqn_xzgsbup',
    'pfqn_conwayms', 'pfqn_egflinearizer', 'pfqn_gflinearizer', 'pfqn_gld_complex',
    'pfqn_gldsingle_complex', 'pfqn_le_hessian', 'pfqn_le_hessianZ', 'pfqn_lldfun',
    'pfqn_mci', 'pfqn_mmint2_gausslegendre', 'pfqn_mmsample2', 'pfqn_mushift',
    'pfqn_cdfun', 'pfqn_nca', 'pfqn_ncld', 'pfqn_pff_delay', 'pfqn_sqni',
    'pfqn_qzgblow', 'pfqn_qzgbup', 'pfqn_nc_sanitize', 'pfqn_comomrm_ld', 'pfqn_mvaldmx_ec',
    'pfqn_nrl', 'pfqn_nrp', 'pfqn_stdf', 'pfqn_stdf_heur', 'pfqn_conwayms_core',
    'pfqn_conwayms_estimate', 'pfqn_conwayms_forwardmva', 'pfqn_mu_ms_gnaux',
    # New PFQN functions (added to fill API gaps)
    'pfqn_ab', 'pfqn_le_fpiz', 'pfqn_le_hessianz', 'pfqn_mom', 'pfqn_procomom2', 'pfqn_schmidt',
    # LCFS queueing network functions
    'pfqn_lcfsqn_ca', 'pfqn_lcfsqn_mva', 'pfqn_lcfsqn_nc',
    'map_pie', 'map_mean', 'map_var', 'map_scv', 'map_skew', 'map_moment', 'map_lambda',
    'map_acf', 'map_acfc', 'map_idc', 'map_gamma', 'map_gamma2', 'map_cdf', 'map_piq',
    'map_embedded', 'map_count_mean', 'map_count_var', 'map_varcount',
    'map2_fit', 'aph_fit', 'aph2_fit', 'aph2_fitall', 'aph2_adjust', 'mmpp2_fit', 'mmpp2_fit1',
    'mmap_mixture_fit', 'mmap_mixture_fit_mmap', 'mamap2m_fit_gamma_fb_mmap', 'mamap2m_fit_gamma_fb',
    'map_exponential', 'map_erlang', 'map_hyperexp', 'map_scale', 'map_normalize',
    'map_timereverse', 'map_mark', 'map_infgen',
    'map_super', 'map_sum', 'map_sumind', 'map_checkfeasible', 'map_isfeasible',
    'map_feastol', 'map_largemap', 'aph2_assemble', 'ph_reindex', 'map_rand', 'map_randn',
    'mmap_lambda', 'mmap_count_mean', 'mmap_count_var', 'mmap_count_idc', 'mmap_idc',
    'mmap_sigma2', 'mmap_exponential', 'mmap_mixture', 'mmap_super', 'mmap_super_safe',
    'mmap_compress', 'mmap_normalize', 'mmap_scale', 'mmap_timereverse', 'mmap_hide',
    'mmap_shorten', 'mmap_maps', 'mmap_pc', 'mmap_forward_moment', 'mmap_backward_moment',
    'mmap_cross_moment', 'mmap_sample', 'mmap_rand',
    'map_sample', 'randp',
    'mmap_count_lambda', 'mmap_isfeasible', 'mmap_mark',
    'aph_bernstein', 'map_jointpdf_derivative', 'map_ccdf_derivative',
    'qbd_R', 'qbd_R_logred', 'qbd_rg',
    'map_pdf', 'map_prob', 'map_joint', 'map_mixture', 'map_max', 'map_renewal',
    'map_stochcomp',
    'qbd_mapmap1', 'qbd_raprap1', 'qbd_bmapbmap1', 'qbd_setupdelayoff',
    'map_kurt', 'mmap_sigma2_cell', 'amap2_adjust_gamma', 'amap2_fitall_gamma',
    'mmpp2_fit_mu00', 'mmpp2_fit_mu11', 'mmpp2_fit_q01', 'mmpp2_fit_q10',
    'assess_compression_quality', 'compress_adaptive', 'compress_autocorrelation',
    'compress_spectral', 'compress_with_quality_control',
    # New MAM functions added
    'amap2_fit_gamma', 'amap2_fit_gamma_map', 'amap2_fit_gamma_trace',
    'aph2_fit_map', 'aph2_fit_trace', 'aph_simplify', 'qbd_r', 'qbd_r_logred',
    # Additional missing MAP/MMAP/MMPP functions from manual
    'map_block', 'map_feasblock', 'map_kpc', 'map_pntiter', 'map_pntquad',
    'mmap_embedded', 'mmap_pie', 'mmap_issym', 'mmap_modulate', 'mmap_mixture_order2', 'mmap_sum',
    'mmpp2_fitc', 'mmpp2_fitc_approx', 'mmpp_rand',
    # MAPQN functions
    'mapqn_bnd_lr', 'mapqn_bnd_lr_mva', 'mapqn_bnd_lr_pf', 'mapqn_bnd_qr', 'mapqn_bnd_qr_delay',
    'mapqn_bnd_qr_ld', 'mapqn_lpmodel', 'mapqn_parameters', 'mapqn_parameters_factory',
    'mapqn_qr_bounds_bas', 'mapqn_qr_bounds_rsrd',
    # M3PP functions
    'm3pp22_fitc_approx_cov', 'm3pp22_fitc_approx_cov_multiclass', 'm3pp22_interleave_fitc',
    'm3pp2m_fitc', 'm3pp2m_fitc_approx', 'm3pp2m_fitc_approx_ag', 'm3pp2m_fitc_approx_ag_multiclass',
    'm3pp2m_interleave', 'm3pp_interleave_fitc', 'm3pp_interleave_fitc_theoretical',
    'm3pp_interleave_fitc_trace', 'm3pp_rand', 'm3pp_superpos_fitc', 'm3pp_superpos_fitc_theoretical',
    # MAMAP functions
    'mamap22_fit_gamma_fs_trace', 'mamap22_fit_multiclass', 'mamap2m_coefficients',
    'mamap2m_fit', 'mamap2m_fit_fb_multiclass', 'mamap2m_fit_mmap', 'mamap2m_fit_trace',
    # MAPH functions
    'maph2m_fit', 'maph2m_fit_mmap', 'maph2m_fit_multiclass', 'maph2m_fit_trace',
    'maph2m_fitc_approx', 'maph2m_fitc_theoretical',
    'npfqn_nonexp_approx', 'npfqn_traffic_merge', 'npfqn_traffic_merge_cs', 'npfqn_traffic_split_cs',
    'qsys_mm1', 'qsys_mmk', 'qsys_gm1', 'qsys_mg1', 'qsys_gig1_approx_lin',
    'qsys_gig1_approx_kk', 'qsys_gig1_approx_whitt', 'qsys_gig1_approx_allencunneen',
    'qsys_gig1_approx_heyman', 'qsys_gig1_approx_kobayashi', 'qsys_gig1_approx_marchal',
    'qsys_gig1_ubnd_kingman', 'qsys_gigk_approx', 'qsys_gigk_approx_kingman',
    # New QSYS functions added
    'qsys_gg1', 'qsys_gig1_approx_gelenbe', 'qsys_gig1_approx_kimura', 'qsys_gig1_approx_klb',
    'qsys_gig1_approx_myskja', 'qsys_gig1_approx_myskja2', 'qsys_gig1_lbnd',
    'qsys_gigk_approx_cosmetatos', 'qsys_gigk_approx_whitt', 'qsys_mg1k_loss',
    'qsys_mmkk', 'qsys_mmm', 'qsys_mminf', 'qsys_mginf', 'qsys_mm1k_loss', 'qsys_mg1k_loss_mgs',
    'lossn_erlangfp',
    'sn_deaggregate_chain_results', 'sn_get_arvr_from_tput', 'sn_get_demands_chain',
    'sn_get_node_arvr_from_tput', 'sn_get_node_tput_from_tput', 'sn_get_product_form_chain_params',
    'sn_get_product_form_params', 'sn_get_residt_from_respt', 'sn_get_state_aggr',
    'sn_is_state_valid', 'sn_refresh_visits', 'sn_has_class_switching', 'sn_has_fork_join',
    'sn_has_load_dependence', 'sn_has_multi_server', 'sn_has_priorities', 'sn_has_product_form',
    'sn_has_closed_classes', 'sn_has_open_classes', 'sn_has_mixed_classes', 'sn_has_multi_chain',
    'sn_is_closed_model', 'sn_is_open_model', 'sn_is_mixed_model',
    'sn_has_product_form_not_het_fcfs', 'sn_print_routing_matrix', 'sn_has_multi_class',
    # New SN functions added
    'sn_has_dps', 'sn_has_dps_prio', 'sn_has_fcfs', 'sn_has_fractional_populations',
    'sn_has_gps', 'sn_has_gps_prio', 'sn_has_hol', 'sn_has_homogeneous_scheduling',
    'sn_has_inf', 'sn_has_lcfs', 'sn_has_lcfspr', 'sn_has_polling', 'sn_has_ps',
    'sn_has_rr', 'sn_has_siro', 'sn_has_sjf', 'sn_has_slc', 'sn_has_snc',
    'sn_has_srpt', 'sn_has_state_dependence', 'sn_print', 'sn_summary', 'sn_validate',
    'sn_get_arv_r_from_tput', 'sn_get_node_arv_r_from_tput', 'sn_has_lcfs_pr',
    'sn_has_multi_class_fcfs', 'sn_has_multi_class_heter_exp_fcfs',
    'sn_has_multi_class_heter_fcfs', 'sn_is_population_model', 'sn_rtnodes_to_rtorig',
    'polling_qsys_1limited', 'polling_qsys_exhaustive', 'polling_qsys_gated',
    'lsn_max_multiplicity',
    'trace_mean', 'trace_var', 'trace_skew', 'mtrace_mean',
    'mtrace_backward_moment', 'mtrace_bootstrap', 'mtrace_count', 'mtrace_cov',
    'mtrace_cross_moment', 'mtrace_forward_moment', 'mtrace_iat2counts',
    'mtrace_joint', 'mtrace_moment', 'mtrace_moment_simple', 'mtrace_pc',
    'mtrace_sigma', 'mtrace_sigma2', 'mtrace_split', 'mtrace_summary',
    # Workflow functions
    'wf_analyzer', 'wf_auto_integration', 'wf_branch_detector',
    'wf_loop_detector', 'wf_parallel_detector', 'wf_pattern_updater',
    'wf_sequence_detector',
    # MS (metrics) functions
    'ms_additivesymmetricchisquared', 'ms_adtest', 'ms_avgl1linfty', 'ms_bhattacharyya',
    'ms_canberra', 'ms_chebyshev', 'ms_chisquared', 'ms_cityblock', 'ms_clark',
    'ms_condentropy', 'ms_cosine', 'ms_cramer_von_mises', 'ms_czekanowski', 'ms_dice',
    'ms_divergence', 'ms_entropy', 'ms_euclidean', 'ms_fidelity', 'ms_gower',
    'ms_harmonicmean', 'ms_hellinger', 'ms_intersection', 'ms_jaccard', 'ms_jeffreys',
    'ms_jensendifference', 'ms_jensenshannon', 'ms_jointentropy', 'ms_kdivergence',
    'ms_kolmogorov_smirnov', 'ms_kuiper', 'ms_kulczynskid', 'ms_kulczynskis',
    'ms_kullbackleibler', 'ms_kumarhassebrook', 'ms_kumarjohnson', 'ms_lorentzian',
    'ms_matusita', 'ms_minkowski', 'ms_motyka', 'ms_mutinfo', 'ms_neymanchisquared',
    'ms_nmi', 'ms_nvi', 'ms_pearsonchisquared', 'ms_probsymmchisquared', 'ms_product',
    'ms_relatentropy', 'ms_ruzicka', 'ms_soergel', 'ms_sorensen', 'ms_squaredchord',
    'ms_squaredeuclidean', 'ms_taneja', 'ms_tanimoto', 'ms_topsoe', 'ms_wasserstein',
    'ms_wavehegdes'
]
