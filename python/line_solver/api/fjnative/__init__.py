"""
Native fork-join support for SolverCTMC/SolverSSA (python native).

Port of matlab/src/api/fj/. The tag-augmentation adapter (fjtag) lives in
io/model_adapter.py; this package holds the validation and metric-foldback
helpers used by the CTMC/SSA fork-join path, together with the analytical
fork-join API of A. Thomasian, "Analysis of Fork/Join and Related Queueing
Systems", ACM Computing Surveys 47(2), Article 17, 2014.

The analytical functions live here rather than under api/fj because that name is
already bound, for backward compatibility, to the third-party fork-join codes in
lib/thirdparty/fj.
"""

from .sn_fj_validate import sn_fj_validate, sn_fj_supports
from .fj_ordstat_exp import fj_ordstat_exp
from .sn_join_quorum import sn_join_quorum
from .sn_join_siblings import sn_join_siblings
from .sn_join_droprate import sn_join_droprate
from .sn_fj_foldback import sn_fj_foldback
from .fj_tail_forktail import (fj_tail_forktail, fj_mg1_respt_moments, ge_fit,
                               forktail_percentiles)
from .fj_tail_ordstat import fj_tail_ordstat
from .thomasian import (
    fj_harmonic,
    fj_qgb,
    fj_amva,
    fj_respt_closed,
    fj_xmax_het,
    fj_lst_max_het,
    fj_xmax_moments_het,
    fj_xmax_hz,
    fj_xmax_hz_het,
    fj_char_max_discrete,
    fj_char_max_blom,
    fj_cox_fit,
    fj_xmax_coxian,
    fj_dispersion,
    fj_delay_opt,
    fj_respt_nosplit,
    fj_respt_bulk,
    fj_ism_green,
    fj_tsm_capacity,
    fj_serialization,
    fj_dag_makespan,
)

__all__ = ['sn_fj_validate', 'sn_fj_supports', 'sn_fj_foldback', 'fj_ordstat_exp', 'sn_join_quorum',
           'sn_join_siblings', 'sn_join_droprate',
           'fj_tail_forktail', 'fj_tail_ordstat', 'fj_mg1_respt_moments', 'ge_fit',
           'forktail_percentiles',
           'fj_harmonic', 'fj_qgb', 'fj_amva', 'fj_respt_closed',
           'fj_xmax_het', 'fj_lst_max_het', 'fj_xmax_moments_het',
           'fj_xmax_hz', 'fj_xmax_hz_het',
           'fj_char_max_discrete', 'fj_char_max_blom',
           'fj_cox_fit', 'fj_xmax_coxian',
           'fj_dispersion', 'fj_delay_opt',
           'fj_respt_nosplit', 'fj_respt_bulk', 'fj_ism_green',
           'fj_tsm_capacity', 'fj_serialization', 'fj_dag_makespan']
