"""Static applicability test for the 'dae' fluid method.

Python twin of the MATLAB ``fluid_dae_applicable`` and the JAR
``FluidDaeApplicable``. Used by SolverFLD to try 'dae' before dropping a
declined 'minnormal' to a first-order method.

WHY THIS EXISTS SEPARATELY FROM ``fluid_minnormal_applicable``. The two methods
state the SAME closure and differ only in how the coupled equations are
discharged, so a model 'minnormal' accepts is almost always one 'dae' accepts
too. Almost: 'dae' carries a finite-difference Jacobian over the whole unknown
vector rather than one Lyapunov solve, so its state cap is lower; it closes on
the per-station variance only, so DPS and GPS are out; and it has no
decomposition route, so a cache model is out. Those three are exactly the
difference set, and naming them here keeps the fallback ladder from entering a
rung that would refuse the model a moment later.

The test is STATIC, for the same reason the minnormal one is: a rung chosen by
trial and rollback would make the reported method depend on a failed run. The
one condition that cannot be static is the NON-HYPERBOLIC fixed point the ladder
exists to route around -- it exists only once the mean is solved -- and 'dae'
fails on it loudly, which is what moves the ladder to its last rung.

Every condition below mirrors a refusal DaeSolver would otherwise raise, so this
module and those refusals must move together.
"""

from typing import Tuple

from line_solver.api.sn import NodeType, SchedStrategy

from .methods.closing import ClosingMethod


def fluid_dae_applicable(sn, options) -> Tuple[bool, str]:
    """Whether the differential-algebraic route can answer this model.

    Parameters
    ----------
    sn : NetworkStruct
        Network structure, after sn_nonmarkov_toph.
    options : SolverFLDOptions
        Solver options.

    Returns
    -------
    (ok, reason)
        ok is True when 'dae' can be selected; reason names the blocking feature
        otherwise and is empty when ok.
    """
    # A cache model is a DECOMPOSITION, and 'dae' has no arm for it: 'minnormal'
    # on a cache model routes through the rmf alternation with the closure inside
    # its network step, while DaeSolver refuses outright, because a decomposition
    # has no single drift for the algebraic constraint to attach to.
    nodetype = getattr(sn, 'nodetype', None)
    if nodetype is None:
        nodetype = []
    for nt in nodetype:
        if nt == NodeType.CACHE:
            return False, ('a cache model is answered by the decomposition analyzer, '
                           'which has no dae route')

    probe = ClosingMethod(sn, options)
    _, _, phases = probe._extract_service_params()
    sched = probe._get_sched()

    # The DPS and GPS shares close on the covariance BETWEEN station coordinates,
    # not on the station total, so their closure state is a matrix block rather
    # than the scalar the Newton vector carries. Mirrors the DaeSolver refusal.
    for i in range(sn.nstations):
        if sched[i] in (SchedStrategy.DPS, SchedStrategy.GPS):
            return False, ('station %d uses %s, whose share closes on the covariance between '
                           'its class coordinates rather than on the station variance'
                           % (i + 1, sched[i]))

    # The simultaneous solve is quartic overall, against the cubic of one Lyapunov
    # solve, so its crossover is lower than the 200 'minnormal' permits and it
    # carries its own cap. Same count the minnormal test forms, so the two limits
    # are read on the same scale.
    maxstate = 100
    cfg = getattr(options, 'config', None)
    if isinstance(cfg, dict) and cfg.get('dae_maxstate'):
        maxstate = int(cfg['dae_maxstate'])
    nstate = int(phases.sum())
    if nstate > maxstate:
        return False, ('the phase-resolved state has %d coordinates, above the dae_maxstate '
                       'limit of %d' % (nstate, maxstate))

    return True, ''
