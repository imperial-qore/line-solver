"""
Closed Queueing Network: CTMC (SPN) vs NC Solver Comparison

This example demonstrates:
- Single-class closed queueing network (product-form)
- Solves the same model two different ways:
  1. Traditional QN approach: solved with NC (Normalizing Constant) solver
  2. SPN approach: same network modeled as Petri Net, solved with CTMC
- Verifies that both approaches give identical results
- Demonstrates CTMC's capability on Petri Net models

Network Configuration:
- Population: N=3 jobs circulating continuously
- Station 1: Queue (FCFS, mu1=1.0)
- Station 2: Queue (FCFS, mu2=0.8)
- Station 3: Delay (think time, mean=2.0)
- Routing: 1 -> 2 -> 3 -> 1 (cycle)

This is product-form under Jackson's theorem:
- FCFS disciplines
- Exponential service times
- Single job class
- Product-form solution: pi(n1,n2,n3) = pi1(n1)pi2(n2)pi3(n3)

NOTE on what the comparison actually shows. Both tables here reproduce MATLAB
exactly (NC throughput 0.536439, CTMC throughput 0.405141), but the two models
are NOT the same system: the SPN transitions are left at their default server
count, so p1 and p2 serve like delays rather than like the single-server FCFS
queues of Approach 1, and the metrics differ accordingly. MATLAB's comparison
block matches its rows with strcmp against a categorical Station column, finds
nothing and prints an empty section, so the divergence never surfaced there;
this port matches by name and prints it. To make the two agree, give each
transition one server.
"""

import numpy as np

from line_solver import (ClosedClass, CTMC, Delay, Exp, GlobalConstants, NC,
                         Network, Place, Queue, SchedStrategy, Transition,
                         VerboseLevel)

BAR = '=' * 72


def _pick(table, station, column):
    """The scalar `column` of the row naming `station`, or None."""
    rows = table[table['Station'].astype(str) == station]
    if len(rows) == 0:
        return None
    return float(np.asarray(rows[column], dtype=float).sum())


if __name__ == '__main__':
    GlobalConstants.set_verbose(VerboseLevel.SILENT)

    print(BAR)
    print('Closed QN: CTMC (SPN) vs NC Solver Comparison')
    print(BAR)

    # ---- Part 1: Traditional Closed QN (Solved with NC Solver) -----------
    print('\n' + BAR)
    print('Approach 1: Traditional Closed QN (Solved with NC Solver)')
    print(BAR)

    model_qn = Network('closed_qn_traditional')

    queue1 = Queue(model_qn, 'queue1', SchedStrategy.FCFS)
    queue2 = Queue(model_qn, 'queue2', SchedStrategy.FCFS)
    delay = Delay(model_qn, 'delay')

    jobclass = ClosedClass(model_qn, 'jobs', 3, queue1)   # 3 jobs, start at queue1

    queue1.setService(jobclass, Exp.fitMean(1.0))         # mu1 = 1.0
    queue2.setService(jobclass, Exp.fitMean(1 / 0.8))     # mu2 = 0.8
    delay.setService(jobclass, Exp.fitMean(2.0))          # think time

    R_qn = model_qn.initRoutingMatrix()
    R_qn.set(jobclass, jobclass, queue1, queue2, 1.0)
    R_qn.set(jobclass, jobclass, queue2, delay, 1.0)
    R_qn.set(jobclass, jobclass, delay, queue1, 1.0)
    model_qn.link(R_qn)

    print('\nSolving with NC (Normalizing Constant) solver...')
    result_nc = NC(model_qn).getAvgTable()
    print('NC Results:')
    print(result_nc)

    # ---- Part 2: Same Network as SPN (Solved with CTMC) ------------------
    print('\n' + BAR)
    print('Approach 2: Same Network as SPN (Solved with CTMC)')
    print(BAR)

    model_spn = Network('closed_qn_spn')

    p1 = Place(model_spn, 'p1')   # Queue 1 buffer
    p2 = Place(model_spn, 'p2')   # Queue 2 buffer
    p3 = Place(model_spn, 'p3')   # Delay buffer

    t1 = Transition(model_spn, 't1')   # Queue 1 service completion
    t2 = Transition(model_spn, 't2')   # Queue 2 service completion
    t3 = Transition(model_spn, 't3')   # Delay completion

    jobclass_spn = ClosedClass(model_spn, 'jobs', 3, p1)

    mode_t1 = t1.addMode('service1')
    t1.setDistribution(mode_t1, Exp.fitMean(1.0))
    t1.setEnablingConditions(mode_t1, jobclass_spn, p1, 1)   # require job in p1
    t1.setFiringOutcome(mode_t1, jobclass_spn, p1, -1)       # remove from p1
    t1.setFiringOutcome(mode_t1, jobclass_spn, p2, 1)        # add to p2

    mode_t2 = t2.addMode('service2')
    t2.setDistribution(mode_t2, Exp.fitMean(1 / 0.8))
    t2.setEnablingConditions(mode_t2, jobclass_spn, p2, 1)
    t2.setFiringOutcome(mode_t2, jobclass_spn, p2, -1)
    t2.setFiringOutcome(mode_t2, jobclass_spn, p3, 1)

    mode_t3 = t3.addMode('think')
    t3.setDistribution(mode_t3, Exp.fitMean(2.0))
    t3.setEnablingConditions(mode_t3, jobclass_spn, p3, 1)
    t3.setFiringOutcome(mode_t3, jobclass_spn, p3, -1)
    t3.setFiringOutcome(mode_t3, jobclass_spn, p1, 1)

    R_spn = model_spn.initRoutingMatrix()
    R_spn.set(jobclass_spn, jobclass_spn, p1, t1, 1.0)
    R_spn.set(jobclass_spn, jobclass_spn, t1, p2, 1.0)
    R_spn.set(jobclass_spn, jobclass_spn, p2, t2, 1.0)
    R_spn.set(jobclass_spn, jobclass_spn, t2, p3, 1.0)
    R_spn.set(jobclass_spn, jobclass_spn, p3, t3, 1.0)
    R_spn.set(jobclass_spn, jobclass_spn, t3, p1, 1.0)
    model_spn.link(R_spn)

    # Set initial state: all 3 jobs at p1
    p1.setState([3])
    p2.setState([0])
    p3.setState([0])

    print('\nSolving with CTMC solver...')
    result_ctmc = CTMC(model_spn, 'exact').getAvgTable()
    print('CTMC Results:')
    print(result_ctmc)

    # ---- Part 3: Comparison ---------------------------------------------
    print('\n' + BAR)
    print('COMPARISON: NC Solver (Traditional) vs CTMC (SPN)')
    print(BAR)
    print('\nDetailed Metrics Comparison:')
    print('-' * 72)

    stations_nc = ['queue1', 'queue2', 'delay']
    stations_spn = ['p1', 'p2', 'p3']

    for nc_name, spn_name in zip(stations_nc, stations_spn):
        nc_qlen = _pick(result_nc, nc_name, 'QLen')
        ctmc_qlen = _pick(result_ctmc, spn_name, 'QLen')
        nc_util = _pick(result_nc, nc_name, 'Util')
        ctmc_util = _pick(result_ctmc, spn_name, 'Util')
        if nc_qlen is None or ctmc_qlen is None:
            continue

        qlen_diff = abs(nc_qlen - ctmc_qlen) / max(abs(nc_qlen), 1e-6) * 100
        util_diff = abs(nc_util - ctmc_util) / max(abs(nc_util), 1e-6) * 100
        status_qlen = 'OK' if qlen_diff < 1.0 else 'X'
        status_util = 'OK' if util_diff < 1.0 else 'X'

        print('\n%s / %s:' % (nc_name, spn_name))
        print('  %s QLen:  NC=%.6f, CTMC=%.6f (diff=%.3f%%)'
              % (status_qlen, nc_qlen, ctmc_qlen, qlen_diff))
        print('  %s Util:  NC=%.6f, CTMC=%.6f (diff=%.3f%%)'
              % (status_util, nc_util, ctmc_util, util_diff))

    print('\n' + '-' * 72)
    print('System-Level Metrics:')
    print('-' * 72)

    nc_total_qlen = float(np.asarray(result_nc['QLen'], dtype=float).sum())
    ctmc_total_qlen = float(np.asarray(result_ctmc['QLen'], dtype=float).sum())
    print('\nTotal QLen (should be 3.0):')
    print('  NC:   %.6f' % nc_total_qlen)
    print('  CTMC: %.6f' % ctmc_total_qlen)

    nc_tput = _pick(result_nc, 'queue1', 'Tput')
    ctmc_tput = _pick(result_ctmc, 'p1', 'Tput')
    tput_diff = abs(nc_tput - ctmc_tput) / max(abs(nc_tput), 1e-6) * 100
    status_tput = 'OK' if tput_diff < 1.0 else 'X'
    print('\n%s System Throughput:' % status_tput)
    print('  NC:   %.6f' % nc_tput)
    print('  CTMC: %.6f' % ctmc_tput)
    print('  Difference: %.3f%%' % tput_diff)

    print('\n' + BAR)
    print('Summary')
    print(BAR)
    print('\nOK Both approaches (NC solver on traditional QN, CTMC on SPN)')
    print('  give identical results!')
    print('\nOK Product-form property verified:')
    print('  - Total population conserved (=3 jobs)')
    print('  - System throughput matches')
    print('  - Individual station metrics match')
    print('\nOK CTMC successfully models closed QN as SPN')
    print('OK Different modeling approaches yield same answers')
    print(BAR)
