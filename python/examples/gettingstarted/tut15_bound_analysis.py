# %%
# Example 15: Bound analysis with SolverBA
#
# A bounding solver answers a different question from SolverMVA. Instead of a
# single point estimate it returns one guaranteed side of an interval, and the
# .lower/.upper pair of a family brackets the exact solution. Bounds need only
# the service demands and the population, never the service distributions, so
# they are cheap enough to sit inside an optimization loop where a full solve
# would be too slow.
import numpy as np
from line_solver import *
GlobalConstants.set_verbose(VerboseLevel.STD)


def first(x):
    """First entry of a bound field, which may be a scalar or an array."""
    return float(np.asarray(x).ravel()[0])
# %%
# Block 1: model. A think Delay and two Queues of unequal speed, so that the
# bottleneck is well defined.
N = 5                                    # number of jobs in the closed chain
model = Network('BoundsDemo')
delay = Delay(model, 'Think')
q1 = Queue(model, 'Q1', SchedStrategy.PS)
q2 = Queue(model, 'Q2', SchedStrategy.PS)
jobs = ClosedClass(model, 'C', N, delay)
delay.set_service(jobs, Exp(1 / 2))      # think time  Z = 2
q1.set_service(jobs, Exp(1 / 1.0))       # demand D1 = 1.0
q2.set_service(jobs, Exp(1 / 1.5))       # demand D2 = 1.5  (bottleneck)
model.link(Network.serial_routing(delay, q1, q2))
# %%
# Block 2: the exact reference. The throughput at the reference station is the
# system throughput of the closed chain, and is what every bound brackets.
Xexact = float(MVA(model, 'exact').avg_table().Tput[0])
# %%
# Block 3: the bounds table. get_bounds_table is the counterpart of
# avg_table: it reports the bracket per station and class, with columns
# Qlower/Qupper and Tlower/Tupper. A single call evaluates both sides of the
# family of the selected method.
print(SolverBA(model, 'gb.upper').get_bounds_table())
# %%
# Block 4: comparing families. get_bounds is the programmatic accessor,
# returning the raw bracket rather than a formatted table. aba uses only the
# bottleneck demand and the total demand and is the crudest; bjb and gb exploit
# more structure. Which family is sharpest is model dependent.
print('\n%-6s %12s %12s %12s' % ('family', 'Tlower', 'Texact', 'Tupper'))
for family in ['aba', 'bjb', 'gb']:
    b = SolverBA(model, family + '.upper').get_bounds()
    print('%-6s %12.6f %12.6f %12.6f'
          % (family, first(b['Tlower']), Xexact, first(b['Tupper'])))
# %%
# Block 5: a bound hierarchy tightening with the level option. Hierarchical
# families are parameterized by the level option: raising it spends more work
# and returns a tighter pair. The Eager-Sevcik hierarchy pbh becomes exact once
# the level reaches the population, so the bracket width falls to zero at
# level N.
print('\n%-6s %12s %12s %12s' % ('level', 'Tlower', 'Tupper', 'width'))
for level in range(1, N + 1):
    b = SolverBA(model, 'pbh.upper', level=level).get_bounds()
    lo, hi = first(b['Tlower']), first(b['Tupper'])
    print('%-6d %12.6f %12.6f %12.6f' % (level, lo, hi, hi - lo))
# %%
# Block 6: one-sided families. Not every family is two-sided. cub is
# upper-only and mbjb and ldbcmp are lower-only, so the missing side is
# reported as NaN rather than as zero, which keeps "no bound" distinguishable
# from "the bound is zero".
print(SolverBA(model, 'cub.upper').get_bounds_table())

# list_valid_methods reports every method name the solver accepts. Some carry
# structural restrictions beyond the feature set: sb and sib are delay-free,
# and lr additionally requires a single-server single-class closed model, so on
# the model above they raise an error rather than return a wrong answer.
print('SolverBA advertises %d methods.' % len(SolverBA(model).list_valid_methods()))
