"""
MMAP-fed small RR cache with two correlated classes.

A marked MMPP2 arrival stream feeds a small Round-Robin (RR) cache. Its two
marks are bound to two open read classes that share the modulating chain, so
the classes are cross-correlated and autocorrelated in time. Each class reads
the cache with a DIFFERENT item-popularity distribution.

Phase 1 (bursty) emits mostly class-1 references at a high rate; phase 2
(calm) emits mostly class-2 references at a low rate. The shared modulating
chain therefore couples "which class arrives" with "how fast requests arrive".

We solve the same system three ways:
  (1) LDES  - discrete-event simulation of the true MMAP-fed cache;
  (2) CTMC  - exact continuous-time Markov chain of the true system;
  (3) ENV   - three random-environment methods that view the MMPP2 phase as
              an environment modulating phase-conditional Poisson arrivals
              (D1 diagonal):
                'avg'   - fast-environment limit: replace the modulation by its
                          time-average per-class Poisson rates and solve ONE
                          rate-averaged cache;
                'dec'   - slow-environment limit: quasi-stationary
                          decomposition, solve each phase independently and
                          average per-phase hit ratios weighted by phase prob;
                'blend' - state-vector coupling: carry the cache-state
                          distribution across phase switches and average each
                          phase's sojourn-weighted distribution. For a
                          Markovian environment this recovers the exact joint
                          (cache x phase) solution, i.e. it matches CTMC.
The 'avg' and 'dec' limits discard the within-phase temporal correlation of
references, so neither is guaranteed to bracket the exact value.
"""

from line_solver import *
import numpy as np


def build_cache_model(n, m, p_access1, p_access2):
    model = Network('MMAPCache')
    source = Source(model, 'Source')
    cache_node = Cache(model, 'Cache', n, m, ReplacementStrategy.RR)
    sink = Sink(model, 'Sink')

    rd1 = OpenClass(model, 'Read1', 0)
    rd2 = OpenClass(model, 'Read2', 0)
    hit1 = OpenClass(model, 'Hit1', 0)
    mis1 = OpenClass(model, 'Miss1', 0)
    hit2 = OpenClass(model, 'Hit2', 0)
    mis2 = OpenClass(model, 'Miss2', 0)

    cache_node.set_read(rd1, p_access1)
    cache_node.set_read(rd2, p_access2)
    cache_node.set_hit_class(rd1, hit1)
    cache_node.set_miss_class(rd1, mis1)
    cache_node.set_hit_class(rd2, hit2)
    cache_node.set_miss_class(rd2, mis2)

    P = model.init_routing_matrix()
    P.set(rd1, rd1, source, cache_node, 1.0)
    P.set(rd2, rd2, source, cache_node, 1.0)
    P.set(hit1, hit1, cache_node, sink, 1.0)
    P.set(mis1, mis1, cache_node, sink, 1.0)
    P.set(hit2, hit2, cache_node, sink, 1.0)
    P.set(mis2, mis2, cache_node, sink, 1.0)
    model.link(P)
    return model


def set_rates(base_model, lambda1, lambda2):
    """One environment stage: phase-conditional Poisson arrivals per read class."""
    model = base_model.copy()
    source = model.get_node_by_name('Source')
    source.set_arrival(model.classes[0], Exp(lambda1))  # Read1 Poisson
    source.set_arrival(model.classes[1], Exp(lambda2))  # Read2 Poisson
    return model


def cache_mmap_rr_env():
    n = 4  # number of items
    m = 2  # cache capacity

    # Per-class item-popularity distributions (deliberately different)
    p_access1 = DiscreteSampler(np.array([8, 4, 2, 1]) / 15.0)  # class 1 favors low-index items
    p_access2 = DiscreteSampler(np.array([1, 2, 4, 8]) / 15.0)  # class 2 favors high-index items

    # Marked MMPP2 (M3A layout D = {D0, D11, D12}); D1 = D11 + D12 diagonal.
    # Phase 1 bursty (rate 4, 90% class1); phase 2 calm (rate 1, 80% class2).
    # Off-diagonal of D0 are the phase-switch rates (both 0.5).
    D0 = np.array([[-4.5, 0.5], [0.5, -1.5]])
    D11 = np.array([[3.6, 0.0], [0.0, 0.2]])  # class-1 arrivals per phase
    D12 = np.array([[0.4, 0.0], [0.0, 0.8]])  # class-2 arrivals per phase
    mmap = MarkedMAP([D0, D11, D12])

    true_model = build_cache_model(n, m, p_access1, p_access2)
    src = true_model.get_node_by_name('Source')
    src.set_marked_arrival(mmap, [true_model.classes[0], true_model.classes[1]])

    env_base = build_cache_model(n, m, p_access1, p_access2)
    return true_model, env_base, D0, D11, D12


if __name__ == "__main__":
    GlobalConstants.set_verbose(VerboseLevel.STD)
    true_model, env_base, D0, D11, D12 = cache_mmap_rr_env()
    cache_true = true_model.get_node_by_name('Cache')

    # (1) LDES - simulation of the true MMAP-fed cache
    print('\nSOLVER: LDES')
    print(LDES(true_model, samples=200000, seed=23000, verbose=True).get_avg_node_table())
    hit_ldes = np.atleast_1d(cache_true.get_hit_ratio())

    # (2) CTMC - exact solution of the true system
    true_model.reset()
    print('\nSOLVER: CTMC')
    print(CTMC(true_model, 'exact', keep=False, cutoff=1).get_avg_node_table())
    hit_ctmc = np.atleast_1d(cache_true.get_hit_ratio())

    # Random environment: the MMPP2 phase modulates phase-conditional Poisson
    # arrivals (D1 diagonal); the environment switches at the MMPP2
    # phase-transition rates (-D0 diagonal minus the total arrival rate).
    env = Environment('MMPPphase', 2)
    env.add_stage(0, 'Phase1', 'bursty', set_rates(env_base, D11[0, 0], D12[0, 0]))  # 3.6/0.4
    env.add_stage(1, 'Phase2', 'calm', set_rates(env_base, D11[1, 1], D12[1, 1]))    # 0.2/0.8
    env.add_transition(0, 1, Exp(-D0[0, 0] - (D11[0, 0] + D12[0, 0])))               # 0.5
    env.add_transition(1, 0, Exp(-D0[1, 1] - (D11[1, 1] + D12[1, 1])))               # 0.5
    env.init()
    print(env.get_stage_table())

    def solver_factory(mdl):
        return CTMC(mdl, 'exact', keep=False, cutoff=1)

    def env_hit(method, factory, **kwargs):
        opt = {'method': method, 'verbose': False}
        opt.update(kwargs)
        solver = ENV(env, factory, opt)
        solver.get_avg()
        return np.atleast_1d(solver.ensemble[0].get_node_by_name('Cache').get_hit_ratio())

    # (3a) 'avg'   - fast-environment limit (rate-averaged single model)
    hit_avg = env_hit('avg', solver_factory)
    # (3b) 'dec'   - slow-environment quasi-stationary decomposition
    hit_dec = env_hit('dec', solver_factory)

    # (3c) 'blend' - state-vector coupling: carries the cache-state distribution
    # across phase switches and averages each phase's sojourn-weighted
    # distribution. Needs a finite-timespan CTMC inner solver.
    def blend_factory(mdl):
        return CTMC(mdl, 'exact', keep=False, cutoff=1, timespan=[0, 1e3])

    hit_blend = env_hit('blend', blend_factory, iter_max=100, iter_tol=1e-4)

    # (3d) ENV default mean-field with an FLD (refined mean-field) inner solver.
    # The cache is analyzed by the RMF drift; the mean occupancy is carried
    # across phase switches (the cache analog of the queue-length handoff) and
    # the hit ratio is the probEnv-weighted, sojourn-averaged (arrival x
    # hit-prob). Needs a finite-timespan fluid inner solver. This is the fully
    # mean-field counterpart of 'blend', trading the exact joint distribution
    # for a fluid approximation.
    def fld_factory(mdl):
        return FLD(mdl, method='rmf', timespan=[0, 50])

    hit_fld = env_hit('default', fld_factory, iter_max=100, iter_tol=1e-4)

    # Summary: per-read-class actual hit ratio (classes 1=Read1, 2=Read2)
    print('\n--- Actual cache hit ratio per read class ---')
    print('                 Read1      Read2')
    for label, h in (('LDES (sim)  ', hit_ldes), ('CTMC (true) ', hit_ctmc),
                     ('ENV (avg)   ', hit_avg), ('ENV (dec)   ', hit_dec),
                     ('ENV (blend) ', hit_blend), ('ENV (mf/FLD)', hit_fld)):
        print('%s : %8.4f  %8.4f' % (label, h[0], h[1]))
