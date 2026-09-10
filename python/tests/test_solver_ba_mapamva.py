"""MAP-AMVA LP bounds (Casale-Smirni, DSN 2009) in SolverBA.

The family is the only one in SolverBA derived FOR a correlated service process,
so these tests check the two things that distinguishes it from its neighbours:
that its bracket actually contains the exact solution of a MAP model, and that
the feature gate lets a MAP reach it and no other family.
"""

import unittest

import numpy as np

from line_solver import (Network, Queue, Delay, ClosedClass, SchedStrategy, Exp,
                         Erlang, MAP, SolverBA, SolverCTMC)


# A two-phase MAP: D0+D1 is a proper generator, D1 has an off-diagonal entry so
# successive services are correlated (it is not a renewal PH).
D0 = np.array([[-3.0, 0.5], [0.2, -0.4]])
D1 = np.array([[2.5, 0.0], [0.2, 0.0]])


def _tandem(map_first=False, njobs=8):
    model = Network('mapamva_tandem')
    q1 = Queue(model, 'Queue1', SchedStrategy.FCFS)
    q2 = Queue(model, 'Queue2', SchedStrategy.FCFS)
    cl = ClosedClass(model, 'Class1', njobs, q1)
    if map_first:
        q1.setService(cl, MAP(D0, D1))
        q2.setService(cl, Exp(2.0))
    else:
        q1.setService(cl, Exp(2.0))
        q2.setService(cl, MAP(D0, D1))
    model.link(Network.serialRouting(q1, q2))
    return model


def _cols(model, method):
    t = (SolverCTMC(model) if method == 'ctmc' else SolverBA(model, method)).getAvgTable()
    return {c: t[c].to_numpy().astype(float) for c in ('QLen', 'Util', 'RespT', 'Tput')}


class TestTheBracketContainsTheExactSolution(unittest.TestCase):

    def _assert_brackets(self, model):
        lo = _cols(model, 'mapamva.lower')
        up = _cols(model, 'mapamva.upper')
        exact = _cols(model, 'ctmc')
        for col in ('QLen', 'Util', 'RespT', 'Tput'):
            for i in range(len(exact[col])):
                self.assertLessEqual(
                    lo[col][i], exact[col][i] + 1e-9,
                    '%s station %d: lower %g above exact %g' % (col, i, lo[col][i], exact[col][i]))
                self.assertLessEqual(
                    exact[col][i], up[col][i] + 1e-9,
                    '%s station %d: exact %g above upper %g' % (col, i, exact[col][i], up[col][i]))

    def test_a_map_at_the_last_station_is_bracketed(self):
        self._assert_brackets(_tandem())

    def test_a_map_at_the_first_station_is_bracketed(self):
        # The LP requires the phase-carrying queue to be index M; the analyzer
        # PERMUTES it there and inverts the permutation before reporting, so the
        # bracket must hold with the stations the other way round too.
        self._assert_brackets(_tandem(map_first=True))

    def test_an_all_exponential_model_degenerates_to_one_level(self):
        # No phase-carrying station: K collapses to 1 and the balances become the
        # product-form ones. Still a bracket, and the exact solution is in it.
        model = Network('mapamva_exp')
        q1 = Queue(model, 'Queue1', SchedStrategy.FCFS)
        q2 = Queue(model, 'Queue2', SchedStrategy.FCFS)
        cl = ClosedClass(model, 'Class1', 5, q1)
        q1.setService(cl, Exp(2.0))
        q2.setService(cl, Exp(1.0))
        model.link(Network.serialRouting(q1, q2))
        self._assert_brackets(model)

    def test_the_two_sides_are_ordered(self):
        model = _tandem()
        lo = _cols(model, 'mapamva.lower')
        up = _cols(model, 'mapamva.upper')
        for col in ('QLen', 'Util', 'RespT', 'Tput'):
            self.assertTrue(np.all(lo[col] <= up[col] + 1e-9), col)


class TestTheFeatureGate(unittest.TestCase):

    def test_a_map_model_reaches_mapamva_and_nothing_else(self):
        # 'MAP' sits in the BASE feature envelope so that this family can accept
        # it, and getMethodFeatureSet takes it back from every other: each of the
        # rest reads the service mean alone and would bracket a DIFFERENT system.
        solver = SolverBA(_tandem())
        supported = [m for m in solver.listValidMethods()
                     if solver.supportsModelMethod(m)[0]]
        self.assertEqual(sorted(['mapamva.lower', 'mapamva.upper']), sorted(supported))

    def test_a_renewal_model_does_not_lose_its_own_families(self):
        # The mirror check: taking MAP away from the others must not take
        # anything else with it.
        model = Network('mapamva_renewal')
        q1 = Queue(model, 'Queue1', SchedStrategy.FCFS)
        q2 = Queue(model, 'Queue2', SchedStrategy.FCFS)
        cl = ClosedClass(model, 'Class1', 5, q1)
        q1.setService(cl, Exp(2.0))
        q2.setService(cl, Erlang.fitMeanAndOrder(1.0, 2))
        model.link(Network.serialRouting(q1, q2))
        solver = SolverBA(model)
        supported = [m for m in solver.listValidMethods()
                     if solver.supportsModelMethod(m)[0]]
        for name in ('gb.upper', 'gb.lower', 'aba.upper', 'mapamva.upper'):
            self.assertIn(name, supported)

    def test_a_delay_station_is_refused(self):
        # The LP is a network of queues, and Casale-Smirni name the delay
        # extension as open work.
        model = Network('mapamva_delay')
        d = Delay(model, 'Think')
        q = Queue(model, 'Queue1', SchedStrategy.FCFS)
        cl = ClosedClass(model, 'Class1', 4, d)
        d.setService(cl, Exp(1.0))
        q.setService(cl, MAP(D0, D1))
        model.link(Network.serialRouting(d, q))
        self.assertNotIn('mapamva.upper', SolverBA(model).listValidMethods())
        with self.assertRaises(Exception):
            SolverBA(model, 'mapamva.upper').getAvgTable()

    def test_two_phase_carrying_stations_are_refused(self):
        # The LP gives queue M the (D0,D1) pair and every other queue a SCALAR
        # rate, so it has nowhere to put a second phase-type station.
        model = Network('mapamva_two_ph')
        q1 = Queue(model, 'Queue1', SchedStrategy.FCFS)
        q2 = Queue(model, 'Queue2', SchedStrategy.FCFS)
        cl = ClosedClass(model, 'Class1', 4, q1)
        q1.setService(cl, Erlang.fitMeanAndOrder(0.5, 2))
        q2.setService(cl, MAP(D0, D1))
        model.link(Network.serialRouting(q1, q2))
        with self.assertRaises(Exception):
            SolverBA(model, 'mapamva.upper').getAvgTable()


class TestTheApiAggregatesOverLevels(unittest.TestCase):

    def test_level_zero_is_tighter_than_the_sum_of_the_levels(self):
        # sum_k max U_i^k is also an upper bound on U_i = sum_k U_i^k, but a
        # looser one: the phases cannot all peak at once. The aggregate
        # objective is what the paper's bounds are stated on.
        from line_solver.api.mapqn import mapqn_bnd_lr_mva
        from line_solver.api.mapqn.parameters import MVAVersionParameters
        v = np.array(D0, copy=True)
        np.fill_diagonal(v, 0.0)
        params = MVAVersionParameters(_M=2, _N=6, K=2, muM=np.array([2.0]),
                                      muMAP=D1, r=np.array([[0.0, 1.0], [1.0, 0.0]]), v=v)
        aggregate = mapqn_bnd_lr_mva(params, 2, 0, 'max', 'UN').objective_value
        per_level = sum(mapqn_bnd_lr_mva(params, 2, k, 'max', 'UN').objective_value
                        for k in (1, 2))
        self.assertLessEqual(aggregate, per_level + 1e-9)

    def test_the_sense_argument_orders_the_two_directions(self):
        from line_solver.api.mapqn import mapqn_bnd_lr_mva
        from line_solver.api.mapqn.parameters import MVAVersionParameters
        v = np.array(D0, copy=True)
        np.fill_diagonal(v, 0.0)
        params = MVAVersionParameters(_M=2, _N=6, K=2, muM=np.array([2.0]),
                                      muMAP=D1, r=np.array([[0.0, 1.0], [1.0, 0.0]]), v=v)
        for var in ('UN', 'QN'):
            lo = mapqn_bnd_lr_mva(params, 1, 0, 'min', var).objective_value
            hi = mapqn_bnd_lr_mva(params, 1, 0, 'max', var).objective_value
            self.assertLessEqual(lo, hi + 1e-9, var)


if __name__ == '__main__':
    unittest.main()
