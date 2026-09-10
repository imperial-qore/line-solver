"""
Native-Python convergence-guard regression for SolverLN.

model_C2_L2_T4_P2_1_doesnotconverge_v3.lqnx exposes an AMVA layer-bistability
(two coexisting fixed points). Cold-init SolverLN diverged on this model before
the init_sol warm-start fix landed in all codebases. This test guards that the
warm start keeps SolverLN(MVA) converging.
"""

import os

from line_solver import LayeredNetwork, SolverLN, SolverMVA, GlobalConstants, VerboseLevel

MODEL = os.path.join(os.path.dirname(os.path.abspath(__file__)), 'models',
                     'model_C2_L2_T4_P2_1_doesnotconverge_v3.lqnx')


def test_ln_doesnotconverge_v3():
    GlobalConstants.set_verbose(VerboseLevel.SILENT)
    model = LayeredNetwork.parse_xml(MODEL)
    solver = SolverLN(model, lambda m: SolverMVA(m))
    avg_table = solver.getAvgTable()
    assert solver.hasconverged, "SolverLN(MVA) should converge on doesnotconverge_v3"
    assert avg_table is not None
