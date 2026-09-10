"""Tests for SolverFLD.exportODEs, the LaTeX export of the mean-field ODE
system. Expected strings are cross-validated against the MATLAB
SolverFLD.exportODEs output, which is numerically verified against the ODE
right-hand sides integrated by the solver."""

import os
import sys

import pytest

sys.path.insert(0, os.path.join(os.path.dirname(os.path.abspath(__file__)), '..'))

from line_solver import (Network, Queue, Delay, Source, Sink, ClosedClass, OpenClass,
                         Exp, Erlang, SchedStrategy)
from line_solver.solvers.solver_fld import SolverFLD


def build_closed_exp_model():
    model = Network('D')
    delay = Delay(model, 'Delay1')
    queue = Queue(model, 'Queue1', SchedStrategy.PS)
    c1 = ClosedClass(model, 'C1', 4, delay)
    delay.setService(c1, Exp(1.0))
    queue.setService(c1, Exp(2.0))
    P = model.init_routing_matrix()
    P[c1] = Network.serial_routing(delay, queue)
    model.link(P)
    return model


def build_dps_model():
    model = Network('B')
    delay = Delay(model, 'Delay1')
    qdps = Queue(model, 'QueueDPS', SchedStrategy.DPS)
    b1 = ClosedClass(model, 'B1', 2, delay)
    b2 = ClosedClass(model, 'B2', 3, delay)
    delay.setService(b1, Exp(1.0))
    delay.setService(b2, Exp(0.5))
    qdps.setService(b1, Erlang.fit_mean_and_order(1.0, 2), 2.0)
    qdps.setService(b2, Exp(1.0), 1.0)
    P = model.init_routing_matrix()
    P[b1] = Network.serial_routing(delay, qdps)
    P[b2] = Network.serial_routing(delay, qdps)
    model.link(P)
    return model


def build_open_model():
    model = Network('C')
    source = Source(model, 'Source')
    queue = Queue(model, 'Queue1', SchedStrategy.PS)
    sink = Sink(model, 'Sink')
    oc = OpenClass(model, 'OC')
    source.setArrival(oc, Exp(0.5))
    queue.setService(oc, Exp(1.0))
    P = model.init_routing_matrix()
    P[oc] = Network.serial_routing(source, queue, sink)
    model.link(P)
    return model


def test_closing_scalar_export_closed_model():
    tex = SolverFLD(build_closed_exp_model(), 'closing').exportODEs('', 'scalar')
    assert '% method: closing' in tex
    assert '% form: dx/dt = J*r(x)' in tex
    assert '% nstates: 2' in tex
    assert '% nevents: 2' in tex
    assert '% STATE 1 station=Delay1 class=C1 phase=1' in tex
    assert '% STATE 2 station=Queue1 class=C1 phase=1' in tex
    assert 'n_{2}(\\mathbf{x}) &= x_{2}' in tex
    assert 'g_{2}(\\mathbf{x}) &= \\frac{\\min(n_{2}(\\mathbf{x}),\\, 1)}{n_{2}(\\mathbf{x})}' in tex
    assert '\\frac{\\mathrm{d}x_{1}}{\\mathrm{d}t} &= -x_{1} + 2\\,x_{2}\\,g_{2}(\\mathbf{x})\\\\' in tex
    assert '\\frac{\\mathrm{d}x_{2}}{\\mathrm{d}t} &= x_{1} - 2\\,x_{2}\\,g_{2}(\\mathbf{x})' in tex
    assert '\\mathbf{x}(0) = \\begin{pmatrix} 4 & 0 \\end{pmatrix}^{\\top}' in tex


def test_closing_scalar_export_dps_model():
    tex = SolverFLD(build_dps_model(), 'closing').exportODEs('', 'scalar')
    # DPS: weights normalized to (2/3, 1/3) and folded into the coefficients;
    # the shares divide the capacity min(n_2, S_2), with no additive seed
    assert '\\tilde{n}_{2}(\\mathbf{x}) &= 0.66666667\\,(x_{3} + x_{4}) + 0.33333333\\,(x_{5})' in tex
    assert 'g_{2}(\\mathbf{x}) &= \\frac{\\min(n_{2}(\\mathbf{x}),\\, 1)}{\\tilde{n}_{2}(\\mathbf{x})}' in tex
    assert '\\frac{\\mathrm{d}x_{1}}{\\mathrm{d}t} &= -x_{1} + 1.3333333\\,x_{4}\\,g_{2}(\\mathbf{x})\\\\' in tex
    assert '\\frac{\\mathrm{d}x_{5}}{\\mathrm{d}t} &= 0.5\\,x_{2} - 0.33333333\\,x_{5}\\,g_{2}(\\mathbf{x})' in tex
    assert '\\mathbf{x}(0) = \\begin{pmatrix} 2 & 3 & 0 & 0 & 0 \\end{pmatrix}^{\\top}' in tex


def test_statedep_and_softmin_exports():
    tex = SolverFLD(build_dps_model(), 'statedep').exportODEs('', 'scalar')
    assert 'g_{2,1}(\\mathbf{x}) &= \\begin{cases} 1 & n_{2}(\\mathbf{x}) \\le 1\\\\' in tex
    assert 'x_{3}\\,g_{2,1}(\\mathbf{x})' in tex
    tex2 = SolverFLD(build_dps_model(), 'softmin').exportODEs('', 'scalar')
    assert '% method: softmin' in tex2


def test_matrix_notation_open_model():
    tex = SolverFLD(build_open_model(), 'matrix').exportODEs('', 'matrix')
    assert '% form: dx/dt = W^T*theta(x) + lambda' in tex
    assert '\\frac{\\mathrm{d}\\mathbf{x}}{\\mathrm{d}t} = W^{\\top}\\,\\theta(\\mathbf{x}) + \\boldsymbol{\\lambda}' in tex
    # Source state theta = 0, arrivals via lambda
    assert '\\theta(\\mathbf{x}) = \\begin{bmatrix} 0 \\\\ x_{2}\\,g_{2}(\\mathbf{x}) \\end{bmatrix}' in tex
    assert '\\boldsymbol{\\lambda} = \\begin{pmatrix} 0 & 0.5 \\end{pmatrix}^{\\top}' in tex


def test_closing_scalar_export_open_model():
    tex = SolverFLD(build_open_model(), 'closing').exportODEs('', 'scalar')
    # pseudo-closed recirculation: source mass conservation + queue departures
    assert '% nevents: 2' in tex
    assert '\\frac{\\mathrm{d}x_{1}}{\\mathrm{d}t} &= -0.5 + x_{2}\\,g_{2}(\\mathbf{x})\\\\' in tex
    assert '\\frac{\\mathrm{d}x_{2}}{\\mathrm{d}t} &= 0.5 - x_{2}\\,g_{2}(\\mathbf{x})' in tex
    assert '\\mathbf{x}(0) = \\begin{pmatrix} 1 & 0 \\end{pmatrix}^{\\top}' in tex


def test_unsupported_method_and_notation_rejected():
    with pytest.raises(ValueError):
        SolverFLD(build_closed_exp_model(), 'mfq').exportODEs('', 'scalar')
    with pytest.raises(ValueError):
        SolverFLD(build_closed_exp_model(), 'closing').exportODEs('', 'vector')
    with pytest.raises(ValueError):
        SolverFLD(build_open_model(), 'statedep').exportODEs('', 'scalar')


def test_scalar_and_matrix_share_structure():
    solver = SolverFLD(build_closed_exp_model(), 'closing')
    scalar = solver.exportODEs('', 'scalar')
    matrix = solver.exportODEs('', 'matrix')
    assert '\\frac{\\mathrm{d}\\mathbf{x}}{\\mathrm{d}t} = J\\,r(\\mathbf{x})' in matrix
    assert 'r_{1}(\\mathbf{x}) &= x_{1}\\\\' in matrix
    assert 'r_{2}(\\mathbf{x}) &= 2\\,x_{2}\\,g_{2}(\\mathbf{x})' in matrix
    assert 'stoichiometry' not in scalar


if __name__ == '__main__':
    sys.exit(pytest.main([__file__, '-v']))
