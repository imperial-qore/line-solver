/*
 * Copyright (c) 2012-2026, QORE Lab, Imperial College London
 * All rights reserved.
 *
 * The joint-dependence evaluation point of solver_amvald, and the fact that
 * solver_amvald applies eta AT ALL.
 *
 * C++ used to assemble its scaling list from `st.cdscaling` alone, so a
 * joint-dependent station reached AMVA (solver_mva.h routes cd||jd there) and
 * was solved UNSCALED -- a silently wrong answer rather than a refusal. MATLAB
 * and the JAR did apply eta, but at the cdterm point, incrementing EVERY class
 * coordinate; that is inert for a beta reading its own marginal and is not
 * inert for an eta reading the whole row.
 *
 * Model: the IS+2xOI model of Casale/Comte/Dorsman, "Notes on Order-Independent
 * and Pass-and-Swap CQNs", under the QD-AMVA closure -- an IS delay with think
 * rates (0.7, 1.1) and two identical PS stations carrying the support-rank eta
 * of a three-server compatibility station, mu = (2, 1, 2) with server 1
 * dedicated to class 1, server 2 shared and server 3 dedicated to class 2. At
 * N = (1, 1): eta = 3 with one class present, 5 with both, bilinear between.
 *
 * The expected values are python's, which has always used this point.
 */
#include <algorithm>
#include <vector>

#include "doctest.h"
#include "line/lang/qn/network_builder.h"
#include "line/solvers/mva/solver_mva_runner.h"

namespace qn = line::qn;
namespace mva = line::mva;
using line::lang::SchedStrategy;
using D = line::lang::Distrib<double>;

namespace {

/** svc(n) = 1/mu(supp n) on the unit box, bilinearly interpolated, as a rate. */
std::vector<double> eta_support(const std::vector<double>& n) {
    const double n1 = std::min(std::max(n[0], 0.0), 1.0);
    const double n2 = std::min(std::max(n[1], 0.0), 1.0);
    const double s10 = 1.0 / 3.0, s01 = 1.0 / 3.0, s11 = 1.0 / 5.0;
    const double svc = n1 * (1 - n2) * s10 + (1 - n1) * n2 * s01 + n1 * n2 * s11;
    return std::vector<double>(1, svc > 0.0 ? 1.0 / svc : 1.0);
}

qn::Network<double> note_model() {
    qn::Network<double> m("IS+2xOI");
    const std::size_t is = m.add_delay("IS");
    const std::size_t q1 = m.add_queue("CD1", SchedStrategy::PS);
    const std::size_t q2 = m.add_queue("CD2", SchedStrategy::PS);
    const std::size_t c1 = m.add_closed_class("C1", 1.0, is);
    const std::size_t c2 = m.add_closed_class("C2", 1.0, is);
    m.set_service(is, c1, D::exp_rate(0.7));
    m.set_service(is, c2, D::exp_rate(1.1));
    const std::vector<double> peak(1, 5.0);
    for (std::size_t q : {q1, q2}) {
        m.set_service(q, c1, D::exp_rate(1.0));
        m.set_service(q, c2, D::exp_rate(1.0));
        m.set_number_of_servers(q, 1);
        m.set_joint_dependence(q, eta_support, peak);
    }
    qn::RoutingMatrix<double> P;
    for (std::size_t c : {c1, c2}) {   // per class: the 3-argument set covers one class only
        P.set(c, c, is, q1, 1.0);
        P.set(c, c, q1, q2, 1.0);
        P.set(c, c, q2, is, 1.0);
    }
    m.link(P);
    return m;
}

}  // namespace

TEST_CASE("jd eval point: AMVA applies eta at the arrival-theorem point") {
    qn::Network<double> m = note_model();
    mva::MvaOptions opt;
    opt.method = "qd";
    line::Matrix<double> init;
    const mva::AvgResult<double> r = mva::solver_mva_run_analyzer(m.get_struct(), opt, init);

    // python line_solver, SolverMVA(method='qd'), same model
    const double expected[3][2] = {{0.6589337427, 0.5565433634},
                                   {0.1705331157, 0.2217281442},
                                   {0.1705331157, 0.2217281442}};
    double tot = 0.0;
    for (std::size_t k = 0; k < 3; ++k)
        for (std::size_t c = 0; c < 2; ++c) {
            CHECK(r.QN(k, c) == doctest::Approx(expected[k][c]).epsilon(1e-5));
            tot += r.QN(k, c);
        }
    CHECK(tot == doctest::Approx(2.0).epsilon(1e-5));
}
