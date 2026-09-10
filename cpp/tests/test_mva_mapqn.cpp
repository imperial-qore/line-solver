/**
 * @file test_mva_mapqn.cpp
 * @brief SolverMVA method 'amva.mapqn': the horizontal-cut mean value analysis
 * (Casale-Smirni, IEEE/IFIP DSN 2009, balances closed by the arrival theorem) for a
 * closed model of one exponential delay and one FCFS single-server queue with a MAP
 * service per class.
 *
 * The reference values are the prototype's, which every codebase reproduces to solver
 * precision; the CTMC bounds the accuracy at the level the method has; the exponential
 * case must be exact multiclass MVA.
 */
#include <string>
#include <vector>

#include "doctest.h"
#include "line/api/mapqn/mapqn_amva.h"
#include "line/lang/qn/network_builder.h"
#include "line/solvers/ctmc/solver_ctmc_analyzer.h"
#include "line/solvers/mva/solver_mva_runner.h"

using namespace line;
using lang::Distrib;
using lang::ProcessType;
using lang::SchedStrategy;

namespace {

using D = Distrib<double>;

Matrix<double> m2(double a, double b, double c, double d) {
    Matrix<double> M(2, 2, 0.0);
    M(0, 0) = a; M(0, 1) = b; M(1, 0) = c; M(1, 1) = d;
    return M;
}

D map_a() { return D::map_dist(m2(-3.0, 0.5, 0.2, -0.4), m2(2.5, 0.0, 0.2, 0.0), ProcessType::MAP); }
D map_b() { return D::map_dist(m2(-1.8, 0.3, 0.6, -0.9), m2(1.5, 0.0, 0.0, 0.3), ProcessType::MAP); }

/** delay first, queue second; every class cycles delay -> queue -> delay */
qn::Network<double> model(const std::vector<double>& think, const std::vector<D>& svc,
                          const std::vector<double>& njobs, SchedStrategy sched = SchedStrategy::FCFS) {
    qn::Network<double> m("mapqn");
    const std::size_t d = m.add_delay("Delay");
    const std::size_t q = m.add_queue("Queue", sched);
    for (std::size_t r = 0; r < njobs.size(); ++r) {
        const std::size_t c = m.add_closed_class("C" + std::to_string(r + 1), njobs[r], d);
        m.set_service(d, c, D::exp_rate(1.0 / think[r]));
        m.set_service(q, c, svc[r]);
    }
    qn::RoutingMatrix<double> P;
    for (std::size_t r = 1; r <= njobs.size(); ++r) {
        P.set(r, r, d, q, 1.0);
        P.set(r, r, q, d, 1.0);
    }
    m.link(P);
    return m;
}

mva::AvgResult<double> run(qn::Network<double>& m, const std::string& method) {
    mva::MvaOptions opt;
    opt.method = method;
    Matrix<double> init;
    return mva::solver_mva_run_analyzer(m.get_struct(), opt, init);
}

bool listed(qn::Network<double>& m, const std::string& name) {
    for (const std::string& s : mva::list_valid_methods(m.get_struct()))
        if (s == name) return true;
    return false;
}

}  // namespace

TEST_CASE("amva.mapqn reproduces the reference values on two classes") {
    qn::Network<double> m = model({2.0, 8.0}, {map_a(), map_b()}, {2, 2});
    const mva::AvgResult<double> r = run(m, "amva.mapqn");
    CHECK(r.actualmethod == "amva.mapqn");
    // station 0 is the delay, station 1 the queue
    CHECK(r.TN(0, 0) == doctest::Approx(0.565971964).epsilon(1e-7));
    CHECK(r.TN(0, 1) == doctest::Approx(0.1934693276).epsilon(1e-7));
    CHECK(r.QN(1, 0) == doctest::Approx(0.8680560719).epsilon(1e-7));
    CHECK(r.QN(1, 1) == doctest::Approx(0.4522453793).epsilon(1e-7));
    CHECK(r.QN(0, 0) == doctest::Approx(1.1319439281).epsilon(1e-7));
    CHECK(r.QN(0, 1) == doctest::Approx(1.5477546207).epsilon(1e-7));
}

TEST_CASE("amva.mapqn reproduces the reference values on one class") {
    qn::Network<double> m = model({2.0}, {map_a()}, {4});
    const mva::AvgResult<double> r = run(m, "amva.mapqn");
    CHECK(r.TN(0, 0) == doctest::Approx(0.9717108597).epsilon(1e-7));
    CHECK(r.QN(1, 0) == doctest::Approx(2.0565782806).epsilon(1e-7));
}

TEST_CASE("mapqn_amva called directly") {
    std::vector<Matrix<double>> D0{m2(-3.0, 0.5, 0.2, -0.4), m2(-1.8, 0.3, 0.6, -0.9)};
    std::vector<Matrix<double>> D1{m2(2.5, 0.0, 0.2, 0.0), m2(1.5, 0.0, 0.0, 0.3)};
    const mapqn::MapqnAmvaResult<double> r = mapqn::mapqn_amva<double>({0.5, 0.125}, D0, D1, {2, 2});
    CHECK(r.X[0] == doctest::Approx(0.565971964).epsilon(1e-7));
    CHECK(r.X[1] == doctest::Approx(0.1934693276).epsilon(1e-7));
    CHECK(r.Qq[0] == doctest::Approx(0.8680560719).epsilon(1e-7));
    CHECK(r.Qq[1] == doctest::Approx(0.4522453793).epsilon(1e-7));
}

TEST_CASE("amva.mapqn stays within the method accuracy of the exact chain") {
    qn::Network<double> m = model({2.0, 8.0}, {map_a(), map_b()}, {2, 2});
    const mva::AvgResult<double> a = run(m, "amva.mapqn");
    ctmc::CtmcOptions copt;
    const mva::AvgResult<double> e = ctmc::solver_ctmc_run_analyzer(m.get_struct(), copt);
    for (std::size_t i = 0; i < 2; ++i)
        for (std::size_t r = 0; r < 2; ++r) {
            CHECK(a.TN(i, r) == doctest::Approx(e.TN(i, r)).epsilon(0.10));
            CHECK(a.QN(i, r) == doctest::Approx(e.QN(i, r)).epsilon(0.10));
            CHECK(a.RN(i, r) == doctest::Approx(e.RN(i, r)).epsilon(0.10));
        }
}

TEST_CASE("amva.mapqn with exponential service is exact multiclass MVA") {
    qn::Network<double> m = model({2.0, 8.0}, {D::exp_rate(1.2), D::exp_rate(1.2)}, {3, 3});
    const mva::AvgResult<double> a = run(m, "amva.mapqn");
    const mva::AvgResult<double> e = run(m, "exact");
    for (std::size_t i = 0; i < 2; ++i)
        for (std::size_t r = 0; r < 2; ++r) {
            CHECK(a.TN(i, r) == doctest::Approx(e.TN(i, r)).epsilon(1e-9));
            CHECK(a.QN(i, r) == doctest::Approx(e.QN(i, r)).epsilon(1e-9));
            CHECK(a.UN(i, r) == doctest::Approx(e.UN(i, r)).epsilon(1e-9));
        }
}

TEST_CASE("amva.mapqn is offered on its shape only and refuses the rest") {
    qn::Network<double> ok = model({2.0, 8.0}, {map_a(), map_b()}, {2, 2});
    CHECK(listed(ok, "amva.mapqn"));
    CHECK(mva::mva_mapqn_reason(ok.get_struct(), "amva.mapqn").empty());
    qn::Network<double> ps = model({2.0}, {map_a()}, {2}, SchedStrategy::PS);
    CHECK_FALSE(listed(ps, "amva.mapqn"));
    CHECK_FALSE(mva::mva_mapqn_reason(ps.get_struct(), "amva.mapqn").empty());
    CHECK_THROWS(run(ps, "amva.mapqn"));
}
