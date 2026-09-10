/*
 * Copyright (c) 2012-2026, QORE Lab, Imperial College London
 * All rights reserved.
 */

/**
 * Mean busy period of order n for a subnetwork (Daduna, J. ACM 35(3), 1988).
 *
 * The oracle is the definition itself, not another implementation: Keilson's
 * mean ergodic sojourn time on the set G = {sum_{i in I} n_i >= n} is
 * pi(G) / h(G -> B), the stationary probability of G over the ergodic flow rate
 * out of it, and the test builds the Gordon-Newell generator directly and reads
 * both quantities off its stationary vector. The formula and the oracle share
 * nothing but the model, so an error in either is visible.
 *
 * Cases: a three-station closed network load-independent and load-dependent
 * (a two-server station and an infinite server), and an open M/M/1 whose busy
 * period of EVERY order is 1/(mu-lambda) in closed form.
 *
 * The MATLAB, Java and Python twins of this test agree with the same oracle to
 * 1e-14, so the tolerance here is a rounding tolerance, not a model tolerance.
 */

#include <cmath>
#include <map>
#include <vector>

#include "doctest.h"
#include "line/api/mc/ctmc_solve.h"
#include "line/api/pfqn/pfqn_busyp.h"
#include "line/api/pfqn/pfqn_busyp_clw.h"
#include "line/api/pfqn/pfqn_busyp_multiclass.h"
#include "line/lang/qn/network_builder.h"
#include "line/solvers/nc/solver_nc_busyp.h"
#include "line/util/matrix.h"

using line::Matrix;
using line::pfqn::pfqn_busyp;

namespace {

using State = std::vector<int>;

/** Every composition of `total` into `J` non-negative parts. */
std::vector<State> compositions(int total, int J) {
    if (J == 1) return {State{total}};
    std::vector<State> out;
    for (int k = 0; k <= total; ++k) {
        for (const State& rest : compositions(total - k, J - 1)) {
            State s{k};
            s.insert(s.end(), rest.begin(), rest.end());
            out.push_back(s);
        }
    }
    return out;
}

/** Gordon-Newell generator over the compositions of N, with rates mu(j, n_j). */
Matrix<double> gn_generator(const std::vector<std::vector<double>>& P,
                            const std::vector<std::vector<double>>& mu,
                            const std::vector<State>& states) {
    const std::size_t J = P.size(), S = states.size();
    std::map<State, std::size_t> index;
    for (std::size_t s = 0; s < S; ++s) index[states[s]] = s;
    Matrix<double> Q(S, S, 0.0);
    for (std::size_t s = 0; s < S; ++s) {
        for (std::size_t j = 0; j < J; ++j) {
            if (states[s][j] == 0) continue;
            const double rate = mu[j][static_cast<std::size_t>(states[s][j]) - 1];
            for (std::size_t i = 0; i < J; ++i) {
                if (P[j][i] <= 0) continue;
                State nst = states[s];
                nst[j] -= 1;
                nst[i] += 1;
                Q(s, index[nst]) = Q(s, index[nst]) + rate * P[j][i];
            }
        }
    }
    for (std::size_t s = 0; s < S; ++s) {
        double row = 0.0;
        for (std::size_t t = 0; t < S; ++t)
            if (t != s) row += Q(s, t);
        Q(s, s) = -row;
    }
    return Q;
}

/** Keilson's mean ergodic sojourn time on {sum_{i in I} n_i >= n}. */
double busy_period_ctmc(const std::vector<double>& pi, const std::vector<State>& states,
                        const std::vector<std::vector<double>>& P,
                        const std::vector<std::vector<double>>& mu,
                        const std::vector<std::size_t>& I, int n) {
    double piG = 0.0, flow = 0.0;
    for (std::size_t s = 0; s < states.size(); ++s) {
        int held = 0;
        for (std::size_t t = 0; t < I.size(); ++t) held += states[s][I[t]];
        if (held < n) continue;
        piG += pi[s];
        if (held != n) continue;
        for (std::size_t t = 0; t < I.size(); ++t) {
            const std::size_t j = I[t];
            if (states[s][j] == 0) continue;
            double out_of_I = 1.0;
            for (std::size_t u = 0; u < I.size(); ++u) out_of_I -= P[j][I[u]];
            flow += pi[s] * mu[j][static_cast<std::size_t>(states[s][j]) - 1] * out_of_I;
        }
    }
    return piG / flow;
}

/** Stochastic solution of x*P = x, the visit ratios of a closed network. */
std::vector<double> visits(const std::vector<std::vector<double>>& P) {
    const std::size_t J = P.size();
    std::vector<double> x(J, 1.0 / static_cast<double>(J));
    for (int it = 0; it < 20000; ++it) {
        std::vector<double> y(J, 0.0);
        for (std::size_t i = 0; i < J; ++i)
            for (std::size_t j = 0; j < J; ++j) y[j] += x[i] * P[i][j];
        double s = 0.0;
        for (std::size_t j = 0; j < J; ++j) s += y[j];
        for (std::size_t j = 0; j < J; ++j) y[j] /= s;
        x.swap(y);
    }
    return x;
}

Matrix<double> as_matrix(const std::vector<std::vector<double>>& v) {
    Matrix<double> m(v.size(), v[0].size(), 0.0);
    for (std::size_t i = 0; i < v.size(); ++i)
        for (std::size_t j = 0; j < v[0].size(); ++j) m(i, j) = v[i][j];
    return m;
}

const std::vector<std::vector<double>> ROUTING = {
    {0.0, 0.6, 0.4}, {0.7, 0.0, 0.3}, {0.5, 0.5, 0.0}};

}  // namespace

TEST_CASE("pfqn_busyp closed load-independent matches the CTMC sojourn time") {
    const int N = 5;
    const std::vector<double> rate = {1.5, 0.9, 2.0};
    std::vector<std::vector<double>> mu(3, std::vector<double>(N, 0.0));
    for (std::size_t j = 0; j < 3; ++j)
        for (int k = 0; k < N; ++k) mu[j][static_cast<std::size_t>(k)] = rate[j];

    const std::vector<State> states = compositions(N, 3);
    const std::vector<double> pi = line::mc::ctmc_solve(gn_generator(ROUTING, mu, states));
    const std::vector<double> alpha = visits(ROUTING);
    const Matrix<double> P = as_matrix(ROUTING);
    const Matrix<double> muM = as_matrix(mu);

    const std::vector<std::vector<std::size_t>> subnets = {{0}, {1}, {2}, {0, 1}, {1, 2}};
    for (const std::vector<std::size_t>& I : subnets) {
        for (int n = 1; n <= N; ++n) {
            const double oracle = busy_period_ctmc(pi, states, ROUTING, mu, I, n);
            const double got =
                pfqn_busyp(alpha, muM, P, static_cast<double>(N), I,
                           std::vector<std::size_t>{static_cast<std::size_t>(n)})
                    .b[0];
            CHECK(got == doctest::Approx(oracle).epsilon(1e-10));
        }
    }
}

TEST_CASE("pfqn_busyp closed load-dependent matches the CTMC sojourn time") {
    const int N = 5;
    // station 0 has two servers, station 1 is an infinite server, station 2 is single
    std::vector<std::vector<double>> mu(3, std::vector<double>(N, 0.0));
    for (int k = 1; k <= N; ++k) {
        mu[0][static_cast<std::size_t>(k - 1)] = 1.5 * std::min(k, 2);
        mu[1][static_cast<std::size_t>(k - 1)] = 0.9 * k;
        mu[2][static_cast<std::size_t>(k - 1)] = 2.0;
    }

    const std::vector<State> states = compositions(N, 3);
    const std::vector<double> pi = line::mc::ctmc_solve(gn_generator(ROUTING, mu, states));
    const std::vector<double> alpha = visits(ROUTING);
    const Matrix<double> P = as_matrix(ROUTING);
    const Matrix<double> muM = as_matrix(mu);

    const std::vector<std::vector<std::size_t>> subnets = {{0}, {1}, {0, 2}};
    for (const std::vector<std::size_t>& I : subnets) {
        for (int n : {1, 3, 5}) {
            const double oracle = busy_period_ctmc(pi, states, ROUTING, mu, I, n);
            const double got =
                pfqn_busyp(alpha, muM, P, static_cast<double>(N), I,
                           std::vector<std::size_t>{static_cast<std::size_t>(n)})
                    .b[0];
            CHECK(got == doctest::Approx(oracle).epsilon(1e-10));
        }
    }
}

TEST_CASE("pfqn_busyp open M/M/1 busy period is 1/(mu-lambda) at every order") {
    const double lambda = 0.7, svc = 1.3;
    // node 0 is the queue, node 1 absorbs its output; the Source is not a node
    // of the Jackson network, its outflow being the external stream gamma
    const std::vector<double> alpha = {lambda, lambda};
    Matrix<double> mu(2, 1, 0.0);
    mu(0, 0) = svc;
    mu(1, 0) = 2.0;
    Matrix<double> P(2, 2, 0.0);
    P(0, 1) = 1.0;
    const std::vector<double> gamma = {lambda, 0.0};

    for (std::size_t n = 1; n <= 4; ++n) {
        const double got = pfqn_busyp(alpha, mu, P, std::numeric_limits<double>::infinity(),
                                      std::vector<std::size_t>{0},
                                      std::vector<std::size_t>{n}, gamma)
                               .b[0];
        CHECK(got == doctest::Approx(1.0 / (svc - lambda)).epsilon(1e-9));
    }
}

TEST_CASE("pfqn_busyp open tandem subnetwork matches the geometric closed form") {
    // both nodes of an open tandem: sum_{m>=1} G(m) = prod 1/(1-rho_j) - 1
    const double lambda = 0.5;
    const std::vector<double> alpha = {lambda, lambda};
    Matrix<double> mu(2, 1, 0.0);
    mu(0, 0) = 1.0;
    mu(1, 0) = 1.2;
    Matrix<double> P(2, 2, 0.0);
    P(0, 1) = 1.0;
    const std::vector<double> gamma = {lambda, 0.0};

    const double rho0 = lambda / 1.0, rho1 = lambda / 1.2;
    const double expected = (1.0 / ((1 - rho0) * (1 - rho1)) - 1.0) / lambda;
    const double got = pfqn_busyp(alpha, mu, P, std::numeric_limits<double>::infinity(),
                                  std::vector<std::size_t>{0, 1},
                                  std::vector<std::size_t>{1}, gamma)
                           .b[0];
    CHECK(got == doctest::Approx(expected).epsilon(1e-10));
}

TEST_CASE("pfqn_busyp rejects a closed subnetwork covering every node") {
    const std::vector<double> alpha = {0.5, 0.5};
    Matrix<double> mu(2, 1, 1.0);
    Matrix<double> P(2, 2, 0.0);
    P(0, 1) = 1.0;
    P(1, 0) = 1.0;
    CHECK_THROWS_AS(pfqn_busyp(alpha, mu, P, 3.0, std::vector<std::size_t>{0, 1},
                               std::vector<std::size_t>{1}),
                    line::InputError);
}

TEST_CASE("solver_nc_busyp reads the same numbers off a NetworkStruct") {
    // the solver-level entry must reproduce the API called by hand, which is
    // what the MATLAB and Python getters are checked against
    using Net = line::qn::Network<double>;
    using D = line::lang::Distrib<double>;
    using Routing = line::qn::RoutingMatrix<double>;
    using line::lang::SchedStrategy;

    const double N = 5.0;
    Net m("busyPeriodModel");
    const std::size_t q1 = m.add_queue("Queue1", SchedStrategy::FCFS);
    const std::size_t q2 = m.add_queue("Queue2", SchedStrategy::FCFS);
    const std::size_t q3 = m.add_queue("Queue3", SchedStrategy::FCFS);
    const std::size_t c1 = m.add_closed_class("Class1", N, q1, 0);
    m.set_service(q1, c1, D::exp_rate(1.5));
    m.set_service(q2, c1, D::exp_rate(0.9));
    m.set_service(q3, c1, D::exp_rate(2.0));
    Routing P;
    P.set(c1, c1, q1, q2, 0.6);
    P.set(c1, c1, q1, q3, 0.4);
    P.set(c1, c1, q2, q1, 0.7);
    P.set(c1, c1, q2, q3, 0.3);
    P.set(c1, c1, q3, q1, 0.5);
    P.set(c1, c1, q3, q2, 0.5);
    m.link(P);

    // the CTMC oracle of the first case, on the same rates and routing
    std::vector<std::vector<double>> mu(3, std::vector<double>(5, 0.0));
    const std::vector<double> rate = {1.5, 0.9, 2.0};
    for (std::size_t j = 0; j < 3; ++j)
        for (int k = 0; k < 5; ++k) mu[j][static_cast<std::size_t>(k)] = rate[j];
    const std::vector<State> states = compositions(5, 3);
    const std::vector<double> pi = line::mc::ctmc_solve(gn_generator(ROUTING, mu, states));

    const std::vector<std::vector<std::size_t>> subnets = {{0}, {1}, {2}, {0, 1}};
    for (const std::vector<std::size_t>& I : subnets) {
        const std::vector<double> b =
            line::nc::solver_nc_busyp(m.get_struct(), I, std::vector<std::size_t>{1, 3, 5});
        const int want[3] = {1, 3, 5};
        for (int t = 0; t < 3; ++t) {
            const double oracle = busy_period_ctmc(pi, states, ROUTING, mu, I, want[t]);
            CHECK(b[static_cast<std::size_t>(t)] ==
                  doctest::Approx(oracle).epsilon(1e-10));
        }
    }
}

TEST_CASE("pfqn_busyp_multiclass collapses to the single-chain formula at R=1") {
    // the R=1 collapse is the regression the generalization must pass: the inner
    // sum holds one term and H(N-m-e_1) = H(N-n), which is Theorem 1 verbatim
    const int N = 5;
    const std::vector<double> rate = {1.5, 0.9, 2.0};
    std::vector<std::vector<double>> mu(3, std::vector<double>(N, 0.0));
    for (std::size_t j = 0; j < 3; ++j)
        for (int k = 0; k < N; ++k) mu[j][static_cast<std::size_t>(k)] = rate[j];
    const std::vector<double> alpha = visits(ROUTING);
    const Matrix<double> P = as_matrix(ROUTING);
    const Matrix<double> muM = as_matrix(mu);

    Matrix<double> alphaCol(3, 1, 0.0), muCol(3, 1, 0.0);
    for (std::size_t i = 0; i < 3; ++i) {
        alphaCol(i, 0) = alpha[i];
        muCol(i, 0) = rate[i];
    }

    const std::vector<std::vector<std::size_t>> subnets = {{0}, {1}, {2}, {0, 1}};
    for (const std::vector<std::size_t>& I : subnets) {
        for (std::size_t n = 1; n <= 5; ++n) {
            const double single =
                pfqn_busyp(alpha, muM, P, static_cast<double>(N), I,
                           std::vector<std::size_t>{n})
                    .b[0];
            const double multi = line::pfqn::pfqn_busyp_multiclass(
                alphaCol, muCol, std::vector<Matrix<double>>{P},
                std::vector<double>{static_cast<double>(N)}, I,
                std::vector<std::size_t>{n})[0];
            CHECK(multi == doctest::Approx(single).epsilon(1e-12));
        }
    }
}

TEST_CASE("pfqn_busyp_multiclass matches a two-chain CTMC with class-dependent rates") {
    // oracle: the same multiclass PS generator the Python twin builds, whose
    // per-class count vector IS Markov (a class-r departure at rate mu_ir n_ir/|n_i|)
    Matrix<double> alpha(3, 2, 0.0), mu(3, 2, 0.0);
    const std::vector<double> a = visits(ROUTING);
    const double muc[3][2] = {{1.5, 2.5}, {0.9, 0.7}, {2.0, 1.1}};
    for (std::size_t i = 0; i < 3; ++i)
        for (std::size_t r = 0; r < 2; ++r) {
            alpha(i, r) = a[i];
            mu(i, r) = muc[i][r];
        }
    const Matrix<double> P = as_matrix(ROUTING);

    // values verified against the directly built CTMC in the Python twin
    const double oracle[4][4] = {
        {0.9126174546, 0.8012355212, 0.6654673220, 0.5000000000},
        {11.1988414467, 5.5361633551, 2.7147367636, 1.2500000000},
        {1.0651084429, 0.9691585895, 0.8326283551, 0.6451612903},
        {81.3803977003, 20.7023047321, 6.4800653776, 1.7771304082}};
    const std::vector<std::vector<std::size_t>> subnets = {{0}, {1}, {2}, {0, 1}};
    for (std::size_t s = 0; s < subnets.size(); ++s) {
        for (std::size_t n = 1; n <= 4; ++n) {
            const double got = line::pfqn::pfqn_busyp_multiclass(
                alpha, mu, std::vector<Matrix<double>>{P}, std::vector<double>{2.0, 2.0},
                subnets[s], std::vector<std::size_t>{n})[0];
            CHECK(got == doctest::Approx(oracle[s][n - 1]).epsilon(1e-9));
        }
    }
}

TEST_CASE("pfqn_busyp_multiclass per class and mixed match the CTMC") {
    // Oracles built by the Python twin against a directly constructed multiclass
    // CTMC: a processor-sharing generator, whose per-class count vector IS Markov.
    Matrix<double> alpha(3, 2, 0.0), mu(3, 2, 0.0);
    const std::vector<double> a = visits(ROUTING);
    const double muc[3][2] = {{1.5, 2.5}, {0.9, 0.7}, {2.0, 1.1}};
    for (std::size_t i = 0; i < 3; ++i)
        for (std::size_t r = 0; r < 2; ++r) {
            alpha(i, r) = a[i];
            mu(i, r) = muc[i][r];
        }
    const Matrix<double> P = as_matrix(ROUTING);

    SUBCASE("per class, closed N=(3,2)") {
        const double oracle[4][3] = {
            {1.1841553940, 1.0400905281, 0.8165028553},
            {0.7554836027, 0.6747297692, 0.0},
            {11.6096377218, 4.4280317562, 1.7178235965},
            {10.5287171490, 2.9841721108, 0.0}};
        const std::size_t Nr[2] = {3, 2};
        std::size_t row = 0;
        for (std::size_t s = 0; s < 2; ++s)
            for (std::size_t r = 0; r < 2; ++r, ++row)
                for (std::size_t n = 1; n <= Nr[r]; ++n) {
                    const double got = line::pfqn::pfqn_busyp_multiclass(
                        alpha, mu, std::vector<Matrix<double>>{P},
                        std::vector<double>{3.0, 2.0}, std::vector<std::size_t>{s},
                        std::vector<std::size_t>{n}, Matrix<double>(), Matrix<double>(),
                        line::pfqn::PFQN_BUSYP_DEFAULT_TOL, static_cast<int>(r))[0];
                    CHECK(got == doctest::Approx(oracle[row][n - 1]).epsilon(1e-9));
                }
    }

    SUBCASE("mixed: one closed chain and one open chain") {
        // the closed chain cycles between the two stations, the open one crosses them
        Matrix<double> alphaM(2, 2, 0.0), muM(2, 2, 0.0), gammaM(2, 2, 0.0);
        Matrix<double> Pcl(2, 2, 0.0), Pop(2, 2, 0.0);
        const double lam = 0.25;
        const double rates[2][2] = {{1.0, 1.3}, {0.8, 1.1}};
        Pcl(0, 1) = 1.0;
        Pcl(1, 0) = 1.0;
        Pop(0, 1) = 1.0;
        for (std::size_t i = 0; i < 2; ++i) {
            alphaM(i, 0) = 1.0;
            alphaM(i, 1) = lam;
            muM(i, 0) = rates[i][0];
            muM(i, 1) = rates[i][1];
        }
        gammaM(0, 1) = lam;
        const double oracle[6][2] = {
            {2.1856936714, 1.2380952381}, {1.8209481724, 1.5991887456},
            {2.2011001800, 1.5639339391}, {3.7312017567, 1.6176470588},
            {2.8593692142, 2.2706104767}, {3.9653801317, 2.2604782254}};
        const int classes[3] = {0, 1, -1};
        std::size_t row = 0;
        for (std::size_t s = 0; s < 2; ++s)
            for (std::size_t c = 0; c < 3; ++c, ++row)
                for (std::size_t n = 1; n <= 2; ++n) {
                    const double got = line::pfqn::pfqn_busyp_multiclass(
                        alphaM, muM, std::vector<Matrix<double>>{Pcl, Pop},
                        std::vector<double>{2.0, std::numeric_limits<double>::infinity()},
                        std::vector<std::size_t>{s}, std::vector<std::size_t>{n}, gammaM,
                        Matrix<double>(), line::pfqn::PFQN_BUSYP_DEFAULT_TOL, classes[c])[0];
                    CHECK(got == doctest::Approx(oracle[row][n - 1]).epsilon(1e-8));
                }
    }
}

TEST_CASE("pfqn_busyp_clw reproduces the ladder from point evaluations") {
    // the CLW form replaces the ladder by G(N) minus the n lowest shells, so the
    // regression is against the exact routines it is meant to replace
    const std::vector<double> rate = {1.5, 0.9, 2.0};
    const std::vector<double> alpha = visits(ROUTING);
    const Matrix<double> P = as_matrix(ROUTING);
    Matrix<double> alphaCol(3, 1, 0.0), muCol(3, 1, 0.0);
    for (std::size_t i = 0; i < 3; ++i) {
        alphaCol(i, 0) = alpha[i];
        muCol(i, 0) = rate[i];
    }

    SUBCASE("closed single chain, against pfqn_busyp") {
        for (int N : {5, 20, 60}) {
            std::vector<std::vector<double>> mu(3, std::vector<double>(N, 0.0));
            for (std::size_t j = 0; j < 3; ++j)
                for (int k = 0; k < N; ++k) mu[j][static_cast<std::size_t>(k)] = rate[j];
            const Matrix<double> muM = as_matrix(mu);
            for (const std::vector<std::size_t>& I :
                 std::vector<std::vector<std::size_t>>{{0}, {0, 1}})
                for (std::size_t n = 1; n <= 3; ++n) {
                    const double exact =
                        pfqn_busyp(alpha, muM, P, static_cast<double>(N), I,
                                   std::vector<std::size_t>{n})
                            .b[0];
                    const double got = line::pfqn::pfqn_busyp_clw(
                        alphaCol, muCol, std::vector<Matrix<double>>{P},
                        std::vector<double>{static_cast<double>(N)}, I,
                        std::vector<std::size_t>{n})[0];
                    CHECK(got == doctest::Approx(exact).epsilon(1e-8));
                }
        }
    }

    SUBCASE("closed two chains, against the lattice routine") {
        Matrix<double> alpha2(3, 2, 0.0), mu2(3, 2, 0.0);
        const double muc[3][2] = {{1.5, 2.5}, {0.9, 0.7}, {2.0, 1.1}};
        for (std::size_t i = 0; i < 3; ++i)
            for (std::size_t r = 0; r < 2; ++r) {
                alpha2(i, r) = alpha[i];
                mu2(i, r) = muc[i][r];
            }
        for (const std::vector<double>& N :
             std::vector<std::vector<double>>{{2.0, 2.0}, {4.0, 3.0}})
            for (const std::vector<std::size_t>& I :
                 std::vector<std::vector<std::size_t>>{{0}, {0, 1}})
                for (std::size_t n = 1; n <= 2; ++n) {
                    const double exact = line::pfqn::pfqn_busyp_multiclass(
                        alpha2, mu2, std::vector<Matrix<double>>{P}, N, I,
                        std::vector<std::size_t>{n})[0];
                    const double got = line::pfqn::pfqn_busyp_clw(
                        alpha2, mu2, std::vector<Matrix<double>>{P}, N, I,
                        std::vector<std::size_t>{n})[0];
                    CHECK(got == doctest::Approx(exact).epsilon(1e-6));
                }
    }

    SUBCASE("open: the tail is exact, not truncated") {
        const double lam = 0.5;
        Matrix<double> ao(2, 1, 0.0), ro(2, 1, 0.0), go(2, 1, 0.0), Po(2, 2, 0.0);
        ao(0, 0) = lam;
        ao(1, 0) = lam;
        ro(0, 0) = 1.0;
        ro(1, 0) = 1.2;
        go(0, 0) = lam;
        Po(0, 1) = 1.0;
        // both nodes: sum_{m>=1} G(m) = prod 1/(1-rho_j) - 1, in closed form
        const double rho0 = lam / 1.0, rho1 = lam / 1.2;
        const double expected = (1.0 / ((1 - rho0) * (1 - rho1)) - 1.0) / lam;
        const double got = line::pfqn::pfqn_busyp_clw(
            ao, ro, std::vector<Matrix<double>>{Po},
            std::vector<double>{std::numeric_limits<double>::infinity()},
            std::vector<std::size_t>{0, 1}, std::vector<std::size_t>{1}, go)[0];
        CHECK(got == doctest::Approx(expected).epsilon(1e-12));
    }
}
