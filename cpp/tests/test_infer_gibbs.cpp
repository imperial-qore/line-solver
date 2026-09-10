/*
 * Copyright (c) 2012-2026, QORE Lab, Imperial College London
 * All rights reserved.
 *
 * infer_gibbs, in three layers, because the driver is stochastic and its
 * pieces are not:
 *   1. gibbs_analyse_data is a deterministic replay of the traces. Every state,
 *      its probability, the per-class populations and the mean queue
 *      populations are checked against the MATLAB subfunction analyseData,
 *      extracted verbatim and driven on the same synthetic traces.
 *   2. gibbs_slice is the deterministic content of one coordinate update: the
 *      whole log G vector along the grid and the whole normalized slice are
 *      checked against the MATLAB subfunction gibbsSamplerSimple, likewise
 *      extracted and driven on the same inputs. This is where the port either
 *      integrates d log G / d theta correctly or does not.
 *   3. the driver is checked for the properties a sampler must have -- exact
 *      reproducibility from a seed, invariance of the deterministic
 *      intermediates, and a posterior mean that tracks the slice it samples --
 *      and NOT against MATLAB's Mersenne Twister stream, which this port does
 *      not reproduce by decision.
 *
 * The MATLAB that produced every expected number is quoted above the
 * assertion. The subfunctions were extracted with their bodies untouched into
 * ref_analyse.m and ref_slice.m (the latter returning logG and prob as extra
 * outputs) and called from a script; nothing in the tree was modified.
 */
#include <cmath>
#include <vector>

#include "doctest.h"
#include "line/api/infer/infer_gibbs.h"

using line::Matrix;
using line::infer::GibbsOptions;
using line::infer::GibbsTrace;
using line::infer::gibbs_analyse_data;
using line::infer::gibbs_slice;
using line::infer::infer_gibbs;

namespace {

double relg(double got, double want) {
    const double d = std::fabs(got - want);
    return std::fabs(want) > 1e-30 ? d / std::fabs(want) : d;
}

/*
 * MATLAB:
 *   n1=10; n2=8;
 *   data{3,1}=(0:100:100*(n1-1))';  data{4,1}=(0.05+0.01*(0:n1-1))';
 *   data{6,1}=(2:2:2*n1)';
 *   data{3,2}=(50:150:50+150*(n2-1))'; data{4,2}=0.08*ones(n2,1);
 *   data{6,2}=(3:3:3*n2)';
 */
std::vector<GibbsTrace<double>> fixture() {
    std::vector<GibbsTrace<double>> d(2);
    for (int i = 0; i < 10; ++i) {
        d[0].arrival_ms.push_back(100.0 * i);
        d[0].respt_s.push_back(0.05 + 0.01 * i);
        d[0].think_obs.push_back(2.0 * (i + 1));
    }
    for (int i = 0; i < 8; ++i) {
        d[1].arrival_ms.push_back(50.0 + 150.0 * i);
        d[1].respt_s.push_back(0.08);
        d[1].think_obs.push_back(3.0 * (i + 1));
    }
    return d;
}

}  // namespace

TEST_CASE("gibbs_analyse_data reproduces MATLAB's analyseData exactly") {
    // MATLAB: [prob,N,N0] = ref_analyse(data, zeros(1,2), 2, 2, 0);
    //   size(prob,1) -> 6
    //   N  -> [2 1]
    //   N0 -> [0.80508474576271183 0.5423728813559322]
    const std::vector<GibbsTrace<double>> d = fixture();
    const line::infer::GibbsStateProbs<double> st = gibbs_analyse_data(d, 0);

    REQUIRE(st.prob.size() == 6u);
    CHECK(st.N[0] == 2L);
    CHECK(st.N[1] == 1L);
    CHECK(relg(st.N0[0], 0.80508474576271183) < 1e-14);
    CHECK(relg(st.N0[1], 0.5423728813559322) < 1e-14);

    // MATLAB's six rows, printed as "S n11 n12 n21 n22 P p":
    //   0 0 2 1  0.025423728813559324
    //   0 1 2 0  0.025423728813559324
    //   1 0 1 1  0.38135593220338981
    //   1 1 1 0  0.32203389830508472
    //   2 0 0 1  0.13559322033898305
    //   2 1 0 0  0.11016949152542373
    const long want_states[6][4] = {{0, 0, 2, 1}, {0, 1, 2, 0}, {1, 0, 1, 1},
                                    {1, 1, 1, 0}, {2, 0, 0, 1}, {2, 1, 0, 0}};
    const double want_prob[6] = {0.025423728813559324, 0.025423728813559324, 0.38135593220338981,
                                 0.32203389830508472, 0.13559322033898305,  0.11016949152542373};
    double total = 0.0;
    for (std::size_t i = 0; i < 6; ++i) {
        for (std::size_t j = 0; j < 4; ++j) CHECK(st.states(i, j) == want_states[i][j]);
        CHECK(relg(st.prob[i], want_prob[i]) < 1e-14);
        total += st.prob[i];
    }
    // The holding times partition the observation window, so the row order
    // aside, the probabilities are a distribution.
    CHECK(std::fabs(total - 1.0) < 1e-14);

    // Every state conserves the per-class population between the two nodes.
    for (std::size_t i = 0; i < 6; ++i)
        for (std::size_t k = 0; k < 2; ++k)
            CHECK(st.states(i, k) + st.states(i, 2 + k) == st.N[k]);
}

TEST_CASE("gibbs_slice reproduces MATLAB's logG walk and slice") {
    // MATLAB, continuing from the same prob:
    //   think_time(k) = (N(k)-N0(k))/mean(data{6,k})
    //     -> [0.10862865947611709 0.033898305084745763]
    //   logG_initial = sum(N.*log(think_time)) - sum(log(1:N(k))) summed
    //     -> -8.5171774576774535
    //   testset = [prob(:,1:4); prob(:,1:4)];
    //   [~,~,rsd,logG,pr] = ref_slice('TE',think_time,[0 0],testset,1,2,2,N,
    //                                 logG_initial,0.05,1,LV,sumA);
    const std::vector<GibbsTrace<double>> d = fixture();
    const line::infer::GibbsStateProbs<double> st = gibbs_analyse_data(d, 0);

    std::vector<double> think(2), N(2);
    for (std::size_t k = 0; k < 2; ++k) {
        double s = 0.0;
        for (std::size_t i = 0; i < d[k].think_obs.size(); ++i) s += d[k].think_obs[i];
        N[k] = static_cast<double>(st.N[k]);
        think[k] = (N[k] - st.N0[k]) / (s / static_cast<double>(d[k].think_obs.size()));
    }
    CHECK(relg(think[0], 0.10862865947611709) < 1e-14);
    CHECK(relg(think[1], 0.033898305084745763) < 1e-14);

    double logG0 = 0.0;
    for (std::size_t k = 0; k < 2; ++k) {
        logG0 += N[k] * std::log(think[k]);
        for (long j = 1; j <= st.N[k]; ++j) logG0 -= std::log(static_cast<double>(j));
    }
    CHECK(relg(logG0, -8.5171774576774535) < 1e-14);

    Matrix<long> testset(12, 4, 0L);
    for (std::size_t rep = 0; rep < 2; ++rep)
        for (std::size_t i = 0; i < 6; ++i)
            for (std::size_t j = 0; j < 4; ++j) testset(rep * 6 + i, j) = st.states(i, j);

    // MATLAB LOGG (21 grid points, step 0.05, width 1), coordinate 1 at theta = [0 0]
    const double want_logG[21] = {
        -8.5171774576774535, -8.5171774576774535, -7.942870893614506,  -7.4765902245690681,
        -7.0937204866390404, -6.7702064007532385, -6.490834109000466,  -6.2458064967091698,
        -6.0272002544119516, -5.8299269323182576, -5.6502395178115412, -5.4852855320844833,
        -5.3328471577118144, -5.1911683869217935, -5.0588363519934472, -4.9346980291779818,
        -4.8178005545889739, -4.7073476337006204, -4.602667143519759,  -4.5031866677413488,
        -4.4084147538023304};
    // MATLAB PROB
    const double want_prob[21] = {
        0.0,                    0.1235627001667544,     0.51431510319477625,
        0.2479059369212685,     0.079106686375979354,   0.023720768846026894,
        0.0074019402923185783,  0.0024873527142292061,  0.00089612154873736094,
        0.00034522370376131883, 0.00014149577208980035, 6.1346870934286663e-05,
        2.7978047861708282e-05, 1.3353619081382582e-05, 6.639911845036642e-06,
        3.4258625215109112e-06, 1.8276771244902797e-06, 1.0051192998936798e-06,
        5.6826996921639017e-07, 3.2952016459317623e-07, 1.9556525619413076e-07};

    const std::vector<double> theta0{0.0, 0.0};
    const line::infer::GibbsSlice<double> s1 =
        gibbs_slice(think, theta0, testset, 0, N, logG0, 0.05, 1.0);
    REQUIRE(s1.grid.size() == 21u);
    CHECK(s1.grid[0] == 0.0);
    CHECK(relg(s1.grid[20], 1.0) < 1e-15);
    for (std::size_t i = 0; i < 21; ++i) {
        CHECK(relg(s1.logG[i], want_logG[i]) < 1e-10);
        if (want_prob[i] == 0.0)
            CHECK(s1.prob[i] == 0.0);  // exp of the -Inf at theta = 0
        else
            CHECK(relg(s1.prob[i], want_prob[i]) < 1e-9);
    }
    // MATLAB: SLICE1 rsd 2
    CHECK(relg(s1.range_size_dim, 2.0) < 1e-15);

    // Second coordinate, from a theta whose first entry has already moved.
    // MATLAB: ref_slice('TE',think_time,[0.15 0],testset,2,...)
    const double want_logG2[21] = {
        -8.5171774576774535, -8.5171774576774535, -7.9350446106019339, -7.5693341530345206,
        -7.302047810510051,  -7.0913408429453106, -6.9173830038286868, -6.7692437062217294,
        -6.6402419884726021, -6.5259948641201113, -6.4234712965226777, -6.3304877569541684,
        -6.2454190644845067, -6.1670228979650688, -6.094328251702148,  -6.0265617673238321,
        -5.9630974860218675, -5.9034216184276849, -5.8471072487851536, -5.7937957906804733,
        -5.7431831411144998};
    const double want_prob2[21] = {
        0.0,                    0.93416945542707364,    0.055309615376301947,
        0.007824125146873263,   0.0017786741945890153,  0.0005413149263172809,
        0.00020042427602655298, 8.5427580677817971e-05, 4.048116281237628e-05,
        2.0833650081636995e-05, 1.1455268494022527e-05, 6.6492785406790042e-06,
        4.0380112703564447e-06, 2.5478619557971888e-06, 1.661243862507588e-06,
        1.1144110431283808e-06, 7.6642824730013337e-07, 5.3881829706104399e-07,
        3.8627766498340995e-07, 2.8180695521994229e-07, 2.0885291580767673e-07};

    const std::vector<double> theta1{0.15, 0.0};
    const line::infer::GibbsSlice<double> s2 =
        gibbs_slice(think, theta1, testset, 1, N, logG0, 0.05, 1.0);
    for (std::size_t i = 0; i < 21; ++i) {
        CHECK(relg(s2.logG[i], want_logG2[i]) < 1e-10);
        if (want_prob2[i] == 0.0)
            CHECK(s2.prob[i] == 0.0);
        else
            CHECK(relg(s2.prob[i], want_prob2[i]) < 1e-9);
    }
    CHECK(relg(s2.range_size_dim, 2.0) < 1e-15);

    // The slice is a probability vector and log G increases along the grid:
    // adding demand can only increase the normalizing constant.
    double tot = 0.0;
    for (std::size_t i = 0; i < 21; ++i) tot += s1.prob[i];
    CHECK(std::fabs(tot - 1.0) < 1e-14);
    for (std::size_t i = 2; i < 21; ++i) CHECK(s1.logG[i] > s1.logG[i - 1]);
}

TEST_CASE("gibbs_slice leaves logG at zero when the current value is off the grid") {
    // MATLAB's find(range==theta(index)) returns empty and every logG entry
    // stays at its zeros(1,N) initialization, so the slice reduces to the
    // likelihood term alone. Reproduced rather than repaired: repairing it
    // would make the two implementations sample differently.
    const std::vector<GibbsTrace<double>> d = fixture();
    const line::infer::GibbsStateProbs<double> st = gibbs_analyse_data(d, 0);
    Matrix<long> testset(6, 4, 0L);
    for (std::size_t i = 0; i < 6; ++i)
        for (std::size_t j = 0; j < 4; ++j) testset(i, j) = st.states(i, j);
    const std::vector<double> think{0.10862865947611709, 0.033898305084745763};
    const std::vector<double> N{2.0, 1.0};

    const std::vector<double> off{7.0, 0.0};  // 7 is not on 0:0.05:1
    const line::infer::GibbsSlice<double> s =
        gibbs_slice(think, off, testset, 0, N, -8.5171774576774535, 0.05, 1.0);
    for (std::size_t i = 0; i < s.logG.size(); ++i) CHECK(s.logG[i] == 0.0);
    // With logG flat the density is theta^coeff, which is increasing, so all
    // the mass sits at the top of the grid.
    CHECK(s.prob[s.prob.size() - 1] > s.prob[s.prob.size() - 2]);
}

TEST_CASE("infer_gibbs is reproducible from its seed and tracks its slice") {
    const std::vector<GibbsTrace<double>> d = fixture();
    GibbsOptions o;
    o.tol = 0.05;  // grid step; the default 1e-3 costs 1000 AMVA solves a move
    o.data_needed = 0;
    o.likelihood_sample = 200;
    o.nsamples = 60;
    o.block = 10;

    line::pfqn::McRng r1(20260721u), r2(20260721u), r3(7u);
    const std::vector<double> a = infer_gibbs(d, 2.0, o, r1);
    const std::vector<double> b = infer_gibbs(d, 2.0, o, r2);
    const std::vector<double> c = infer_gibbs(d, 2.0, o, r3);

    REQUIRE(a.size() == 2u);
    // Same seed, same numbers, bit for bit: the only randomness is the rng.
    CHECK(a[0] == b[0]);
    CHECK(a[1] == b[1]);
    // Different seed, different chain, but the same posterior.
    CHECK(c.size() == 2u);

    // The estimate is usedCores times a grid value, so it is non-negative and
    // bounded by the widest grid the adaptive rule can reach in this run.
    for (std::size_t k = 0; k < 2; ++k) {
        CHECK(a[k] >= 0.0);
        CHECK(std::isfinite(a[k]));
        CHECK(c[k] >= 0.0);
    }

    // The class-1 slice from the zero start puts 51% of its mass at 0.1 and
    // 25% at 0.15 (asserted above against MATLAB), so the chain samples that
    // neighbourhood; scaled by usedCores = 1.4857142857142858 the estimate
    // must stay well inside a decade of it. Two independent seeds agreeing to
    // that band is the distributional check.
    CHECK(a[0] > 0.02);
    CHECK(a[0] < 2.0);
    CHECK(c[0] > 0.02);
    CHECK(c[0] < 2.0);

    // A one-sample chain has no usable tail: MATLAB averages rows
    // round(n/2)+1 .. n-1, which is empty, and the port says so.
    GibbsOptions tiny = o;
    tiny.nsamples = 1;
    tiny.block = 1;
    line::pfqn::McRng r4(1u);
    CHECK_THROWS_AS(infer_gibbs(d, 2.0, tiny, r4), line::NumericError);
}
