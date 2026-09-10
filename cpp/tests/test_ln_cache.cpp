/*
 * Copyright (c) 2012-2026, QORE Lab, Imperial College London
 * All rights reserved.
 */

/**
 * A cache task: an item lookup that switches the job onto a hit or a miss.
 *
 * A closed client calls the item entry of a CacheTask holding 2 of 4 items
 * under LRU. The read itself is Immediate; the work is on the two branches,
 * a fast hit (0.2) and a slow miss (1.0), so the entry's service time is the
 * hit-ratio mixture of the two and the whole model turns on that ratio.
 *
 * The Cache NODE belongs to the cache task's HOST layer, not to the task's own
 * layer (`iscachelayer` is a host-layer test, buildLayersRecursive.m:88): the
 * lookup is what the task does on its processor, so hit and miss have to queue
 * at the server the read arrived at.
 *
 * Reference numbers are MATLAB SolverLN.getAvgTable() on the identical model
 * (defaultOptions, 2026-07-29), where the popularity is Zipf(1.0, 4).
 */

#include <limits>
#include <string>
#include <vector>

#include "doctest.h"
#include "line/lang/lqn/lqn_builder.h"
#include "line/solvers/ln/solver_ln.h"

using namespace line;
using lang::Distrib;
using lang::ReplacementStrategy;
using lang::SchedStrategy;

namespace {

using D = Distrib<double>;
D E(double m) { return D::exp_mean(m); }

/** Zipf(1.0, n) as an explicit pmf: p_k proportional to 1/k. */
std::vector<double> zipf(std::size_t n) {
    std::vector<double> p(n, 0.0);
    double norm = 0.0;
    for (std::size_t k = 0; k < n; ++k) norm += 1.0 / double(k + 1);
    for (std::size_t k = 0; k < n; ++k) p[k] = (1.0 / double(k + 1)) / norm;
    return p;
}

lqn::LqnStruct<double> build_cache() {
    lqn::LqnBuilder<double> b;
    b.processor("P1", std::numeric_limits<double>::infinity(), SchedStrategy::INF);
    b.processor("P2", 1, SchedStrategy::PS);

    b.task("T1", 2, SchedStrategy::REF, "P1");
    b.think_time("T1", E(1.0));
    b.cache_task("CT", 1, SchedStrategy::FCFS, "P2", 4, {2}, ReplacementStrategy::LRU);

    b.entry("E1", "T1");
    b.item_entry("IE", "CT", 4, zipf(4));

    b.activity("A1", E(0.5), "T1");
    b.bound_to("A1", "E1");
    b.sync_call("A1", "IE", 1.0);

    b.activity("AC", D::immediate(), "CT");
    b.bound_to("AC", "IE");
    b.activity("AH", E(0.2), "CT");
    b.activity("AM", E(1.0), "CT");
    b.cache_access("AC", "AH", "AM");
    b.replies_to("AH", "IE");
    b.replies_to("AM", "IE");

    return b.build();
}

std::size_t idx_of(const lqn::LqnStruct<double>& l, const std::string& hashname) {
    for (std::size_t i = 1; i <= l.nidx; ++i)
        if (l.hashnames[i] == hashname) return i;
    return 0;
}

}  // namespace

TEST_CASE("cache: the task and the item entry reach the struct") {
    const lqn::LqnStruct<double> l = build_cache();
    // MATLAB hashnames the cache task C: and the item entry I:
    const std::size_t ct = idx_of(l, "C:CT"), ie = idx_of(l, "I:IE");
    REQUIRE(ct > 0);
    REQUIRE(ie > 0);
    CHECK(l.iscache[ct]);
    // nitems is carried on BOTH, as getStruct.m writes it
    CHECK(l.nitems[ct] == 4);
    CHECK(l.nitems[ie] == 4);
    REQUIRE(l.itemcap[ct].size() == 1);
    CHECK(l.itemcap[ct][0] == 2);
    CHECK(l.replacestrat[ct] == ReplacementStrategy::LRU);
    REQUIRE(l.itemproc[ie].size() == 4);
    CHECK(l.itemproc[ie][0] == doctest::Approx(0.48).epsilon(1e-2));
    // POST_CACHE marks the BRANCHES, not the read
    CHECK(l.actposttype[idx_of(l, "A:AH")] == lang::PrecedenceType::POST_CACHE);
    CHECK(l.actposttype[idx_of(l, "A:AM")] == lang::PrecedenceType::POST_CACHE);
    CHECK(l.actposttype[idx_of(l, "A:AC")] != lang::PrecedenceType::POST_CACHE);
}

TEST_CASE("cache: the Cache node is built in the HOST layer") {
    const lqn::LqnStruct<double> l = build_cache();
    ln::LnOptions opt;
    ln::SolverLN<double> s(l, opt);
    CHECK(s.nlayers() == 3);

    std::size_t ncachenodes = 0;
    for (const qn::Layer<double>& L : s.layers()) {
        CAPTURE(L.name);
        std::size_t here = 0;
        for (std::size_t n = 0; n < L.nodes.size(); ++n)
            if (L.nodes[n].nodetype == qn::NodeType::Cache) ++here;
        if (here > 0) {
            // P:P2 hosts CT, so the Cache belongs to it -- NOT to layer C:CT
            CHECK(L.name == "P:P2");
            CHECK(L.nclasses == 5);  // C:CT, I:IE, A:AC, A:AH, A:AM
        }
        ncachenodes += here;
    }
    CHECK(ncachenodes == 1);
}

TEST_CASE("cache: the model matches MATLAB SolverLN") {
    const lqn::LqnStruct<double> l = build_cache();
    ln::LnOptions opt;
    ln::SolverLN<double> s(l, opt);
    const ln::LnSolution<double> sol = s.get_ensemble_avg();
    CHECK(sol.converged);

    // MATLAB AvgTable: P1 Util 0.45654, P2 Util 0.49772
    //   T1 QLen 1.0869 Tput 0.91308; CT QLen 0.49772 ResidT 0.5451
    //   E1 RespT 1.1904; IE RespT 0.5451 Tput 0.91308
    //   AH RespT 0.2 Tput 0.5192; AM RespT 1.0 Tput 0.39388
    struct Row { const char* hn; double util; double tput; double respt; };
    const Row rows[] = {
        {"P:P1", 0.45654, 0.0, 0.0},
        {"P:P2", 0.49772, 0.0, 0.0},
        {"R:T1", 0.45654, 0.91308, 0.0},
        {"C:CT", 0.49772, 0.91308, 0.0},
        {"E:E1", 0.45654, 0.91308, 1.1904},
        {"I:IE", 0.49772, 0.91308, 0.5451},
        {"A:AH", 0.10384, 0.5192, 0.2},
        {"A:AM", 0.39388, 0.39388, 1.0},
    };
    for (const Row& row : rows) {
        const std::size_t i = idx_of(l, row.hn);
        REQUIRE(i > 0);
        CAPTURE(std::string(row.hn));
        if (row.util > 1e-7) CHECK(sol.UN[i] == doctest::Approx(row.util).epsilon(1e-3));
        if (row.tput > 1e-7) CHECK(sol.TN[i] == doctest::Approx(row.tput).epsilon(1e-3));
        if (row.respt > 1e-7) CHECK(sol.RN[i] == doctest::Approx(row.respt).epsilon(1e-3));
    }

    // Every read ends on exactly one branch, so the two throughputs must add up
    // to the item entry's -- this is what a wrong hit/miss class wiring breaks.
    CHECK(sol.TN[idx_of(l, "A:AH")] + sol.TN[idx_of(l, "A:AM")] ==
          doctest::Approx(sol.TN[idx_of(l, "I:IE")]).epsilon(1e-3));
    // 2 of 4 items under LRU with a Zipf read: the hit ratio is well inside (0,1)
    const double hitratio = sol.TN[idx_of(l, "A:AH")] / sol.TN[idx_of(l, "I:IE")];
    CHECK(hitratio > 0.0);
    CHECK(hitratio < 1.0);
    CHECK(hitratio == doctest::Approx(0.5686).epsilon(1e-2));
}
