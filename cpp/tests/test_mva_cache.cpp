/*
 * Copyright (c) 2012-2026, QORE Lab, Imperial College London
 * All rights reserved.
 */

/**
 * The non-reentrant cache analyzer (ladder branch 8), a Source-Cache-Sink model.
 *
 * The hit ratio is what the analyzer exists to compute, and it is checkable
 * without either codebase for the top item: with list capacity 2 over 4 items
 * and a Zipf(1,4) read law, item 1 (popularity 0.48) is almost always cached.
 * The two numbers are MATLAB's `Cache.getHitRatio()`: 0.588571428571 from the
 * exact RR recursion (cache_mva) and 0.568623668307 from the fixed-point
 * approximation (cache_prob_fpi).
 */

#include <vector>

#include "doctest.h"
#include "line/lang/qn/network_builder.h"
#include "line/solvers/mva/solver_mva_cache.h"

using namespace line;
using lang::Distrib;
using lang::ReplacementStrategy;
using lang::SchedStrategy;

namespace {

using D = Distrib<double>;

qn::Network<double> cache_model() {
    qn::Network<double> m("cache");
    const std::size_t src = m.add_source("Source");
    qn::CacheParam<double> ch;
    ch.nitems = 4;
    ch.itemcap = std::vector<int>{2};  // one list of capacity 2
    ch.replacestrat = ReplacementStrategy::RR;
    // Zipf(1, 4): 1/(k H), H = 1 + 1/2 + 1/3 + 1/4 = 25/12
    ch.pread = std::vector<std::vector<double>>{{0.48, 0.24, 0.16, 0.12}};
    ch.hitclass = std::vector<std::size_t>{0};   // single class: no class switch
    ch.missclass = std::vector<std::size_t>{0};
    const std::size_t cnode = m.add_cache("Cache", ch);
    const std::size_t snk = m.add_sink("Sink");
    const std::size_t rd = m.add_open_class("Read");
    m.set_arrival(src, rd, D::exp_rate(1.0));
    qn::RoutingMatrix<double> P;
    P.set(src, cnode, 1.0);
    P.set(cnode, snk, 1.0);
    m.link(P);
    return m;
}

TEST_CASE("the exact RR cache recursion reproduces the reference hit ratio") {
    qn::Network<double> m = cache_model();
    mva::MvaOptions opt;
    opt.method = "exact";
    const mva::CacheResult<double> r = mva::solver_mva_cache_analyzer(m.get_struct(), opt);
    CHECK(r.actualmethod == "exact");
    CHECK(r.hitprob[0] == doctest::Approx(0.588571428571).epsilon(1e-9));
    CHECK(r.missprob[0] == doctest::Approx(1.0 - 0.588571428571).epsilon(1e-9));
    // the item-1 popularity dominates, so more than half the reads hit even at
    // capacity 2 out of 4
    CHECK(r.hitprob[0] > 0.5);
}

TEST_CASE("the fixed-point approximation reproduces the reference hit ratio") {
    qn::Network<double> m = cache_model();
    mva::MvaOptions opt;
    opt.method = "default";
    const mva::CacheResult<double> r = mva::solver_mva_cache_analyzer(m.get_struct(), opt);
    CHECK(r.actualmethod == "fpi");
    CHECK(r.hitprob[0] == doctest::Approx(0.568623668307).epsilon(1e-9));
    CHECK(r.missprob[0] == doctest::Approx(1.0 - 0.568623668307).epsilon(1e-9));
}

TEST_CASE("the LRU characteristic-time approximation reproduces the reference hit ratio") {
    qn::Network<double> m("cache-lru");
    const std::size_t src = m.add_source("Source");
    qn::CacheParam<double> ch;
    ch.nitems = 4;
    ch.itemcap = std::vector<int>{2};
    ch.replacestrat = ReplacementStrategy::LRU;
    ch.pread = std::vector<std::vector<double>>{{0.48, 0.24, 0.16, 0.12}};
    ch.hitclass = std::vector<std::size_t>{0};
    ch.missclass = std::vector<std::size_t>{0};
    const std::size_t cnode = m.add_cache("Cache", ch);
    const std::size_t snk = m.add_sink("Sink");
    const std::size_t rd = m.add_open_class("Read");
    m.set_arrival(src, rd, D::exp_rate(1.0));
    qn::RoutingMatrix<double> P;
    P.set(src, cnode, 1.0);
    P.set(cnode, snk, 1.0);
    m.link(P);
    mva::MvaOptions opt;
    const mva::CacheResult<double> r = mva::solver_mva_cache_analyzer(m.get_struct(), opt);
    CHECK(r.actualmethod == "ttl");
    // The LRU-A occupancy is a fixed point stopped on a tolerance in BOTH
    // codebases (the port at FineTol, the reference at its own random-seeded
    // one), so the two converge to points ~3e-7 apart; 1e-6 is the real
    // agreement bound for this doubly-approximate quantity, not 1e-9.
    CHECK(r.hitprob[0] == doctest::Approx(0.597228013344).epsilon(1e-6));
}

TEST_CASE("the h-LRU characteristic-time approximation reproduces the reference hit ratio") {
    qn::Network<double> m("cache-hlru");
    const std::size_t src = m.add_source("Source");
    qn::CacheParam<double> ch;
    ch.nitems = 4;
    ch.itemcap = std::vector<int>{2, 1};  // two lists
    ch.replacestrat = ReplacementStrategy::HLRU;
    ch.pread = std::vector<std::vector<double>>{{0.48, 0.24, 0.16, 0.12}};
    ch.hitclass = std::vector<std::size_t>{0};
    ch.missclass = std::vector<std::size_t>{0};
    const std::size_t cnode = m.add_cache("Cache", ch);
    const std::size_t snk = m.add_sink("Sink");
    const std::size_t rd = m.add_open_class("Read");
    m.set_arrival(src, rd, D::exp_rate(1.0));
    qn::RoutingMatrix<double> P;
    P.set(src, cnode, 1.0);
    P.set(cnode, snk, 1.0);
    m.link(P);
    mva::MvaOptions opt;
    const mva::CacheResult<double> r = mva::solver_mva_cache_analyzer(m.get_struct(), opt);
    CHECK(r.actualmethod == "ttl");
    // a second, smaller list catches more items, so h-LRU hits more than LRU
    CHECK(r.hitprob[0] == doctest::Approx(0.833725168654).epsilon(1e-9));
    CHECK(r.hitprob[0] > 0.597228013344);
}

TEST_CASE("a hit/miss class switch splits the source throughput") {
    // read -> hit / miss, checked through the analyzer's X vector: the hit-class
    // throughput is the hit fraction of the source rate, the miss-class the rest
    qn::Network<double> m("cache-split");
    const std::size_t src = m.add_source("Source");
    qn::CacheParam<double> ch;
    ch.nitems = 4;
    ch.itemcap = std::vector<int>{2};
    ch.replacestrat = ReplacementStrategy::RR;
    ch.pread = std::vector<std::vector<double>>{{0.48, 0.24, 0.16, 0.12}, {}, {}};
    ch.hitclass = std::vector<std::size_t>{2, 0, 0};   // read (class 1) -> hit (class 2)
    ch.missclass = std::vector<std::size_t>{3, 0, 0};  //                 -> miss (class 3)
    const std::size_t cnode = m.add_cache("Cache", ch);
    const std::size_t snk = m.add_sink("Sink");
    const std::size_t rd = m.add_open_class("Read");
    const std::size_t hit = m.add_open_class("Hit");
    const std::size_t miss = m.add_open_class("Miss");
    (void)hit;
    (void)miss;
    m.set_arrival(src, rd, D::exp_rate(1.0));
    qn::RoutingMatrix<double> P;
    P.set(rd, rd, src, cnode, 1.0);
    P.set(hit, hit, cnode, snk, 1.0);
    P.set(miss, miss, cnode, snk, 1.0);
    m.link(P);

    mva::MvaOptions opt;
    opt.method = "exact";
    const mva::CacheResult<double> r = mva::solver_mva_cache_analyzer(m.get_struct(), opt);
    // X is indexed by class: hit class 2, miss class 3 (0-based 1 and 2)
    CHECK(r.sol.X[1] == doctest::Approx(0.588571428571).epsilon(1e-9));
    CHECK(r.sol.X[2] == doctest::Approx(1.0 - 0.588571428571).epsilon(1e-9));
    // the two shares sum to the source rate, which is conservation of reads
    CHECK(r.sol.X[1] + r.sol.X[2] == doctest::Approx(1.0).epsilon(1e-9));
}

}  // namespace
