/*
 * Copyright (c) 2012-2026, QORE Lab, Imperial College London
 * All rights reserved.
 */

/**
 * `python/examples/basic/cacheModel`, in C++: the nine replacement-policy
 * models, the cost-capped cache, the refined-mean-field transient, and the five
 * delayed-hit retrieval models.
 *
 * WHERE THE HIT AND MISS RATIOS COME FROM. The reference reads them back off the
 * Cache OBJECT (`cache_node.get_hit_ratio()`), which every solver overwrites as
 * it runs, so the printed vector belongs to the LAST solver of the block. This
 * port has no model object to write back into: each analyzer returns its own
 * split beside its metric table (`FluidCacheqnSolution::hitprob`,
 * `NcCacheSolution::missrate`, `SsaSerialSolution::cache`, ...), so the ratios
 * below are read from the analyzer the reference would have read them from and
 * from no other. The three analyzers that return no split at all -- the C++
 * CTMC on any cache, and the OPEN MVA delayed-hit analyzer -- say so by name
 * instead of deriving a number the solver never produced.
 *
 * THE FEATURE GATE IS PART OF THE REFERENCE'S OUTPUT. `runAnalyzerChecks` runs
 * before every Python solve and raises naming the offending feature; that
 * refusal is what `cache_replc_lru`'s Fluid block and `cache_compare_replc`'s
 * NC column actually print, so `gate()` below reproduces it rather than letting
 * an unguarded solve answer a model the solver does not declare.
 */

#include <cmath>
#include <cstdio>
#include <exception>
#include <limits>
#include <map>
#include <string>
#include <vector>

#include "example_node_table.h"
#include "examples_common.h"
#include "line/api/cache/cache_miss_rmf.h"
#include "line/lang/qn/feature_set.h"
#include "line/lang/qn/environment.h"
#include "line/lang/qn/solver_feature_sets.h"
#include "line/solvers/ctmc/solver_ctmc_analyzer.h"
#include "line/solvers/fluid/fluid_cacheqn.h"
#include "line/solvers/mva/solver_mva_cache.h"
#include "line/solvers/mva/solver_mva_retrieval.h"
#include "line/solvers/mva/solver_mva_runner.h"
#include "line/solvers/nc/solver_nc_cache.h"
#include "line/solvers/nc/solver_nc_retrieval.h"
#include "line/solvers/nc/solver_nc_runner.h"
#include "line/solvers/ssa/solver_ssa_serial.h"
#include "line/solvers/ssa/ssa_dispatch.h"
#include "line/solvers/env/env_dispatch.h"
#include "line/solvers/env/solver_env.h"
#include "line/solvers/wrappers/ldes/solver_ldes.h"

namespace line {
namespace examples {

namespace {

using lang::ReplacementStrategy;

const double kNaN = std::numeric_limits<double>::quiet_NaN();

// ---------------------------------------------------------------------------
// Popularity vectors
// ---------------------------------------------------------------------------

/** `Zipf(alpha, n)` as the pmf `setRead` stores: k^-alpha over its harmonic sum. */
std::vector<double> zipf_pmf(double alpha, std::size_t n) {
    std::vector<double> p(n, 0.0);
    double h = 0.0;
    for (std::size_t k = 1; k <= n; ++k) h += std::pow(static_cast<double>(k), -alpha);
    for (std::size_t k = 1; k <= n; ++k)
        p[k - 1] = std::pow(static_cast<double>(k), -alpha) / h;
    return p;
}

/** `DiscreteSampler([1/n] * n)`, the uniform reference stream. */
std::vector<double> uniform_pmf(std::size_t n) {
    return std::vector<double>(n, 1.0 / static_cast<double>(n));
}

// ---------------------------------------------------------------------------
// The gate and the cache-specific printers
// ---------------------------------------------------------------------------

/**
 * `NetworkSolver.runAnalyzerChecks`: does the solver declare what the model
 * uses? A refusal is printed by name and the block is skipped, which is exactly
 * what the reference does before it computes anything.
 */
bool gate(const std::string& solver, const qn::FeatureSet& declared, const Sn& sn) {
    const qn::SupportResult r =
        qn::feature_set_supports(solver, declared, qn::used_lang_features(sn));
    if (!r.ok) na(solver, r.reason);
    return r.ok;
}

/** `print(f'Hit Ratio: {...}')`: the per-class vector on one line. */
void print_ratio(const std::string& label, const std::vector<double>& v) {
    std::string s;
    char buf[32];
    for (std::size_t i = 0; i < v.size(); ++i) {
        std::snprintf(buf, sizeof(buf), i ? " %.6g" : "%.6g", v[i]);
        s += buf;
    }
    kv(label, s);
}

/** The same from one row of an (ncaches x nclasses) split matrix. */
std::vector<double> matrix_row(const Matrix<double>& m, std::size_t row) {
    std::vector<double> v;
    for (std::size_t c = 0; c < m.cols(); ++c) v.push_back(m(row, c));
    return v;
}

/** The getAvgCacheTable columns this port's analyzers actually return. */
/**
 * The `getAvgCacheTable` header, INCLUDING THE TWO COLUMNS THE GOLDENS KEY ON.
 *
 * The reference's table carries HitProb, DelayedHitProb, MissProb, the three
 * matching rates, ArvR, ResidT and ListCost. This twin printed the three
 * probabilities alone, so `baselines/retrieval_*.json` -- which asserts ArvR and
 * ResidT -- matched NOTHING and the row failed with "solver X missing from
 * output" having compared zero cells. A row that measures nothing is the failure
 * mode the harness exists to prevent, so the columns are printed rather than the
 * golden narrowed.
 */
void cache_header() {
    std::printf("%-10s %-14s %5s %8s %6s %12s %12s %12s %12s %12s\n", "Node", "JobClass", "List",
                "ListCap", "Items", "HitProb", "DelayedHitP", "MissProb", "ArvR", "ResidT");
    // THE CACHE TABLE IS A RESULT TABLE and needs recording like any other: it
    // does not go through `avg_rows`, so without this the five `retrieval_*`
    // rows are printed, attributed to nothing, and read downstream as four
    // solvers that produced no output. One `begin_table` per solver block, which
    // is what one `cache_header` per block already is.
    namespace parity = line::examples::parity;
    if (parity::enabled()) parity::begin_table("cache", "Station", "JobClass");
}

void cache_row(const std::string& node, const std::string& cls, int list, double listcap,
               std::size_t items, double hit, double dhit, double miss, double arvr = kNaN,
               double residt = kNaN) {
    std::printf("%-10s %-14s %5d %8g %6zu %12.6g %12.6g %12.6g %12.6g %12.6g\n", node.c_str(),
                cls.c_str(), list, listcap, items, hit, dhit, miss, arvr, residt);
    namespace parity = line::examples::parity;
    if (!parity::enabled()) return;
    // ArvR AND ResidT ARE THE TWO THE GOLDENS KEY ON: they are the retrieval
    // system's arrival rate and latency, the pair that reports how many requests
    // are in flight. The three probabilities beside them are not among the
    // metric names a golden carries, so they are printed and not recorded.
    std::vector<parity::Cell> cells;
    cells.push_back(parity::Cell{"ArvR", arvr});
    cells.push_back(parity::Cell{"ResidT", residt});
    parity::add_row(node, cls, cells);
}

/**
 * `ArvR` of the cache table: the RETRIEVAL-SYSTEM throughput,
 * `lambda_read * (missprob + delayedprob)`, i.e. the rate of requests that enter
 * the retrieval system. It is the arrival rate Little-consistent with
 * `ResidT = latency`, so the pair reports the mean number of requests in the
 * retrieval system. `lambda_read` is the read class's own source rate, which is
 * what every solver agrees on, simulators included.
 */
double cache_retrieval_arvr(const Sn& sn, std::size_t cls, double miss, double dhit) {
    const std::size_t src = sn.sourceIdx;
    if (src == 0 || cls == 0 || cls > sn.nclasses) return kNaN;
    if (sn.disabled[src - 1][cls - 1]) return 0.0;
    const double lam = sn.rates(src - 1, cls - 1);
    if (!std::isfinite(lam)) return kNaN;
    const double m = std::isnan(miss) ? 0.0 : miss;
    const double d = std::isnan(dhit) ? 0.0 : dhit;
    return lam * (m + d);
}

/** The getAvgItemTable rows: per item, the miss column then one per cache list. */
void print_item_table(const Matrix<double>& itemprob) {
    if (itemprob.rows() == 0 || itemprob.cols() == 0) {
        note("(no per-item table: the analyzer returned none)");
        return;
    }
    std::printf("%-6s %12s", "Item", "Miss");
    for (std::size_t l = 1; l < itemprob.cols(); ++l) std::printf(" %11s%zu", "List", l);
    std::printf("\n");
    for (std::size_t i = 0; i < itemprob.rows(); ++i) {
        std::printf("%-6zu", i + 1);
        for (std::size_t l = 0; l < itemprob.cols(); ++l) std::printf(" %12.6g", itemprob(i, l));
        std::printf("\n");
    }
}

/** The Cache node's name and item count, for the cache-table rows. */
struct CacheInfo {
    std::string name;
    std::size_t node = 0;
    std::size_t nitems = 0;
    double totcap = 0.0;
    std::size_t nlists = 0;
};

CacheInfo cache_info(const Sn& sn) {
    CacheInfo ci;
    for (const auto& kv2 : sn.nodeparam) {
        if (sn.nodes[kv2.first - 1].nodetype != qn::NodeType::Cache) continue;
        ci.node = kv2.first;
        ci.name = sn.nodes[kv2.first - 1].name;
        ci.nitems = kv2.second.nitems;
        ci.nlists = kv2.second.itemcap.size();
        for (int c : kv2.second.itemcap) ci.totcap += c;
        break;
    }
    return ci;
}

// ---------------------------------------------------------------------------
// Model factories shared by the replacement-policy examples
// ---------------------------------------------------------------------------

/**
 * The CLOSED shape of `cache_replc_lru` and its siblings: one Delay holding the
 * single job, one Cache, and the hit / miss classes switching back to the read
 * class on the way home.
 */
Net closed_cache_model(ReplacementStrategy strat, const std::vector<double>& pread,
                       const std::vector<int>& itemcap, double qlru = 1.0) {
    Net net("model");
    Delay delay(net, "Delay");
    qn::CacheParam<double> ch;
    ch.nitems = pread.size();
    ch.itemcap = itemcap;
    ch.replacestrat = strat;
    ch.qlru = qlru;
    ch.pread = std::vector<std::vector<double> >{pread, {}, {}};
    ch.hitclass = std::vector<std::size_t>{2, 0, 0};
    ch.missclass = std::vector<std::size_t>{3, 0, 0};
    const std::size_t cache = net.add_cache("Cache", ch);
    ClosedClass job(net, "JobClass", 1.0, delay);
    ClosedClass hit(net, "HitClass", 0.0, delay);
    ClosedClass miss(net, "MissClass", 0.0, delay);
    delay.set_service(job, Exp(1.0));
    Routing P;
    P.set(job, job, delay, cache, 1.0);
    P.set(hit, job, cache, delay, 1.0);
    P.set(miss, job, cache, delay, 1.0);
    net.link(P);
    return net;
}

/** The OPEN Source-Cache-Sink shape of `cache_replc_rr` and its siblings. */
Net open_cache_model(ReplacementStrategy strat, const std::vector<double>& pread,
                     const std::vector<int>& itemcap, double arrival_rate,
                     const std::vector<int>& itemsize = std::vector<int>(),
                     const std::vector<int>& costcap = std::vector<int>()) {
    Net net("model");
    Source source(net, "Source");
    qn::CacheParam<double> ch;
    ch.nitems = pread.size();
    ch.itemcap = itemcap;
    ch.replacestrat = strat;
    ch.itemsize = itemsize;
    ch.costcap = costcap;
    ch.pread = std::vector<std::vector<double> >{pread, {}, {}};
    ch.hitclass = std::vector<std::size_t>{2, 0, 0};
    ch.missclass = std::vector<std::size_t>{3, 0, 0};
    const std::size_t cache = net.add_cache("Cache", ch);
    Sink sink(net, "Sink");
    OpenClass init(net, "InitClass");
    OpenClass hit(net, "HitClass");
    OpenClass miss(net, "MissClass");
    source.set_arrival(init, Exp(arrival_rate));
    Routing P;
    P.set(init, init, source, cache, 1.0);
    P.set(hit, hit, cache, sink, 1.0);
    P.set(miss, miss, cache, sink, 1.0);
    net.link(P);
    return net;
}

// ---------------------------------------------------------------------------
// Solver blocks the replacement-policy examples share
//
// EVERY ONE OF THESE REFERENCES CALLS `avg_node_table()`, NOT `getAvgTable()`,
// and on a cache model that is the whole answer rather than a presentation
// choice: the hit and miss rates live at the Cache node and at the ClassSwitch
// under it, and neither is a station. Printing the station table instead left
// `cache_replc_fifo` with one row of its seven and dropped 71 golden cells
// across the four goldened models. See example_node_table.h.
// ---------------------------------------------------------------------------

/** `CTMC(model, keep=False, cutoff=...)` then its node table. */
void run_ctmc(const Sn& sn, double cutoff, const std::string& method = "default") {
    if (!gate("SolverCTMC", qn::ctmc_feature_set(method), sn)) return;
    section("CTMC");
    ctmc::CtmcOptions opt;
    opt.method = method;
    opt.cutoff = cutoff;
    print_avg_node(sn, ctmc::solver_ctmc_run_analyzer(sn, opt));
}

/** `SSA(model, samples=..., seed=...)` then its node table. */
void run_ssa(const Sn& sn, std::size_t samples, unsigned long seed) {
    if (!gate("SolverSSA", qn::ssa_feature_set("default"), sn)) return;
    section("SSA");
    ssa::SsaOptions opt;
    opt.samples = samples;
    opt.seed = seed;
    // The MEASURED hit/miss split, not `link()`'s offered one, in BOTH places it
    // is read: the derived ResidT and ArvR columns divide by visit ratios, and
    // the node table splits a cache's flow by the share the run converged to.
    // Same collection the CLI's `-s ssa` arm makes.
    std::vector<ssa::SsaCacheRatio> cacheratio;
    const ssa::SsaSolution r = ssa::solver_ssa(sn, opt, &cacheratio);
    mva::AvgResult<double> a = solvers::avg_result_from_sim<double>(
        ssa::sn_with_ssa_cache_split<double>(sn, cacheratio), r.QN, r.UN, r.RN, r.TN, r.CN, r.XN,
        r.method);
    a.cache = ssa::cache_metrics_of_ssa<double>(sn, cacheratio);
    print_avg_node(sn, a);
}

/** `MVA(model)` then its node table. */
void run_mva(const Sn& sn) {
    if (!gate("SolverMVA", qn::mva_feature_set("default"), sn)) return;
    section("MVA");
    mva::MvaOptions opt;
    const Matrix<double> init;
    print_avg_node(sn, mva::solver_mva_run_analyzer(sn, opt, init));
}

/** `NC(model)` then its node table. */
void run_nc(const Sn& sn, const std::string& method = "default") {
    if (!gate("SolverNC", qn::nc_feature_set(method), sn)) return;
    section("NC");
    nc::NcSolverOptions opt;
    opt.method = method;
    print_avg_node(sn, nc::solver_nc_run_analyzer(sn, opt));
}

/**
 * `FLD(model, method='rmf')` then its node table AND the hit / miss ratios the
 * reference prints afterwards, which on these models is the split the fluid
 * decomposition converged to. One solve, not two: `solver_fluid_run_analyzer`'s cache arm
 * IS this analyzer, and its return carries the split, the renormalized struct
 * and the solution together.
 */
void run_fluid_rmf(const Sn& sn, bool print_ratios) {
    if (!gate("SolverFluid", qn::fluid_feature_set("rmf"), sn)) {
        if (print_ratios)
            na("Hit Ratio / Miss Ratio",
               "the reference reads them back off the Cache object after the Fluid solve above, "
               "and that solve was refused");
        return;
    }
    // THE GOLDEN'S KEY IS `FLD`, which is what the other seventeen fluid blocks
    // in the corpus declare. `Fluid` is a spelling no golden carries, so the
    // table was recorded under a name nothing could align and the row read
    // "solver FLD missing from the recorded results".
    section("FLD");
    fluid::FluidOptions opt;
    opt.method = "rmf";
    const fluid::FluidCacheqnSolution<double> r = fluid::solver_fld_cacheqn_analyzer(sn, opt);
    // `r.refreshed` carries the CONVERGED hit/miss split; `sn` still carries
    // `link()`'s offered one. Both derived columns and the node-level flow split
    // are read off the routing, so the offered split reported an even 0.5/0.5
    // where `cache_replc_routing` converged to 0.4/0.6, and put zero arrivals on
    // stations that do take them.
    const bool has_ref = !r.refreshed.nodes.empty();
    mva::AvgResult<double> a = solvers::avg_result_from_sim<double>(
        has_ref ? r.refreshed : sn, r.sol.QN, r.sol.UN, r.sol.RN, r.sol.TN, r.sol.CN, r.sol.XN,
        r.sol.method);
    if (has_ref) a.refreshed_struct.reset(new Sn(r.refreshed));
    // AND THE SPLIT ITSELF, not only the struct it renormalized: the node table
    // prefers a stated share over the visit ratios, which is the difference
    // between 0.4/0.6 and the 1/2-1/2 `link()` offers before anything is solved.
    // `cache_metrics_of` is the one assembler every cache analyzer answers
    // through, so what is stated here is what `-a cache` would state.
    a.cache = solvers::cache_metrics_of<double>(
        sn, matrix_row(r.hitprob, 0), matrix_row(r.missprob, 0), std::vector<double>(),
        std::vector<double>(), Matrix<double>(), Matrix<double>(), std::vector<double>());
    print_avg_node(sn, a);
    if (!print_ratios) return;
    print_ratio("Hit Ratio", matrix_row(r.hitprob, 0));
    print_ratio("Miss Ratio", matrix_row(r.missprob, 0));
}

// ---------------------------------------------------------------------------
// The delayed-hit retrieval models
// ---------------------------------------------------------------------------

/** The retrieval classes `set_retrieval_system` mints, one per item. */
std::vector<std::size_t> retrieval_classes(std::size_t last_class_before, std::size_t nitems) {
    std::vector<std::size_t> rc;
    for (std::size_t i = 0; i < nitems; ++i) rc.push_back(last_class_before + 1 + i);
    return rc;
}

/** `SSA(model, samples=100000, method='serial', seed=1).get_avg_cache_table()`. */
void retrieval_ssa_cache_table(const Sn& sn) {
    if (!gate("SolverSSA", qn::ssa_feature_set("serial"), sn)) return;
    section("SSA");
    ssa::SsaSerialOptions opt;
    opt.method = "serial";
    opt.samples = 100000;
    opt.seed = 1;
    const ssa::SsaSerialSolution<double> r = ssa::solver_ssa_serial(sn, opt);
    const CacheInfo ci = cache_info(sn);
    cache_header();
    for (std::size_t k = 0; k < r.cache.size(); ++k)
        for (std::size_t c = 0; c < sn.nclasses; ++c) {
            const double h = r.cache[k].hitprob[c], mi = r.cache[k].missprob[c];
            if (std::isnan(h) && std::isnan(mi)) continue;
            // The delayed share is a MEASURED column here, not a guess: the
            // engine counts the merge transitions, so a model without a
            // retrieval system leaves the field empty and reports NaN.
            const double d = r.cache[k].delayedprob.empty() ? kNaN : r.cache[k].delayedprob[c];
            // ResidT is NaN and that IS parity: `SsaCacheRatio` records that the
            // reference warns "Retrieval-system expected latency is not
            // currently implemented" and reports NaN in every codebase.
            cache_row(sn.nodes[r.cache[k].node - 1].name, sn.classes[c].name, 0, ci.totcap,
                      ci.nitems, h, d, mi, cache_retrieval_arvr(sn, c + 1, mi, d), kNaN);
        }
}

/**
 * `MVA(model).get_avg_cache_table()` on an OPEN delayed-hit model.
 *
 * `solver_mva_retrieval_analyzer` returns the metric table only; the aggregate
 * it does publish is the hit-class throughput, which carries TRUE hits and
 * delayed hits SUMMED, so the two are reported together and the columns the
 * analyzer never forms are refused rather than split by guess.
 */
void retrieval_mva_cache_table(const Sn& sn, std::size_t read_class) {
    if (!gate("SolverMVA", qn::mva_feature_set("default"), sn)) return;
    section("MVA");
    mva::MvaOptions opt;
    const mva::MvaSolution<double> s = mva::solver_mva_retrieval_analyzer(sn, opt);
    const CacheInfo ci = cache_info(sn);
    const std::size_t src = sn.sourceIdx;
    const double lam = sn.disabled[src - 1][read_class - 1] ? 0.0 : sn.rates(src - 1, read_class - 1);
    std::size_t hc = 0, mc = 0;
    for (const auto& kv2 : sn.nodeparam)
        if (sn.nodes[kv2.first - 1].nodetype == qn::NodeType::Cache) {
            hc = kv2.second.hitclass[read_class - 1];
            mc = kv2.second.missclass[read_class - 1];
            break;
        }
    std::printf("%-10s %-14s %8s %6s %18s %12s\n", "Node", "JobClass", "ListCap", "Items",
                "Hit+DelayedHitProb", "MissProb");
    const double hitshare = lam > 0.0 ? s.X[hc - 1] / lam : kNaN;
    const double missshare = lam > 0.0 ? s.X[mc - 1] / lam : kNaN;
    std::printf("%-10s %-14s %8g %6zu %18.6g %12.6g\n", ci.name.c_str(),
                sn.classes[read_class - 1].name.c_str(), ci.totcap, ci.nitems, hitshare,
                missshare);
    // ONE ROW, WITH NO CELL THE GOLDEN NAMES. `baselines/retrieval_*.json` asks
    // MVA for a single ('Cache', 'Jobs') row holding a NaN QLen: it asserts that
    // the solver ANSWERED, not what it answered, because this analyzer returns
    // the two shares above and no rate. The row must therefore exist and carry
    // the NaN, which is what the golden holds; dropping it reads as a solver
    // that produced nothing.
    namespace parity = line::examples::parity;
    if (parity::enabled()) {
        parity::begin_table("cache", "Station", "JobClass");
        std::vector<parity::Cell> cells;
        cells.push_back(parity::Cell{"QLen", kNaN});
        parity::add_row(ci.name, sn.classes[read_class - 1].name, cells);
    }
    note("N/A: HitProb and DelayedHitProb separately -- the C++ open delayed-hit MVA analyzer "
         "returns only the metric table, whose hit-class throughput sums the two");
}

/** `NC(model).get_avg_cache_table()` and `.get_avg_item_table()`. */
void retrieval_nc_tables(const Sn& sn) {
    if (!gate("SolverNC", qn::nc_feature_set("default"), sn)) return;
    section("NC");
    nc::NcSolverOptions opt;
    const nc::NcRetrievalSolution<double> r = nc::solver_nc_retrieval_analyzer(sn, opt);
    const CacheInfo ci = cache_info(sn);
    cache_header();
    for (std::size_t c = 0; c < sn.nclasses; ++c) {
        if (std::isnan(r.hitprob[c]) && std::isnan(r.missprob[c])) continue;
        cache_row(ci.name, sn.classes[c].name, 0, ci.totcap, ci.nitems, r.hitprob[c],
                  r.delayedprob[c], r.missprob[c],
                  cache_retrieval_arvr(sn, c + 1, r.missprob[c], r.delayedprob[c]),
                  c < r.latency.size() ? r.latency[c] : kNaN);
        for (std::size_t l = 0; l < r.hitproblist.cols(); ++l)
            cache_row(ci.name, sn.classes[c].name, static_cast<int>(l + 1), kNaN, ci.nitems,
                      r.hitproblist(c, l), kNaN, kNaN);
    }
    note("Item table (NC):");
    print_item_table(r.itemprob);
}

/**
 * `LDES(model, samples=1000000, seed=1).getAvgCacheTable()`.
 *
 * The delayed-hit fractions are sample-path events in the engine, not an
 * estimator this port reimplements: `ldes::solver_ldes` runs the same
 * `common/ldes` binary the MATLAB and Python clients drive and reads the
 * `cacheMetrics` block back.
 */
void retrieval_ldes_cache_table(const Sn& sn) {
    if (!gate("SolverLDES", qn::ldes_feature_set("default"), sn)) return;
    section("LDES");
    ldes::LdesOptions opt;
    opt.samples = 1000000;
    opt.seed = 1;
    ldes::LdesResult r;
    try {
        r = ldes::solver_ldes(sn, opt);
    } catch (const std::exception& e) {
        na("LDES", e.what());
        return;
    }
    const CacheInfo ci = cache_info(sn);
    cache_header();
    for (std::map<std::string, ldes::LdesCacheMetrics>::const_iterator it = r.cache_metrics.begin();
         it != r.cache_metrics.end(); ++it) {
        const ldes::LdesCacheMetrics& cm = it->second;
        // Only the READ classes get a row, which is how the reference selects
        // them: a class with no hit class configured never reads this cache.
        const std::vector<std::size_t>& hitcls = sn.nodeparam.at(ci.node).hitclass;
        for (std::size_t c = 0; c < sn.nclasses; ++c) {
            if (c >= hitcls.size() || hitcls[c] == 0) continue;
            const double h = c < cm.hit.cols() ? cm.hit(0, c) : kNaN;
            const double d = c < cm.delayed.cols() ? cm.delayed(0, c) : kNaN;
            const double mi = c < cm.miss.cols() ? cm.miss(0, c) : kNaN;
            if (std::isnan(h) && std::isnan(d) && std::isnan(mi)) continue;
            const double lat = (cm.latency.rows() > 0 && c < cm.latency.cols())
                                   ? cm.latency(0, c)
                                   : kNaN;
            cache_row(it->first, sn.classes[c].name, 0, ci.totcap, ci.nitems, h, d, mi,
                      cache_retrieval_arvr(sn, c + 1, mi, d), lat);
        }
        for (std::size_t c = 0; c < cm.hitList.rows(); ++c)
            for (std::size_t l = 0; l < cm.hitList.cols(); ++l)
                cache_row(it->first, sn.classes[c].name, static_cast<int>(l + 1), kNaN, ci.nitems,
                          cm.hitList(c, l), kNaN, kNaN);
    }
    if (r.cache_metrics.empty())
        note("(the engine reported no cacheMetrics block for this model)");
}

}  // namespace

// ---------------------------------------------------------------------------
// Replacement policies
// ---------------------------------------------------------------------------

/**
 * Cache with LRU replacement, closed model, uniform references over 5 items.
 *
 * No Fluid row: LRU has no drift-based fluid model, so the refined mean field
 * carries RANDOM(m)/FIFO(m) and strict FIFO(m) only, and the reference dropped
 * its SolverFLD call for the same reason.
 */
void cache_replc_lru() {
    const std::size_t n = 5;
    Net net =
        closed_cache_model(ReplacementStrategy::LRU, uniform_pmf(n), std::vector<int>{2});
    const Sn& sn = net.get_struct();
    run_ctmc(sn, -1.0);
    run_ssa(sn, 100000, 23000);
    run_mva(sn);
}
LINE_EXAMPLE("basic/cacheModel", cache_replc_lru);

/** Cache with FIFO replacement, closed model, uniform references over 5 items. */
void cache_replc_fifo() {
    const std::size_t n = 5;
    Net net =
        closed_cache_model(ReplacementStrategy::FIFO, uniform_pmf(n), std::vector<int>{2});
    const Sn& sn = net.get_struct();
    run_ctmc(sn, -1.0);
    run_ssa(sn, 100000, 23000);
    run_mva(sn);
    run_fluid_rmf(sn, true);
}
LINE_EXAMPLE("basic/cacheModel", cache_replc_fifo);

/** Cache with random replacement, open model, Zipf(1.4) references. */
void cache_replc_rr() {
    const std::size_t n = 5;
    Net net = open_cache_model(ReplacementStrategy::RR, zipf_pmf(1.4, n),
                                     std::vector<int>{2}, 2.0);
    const Sn& sn = net.get_struct();
    run_ctmc(sn, 1.0);
    run_ssa(sn, 10000, 23000);
    run_mva(sn);
    run_nc(sn);
    run_fluid_rmf(sn, true);
}
LINE_EXAMPLE("basic/cacheModel", cache_replc_rr);

/**
 * Cache with h-LRU / LRU(m): two lists of capacities 2 and 1 over 6 items, with
 * MVA taking the characteristic-time (TTL) approximation of Gast and Van Houdt.
 */
void cache_replc_hlru() {
    const std::size_t n = 6;
    Net net = open_cache_model(ReplacementStrategy::HLRU, zipf_pmf(1.2, n),
                                     std::vector<int>{2, 1}, 1.0);
    const Sn& sn = net.get_struct();
    run_ctmc(sn, 1.0);
    run_mva(sn);
    run_ssa(sn, 10000, 23000);
}
LINE_EXAMPLE("basic/cacheModel", cache_replc_hlru);

/** Cache with the CLIMB (transposition) rule: exact in CTMC, refused elsewhere. */
void cache_replc_climb() {
    const std::size_t n = 5;
    Net net =
        closed_cache_model(ReplacementStrategy::CLIMB, zipf_pmf(1.2, n), std::vector<int>{2});
    const Sn& sn = net.get_struct();
    run_ctmc(sn, -1.0, "exact");
}
LINE_EXAMPLE("basic/cacheModel", cache_replc_climb);

/** Cache with q-LRU: a miss is admitted with probability q = 0.5. */
void cache_replc_qlru() {
    const std::size_t n = 5;
    Net net = closed_cache_model(ReplacementStrategy::QLRU, zipf_pmf(1.2, n),
                                       std::vector<int>{2}, 0.5);
    const Sn& sn = net.get_struct();
    run_ctmc(sn, -1.0, "exact");
}
LINE_EXAMPLE("basic/cacheModel", cache_replc_qlru);

/**
 * Cache feeding a Router that splits hits and misses uniformly over two delays
 * with different service rates.
 */
void cache_replc_routing() {
    const std::size_t n = 5;
    Net net("model");
    Source source(net, "Source");
    qn::CacheParam<double> ch;
    ch.nitems = n;
    ch.itemcap = std::vector<int>{2};
    ch.replacestrat = ReplacementStrategy::FIFO;
    ch.pread = std::vector<std::vector<double> >{uniform_pmf(n), {}, {}};
    ch.hitclass = std::vector<std::size_t>{2, 0, 0};
    ch.missclass = std::vector<std::size_t>{3, 0, 0};
    const std::size_t cache = net.add_cache("Cache", ch);
    Router router(net, "Router");
    Delay delay1(net, "Delay1");
    Delay delay2(net, "Delay2");
    Sink sink(net, "Sink");

    OpenClass init(net, "InitClass");
    OpenClass hit(net, "HitClass");
    OpenClass miss(net, "MissClass");

    source.set_arrival(init, Exp(2.0));
    delay1.set_service(hit, Exp(10.0));
    delay1.set_service(miss, Exp(1.0));
    delay2.set_service(hit, Exp(20.0));
    delay2.set_service(miss, Exp(2.0));

    // A RAND dispatcher reads only the CONNECTIONS; the refresh replaces the
    // positive entries below by the uniform split.
    net.set_routing(router, hit, RoutingStrategy::RAND);
    net.set_routing(router, miss, RoutingStrategy::RAND);

    Routing P;
    P.set(init, init, source, cache, 1.0);
    P.set(hit, hit, cache, router, 1.0);
    P.set(miss, miss, cache, router, 1.0);
    P.set(hit, hit, router, delay1, 1.0);
    P.set(hit, hit, router, delay2, 1.0);
    P.set(miss, miss, router, delay1, 1.0);
    P.set(miss, miss, router, delay2, 1.0);
    P.set(hit, hit, delay1, sink, 1.0);
    P.set(hit, hit, delay2, sink, 1.0);
    P.set(miss, miss, delay1, sink, 1.0);
    P.set(miss, miss, delay2, sink, 1.0);
    net.link(P);

    const Sn& sn = net.get_struct();
    run_ctmc(sn, 1.0);
    run_ssa(sn, 10000, 23000);
    run_mva(sn);
    run_nc(sn);
    run_fluid_rmf(sn, true);
}
LINE_EXAMPLE("basic/cacheModel", cache_replc_routing);

/**
 * RR, FIFO and LRU compared on one open two-level cache, one hit-ratio line per
 * policy. The CTMC column is refused: this port's CTMC returns no cache split.
 */
void cache_compare_replc() {
    const std::size_t n = 5;
    const double alpha = 1.0;
    const std::vector<int> m{2, 1};
    const ReplacementStrategy strat[3] = {ReplacementStrategy::RR, ReplacementStrategy::FIFO,
                                          ReplacementStrategy::LRU};
    const char* names[3] = {"rr", "fifo", "lru"};

    for (int s = 0; s < 3; ++s) {
        Net net = open_cache_model(strat[s], zipf_pmf(alpha, n), m, 2.0);
        const Sn& sn = net.get_struct();

        // `cacheNode.getHitRatio` AFTER the CTMC solve, which is the hit share
        // the analyzer writes back into `actualhitprob`: the hit-class
        // departure rate at the cache over the hit- plus miss-class one. It is
        // read off the solved chain, not off the routing matrix, whose cache
        // entries `refresh_routing` resolved to a uniform split.
        double ctmc_hr = kNaN;
        if (gate("SolverCTMC", qn::ctmc_feature_set("exact"), sn)) {
            ctmc::CtmcOptions copt;
            copt.method = "exact";
            copt.cutoff = 1.0;
            const ctmc::CtmcSolution<double> cr = ctmc::solver_ctmc_analyzer(sn, copt);
            if (!cr.cache.caches.empty() && !cr.cache.caches[0].hitprob.empty())
                ctmc_hr = cr.cache.caches[0].hitprob[0];
        }

        double mva_hr = kNaN, nc_hr = kNaN;
        if (gate("SolverMVA", qn::mva_feature_set("default"), sn)) {
            mva::MvaOptions opt;
            mva_hr = mva::solver_mva_cache_analyzer(sn, opt).hitprob[0];
        }
        if (gate("SolverNC", qn::nc_feature_set("default"), sn)) {
            nc::NcSolverOptions opt;
            const nc::NcCacheSolution<double> r = nc::solver_nc_cache_analyzer(sn, opt);
            const double lam = sn.rates(sn.sourceIdx - 1, 0);
            if (lam > 0.0) nc_hr = 1.0 - r.missrate[0] / lam;
        }
        std::printf("%s: %.8f, %.8f, %.8f\n", names[s], ctmc_hr, mva_hr, nc_hr);
        // THE GOLDEN'S `item0`, `item1` AND `item2` ARE NOT ITEMS. They are the
        // CTMC, MVA and NC columns of the line above, and the reader that
        // generated the golden named them by position. Reproduced positionally
        // and under the same names so the golden still matches; regenerating it
        // with honest column names is a separate change. The shape key is
        // `CACHE`: a hit ratio compared across three solvers is not any one
        // solver's row.
        const double column[3] = {ctmc_hr, mva_hr, nc_hr};
        for (int k = 0; k < 3; ++k)
            derived("CACHE", names[s], "item" + std::to_string(k), column[k]);
    }
}
LINE_EXAMPLE("basic/cacheModel", cache_compare_replc);

/**
 * Per-item storage costs with per-list caps: list 2 admits small items only.
 * `SolverNC` evaluates the constrained normalizing constant E(m,k) of
 * Casale-Gast, IEEE/ACM ToN 29(2), 2021, Sec. IX.
 */
void cache_itemsize_costcap() {
    const std::size_t n = 6;
    Net net = open_cache_model(ReplacementStrategy::RR, uniform_pmf(n),
                                     std::vector<int>{1, 1}, 2.0,
                                     std::vector<int>{1, 1, 1, 2, 2, 2}, std::vector<int>{2, 1});
    const Sn& sn = net.get_struct();
    if (!gate("SolverNC", qn::nc_feature_set("exact"), sn)) return;

    section("NC");
    nc::NcSolverOptions opt;
    opt.method = "exact";
    print_avg(sn, nc::solver_nc_run_analyzer(sn, opt));

    const nc::NcCacheSolution<double> r = nc::solver_nc_cache_analyzer(sn, opt);
    const CacheInfo ci = cache_info(sn);
    const double lam = sn.rates(sn.sourceIdx - 1, 0);

    // The reference's per-class hit ratio: one over each read stream, zero where
    // the class never reads, which is what the miss RATE over the read rate is.
    std::vector<double> hit(sn.nclasses, 0.0);
    for (std::size_t c = 0; c < sn.nclasses && c < r.missrate.size(); ++c)
        if (c == 0 && lam > 0.0) hit[c] = 1.0 - r.missrate[c] / lam;
    const double miss = 1.0 - hit[0];

    note("Cache table (NC):");
    cache_header();
    cache_row(ci.name, sn.classes[0].name, 0, ci.totcap, ci.nitems, hit[0], kNaN, miss);
    for (std::size_t l = 0; l < r.hitproblist.cols(); ++l)
        cache_row(ci.name, sn.classes[0].name, static_cast<int>(l + 1), kNaN, ci.nitems,
                  r.hitproblist(0, l), kNaN, kNaN);

    note("Item table (NC):");
    print_item_table(r.itemprob);

    print_ratio("Hit Ratio", hit);
    print_ratio("Mean per-list storage cost", r.listcost);
    if (r.costcap_method_switched)
        note("(the requested method had no cost-capped counterpart and was switched)");
    if (!r.costcap_blocked.empty())
        note("(a cost cap blocks a promotion path: cross-check with SolverLDES)");
}
LINE_EXAMPLE("basic/cacheModel", cache_itemsize_costcap);

/**
 * Refined mean field of a two-list RANDOM(m) cache: the steady-state hit rates
 * with the 1/N correction, then the transient the reference also prints.
 */
void cache_rmf_transient() {
    const std::size_t n = 10;
    const std::vector<int> m{3, 2};
    const double alpha = 0.8;
    const std::vector<double> p = zipf_pmf(alpha, n);

    std::printf("Cache parameters: n=%zu, m=[%d, %d], Zipf(%.1f)\n", n, m[0], m[1], alpha);

    Matrix<double> lambda(1, n, 0.0);
    for (std::size_t i = 0; i < n; ++i) lambda(0, i) = p[i];
    const cache::CacheMissRmfResult<double> r =
        cache::cache_miss_rmf(std::vector<double>(), m, lambda, 10000.0);

    std::printf("\nSteady-state results (refined mean field):\n");
    double total_hit = 0.0;
    for (std::size_t k = 1; k <= m.size(); ++k) {
        const double hr = cache::rmf_detail::hit_rate(r.xss, p, k, n);
        total_hit += hr;
        std::printf("  Hit rate (list %zu): %.6f\n", k, hr);
        // The example SUMS `hit_rate` over the lists and forms the transient
        // state as `X[t] + V[t]/n` itself, so these totals exist nowhere else --
        // no result table carries them and no getter returns them. The shape key
        // is `CACHE`; no solver produced these.
        derived("CACHE", "HitRate_L" + std::to_string(k), "Steady", hr);
    }
    const double miss_rate = cache::rmf_detail::hit_rate(r.xss, p, std::size_t(0), n);
    std::printf("  Miss rate:         %.6f\n", miss_rate);
    std::printf("  Total hit prob:    %.6f\n", total_hit);
    std::printf("  Total miss prob:   %.6f\n", miss_rate);
    derived("CACHE", "MissRate", "Steady", miss_rate);
    derived("CACHE", "TotalHitProb", "Steady", total_hit);
    derived("CACHE", "TotalMissProb", "Steady", miss_rate);
    if (!r.refined) note("(the 1/N correction was declined; the plain mean field is reported)");

    // The reference's `meanFieldExpansionTransient(50, 200, 1)`, read at the
    // same five indices of its 200-point grid: t = 0, ~5, ~12.5, ~25 and 50.
    const cache::CacheRmfExpansionTransient<double> tr =
        cache::cache_miss_rmf_expansion_transient(m, lambda, 50.0, 200, 1);
    std::printf("\nTransient hit rates (refined, N=%zu):\n", n);
    const std::size_t idx[5] = {0, 20, 50, 100, 199};
    for (std::size_t s = 0; s < 5; ++s) {
        const std::size_t j = idx[s];
        std::vector<double> xt(tr.X.cols(), 0.0);
        for (std::size_t i = 0; i < xt.size(); ++i)
            xt[i] = tr.X(j, i) + tr.V(j, i) / static_cast<double>(n);
        double hr = 0.0;
        for (std::size_t k = 1; k <= m.size(); ++k) hr += cache::rmf_detail::hit_rate(xt, p, k, n);
        std::printf("  t=%7.3f: hit_rate=%.6f\n", tr.t[j], hr);
        // Keyed by the time as the golden spells it, `t=%.3f`: the grid point
        // is the row label, and there is no table these belong to.
        char label[32];
        std::snprintf(label, sizeof(label), "t=%.3f", tr.t[j]);
        derived("CACHE", label, "Transient", hr);
    }
}
LINE_EXAMPLE("basic/cacheModel", cache_rmf_transient);

// ---------------------------------------------------------------------------
// Cache networks
// ---------------------------------------------------------------------------

/**
 * Tandem of two caches sharing one item set, the C++ twin of
 * `matlab/examples/basic/cacheModel/cachenet_tandem.m`.
 *
 * A request reads Cache1; on a miss it reads Cache2 for THE SAME item; on a
 * second miss it leaves the network. Item identity is carried by giving each
 * item its OWN job class for the whole cycle, so no class switching happens on
 * an arc into a cache and the miss class of item i at Cache1 IS the read class
 * of item i at Cache2.
 *
 * Capacities: with LRU(1) at BOTH levels Cache2 can never hit, because a miss
 * at Cache1 inserts the item there and the next Cache1 miss is necessarily a
 * different item. Cache2 therefore gets room for two items.
 */
void cachenet_tandem() {
    const std::size_t n = 3;
    const double pAccess[3] = {0.5, 0.3, 0.2};

    Net net("CacheTandem");
    Delay think(net, "Think");

    qn::CacheParam<double> c1p;
    c1p.nitems = n;
    c1p.itemcap = std::vector<int>{1};
    c1p.replacestrat = ReplacementStrategy::LRU;
    const std::size_t cache1 = net.add_cache("Cache1", c1p);

    qn::CacheParam<double> c2p;
    c2p.nitems = n;
    c2p.itemcap = std::vector<int>{2};
    c2p.replacestrat = ReplacementStrategy::LRU;
    const std::size_t cache2 = net.add_cache("Cache2", c2p);

    // one dedicated class per item for each role
    std::vector<std::size_t> read(n), hit1(n), hit2(n), miss(n);
    for (std::size_t f = 0; f < n; ++f) {
        read[f] = net.add_closed_class("Read" + std::to_string(f + 1), 1.0, think);
        hit1[f] = net.add_closed_class("Hit1_" + std::to_string(f + 1), 0.0, think);
        hit2[f] = net.add_closed_class("Hit2_" + std::to_string(f + 1), 0.0, think);
        miss[f] = net.add_closed_class("Miss_" + std::to_string(f + 1), 0.0, think);
        think.set_service(read[f], Exp(pAccess[f]));
    }

    net.set_item_read_classes(cache1, read, hit1);
    net.set_miss_cache(cache1, cache2, hit2);
    net.set_item_miss_class(cache2, miss);

    Routing P;
    for (std::size_t f = 0; f < n; ++f) {
        P.set(read[f], read[f], think, cache1, 1.0);
        P.set(hit1[f], read[f], cache1, think, 1.0);
        P.set(hit2[f], read[f], cache2, think, 1.0);
        P.set(miss[f], read[f], cache2, think, 1.0);
        // the Cache1 -> Cache2 miss hop for c2read[f] is registered by
        // set_miss_cache and injected by link(), so it is NOT drawn here
    }
    net.link(P);

    const Sn& sn = net.get_struct();
    if (!gate("SolverSSA", qn::ssa_feature_set("serial"), sn)) return;
    section("SSA");
    ssa::SsaSerialOptions opt;
    opt.method = "serial";
    opt.samples = 200000;
    opt.seed = 1;
    const ssa::SsaSerialSolution<double> r = ssa::solver_ssa_serial(sn, opt);
    // Per (cache, read class) hit and miss. `cache_info` stops at the FIRST
    // cache, so the per-node capacity is read from `sn.nodeparam` here instead.
    cache_header();
    for (std::size_t k = 0; k < r.cache.size(); ++k) {
        const std::size_t nd = r.cache[k].node;
        const qn::CacheParam<double>& cp = sn.nodeparam.at(nd);
        double cap = 0.0;
        for (int c : cp.itemcap) cap += c;
        for (std::size_t c = 0; c < sn.nclasses; ++c) {
            const double h = r.cache[k].hitprob[c], mi = r.cache[k].missprob[c];
            if (std::isnan(h) && std::isnan(mi)) continue;
            cache_row(sn.nodes[nd - 1].name, sn.classes[c].name, 0, cap, cp.nitems, h, kNaN, mi);
        }
    }
}
LINE_EXAMPLE("basic/cacheModel", cachenet_tandem);

/**
 * MMAP-fed small RR cache with two correlated classes.
 *
 * A marked MMPP2 arrival stream feeds a small Round-Robin cache. Its two marks
 * are bound to two open read classes that share the modulating chain, so the
 * classes are cross-correlated and autocorrelated in time, and each reads the
 * cache with a DIFFERENT item popularity. Phase 1 (bursty) emits mostly class-1
 * references at a high rate, phase 2 (calm) mostly class-2 at a low rate, so
 * the shared chain couples "which class arrives" with "how fast".
 *
 * The same system is solved three ways: LDES simulates the true MMAP-fed cache,
 * CTMC solves it exactly, and SolverENV views the MMPP2 phase as an environment
 * modulating phase-conditional Poisson arrivals (D1 diagonal). Of the four ENV
 * methods the reference runs, only `blend` is answerable here; see the refusals
 * below.
 */
void cache_mmap_rr_env() {
    const std::size_t n = 4;                  // number of items
    const std::vector<int> m{2};              // cache capacity
    // Per-class item-popularity distributions (deliberately different)
    const std::vector<double> p1{8.0 / 15, 4.0 / 15, 2.0 / 15, 1.0 / 15};
    const std::vector<double> p2{1.0 / 15, 2.0 / 15, 4.0 / 15, 8.0 / 15};

    // Marked MMPP2 (M3A layout D = {D0, D11, D12}); D1 = D11 + D12 diagonal.
    // Phase 1 bursty (rate 4, 90% class1); phase 2 calm (rate 1, 80% class2).
    // Off-diagonal of D0 are the phase-switch rates (both 0.5).
    Matrix<double> D0(2, 2), D11(2, 2, 0.0), D12(2, 2, 0.0);
    D0(0, 0) = -4.5; D0(0, 1) = 0.5; D0(1, 0) = 0.5; D0(1, 1) = -1.5;
    D11(0, 0) = 3.6; D11(1, 1) = 0.2;   // class-1 arrivals per phase
    D12(0, 0) = 0.4; D12(1, 1) = 0.8;   // class-2 arrivals per phase

    // `buildCacheModel` with the arrival law of one stage, or the MMAP itself.
    struct Build {
        static Net run(std::size_t n, const std::vector<int>& m, const std::vector<double>& p1,
                       const std::vector<double>& p2, std::size_t* cache_out) {
            Net net("MMAPCache");
            Source source(net, "Source");
            qn::CacheParam<double> ch;
            ch.nitems = n;
            ch.itemcap = m;
            ch.replacestrat = ReplacementStrategy::RR;
            ch.pread = std::vector<std::vector<double> >{p1, p2, {}, {}, {}, {}};
            ch.hitclass = std::vector<std::size_t>{3, 5, 0, 0, 0, 0};
            ch.missclass = std::vector<std::size_t>{4, 6, 0, 0, 0, 0};
            const std::size_t cache = net.add_cache("Cache", ch);
            Sink sink(net, "Sink");
            OpenClass rd1(net, "Read1");
            OpenClass rd2(net, "Read2");
            OpenClass hit1(net, "Hit1");
            OpenClass mis1(net, "Miss1");
            OpenClass hit2(net, "Hit2");
            OpenClass mis2(net, "Miss2");
            Routing P;
            P.set(rd1, rd1, source, cache, 1.0);
            P.set(rd2, rd2, source, cache, 1.0);
            P.set(hit1, hit1, cache, sink, 1.0);
            P.set(mis1, mis1, cache, sink, 1.0);
            P.set(hit2, hit2, cache, sink, 1.0);
            P.set(mis2, mis2, cache, sink, 1.0);
            net.link(P);
            if (cache_out) *cache_out = cache;
            (void)source;
            return net;
        }
    };

    // The true system: one MMAP whose two marks are the two read classes.
    std::size_t cache_node = 0;
    Net true_model = Build::run(n, m, p1, p2, &cache_node);
    std::vector<Matrix<double> > D1k;
    D1k.push_back(D11);
    D1k.push_back(D12);
    const D mmap = D::mmap(D0, D1k);
    true_model.set_arrival(0, 0, mmap);
    true_model.set_arrival(0, 1, mmap);
    true_model.set_marked_classes(0, std::vector<std::size_t>{1, 2});
    const Sn& sn_true = true_model.get_struct();

    // (1) LDES - simulation of the true MMAP-fed cache
    if (gate("SolverLDES", qn::ldes_feature_set("default"), sn_true)) {
        section("LDES");
        ldes::LdesOptions opt;
        opt.samples = 200000;
        opt.seed = 23000;
        try {
            print_avg_sim(sn_true, ldes::solver_ldes(sn_true, opt));
        } catch (const std::exception& e) {
            na("LDES", e.what());
        }
    }

    // (2) CTMC - exact solution of the true system
    if (gate("SolverCTMC", qn::ctmc_feature_set("exact"), sn_true)) {
        section("CTMC");
        ctmc::CtmcOptions copt;
        copt.method = "exact";
        copt.cutoff = 1;
        print_avg(sn_true, ctmc::solver_ctmc_run_analyzer(sn_true, copt));
    }

    // (3) The random environment: the MMPP2 phase modulates phase-conditional
    // Poisson arrivals, switching at the MMPP2 phase-transition rates, i.e. the
    // -D0 diagonal minus that phase's total arrival rate.
    Net s1 = Build::run(n, m, p1, p2, NULL);
    s1.set_arrival(0, 0, Exp(D11(0, 0)));
    s1.set_arrival(0, 1, Exp(D12(0, 0)));
    Net s2 = Build::run(n, m, p1, p2, NULL);
    s2.set_arrival(0, 0, Exp(D11(1, 1)));
    s2.set_arrival(0, 1, Exp(D12(1, 1)));

    env::Environment<double> e("MMPPphase", 2);
    e.set_stage(0, "Phase1", "bursty", s1.get_struct());
    e.set_stage(1, "Phase2", "calm", s2.get_struct());
    e.add_transition(0, 1, Exp(-D0(0, 0) - (D11(0, 0) + D12(0, 0))));  // 0.5
    e.add_transition(1, 0, Exp(-D0(1, 1) - (D11(1, 1) + D12(1, 1))));  // 0.5

    // 'avg' and 'dec' are the fast- and slow-environment limits; both refuse a
    // Cache stage here by name, because the blend they would report drops the
    // hit and miss ratios the model is asked for.
    for (std::size_t k = 0; k < 2; ++k) {
        const std::string method = k == 0 ? "avg" : "dec";
        env::EnvOptions o;
        o.method = method;
        try {
            const env::EnvAnalyzerSolution<double> r = env::solver_env(e, o);
            section("ENV (" + method + ")");
            (void)r;
        } catch (const std::exception& ex) {
            na("ENV (" + method + ")", ex.what());
        }
    }

    // 'blend' - state-vector coupling: carries the cache-state distribution
    // across phase switches and averages each phase's sojourn-weighted
    // distribution. For a Markovian environment this recovers the exact joint
    // (cache x phase) solution, i.e. it matches CTMC.
    {
        env::EnvOptions o;
        o.method = "blend";
        o.iter_max = 100;
        o.iter_tol = 1e-4;
        o.timespan_end = 1e3;
        try {
            const env::EnvAnalyzerSolution<double> r = env::solver_env(e, o);
            section("ENV (blend)");
            std::map<std::size_t, std::vector<double> >::const_iterator it =
                r.statevec.hit_prob.find(cache_node);
            if (it == r.statevec.hit_prob.end()) {
                note("(the state-vector blend reported no cache hit probability)");
            } else {
                cache_header();
                const CacheInfo ci = cache_info(e.stage(0).model);
                for (std::size_t c = 0; c < 2 && c < it->second.size(); ++c)
                    cache_row(ci.name, e.stage(0).model.classes[c].name, 0, ci.totcap, ci.nitems,
                              it->second[c], kNaN,
                              r.statevec.miss_prob.at(cache_node)[c]);
            }
        } catch (const std::exception& ex) {
            na("ENV (blend)", ex.what());
        }
    }

    na("ENV (mean field over FLD stages)",
       "the reference's default ENV method carries the cache through "
       "aggregateCacheMeanfield_, whose every sweep is a solver_fld_cacheqn_tran call; that "
       "transient cache fluid entry point is not ported, so a Cache in an environment stage "
       "has no mean-field blend here");
}
LINE_EXAMPLE("basic/cacheModel", cache_mmap_rr_env);

// ---------------------------------------------------------------------------
// Delayed-hit retrieval systems
// ---------------------------------------------------------------------------

/** Cache whose misses are fetched by a single infinite-server retrieval station. */
void retrieval_simple() {
    const std::vector<double> access{0.6, 0.3, 0.1};
    const std::size_t n = access.size();
    Net net("Simple Model");
    Source source(net, "Source");
    qn::CacheParam<double> ch;
    ch.nitems = n;
    ch.itemcap = std::vector<int>{1};
    ch.replacestrat = ReplacementStrategy::FIFO;
    ch.pread = std::vector<std::vector<double> >{access, {}, {}};
    ch.hitclass = std::vector<std::size_t>{2, 0, 0};
    ch.missclass = std::vector<std::size_t>{3, 0, 0};
    const std::size_t cache = net.add_cache("Cache", ch);
    Queue queue(net, "Queue", SchedStrategy::INF);
    Sink sink(net, "Sink");

    OpenClass init(net, "InitClass");
    OpenClass hit(net, "HitClass");
    OpenClass miss(net, "MissClass");
    source.set_arrival(init, Exp(1.0));
    queue.set_service(init, Exp(2.0));
    net.set_retrieval_system(cache, init, miss, std::vector<std::size_t>{queue});

    Routing P;
    P.set(init, init, source, cache, 1.0);
    P.set(init, init, cache, queue, 1.0);
    P.set(init, init, queue, cache, 1.0);
    P.set(hit, hit, cache, sink, 1.0);
    P.set(miss, miss, cache, sink, 1.0);
    net.link(P);

    const Sn& sn = net.get_struct();
    retrieval_ssa_cache_table(sn);
    retrieval_ldes_cache_table(sn);
    retrieval_mva_cache_table(sn, init);
    retrieval_nc_tables(sn);
}
LINE_EXAMPLE("basic/cacheModel", retrieval_simple);

/** The processor-sharing variant: one PS retrieval station over seven items. */
void retrieval_ps() {
    const double w[7] = {49, 49, 49, 49, 7, 1, 1};
    std::vector<double> access(7, 0.0);
    for (std::size_t i = 0; i < 7; ++i) access[i] = w[i] / 205.0;
    const std::size_t n = access.size();

    Net net("PS Model");
    Source source(net, "Source");
    qn::CacheParam<double> ch;
    ch.nitems = n;
    ch.itemcap = std::vector<int>{6};
    ch.replacestrat = ReplacementStrategy::RR;
    ch.pread = std::vector<std::vector<double> >{access, {}, {}};
    ch.hitclass = std::vector<std::size_t>{2, 0, 0};
    ch.missclass = std::vector<std::size_t>{3, 0, 0};
    const std::size_t cache = net.add_cache("Cache", ch);
    Queue queue(net, "Queue", SchedStrategy::PS);
    Sink sink(net, "Sink");

    OpenClass init(net, "InitClass");
    OpenClass hit(net, "HitClass");
    OpenClass miss(net, "MissClass");
    source.set_arrival(init, Exp(1.0));
    queue.set_service(init, Exp(1.0));
    net.set_retrieval_system(cache, init, miss, std::vector<std::size_t>{queue});

    Routing P;
    P.set(init, init, source, cache, 1.0);
    P.set(init, init, cache, queue, 1.0);
    P.set(init, init, queue, cache, 1.0);
    P.set(hit, hit, cache, sink, 1.0);
    P.set(miss, miss, cache, sink, 1.0);
    net.link(P);

    const Sn& sn = net.get_struct();
    retrieval_ssa_cache_table(sn);
    retrieval_ldes_cache_table(sn);
    retrieval_mva_cache_table(sn, init);
    retrieval_nc_tables(sn);
}
LINE_EXAMPLE("basic/cacheModel", retrieval_ps);

/** A retrieval chain: Cache -> Queue_1 -> Queue_2 -> Cache, both FCFS. */
void retrieval_chain() {
    const std::vector<double> access{0.6, 0.3, 0.1};
    const std::size_t n = access.size();

    Net net("Chain Model");
    Source source(net, "Source");
    qn::CacheParam<double> ch;
    ch.nitems = n;
    ch.itemcap = std::vector<int>{1};
    ch.replacestrat = ReplacementStrategy::FIFO;
    ch.pread = std::vector<std::vector<double> >{access, {}, {}};
    ch.hitclass = std::vector<std::size_t>{2, 0, 0};
    ch.missclass = std::vector<std::size_t>{3, 0, 0};
    const std::size_t cache = net.add_cache("Cache", ch);
    Queue q1(net, "Queue_1", SchedStrategy::FCFS);
    Queue q2(net, "Queue_2", SchedStrategy::FCFS);
    Sink sink(net, "Sink");

    OpenClass init(net, "InitClass");
    OpenClass hit(net, "HitClass");
    OpenClass miss(net, "MissClass");
    source.set_arrival(init, Exp(1.0));
    q1.set_service(init, Exp(2.0));
    q2.set_service(init, Exp(3.0));
    net.set_retrieval_system(cache, init, miss, std::vector<std::size_t>{q1, q2});

    Routing P;
    P.set(init, init, source, cache, 1.0);
    P.set(init, init, cache, q1, 1.0);
    P.set(init, init, q1, q2, 1.0);
    P.set(init, init, q2, cache, 1.0);
    P.set(hit, hit, cache, sink, 1.0);
    P.set(miss, miss, cache, sink, 1.0);
    net.link(P);

    const Sn& sn = net.get_struct();
    retrieval_ssa_cache_table(sn);
    retrieval_ldes_cache_table(sn);
    retrieval_mva_cache_table(sn, init);
    retrieval_nc_tables(sn);
}
LINE_EXAMPLE("basic/cacheModel", retrieval_chain);

/**
 * The default retrieval topology with per-item overrides: item 1 skips Queue_2
 * and is fetched faster at Queue_1.
 *
 * THE OVERRIDES ARE WRITTEN ON THE RETRIEVAL CLASSES THEMSELVES, not as a
 * template plus a patch. `link()` copies the read class's edges over the queue
 * set onto every retrieval class WHEREVER THEY ARE POSITIVE, so a template edge
 * cannot be deleted for one item afterwards; the reference's
 * `set_item_routing_prob` writes the per-item matrix directly, and giving each
 * retrieval class its own routing here is the same model without a template to
 * contradict.
 */
void retrieval_default() {
    const std::vector<double> access{0.6, 0.3, 0.1};
    const std::size_t n = access.size();

    Net net("DelayedHits");
    Source source(net, "Source");
    qn::CacheParam<double> ch;
    ch.nitems = n;
    ch.itemcap = std::vector<int>{1};
    ch.replacestrat = ReplacementStrategy::FIFO;
    ch.pread = std::vector<std::vector<double> >{access, {}, {}};
    ch.hitclass = std::vector<std::size_t>{2, 0, 0};
    ch.missclass = std::vector<std::size_t>{3, 0, 0};
    const std::size_t cache = net.add_cache("Cache", ch);
    Queue q1(net, "Queue_1", SchedStrategy::PS);
    Queue q2(net, "Queue_2", SchedStrategy::PS);
    Sink sink(net, "Sink");

    OpenClass init(net, "InitClass");
    OpenClass hit(net, "HitClass");
    OpenClass miss(net, "MissClass");
    source.set_arrival(init, Exp(1.0));
    q1.set_service(init, Exp(2.0));
    q2.set_service(init, Exp(3.0));
    net.set_retrieval_system(cache, init, miss, std::vector<std::size_t>{q1, q2});
    const std::vector<std::size_t> rc = retrieval_classes(miss, n);

    q1.set_service(rc[0], Exp(5.0));  // faster item-1 fetch

    Routing P;
    P.set(init, init, source, cache, 1.0);
    P.set(hit, hit, cache, sink, 1.0);
    P.set(miss, miss, cache, sink, 1.0);
    P.set(rc[0], rc[0], cache, q1, 1.0);
    P.set(rc[0], rc[0], q1, cache, 1.0);
    for (std::size_t i = 1; i < n; ++i) {
        P.set(rc[i], rc[i], cache, q1, 1.0);
        P.set(rc[i], rc[i], q1, q2, 1.0);
        P.set(rc[i], rc[i], q2, cache, 1.0);
    }
    net.link(P);

    const Sn& sn = net.get_struct();
    retrieval_ssa_cache_table(sn);
    retrieval_ldes_cache_table(sn);
    retrieval_mva_cache_table(sn, init);
    retrieval_nc_tables(sn);
}
LINE_EXAMPLE("basic/cacheModel", retrieval_default);

/** Per-item probabilistic routing over an IS station and two FCFS queues. */
void retrieval_routing() {
    const std::vector<double> access{0.6, 0.3, 0.1};
    const std::size_t n = access.size();

    Net net("Probabilistic Routing");
    Source source(net, "Source");
    qn::CacheParam<double> ch;
    ch.nitems = n;
    ch.itemcap = std::vector<int>{2};
    ch.replacestrat = ReplacementStrategy::FIFO;
    ch.pread = std::vector<std::vector<double> >{access, {}, {}};
    ch.hitclass = std::vector<std::size_t>{2, 0, 0};
    ch.missclass = std::vector<std::size_t>{3, 0, 0};
    const std::size_t cache = net.add_cache("Cache", ch);
    Queue isq(net, "IS Queue", SchedStrategy::INF);
    Queue q1(net, "Queue 1", SchedStrategy::FCFS);
    Queue q2(net, "Queue 2", SchedStrategy::FCFS);
    Sink sink(net, "Sink");
    const std::vector<std::size_t> queues{isq, q1, q2};

    OpenClass init(net, "InitClass");
    OpenClass hit(net, "HitClass");
    OpenClass miss(net, "MissClass");
    source.set_arrival(init, Exp(1.0));
    isq.set_service(init, Exp(2.0));
    q1.set_service(init, Exp(3.0));
    q2.set_service(init, Exp(3.0));
    net.set_retrieval_system(cache, init, miss, queues);
    const std::vector<std::size_t> rc = retrieval_classes(miss, n);

    // Per-item routing over [IS(0), Queue1(1), Queue2(2), Cache(3)]: row = from,
    // col = to, index 3 being the cache.
    const double R[3][4][4] = {{{0.00, 0.50, 0.00, 0.50},
                                {0.00, 0.00, 0.70, 0.30},
                                {0.00, 0.00, 0.00, 1.00},
                                {0.70, 0.30, 0.00, 0.00}},
                               {{0.00, 0.30, 0.00, 0.70},
                                {0.00, 0.00, 0.50, 0.50},
                                {0.00, 0.00, 0.00, 1.00},
                                {0.20, 0.80, 0.00, 0.00}},
                               {{0.00, 0.60, 0.00, 0.40},
                                {0.00, 0.00, 0.40, 0.60},
                                {0.00, 0.00, 0.00, 1.00},
                                {0.50, 0.50, 0.00, 0.00}}};

    Routing P;
    P.set(init, init, source, cache, 1.0);
    P.set(hit, hit, cache, sink, 1.0);
    P.set(miss, miss, cache, sink, 1.0);
    const std::size_t nq = queues.size();
    for (std::size_t i = 0; i < n; ++i)
        for (std::size_t a = 0; a < nq; ++a) {
            if (R[i][nq][a] > 0.0) P.set(rc[i], rc[i], cache, queues[a], R[i][nq][a]);
            if (R[i][a][nq] > 0.0) P.set(rc[i], rc[i], queues[a], cache, R[i][a][nq]);
            for (std::size_t b = 0; b < nq; ++b)
                if (R[i][a][b] > 0.0) P.set(rc[i], rc[i], queues[a], queues[b], R[i][a][b]);
        }
    net.link(P);

    const Sn& sn = net.get_struct();
    retrieval_ssa_cache_table(sn);
    retrieval_ldes_cache_table(sn);
    retrieval_mva_cache_table(sn, init);
    retrieval_nc_tables(sn);
}
LINE_EXAMPLE("basic/cacheModel", retrieval_routing);

}  // namespace examples
}  // namespace line
