/*
 * Copyright (c) 2012-2026, QORE Lab, Imperial College London
 * All rights reserved.
 */

/**
 * `advanced/passAndSwap`: pass-and-swap (PAS) and order-independent (OI)
 * stations.
 *
 * A PAS station is parameterized by the TOTAL service rate mu(c) of the ordered
 * microstate c and by a swapping graph G saying which class may take another's
 * place when a job completes. An empty G is a plain order-independent queue.
 * Both live on the station: `set_service_rate_function` is the model-level
 * declaration (it also installs the representative per-class Exp(mu([r])) the
 * rate machinery reads), and `set_pas` is what the CTMC state walk reads --
 * `after_event_station_pas` requires `sn.pasparam`, so a model meant for CTMC
 * declares both.
 *
 * INDEX BASE. The reference's mu and G are 0-based Python class indices; the
 * C++ rate function receives the ordered list as 1-BASED class indices, so
 * every mu below subtracts one before indexing its own tables.
 *
 * WHAT IS REFUSED. LDES is not carried by this port, so every simulation
 * comparison is refused by name rather than answered by CTMC under LDES's
 * label. The `save_model` inspection of `pas_saturation` IS ported, on
 * `network_json_envelope`, but its table is materialized over `lattice_cutoffs`
 * -- 10 per open class -- rather than over the saturation cutoffs the reference
 * derives, so the table is constant in the buffer (which is the claim) at a
 * different size (which is a divergence, and said so at the call site).
 *
 * NONE OF THE CTMC NUMBERS BELOW MATCH THE REFERENCE, and the examples say so
 * rather than being trimmed to what works. `state.h`'s `from_marginal_node` has
 * no PAS/OI branch -- `state_detail::buffer_is_class_tag` (state.h:75-85) names
 * FCFS/HOL/LCFS/LCFSPRIO only -- so a PAS station takes the COUNT-shaped
 * default: one column per class, where `after_event_station_pas` reads an
 * ORDERED LIST of width `cap`. Measured:
 *
 *   pas_mmk (open, cap 4, 2 classes)  9 states, not the 31 of the ordered
 *                                     space; MATLAB and native Python both
 *                                     give QLen 0.80328 / 0.57377 against this
 *                                     port's 0.93151 / 0.34932
 *   pas_cyclic (closed)               the initial row is (2, 2) instead of
 *                                     (1, 1, 2, 2), 4 states, nothing departs
 *
 * `pas_closed_tandem_fig6` needs one thing more: a closed PAS network with a
 * non-empty swapping graph is REDUCIBLE, so the initial PLACEMENT selects the
 * recurrent component. The reference says so with `q1.setState([1..6])`, and
 * this port exposes no PAS setState. The two closed examples therefore print
 * their own product-form reference beside the CTMC answer and state the
 * deviation instead of asserting it away.
 */

#include <algorithm>
#include <cmath>
#include <cstdio>
#include <map>
#include <string>
#include <vector>

#include "examples_common.h"
#include "line/api/mc/ctmc_solve.h"
#include "line/io/network_writer.h"
#include "line/solvers/ctmc/solver_ctmc_analyzer.h"
#include "line/solvers/nc/solver_nc_runner.h"

namespace line {
namespace examples {

namespace {

using PasRate = std::function<double(const std::vector<std::size_t>&)>;
using BoolGraph = std::vector<std::vector<bool>>;

/** The swapping graph as an undirected edge list over 0-based classes. */
struct Edge {
    std::size_t a, b;
};

BoolGraph swap_bool(std::size_t R, const std::vector<Edge>& edges, bool symmetric = true) {
    BoolGraph g(R, std::vector<bool>(R, false));
    for (std::size_t i = 0; i < edges.size(); ++i) {
        g[edges[i].a][edges[i].b] = true;
        if (symmetric) g[edges[i].b][edges[i].a] = true;
    }
    return g;
}

Matrix<double> swap_matrix(const BoolGraph& g) {
    Matrix<double> G(g.size(), g.size(), 0.0);
    for (std::size_t i = 0; i < g.size(); ++i)
        for (std::size_t j = 0; j < g.size(); ++j) G(i, j) = g[i][j] ? 1.0 : 0.0;
    return G;
}

/** Declare mu(c) and G on a station, for both the model and the CTMC walk. */
void declare_pas(Net& m, std::size_t node, const PasRate& mu, const BoolGraph& g) {
    m.set_service_rate_function(node, mu, swap_matrix(g));
    m.set_pas(node, mu, g);
}

/** The AvgTable of `CTMC(model, cutoff=k)`, with the PAS state-width caveat. */
void ctmc_avg(const Sn& sn, double cutoff) {
    note("NOTE: a PAS station's state space is built count-shaped here, not as the ordered "
         "list;");
    note("      these numbers do not match MATLAB (see the file header).");
    ctmc::CtmcOptions o;
    o.cutoff = cutoff;
    print_avg(sn, ctmc::solver_ctmc_run_analyzer(sn, o));
}

/** An open Source -> PASQueue -> Sink model with `R` classes at rates `lam`. */
Net pas_open_model(const std::string& nm, const std::vector<double>& lam, const PasRate& mu,
                   const BoolGraph& g, double nservers, double cap) {
    Net m(nm);
    Source src(m, "Source");
    Queue q(m, "PASQueue", SchedStrategy::PAS);
    Sink snk(m, "Sink");
    std::vector<std::size_t> cls;
    for (std::size_t r = 0; r < lam.size(); ++r)
        cls.push_back(m.add_open_class("Class" + std::to_string(r + 1)));
    for (std::size_t r = 0; r < lam.size(); ++r)
        src.set_arrival(cls[r], Exp(lam[r]));
    declare_pas(m, q, mu, g);
    q.set_number_of_servers(nservers);
    q.set_capacity(cap);
    Routing P;
    for (std::size_t r = 0; r < lam.size(); ++r) serial(P, cls[r], {src, q, snk});
    m.link(P);
    return m;
}

/** The QLen column of an AvgResult at one station. */
std::vector<double> qlen_at(const Sn& sn, const mva::AvgResult<double>& r, const std::string& st) {
    std::vector<double> q;
    for (std::size_t i = 0; i < sn.nstations; ++i)
        if (sn.stations[i].name == st)
            for (std::size_t c = 0; c < sn.nclasses; ++c) q.push_back(r.QN(i, c));
    return q;
}

}  // namespace

// ---------------------------------------------------------------------------
// pas_mmk
// ---------------------------------------------------------------------------

/**
 * A PAS station with an empty swapping graph is a plain OI queue, and with K
 * unit-rate servers mu(c) = min(|c|, K) is the M/M/K rate.
 */
void pas_mmk() {
    const std::size_t K = 2;
    const PasRate mu = [K](const std::vector<std::size_t>& c) {
        return static_cast<double>(std::min(c.size(), K));
    };
    Net m = pas_open_model("PASmmk", {0.7, 0.5}, mu, swap_bool(2, {}), 2.0, 4.0);
    section("CTMC");
    ctmc_avg(m.get_struct(), 4.0);
}

// ---------------------------------------------------------------------------
// pas_compatibility_5class
// ---------------------------------------------------------------------------

/**
 * The five-class, three-server compatibility queue of Dorsman & Gardner (2024),
 * Figs 1-2: mu(c) is the number of servers compatible with a class present in
 * c, and the product form is invariant to the swapping graph.
 */
void pas_compatibility_5class() {
    // comp[s][r] = 1 iff server s serves class r (0-based classes).
    static const int comp[3][5] = {{1, 0, 0, 1, 0}, {0, 1, 0, 1, 0}, {0, 0, 1, 0, 1}};
    const PasRate mu = [](const std::vector<std::size_t>& c) {
        double n = 0.0;
        for (std::size_t s = 0; s < 3; ++s) {
            bool any = false;
            for (std::size_t i = 0; i < c.size(); ++i)
                if (comp[s][c[i] - 1]) any = true;
            if (any) n += 1.0;
        }
        return n;
    };
    const std::vector<Edge> edges{{0, 2}, {0, 4}, {1, 3}, {2, 3}, {3, 4}};
    Net m = pas_open_model("PAScompatibility", {0.5, 0.4, 0.3, 0.2, 0.1}, mu,
                           swap_bool(5, edges), 3.0, 3.0);
    section("CTMC");
    ctmc_avg(m.get_struct(), 3.0);
}

// ---------------------------------------------------------------------------
// pas_selfloop
// ---------------------------------------------------------------------------

/**
 * A swapping graph with a self-loop, which Dorsman & Gardner (2024) Sect. 2.3
 * permits: a completing class-i job may take the place of another class-i job
 * further back in the queue.
 */
void pas_selfloop() {
    // Three unit-rate servers, one per class: mu(c) = number of distinct classes.
    const PasRate mu = [](const std::vector<std::size_t>& c) {
        bool seen[3] = {false, false, false};
        for (std::size_t i = 0; i < c.size(); ++i) seen[c[i] - 1] = true;
        return static_cast<double>((seen[0] ? 1 : 0) + (seen[1] ? 1 : 0) + (seen[2] ? 1 : 0));
    };
    BoolGraph g = swap_bool(3, {});
    g[0][0] = true;  // the self-loop
    g[1][2] = true;
    g[2][1] = true;
    Net m = pas_open_model("PASselfloop", {0.6, 0.4, 0.3}, mu, g, 3.0, 3.0);
    section("CTMC");
    ctmc_avg(m.get_struct(), 3.0);
}

// ---------------------------------------------------------------------------
// pas_saturation
// ---------------------------------------------------------------------------

/**
 * Saturation of an order-independent rate function.
 *
 * mu saturates once every present class has one job, so the rate table an
 * engine has to carry stays compact whatever the buffer is. The reference shows
 * that by SERIALIZING the model and counting the table entries, and by matching
 * a large-buffer LDES run against CTMC; neither the writer nor LDES exists in
 * this port, so what is ported is the CTMC answer at the same buffer and the
 * saturation itself, exhibited directly on mu.
 */
void pas_saturation() {
    // comp[s][r] = 1 iff server s serves class r; cap_s is its capacity.
    static const int comp[3][2] = {{1, 0}, {1, 1}, {0, 1}};
    static const double cap_s[3] = {2.0, 1.0, 2.0};
    const PasRate mu_rank = [](const std::vector<std::size_t>& c) {
        double s = 0.0;
        for (std::size_t k = 0; k < 3; ++k) {
            bool any = false;
            for (std::size_t i = 0; i < c.size(); ++i)
                if (comp[k][c[i] - 1]) any = true;
            if (any) s += cap_s[k];
        }
        return s;
    };

    // save_model: mu(c) crosses the wire as a rate table over a box lattice,
    // and the table size is constant in the buffer, which is the reference's
    // point. It is NOT the reference's SIZE: `lattice_cutoffs` materializes an
    // open class over a fixed box of 10, so the cutoffs are [10, 10] and the
    // table 120 entries where the reference's saturation cutoffs give [1, 1]
    // and 3.
    note("the order-independent rate table crosses the wire at a size independent of the buffer:");
    const double caps[3] = {8.0, 100.0, 1000.0};
    for (std::size_t k = 0; k < 3; ++k) {
        const double cap = caps[k];
        Net w = pas_open_model("PASsaturation", {1.5, 1.0}, mu_rank, swap_bool(2, {}), 3.0, cap);
        const io::detail::json j = io::network_json_envelope(w.get_struct());
        for (io::detail::json::const_iterator n = j.at("model").at("nodes").begin();
             n != j.at("model").at("nodes").end(); ++n)
            if (n->contains("oiServiceRate"))
                std::printf("  cap=%-6g table entries = %zu  cutoffs = %s\n", cap,
                            n->at("oiServiceRate").size(), n->at("oiCutoffs").dump().c_str());
    }

    const double CAP = 8.0;
    Net m = pas_open_model("PASsaturation", {1.5, 1.0}, mu_rank, swap_bool(2, {}), 3.0, CAP);
    const Sn& sn = m.get_struct();
    ctmc::CtmcOptions o;
    o.cutoff = CAP;
    const mva::AvgResult<double> r = ctmc::solver_ctmc_run_analyzer(sn, o);
    std::printf("\nopen compatibility PAS, buffer=%g:\n", CAP);
    section("CTMC");
    note("NOTE: the PAS state space is count-shaped here; see the file header.");
    print_avg(sn, r);
    const std::vector<double> q = qlen_at(sn, r, "PASQueue");
    for (std::size_t i = 0; i < q.size(); ++i)
        std::printf("  CTMC QLen Class%zu = %.5f\n", i + 1, q[i]);

    // TODO(cpp): Ql = qlen(LDES(build(mu_rank, 8), samples=300000, seed=23000).getAvgTable())
    // TODO(cpp): assert max(abs(Qc - Ql) / Qc) <= 0.03
    // TODO(cpp): for mu in (mu_rank, mu_additive): LDES(build(mu, None), samples=1000, seed=1)
    // TODO(cpp):     -> a clean, actionable finite-buffer error naming the PAS station
    na("LDES",
       "the reference matches a 3e5-sample LDES run against CTMC at buffer 8 and then checks the "
       "actionable error an unset (infinite) buffer raises; this port carries no LDES engine, and "
       "answering a simulation block with the CTMC number would be a silent solver substitution");
}

// ---------------------------------------------------------------------------
// pas_cyclic_ctmc_vs_bruteforce
// ---------------------------------------------------------------------------

namespace {

const double kCycS1 = 1.0;
const std::size_t kCycK1 = 2;
const double kCycBeta2[2] = {1.5, 1.0};
const std::size_t kCycPop[2] = {2, 2};

double cyc_mu1(const std::vector<std::size_t>& c) {
    return static_cast<double>(std::min(c.size(), kCycK1)) * kCycS1;
}

double cyc_mu2(const std::vector<std::size_t>& c) {
    double s = 0.0;
    for (std::size_t i = 0; i < c.size(); ++i) s += kCycBeta2[c[i] - 1];
    return s;
}

/**
 * The OI balance function summed over every distinct ORDERED arrangement of a
 * class-count vector: sum over orders of prod_j e_{c_j} / mu(c_1..c_j).
 *
 * The prefix sums depend on the order at station 2, so the enumeration is
 * genuinely over orders and does not collapse to class counts.
 */
double oi_dfs(std::vector<std::size_t>& counts, std::vector<std::size_t>& prefix, double w,
              const std::vector<double>& e, double (*mu)(const std::vector<std::size_t>&)) {
    bool empty = true;
    for (std::size_t r = 0; r < counts.size(); ++r)
        if (counts[r] > 0) empty = false;
    if (empty) return w;
    double s = 0.0;
    for (std::size_t r = 0; r < counts.size(); ++r) {
        if (counts[r] == 0) continue;
        counts[r] -= 1;
        prefix.push_back(r + 1);
        s += oi_dfs(counts, prefix, w * e[r] / mu(prefix), e, mu);
        prefix.pop_back();
        counts[r] += 1;
    }
    return s;
}

double oi_sum(const std::vector<std::size_t>& counts, const std::vector<double>& e,
              double (*mu)(const std::vector<std::size_t>&)) {
    std::vector<std::size_t> c = counts;
    std::vector<std::size_t> prefix;
    return oi_dfs(c, prefix, 1.0, e, mu);
}

std::vector<std::vector<std::size_t>> enum_splits(const std::vector<std::size_t>& K) {
    std::vector<std::vector<std::size_t>> out(1, std::vector<std::size_t>(K.size(), 0));
    for (std::size_t r = 0; r < K.size(); ++r) {
        std::vector<std::vector<std::size_t>> next;
        for (std::size_t i = 0; i < out.size(); ++i)
            for (std::size_t v = 0; v <= K[r]; ++v) {
                std::vector<std::size_t> row = out[i];
                row[r] = v;
                next.push_back(row);
            }
        out = next;
    }
    return out;
}

double cyc_norm_const(const std::vector<std::size_t>& K, const std::vector<double>& e1,
                      const std::vector<double>& e2) {
    double G = 0.0;
    const std::vector<std::vector<std::size_t>> splits = enum_splits(K);
    for (std::size_t i = 0; i < splits.size(); ++i) {
        std::vector<std::size_t> comp(K.size(), 0);
        for (std::size_t r = 0; r < K.size(); ++r) comp[r] = K[r] - splits[i][r];
        G += oi_sum(splits[i], e1, cyc_mu1) * oi_sum(comp, e2, cyc_mu2);
    }
    return G;
}

}  // namespace

/**
 * A closed cycle of two OI stations against the exact product-form normalizing
 * constant, enumerated by brute force over the whole ORDERED state space.
 */
void pas_cyclic_ctmc_vs_bruteforce() {
    const std::size_t R = 2;
    const std::vector<std::size_t> K(kCycPop, kCycPop + R);
    const std::vector<double> e1(R, 1.0), e2(R, 1.0);

    const double Gbf = cyc_norm_const(K, e1, e2);
    std::printf("Brute-force normalizing constant G = %.12g\n\n", Gbf);

    std::vector<double> Xbf(R, 0.0);
    for (std::size_t r = 0; r < R; ++r) {
        std::vector<std::size_t> Km = K;
        Km[r] -= 1;
        Xbf[r] = cyc_norm_const(Km, e1, e2) / Gbf;
    }
    std::vector<std::vector<double>> Qbf(2, std::vector<double>(R, 0.0));
    const std::vector<std::vector<std::size_t>> splits = enum_splits(K);
    for (std::size_t i = 0; i < splits.size(); ++i) {
        std::vector<std::size_t> comp(R, 0);
        for (std::size_t r = 0; r < R; ++r) comp[r] = K[r] - splits[i][r];
        const double w = oi_sum(splits[i], e1, cyc_mu1) * oi_sum(comp, e2, cyc_mu2);
        for (std::size_t r = 0; r < R; ++r) Qbf[0][r] += static_cast<double>(splits[i][r]) * w;
    }
    for (std::size_t r = 0; r < R; ++r) {
        Qbf[0][r] /= Gbf;
        Qbf[1][r] = static_cast<double>(K[r]) - Qbf[0][r];
    }

    double total = 0.0;
    for (std::size_t r = 0; r < R; ++r) total += static_cast<double>(K[r]);
    Net m("PAScyclic");
    Queue q1(m, "PASQueue1", SchedStrategy::PAS);
    Queue q2(m, "PASQueue2", SchedStrategy::PAS);
    std::vector<std::size_t> cls;
    for (std::size_t r = 0; r < R; ++r)
        cls.push_back(m.add_closed_class("Class" + std::to_string(r + 1),
                                         static_cast<double>(K[r]), q1));
    declare_pas(m, q1, cyc_mu1, swap_bool(R, {}));
    declare_pas(m, q2, cyc_mu2, swap_bool(R, {}));
    q1.set_number_of_servers(static_cast<double>(kCycK1));
    q1.set_capacity(total);
    q2.set_number_of_servers(total);
    q2.set_capacity(total);
    Routing P;
    for (std::size_t r = 0; r < R; ++r) cyclic(P, cls[r], {q1, q2});
    m.link(P);

    const Sn& sn = m.get_struct();
    ctmc::CtmcOptions o;
    o.cutoff = total;
    const mva::AvgResult<double> r = ctmc::solver_ctmc_run_analyzer(sn, o);
    section("CTMC");
    print_avg(sn, r);

    const std::vector<double> Qc1 = qlen_at(sn, r, "PASQueue1");
    const std::vector<double> Qc2 = qlen_at(sn, r, "PASQueue2");
    std::vector<double> Xc(R, 0.0);
    for (std::size_t i = 0; i < sn.nstations; ++i)
        if (sn.stations[i].name == "PASQueue1")
            for (std::size_t c = 0; c < R; ++c) Xc[c] = r.TN(i, c);

    note("\n--- CTMC vs brute-force ---");
    std::printf("%-8s %-8s %14s %14s %10s\n", "metric", "class", "brute-force", "CTMC", "abs.err");
    double errX = 0.0, errQ = 0.0;
    for (std::size_t c = 0; c < R; ++c) {
        std::printf("Tput     Class%-2zu %14.10f %14.10f %10.2e\n", c + 1, Xbf[c], Xc[c],
                    std::fabs(Xbf[c] - Xc[c]));
        errX = std::max(errX, std::fabs(Xbf[c] - Xc[c]));
    }
    for (std::size_t c = 0; c < R; ++c) {
        const double v = c < Qc1.size() ? Qc1[c] : 0.0;
        std::printf("QLen q1   Class%-2zu %14.10f %14.10f %10.2e\n", c + 1, Qbf[0][c], v,
                    std::fabs(Qbf[0][c] - v));
        errQ = std::max(errQ, std::fabs(Qbf[0][c] - v));
    }
    for (std::size_t c = 0; c < R; ++c) {
        const double v = c < Qc2.size() ? Qc2[c] : 0.0;
        std::printf("QLen q2   Class%-2zu %14.10f %14.10f %10.2e\n", c + 1, Qbf[1][c], v,
                    std::fabs(Qbf[1][c] - v));
        errQ = std::max(errQ, std::fabs(Qbf[1][c] - v));
    }
    std::printf("\nCTMC:  max|dX| = %.3e   max|dQ| = %.3e\n", errX, errQ);
    note(errX <= 1e-9 && errQ <= 1e-9
             ? "PASS: CTMC matches brute-force product form within 1.0e-09."
             : "MISMATCH: the closed PAS walk starts from a COUNT-shaped initial row, so it "
               "never leaves it (see the file header).");

    // TODO(cpp): ldes = LDES(build_model(), samples=200000, seed=23000)
    // TODO(cpp): Xl, Ql = extract(ldes.getAvgTable())
    // TODO(cpp): assert max rel|dX| <= 0.02 and max rel|dQ| <= 0.02 against the brute force
    na("LDES",
       "the reference then repeats the comparison with a 2e5-sample simulation at 2% tolerance; "
       "this port carries no LDES engine");
}

// ---------------------------------------------------------------------------
// pas_cyclic_normconst_bruteforce
// ---------------------------------------------------------------------------

namespace {

const double kGnMu1 = 3.0;
const double kGnMu2 = 2.0;

/** `pas_rate(c, mu1, 1)`: a single-server PAS station serves only the head. */
double gn_mu1(const std::vector<std::size_t>& c) { return c.empty() ? 0.0 : kGnMu1; }

double gn_mu2(const std::vector<std::size_t>& c) { return c.empty() ? 0.0 : kGnMu2; }

/** `norm_const` over an arbitrary pair of OI rate functions. */
double cyc_norm_const_mu(const std::vector<std::size_t>& K, const std::vector<double>& e1,
                         const std::vector<double>& e2,
                         double (*mu1)(const std::vector<std::size_t>&),
                         double (*mu2)(const std::vector<std::size_t>&)) {
    double G = 0.0;
    const std::vector<std::vector<std::size_t>> splits = enum_splits(K);
    for (std::size_t i = 0; i < splits.size(); ++i) {
        std::vector<std::size_t> comp(K.size(), 0);
        for (std::size_t r = 0; r < K.size(); ++r) comp[r] = K[r] - splits[i][r];
        G += oi_sum(splits[i], e1, mu1) * oi_sum(comp, e2, mu2);
    }
    return G;
}

}  // namespace

/**
 * The brute-force normalizing constant of a closed cycle of two PAS stations,
 * on its own: G by exact enumeration of the ordered state space, the per-class
 * throughputs X_r = G(K - 1_r)/G(K) it implies, and the single-class reduction
 * to the Gordon-Newell geometric sum as a self-test.
 *
 * No solver runs here. The reference's point is that the OI product form is
 * INVARIANT to the swap graph, so the constant is computable from mu and the
 * visit ratios alone, which is what `pas_nc_sampling` then estimates.
 */
void pas_cyclic_normconst_bruteforce() {
    const std::size_t R = 2;
    const std::vector<std::size_t> K(kCycPop, kCycPop + R);
    const std::vector<double> e1(R, 1.0), e2(R, 1.0);

    const double G = cyc_norm_const(K, e1, e2);
    std::printf("Closed cyclic PAS network: R=%zu classes, population K=[", R);
    for (std::size_t r = 0; r < R; ++r) std::printf("%s%zu", r ? " " : "", K[r]);
    std::printf("]\n");
    std::printf("Brute-force normalizing constant G = %.15g\n", G);

    std::vector<double> X(R, 0.0);
    for (std::size_t r = 0; r < R; ++r) {
        if (K[r] == 0) continue;
        std::vector<std::size_t> Km = K;
        Km[r] -= 1;
        X[r] = cyc_norm_const(Km, e1, e2) / G;
    }
    std::printf("Per-class throughput X = [");
    for (std::size_t r = 0; r < R; ++r) std::printf("%s%.6g", r ? " " : "", X[r]);
    std::printf("]\n");

    // Single class, two single-server PAS stations: the closed form is the
    // geometric sum G = sum_n (1/mu1)^n (1/mu2)^(N-n).
    const std::size_t N = 5;
    const std::vector<std::size_t> Kn(1, N);
    const std::vector<double> en(1, 1.0);
    const double Gbf = cyc_norm_const_mu(Kn, en, en, gn_mu1, gn_mu2);
    double Gcf = 0.0;
    for (std::size_t n = 0; n <= N; ++n)
        Gcf += std::pow(1.0 / kGnMu1, static_cast<double>(n)) *
               std::pow(1.0 / kGnMu2, static_cast<double>(N - n));
    if (std::fabs(Gbf - Gcf) > 1e-12 * Gcf)
        throw NumericError("self-test failed: brute-force " + std::to_string(Gbf) +
                           " vs closed-form " + std::to_string(Gcf));
    std::printf("Self-test (single-class Gordon-Newell): PASS (G=%.12g)\n", Gbf);
}

// ---------------------------------------------------------------------------
// pas_nc_sampling
// ---------------------------------------------------------------------------

namespace {

const double kIsMu1 = 1.0;
const double kIsMu2 = 1.3;

/** The head-only single server of the Fig. 6 tandem: mu(c) = mu, independent of c. */
double is_mu1(const std::vector<std::size_t>& c) { return c.empty() ? 0.0 : kIsMu1; }

double is_mu2(const std::vector<std::size_t>& c) { return c.empty() ? 0.0 : kIsMu2; }

/** The Fig. 6a swap graph: the 7 undirected edges over the six classes. */
std::vector<std::vector<bool>> fig6_swap_graph() {
    static const std::size_t edges[7][2] = {{1, 3}, {1, 4}, {2, 4}, {2, 5}, {3, 6}, {4, 6}, {5, 6}};
    std::vector<std::vector<bool>> G(6, std::vector<bool>(6, false));
    for (std::size_t k = 0; k < 7; ++k) {
        G[edges[k][0] - 1][edges[k][1] - 1] = true;
        G[edges[k][1] - 1][edges[k][0] - 1] = true;
    }
    return G;
}

}  // namespace

/**
 * The closed pass-and-swap tandem of Comte and Dorsman (2021, Fig. 6) solved by
 * SolverNC's importance sampling (`pfqn_pas_is`) against the exact CTMC.
 *
 * A non-empty swap graph makes the ordered-state chain reducible; the recurrent
 * class carries the per-class product form pi(c) = Phi_1 Phi_2 / G_C, which
 * auto-normalized importance sampling estimates without enumerating it.
 */
void pas_nc_sampling() {
    const std::size_t R = 6;
    const std::vector<std::vector<bool>> G = fig6_swap_graph();

    Net m("PASsampling");
    Queue q1(m, "PASQueue1", SchedStrategy::PAS);
    Queue q2(m, "PASQueue2", SchedStrategy::PAS);
    std::vector<std::size_t> cls;
    for (std::size_t r = 0; r < R; ++r)
        cls.push_back(m.add_closed_class("Class" + std::to_string(r + 1), 1.0, q1));
    declare_pas(m, q1, is_mu1, G);
    declare_pas(m, q2, is_mu2, G);
    q1.set_number_of_servers(1.0);
    q1.set_capacity(static_cast<double>(R));
    q2.set_number_of_servers(1.0);
    q2.set_capacity(static_cast<double>(R));
    Routing P;
    for (std::size_t r = 0; r < R; ++r) cyclic(P, cls[r], {q1, q2});
    m.link(P);

    const Sn& sn = m.get_struct();

    std::printf("=== CTMC (exact) ===\n");
    ctmc::CtmcOptions co;
    co.cutoff = static_cast<double>(R);
    const mva::AvgResult<double> rc = ctmc::solver_ctmc_run_analyzer(sn, co);
    print_avg(sn, rc);

    std::printf("=== SolverNC method='sampling' (importance sampling, 5e5 samples) ===\n");
    nc::NcSolverOptions no;
    no.method = "sampling";
    no.samples = 500000;
    no.seed = 777;
    const mva::AvgResult<double> rn = nc::solver_nc_run_analyzer(sn, no);
    print_avg(sn, rn);

    std::printf("\n  station     class    CTMC      NC-samp\n");
    double err = 0.0;
    for (std::size_t i = 0; i < sn.nstations; ++i)
        for (std::size_t c = 0; c < sn.nclasses; ++c) {
            const double qc = rc.QN(i, c), qn = rn.QN(i, c);
            std::printf("  %-9s  %-6s %9.5f %9.5f\n", sn.stations[i].name.c_str(),
                        sn.classes[c].name.c_str(), qc, qn);
            err = std::max(err, std::fabs(qn - qc));
        }
    std::printf("\nNC-samp vs CTMC:  max|dQ| = %.3e (importance-sampling noise)\n", err);
    note(err <= 5e-2
             ? "PASS: SolverNC 'sampling' matches the exact CTMC within importance-sampling noise."
             : "MISMATCH: the closed PAS CTMC walk starts from a COUNT-shaped initial row and "
               "never leaves it (see the file header), so the exact side of this comparison is "
               "not the reference's.");
}

// ---------------------------------------------------------------------------
// pas_closed_tandem_fig6
// ---------------------------------------------------------------------------

namespace {

const double kFigMu1 = 1.0, kFigMu2 = 1.3;

BoolGraph fig5_graph() {
    const std::vector<Edge> edges{{0, 2}, {0, 3}, {1, 3}, {1, 4}, {2, 5}, {3, 5}, {4, 5}};
    return swap_bool(6, edges);
}

/**
 * The pass-and-swap scan, written out independently of the library's own:
 * the job at position p completes, the classes shift one step along the swap
 * chain, and the LAST link of the chain is what actually departs.
 */
std::pair<std::vector<std::size_t>, std::size_t> ps_algorithm(const std::vector<std::size_t>& lst,
                                                              std::size_t p, const BoolGraph& G) {
    std::vector<std::size_t> chain(1, p);
    std::size_t cur = p;
    while (true) {
        std::size_t nxt = lst.size();
        for (std::size_t j = cur + 1; j < lst.size(); ++j)
            if (G[lst[cur] - 1][lst[j] - 1]) {
                nxt = j;
                break;
            }
        if (nxt == lst.size()) break;
        chain.push_back(nxt);
        cur = nxt;
    }
    const std::size_t dep = lst[chain.back()];
    std::vector<std::size_t> nw = lst;
    for (std::size_t i = 0; i + 1 < chain.size(); ++i) nw[chain[i + 1]] = lst[chain[i]];
    nw.erase(nw.begin() + static_cast<long>(chain[0]));
    return std::make_pair(nw, dep);
}

using FigState = std::pair<std::vector<std::size_t>, std::vector<std::size_t>>;

/** The product-form reference: the Markov chain built straight from the scan. */
void fig6_reference(std::vector<double>& Q1, std::vector<double>& Q2) {
    const BoolGraph G = fig5_graph();
    FigState init;
    for (std::size_t r = 1; r <= 6; ++r) init.first.push_back(r);
    std::vector<FigState> states(1, init);
    std::map<FigState, std::size_t> index;
    index[init] = 0;
    std::vector<std::size_t> from, to;
    std::vector<double> rate;
    for (std::size_t f = 0; f < states.size(); ++f) {
        const std::vector<std::size_t> l1 = states[f].first, l2 = states[f].second;
        if (!l1.empty()) {
            const std::pair<std::vector<std::size_t>, std::size_t> s = ps_algorithm(l1, 0, G);
            FigState t(s.first, l2);
            t.second.push_back(s.second);
            std::map<FigState, std::size_t>::iterator it = index.find(t);
            if (it == index.end()) {
                it = index.insert(std::make_pair(t, states.size())).first;
                states.push_back(t);
            }
            from.push_back(f);
            to.push_back(it->second);
            rate.push_back(kFigMu1);
        }
        if (!l2.empty()) {
            const std::pair<std::vector<std::size_t>, std::size_t> s = ps_algorithm(l2, 0, G);
            FigState t(l1, s.first);
            t.first.push_back(s.second);
            std::map<FigState, std::size_t>::iterator it = index.find(t);
            if (it == index.end()) {
                it = index.insert(std::make_pair(t, states.size())).first;
                states.push_back(t);
            }
            from.push_back(f);
            to.push_back(it->second);
            rate.push_back(kFigMu2);
        }
    }
    const std::size_t ns = states.size();
    Matrix<double> Q(ns, ns, 0.0);
    for (std::size_t k = 0; k < from.size(); ++k) Q(from[k], to[k]) += rate[k];
    for (std::size_t i = 0; i < ns; ++i) {
        double s = 0.0;
        for (std::size_t j = 0; j < ns; ++j)
            if (j != i) s += Q(i, j);
        Q(i, i) = -s;
    }
    const std::vector<double> pi = mc::ctmc_solve(Q);
    Q1.assign(6, 0.0);
    Q2.assign(6, 0.0);
    for (std::size_t s = 0; s < ns; ++s)
        for (std::size_t r = 1; r <= 6; ++r) {
            std::size_t n1 = 0, n2 = 0;
            for (std::size_t i = 0; i < states[s].first.size(); ++i)
                if (states[s].first[i] == r) ++n1;
            for (std::size_t i = 0; i < states[s].second.size(); ++i)
                if (states[s].second[i] == r) ++n2;
            Q1[r - 1] += pi[s] * static_cast<double>(n1);
            Q2[r - 1] += pi[s] * static_cast<double>(n2);
        }
}

}  // namespace

/**
 * The closed tandem of two pass-and-swap queues of Comte & Dorsman (2021),
 * Figs 5-6: six classes of one job each, both queues carrying the Fig. 5
 * swapping graph and head-only service.
 */
void pas_closed_tandem_fig6() {
    const BoolGraph G = fig5_graph();
    const PasRate mu1 = [](const std::vector<std::size_t>&) { return kFigMu1; };
    const PasRate mu2 = [](const std::vector<std::size_t>&) { return kFigMu2; };

    Net m("PASclosedTandem");
    Queue q1(m, "PASQueue1", SchedStrategy::PAS);
    Queue q2(m, "PASQueue2", SchedStrategy::PAS);
    std::vector<std::size_t> cls;
    for (std::size_t r = 0; r < 6; ++r)
        cls.push_back(m.add_closed_class("Class" + std::to_string(r + 1), 1.0, q1));
    declare_pas(m, q1, mu1, G);
    declare_pas(m, q2, mu2, G);
    q1.set_number_of_servers(1.0);
    q1.set_capacity(6.0);
    q2.set_number_of_servers(1.0);
    q2.set_capacity(6.0);
    Routing P;
    for (std::size_t r = 0; r < 6; ++r) cyclic(P, cls[r], {q1, q2});
    m.link(P);

    note("NOTE: the reference selects the recurrent component with "
         "q1.setState([1,2,3,4,5,6]);");
    note("      this port has no PAS setState, and a PAS station's state space is built");
    note("      count-shaped rather than as the ordered list the walk reads.");

    const Sn& sn = m.get_struct();
    ctmc::CtmcOptions o;
    o.cutoff = 6.0;
    const mva::AvgResult<double> r = ctmc::solver_ctmc_run_analyzer(sn, o);
    section("CTMC");
    print_avg(sn, r);

    std::vector<double> Qr1, Qr2;
    fig6_reference(Qr1, Qr2);
    const std::vector<double> Qc1 = qlen_at(sn, r, "PASQueue1");
    const std::vector<double> Qc2 = qlen_at(sn, r, "PASQueue2");

    note("\n  station    class   reference     CTMC");
    double err = 0.0;
    for (std::size_t i = 0; i < 6; ++i) {
        const double v = i < Qc1.size() ? Qc1[i] : 0.0;
        std::printf("  %-9s  Class%-2zu %9.5f %9.5f\n", "PASQueue1", i + 1, Qr1[i], v);
        err = std::max(err, std::fabs(Qr1[i] - v));
    }
    for (std::size_t i = 0; i < 6; ++i) {
        const double v = i < Qc2.size() ? Qc2[i] : 0.0;
        std::printf("  %-9s  Class%-2zu %9.5f %9.5f\n", "PASQueue2", i + 1, Qr2[i], v);
        err = std::max(err, std::fabs(Qr2[i] - v));
    }
    std::printf("\nCTMC vs reference: max|dQ| = %.3e\n", err);
    note(err <= 1e-6 ? "PASS: CTMC matches the pass-and-swap reference."
                     : "MISMATCH: the CTMC walk never left its initial placement (see the note "
                       "above).");
}

LINE_EXAMPLE("advanced/passAndSwap", pas_mmk);
LINE_EXAMPLE("advanced/passAndSwap", pas_compatibility_5class);
LINE_EXAMPLE("advanced/passAndSwap", pas_selfloop);
LINE_EXAMPLE("advanced/passAndSwap", pas_saturation);
LINE_EXAMPLE("advanced/passAndSwap", pas_cyclic_ctmc_vs_bruteforce);
LINE_EXAMPLE("advanced/passAndSwap", pas_cyclic_normconst_bruteforce);
LINE_EXAMPLE("advanced/passAndSwap", pas_nc_sampling);
LINE_EXAMPLE("advanced/passAndSwap", pas_closed_tandem_fig6);

}  // namespace examples
}  // namespace line
