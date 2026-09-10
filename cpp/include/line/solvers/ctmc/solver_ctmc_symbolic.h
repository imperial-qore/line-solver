/*
 * Copyright (c) 2012-2026, QORE Lab, Imperial College London
 * All rights reserved.
 */
#ifndef LINE_SOLVERS_CTMC_SOLVER_CTMC_SYMBOLIC_H
#define LINE_SOLVERS_CTMC_SOLVER_CTMC_SYMBOLIC_H

/**
 * Port of `@@SolverCTMC/getSymbolicGenerator` and `getSymbolicSolution`.
 *
 * THE TWO HALVES ARE NOT THE SAME KIND OF PROBLEM, which is why they live in one
 * header but only one of them needs a backend. The generator is LINEAR in the
 * event symbols: each synchronization contributes one numeric filtration,
 * normalized by its own minimum positive rate, scaled by its symbol x1..xE. So
 * the symbolic generator is a sum of numeric matrices with symbolic
 * coefficients, assembled here with no algebra system in sight and printed as
 * expression strings at the end. Solving pi Q = 0 with it is a linear solve
 * whose pivots are multivariate polynomials, which needs exact arithmetic with
 * cancellation in a rational function field; that half is delegated to the
 * backend `api/sym` resolves, exactly as the reference delegates it to the
 * Symbolic Math Toolbox or to line-sage-rest.
 *
 * WHY THE FILTRATION IS NORMALIZED. Dividing event e's filtration by its
 * smallest positive rate makes the nominal value of x_e that rate, so the
 * printed coefficients are RATIOS within one event and are 1 wherever the event
 * fires at its base rate. Substituting the minimum positive rates therefore
 * reproduces the numeric generator; `ctmc_symbolic_eval_infgen` does exactly
 * that and is the check to run before trusting a printed normal form.
 *
 * COEFFICIENT TEXT IS READ AS AN EXACT RATIONAL SERVER SIDE, so how a
 * non-integer coefficient is printed changes the answer, not merely its
 * appearance. This port prints the SHORTEST decimal that round-trips to the
 * double, which is what the JAR's `Double.toString` emits; MATLAB prints
 * `%.17g`, so on a coefficient such as 1/10 the two send different exact
 * rationals to the same service. The divergence is deliberate here: the shortest
 * form is the one whose exact value is the coefficient's intended decimal. The
 * expressions were never text-comparable across codebases anyway -- symbol
 * numbering follows event enumeration order and normal forms depend on the
 * engine -- so compare by substituting rates and comparing numbers.
 *
 * AN EVENT WITH NO POSITIVE RATE CONTRIBUTES NO SYMBOL. Its slot in `symbols`
 * is the empty string and its filtration and term are empty matrices, mirroring
 * MATLAB's empty cell and the JAR's null. Dropping it from the vectors instead
 * would renumber every later symbol and silently change which rate x_e means.
 */

#include <cmath>
#include <cstddef>
#include <cstdio>
#include <memory>
#include <string>
#include <vector>

#include "line/api/mc/ctmc_solve.h"
#include "line/api/sym/sage_rest_engine.h"
#include "line/api/sym/sym_engine.h"
#include "line/api/sym/sym_engines.h"
#include "line/lang/qn/network_struct.h"
#include "line/lang/qn/state.h"
#include "line/lang/qn/state_events.h"
#include "line/solvers/ctmc/solver_ctmc_analyzer.h"
#include "line/solvers/ctmc/solver_ctmc_getters.h"
#include "line/util/error.h"
#include "line/util/matrix.h"

namespace line {
namespace ctmc {

/** Backend selection, mirroring `options.config.symbolic` and its timeout. */
struct CtmcSymbolicOptions {
    /** `auto` to search, a URL, an image name, or `none` to stay local. */
    std::string backend = "auto";
    /** `options.config.symbolic_timeout`, seconds; the reference defaults to 300. */
    int timeout_s = 300;
};

/**
 * The outputs of `@@SolverCTMC/getSymbolicGenerator`.
 *
 * `filt`, `terms`, `rate0` and `symbols` are all indexed by SYNCHRONIZATION, in
 * the order `sync` lists them, and an inactive event keeps its slot.
 */
template <class T>
struct CtmcSymbolicGenerator {
    /** The generator as expression strings, row major; "0" where the entry is zero. */
    std::vector<std::vector<std::string> > Q;
    /** `x1..xE`, empty for an event with no positive rate. */
    std::vector<std::string> symbols;
    /** Event filtration divided by its minimum positive rate; empty if inactive. */
    std::vector<Matrix<T> > filt;
    /** `ctmc_makeinfgen(filt[e])`, the numeric term symbol e scales; empty if inactive. */
    std::vector<Matrix<T> > terms;
    /** Minimum positive rate of each event, i.e. the nominal value of its symbol. */
    std::vector<T> rate0;
    std::vector<NetState<T> > space;  ///< row i of Q is space[i]
    std::vector<Sync<T> > sync;       ///< what `filt` is indexed by
    bool invert_symbol = false;       ///< entries carry `c/x_e` instead of `c*x_e`

    /** The symbols that actually occur, i.e. the non-empty ones. */
    std::vector<std::string> active_symbols() const {
        std::vector<std::string> out;
        for (std::size_t e = 0; e < symbols.size(); ++e)
            if (!symbols[e].empty()) out.push_back(symbols[e]);
        return out;
    }
};

namespace symbolic_detail {

/**
 * A coefficient as text the backend reads exactly: an integer when it is one,
 * the shortest round-tripping decimal otherwise. See the header note on why the
 * shortest form and not MATLAB's `%.17g`.
 */
inline std::string coeff_string(double c) {
    if (c == std::floor(c) && std::fabs(c) < 1e15) {
        char buf[32];
        std::snprintf(buf, sizeof(buf), "%lld", static_cast<long long>(c));
        return std::string(buf);
    }
    return sym::detail::decimal_string(c);
}

/**
 * Entry (i,j) of the symbolic generator, e.g. "2*x1 - 3*x2".
 *
 * The format is the JAR's `getSymbolicEntry` verbatim -- leading unary minus,
 * " + " / " - " between terms, the coefficient omitted when it is 1 -- because
 * that is the form the service parses and the only thing that makes two
 * codebases' output comparable at all.
 */
template <class T>
std::string symbolic_entry(const CtmcSymbolicGenerator<T>& g, std::size_t i, std::size_t j) {
    std::string s;
    for (std::size_t e = 0; e < g.symbols.size(); ++e) {
        if (g.symbols[e].empty()) continue;
        const double c = num_traits<T>::to_double(g.terms[e](i, j));
        if (c == 0.0) continue;
        if (s.empty()) {
            if (c < 0.0) s += "-";
        } else {
            s += c < 0.0 ? " - " : " + ";
        }
        const double a = std::fabs(c);
        if (g.invert_symbol) {
            s += coeff_string(a) + "/" + g.symbols[e];
        } else {
            if (a != 1.0) s += coeff_string(a) + "*";
            s += g.symbols[e];
        }
    }
    return s.empty() ? std::string("0") : s;
}

/** Minimum positive entry of a filtration, or `has = false` when it has none. */
template <class T>
struct MinPositive {
    T value = num_traits<T>::from_int(0);
    bool has = false;
};

template <class T>
MinPositive<T> min_positive(const Matrix<T>& F) {
    MinPositive<T> out;
    const T zero = num_traits<T>::from_int(0);
    for (std::size_t i = 0; i < F.rows(); ++i)
        for (std::size_t j = 0; j < F.cols(); ++j) {
            if (!(F(i, j) > zero)) continue;
            if (!out.has || F(i, j) < out.value) {
                out.value = F(i, j);
                out.has = true;
            }
        }
    return out;
}

/** Sets the timeout on the engine, which only the REST one carries. */
inline void apply_timeout(const std::shared_ptr<sym::SymEngine>& engine, int timeout_s) {
    if (timeout_s <= 0) return;
    sym::SageRestEngine* rest = dynamic_cast<sym::SageRestEngine*>(engine.get());
    if (rest != nullptr) rest->setTimeoutSeconds(timeout_s);
}

/** The reference's `SAGE.require()` message, verbatim in substance. */
inline std::string backend_missing(const std::string& what) {
    return std::string("SolverCTMC: ") + what +
           " needs a computer-algebra backend, and none is available. Start one with\n  docker "
           "run -d -p 8080:8080 " +
           sym::SYM_DOCKER_IMAGE + "\npoint the " + sym::SYM_URL_ENV +
           " environment variable at a running service, or set the backend to its URL";
}

/**
 * Resolves a backend or throws with the guidance above.
 *
 * Refusing beats returning nothing: a caller that got an empty solution back
 * would read it as "this chain has no symbolic stationary distribution", which
 * is a different and false statement.
 */
inline std::shared_ptr<sym::SymEngine> require_engine(const CtmcSymbolicOptions& opt,
                                                      const std::string& what) {
    std::shared_ptr<sym::SymEngine> engine = sym::sym_resolve(opt.backend);
    if (!engine) throw sym::SymEngineError(backend_missing(what));
    apply_timeout(engine, opt.timeout_s);
    return engine;
}

}  // namespace symbolic_detail

/**
 * Port of `@@SolverCTMC/getSymbolicGenerator.m`.
 *
 * No backend is contacted: the generator is linear in the symbols, so it is
 * assembled from the numeric filtration `ctmc_get_generator` already returns.
 *
 * @param sn the refreshed network struct
 * @param opt CTMC options; `keep_filtration` is forced on, the filtration being
 *            the whole content of the answer
 * @param invert_symbol divide each filtration by its symbol instead of
 *                      multiplying, i.e. parameterize by mean times not rates
 * @return the expression matrix, the symbols, and the numeric terms they scale
 */
template <class T>
CtmcSymbolicGenerator<T> ctmc_symbolic_generator(const NetworkStruct<T>& sn,
                                                 const CtmcOptions& opt,
                                                 bool invert_symbol = false) {
    const CtmcGenerator<T> gen = ctmc_get_generator(sn, opt);

    CtmcSymbolicGenerator<T> g;
    g.invert_symbol = invert_symbol;
    g.space = gen.space;
    g.sync = gen.sync;
    const std::size_t n = gen.Q.rows();
    const std::size_t ne = gen.filt.size();
    g.symbols.resize(ne);
    g.filt.resize(ne);
    g.terms.resize(ne);
    g.rate0.assign(ne, num_traits<T>::from_int(0));

    for (std::size_t e = 0; e < ne; ++e) {
        const symbolic_detail::MinPositive<T> m = symbolic_detail::min_positive(gen.filt[e]);
        if (!m.has) continue;  // inactive event: no symbol, empty term, slot kept
        g.rate0[e] = m.value;
        Matrix<T> F(gen.filt[e].rows(), gen.filt[e].cols(), num_traits<T>::from_int(0));
        for (std::size_t i = 0; i < F.rows(); ++i)
            for (std::size_t j = 0; j < F.cols(); ++j)
                if (!(gen.filt[e](i, j) == num_traits<T>::from_int(0)))
                    F(i, j) = T(gen.filt[e](i, j) / m.value);
        g.filt[e] = F;
        // ctmc_makeinfgen is linear, so the symbolic generator is the sum of the
        // per-event terms scaled by their symbols, diagonal included.
        g.terms[e] = mc::ctmc_makeinfgen(F);
        g.symbols[e] = "x" + std::to_string(e + 1);
    }

    g.Q.assign(n, std::vector<std::string>(n, "0"));
    for (std::size_t i = 0; i < n; ++i)
        for (std::size_t j = 0; j < n; ++j) g.Q[i][j] = symbolic_detail::symbolic_entry(g, i, j);
    return g;
}

/**
 * Evaluates the symbolic generator at a symbol assignment, the twin of the
 * JAR's `evalInfGen`.
 *
 * With `x[e] = rate0[e]` this reproduces the numeric generator, which is the
 * cheapest way to check a symbolic build against `ctmc_get_generator`.
 *
 * @param g the symbolic generator
 * @param x one value per event, in the same order; an inactive event's value is
 *          ignored
 */
template <class T>
Matrix<T> ctmc_symbolic_eval_infgen(const CtmcSymbolicGenerator<T>& g, const std::vector<T>& x) {
    if (x.size() != g.symbols.size())
        throw InputError("ctmc_symbolic_eval_infgen: expected " +
                         std::to_string(g.symbols.size()) + " symbol values, got " +
                         std::to_string(x.size()));
    const std::size_t n = g.Q.size();
    Matrix<T> Q(n, n, num_traits<T>::from_int(0));
    for (std::size_t e = 0; e < g.symbols.size(); ++e) {
        if (g.symbols[e].empty()) continue;
        if (g.invert_symbol && x[e] == num_traits<T>::from_int(0))
            throw InputError("ctmc_symbolic_eval_infgen: symbol " + g.symbols[e] +
                             " is inverted in the generator and cannot be zero");
        const T c = g.invert_symbol ? T(num_traits<T>::from_int(1) / x[e]) : x[e];
        for (std::size_t i = 0; i < n; ++i)
            for (std::size_t j = 0; j < n; ++j) Q(i, j) += T(c * g.terms[e](i, j));
    }
    return Q;
}

/** What `@@SolverCTMC/getSymbolicSolution.m` returns, plus the engine that answered. */
template <class T>
struct CtmcSymbolicSolution {
    std::vector<std::string> pi;   ///< stationary probability of each state
    std::vector<std::string> num;  ///< numerator of each entry over `den`
    std::string den = "1";         ///< common denominator of the vector
    int nconncomp = 1;             ///< weakly connected components of the generator
    std::vector<int> conncomp;     ///< component index of each state, one based
    std::vector<NetState<T> > space;  ///< entry i of `pi` is `space[i]`
    std::string engine;               ///< backend that solved it, e.g. `sage`
};

/**
 * Port of `@@SolverCTMC/getSymbolicSolution.m`: pi Q = 0 with sum(pi) = 1 over
 * the field of rational functions in x1..xE.
 *
 * @param sn the refreshed network struct
 * @param opt CTMC options
 * @param symopt backend selection and per-request timeout
 * @return the distribution, also split over one common denominator
 */
template <class T>
CtmcSymbolicSolution<T> ctmc_symbolic_solution(
    const NetworkStruct<T>& sn, const CtmcOptions& opt,
    const CtmcSymbolicOptions& symopt = CtmcSymbolicOptions()) {
    const CtmcSymbolicGenerator<T> g = ctmc_symbolic_generator(sn, opt);
    const std::shared_ptr<sym::SymEngine> engine =
        symbolic_detail::require_engine(symopt, "the symbolic stationary distribution");
    const sym::CtmcSolution sol = engine->solveCTMC(g.Q, g.active_symbols());

    CtmcSymbolicSolution<T> out;
    out.pi = sol.pi;
    out.num = sol.num;
    out.den = sol.den;
    out.nconncomp = sol.nConnComp;
    out.conncomp = sol.connComp;
    out.space = g.space;
    out.engine = engine->name();
    if (out.pi.size() != g.space.size())
        throw sym::SymEngineError("ctmc_symbolic_solution: backend '" + out.engine +
                                  "' returned " + std::to_string(out.pi.size()) +
                                  " entries for a " + std::to_string(g.space.size()) +
                                  " state chain");
    return out;
}

}  // namespace ctmc
}  // namespace line

#endif  // LINE_SOLVERS_CTMC_SOLVER_CTMC_SYMBOLIC_H
