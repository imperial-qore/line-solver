/*
 * Copyright (c) 2012-2026, QORE Lab, Imperial College London
 * All rights reserved.
 */
#ifndef LINE_API_SYM_SYM_ENGINE_H
#define LINE_API_SYM_SYM_ENGINE_H

/**
 * @file
 * @ingroup api_sym
 * Computer algebra operations LINE needs, as seen by this port.
 *
 * Port of jline.api.sym.SymEngine. C++ has no computer algebra system any more
 * than Java does, so every operation here is delegated to an external engine;
 * SageRestEngine is the SageMath implementation and is resolved by
 * sym_engines.h. Expressions cross the interface as plain ASCII infix strings,
 * e.g. "2*x1 - 3*x2".
 *
 * EXPRESSION STRINGS ARE NOT COMPARABLE ACROSS CODEBASES, and not even across
 * engine versions: the symbol numbering x1..xE follows event enumeration order
 * and printed normal forms depend on the engine. Compare by substituting values
 * with eval() and comparing numbers; a string diff against the JAR's output is
 * a false failure waiting to happen.
 *
 * ABSENT VALUES. Java returns null for a field the request did not ask for.
 * Here an absent scalar is the empty string and an absent list or matrix is
 * empty, with an explicit `has*` flag where emptiness would otherwise be
 * ambiguous.
 */

#include <map>
#include <string>
#include <vector>

#include "line/util/error.h"

namespace line {
namespace sym {

/** The symbolic backend is unreachable, or rejected the request. */
class SymEngineError : public Error {
public:
    explicit SymEngineError(const std::string& what) : Error(what) {}
};

/** Symbolic stationary distribution of a CTMC. */
/** One weighted sum of the stationary distribution, or a ratio of two. */
struct CtmcMeasure {
    std::string name;    ///< Caller supplied name
    std::string expr;    ///< The value as a rational function, or "undefined"
    std::string num;     ///< Numerator of `expr`, empty when undefined
    std::string den;     ///< Denominator of `expr`, empty when undefined
    /**
     * Why the value is undefined, empty otherwise.
     *
     * A ratio whose denominator measure is identically zero over the whole rate
     * space has no value. The numeric arms sweep such a 0/0 to zero with an
     * isnan pass, which is a defensible cleanup of floating point dust and an
     * indefensible answer for an exact one.
     */
    std::string reason;

    /** True if this measure has a value. */
    bool isDefined() const { return reason.empty(); }
};

/** A ratio of two measures, by their index into the weight block. */
struct RatioSpec {
    std::string name;      ///< Name of the resulting measure
    std::size_t num = 0;   ///< Index of the numerator measure
    std::size_t den = 0;   ///< Index of the denominator measure
};

/** Stationary distribution together with the measures taken from it. */
struct CtmcMeasures {
    std::vector<std::string> pi;   ///< Stationary probability of each state
    std::vector<std::string> num;  ///< Numerator of each entry over `den`
    std::string den = "1";         ///< Common denominator of the whole vector
    int nConnComp = 1;             ///< Weakly connected components
    std::vector<int> connComp;     ///< Component index of each state, one based
    std::vector<CtmcMeasure> measures;  ///< One per weight vector, in order
    std::vector<CtmcMeasure> ratios;    ///< One per requested ratio, in order
};

/** First passage time transform and moments. */
struct Passage {
    std::string lst;      ///< The transform L(s), empty if not requested
    std::string lstNum;   ///< Numerator of `lst`
    std::string lstDen;   ///< Denominator of `lst`
    std::vector<std::string> lstAll;               ///< Per start state transform
    std::vector<std::string> moments;              ///< Moments 1..nmax for the initial law
    std::vector<std::vector<std::string>> momAll;  ///< Per start state moments
    std::vector<int> unreachable;                  ///< States that cannot reach the target
};

struct CtmcSolution {
    std::vector<std::string> pi;   ///< Stationary probability of each state, as an expression
    std::vector<std::string> num;  ///< Numerator of each entry over the common denominator
    std::string den = "1";         ///< Common denominator of the whole vector
    int nConnComp = 1;             ///< Weakly connected components of the generator
    std::vector<int> connComp;     ///< Component index of each state, one based
};

/**
 * Exact parametric sensitivity, following Trivedi and Bobbio (2017), Sec. 9.7.
 *
 * As in `@SolverCTMC/getSensitivity.m`, dr/dtheta is taken to be zero: a reward
 * whose rates themselves depend on theta needs the second term of Eq. (9.83)
 * and is not covered here.
 */
struct SymSensitivity {
    std::vector<std::string> pi;   ///< Stationary distribution
    std::vector<std::string> dpi;  ///< Derivative of the distribution with respect to theta
    std::string Er;                ///< Mean reward, empty if no reward was given
    std::string S;                 ///< Unscaled sensitivity d(E[r])/dtheta, Eq. (9.79)
    std::string SS;                ///< Scaled sensitivity (theta/E[r]) d(E[r])/dtheta, Eq. (9.80)
    bool hasReward = false;        ///< Whether Er, S and SS were computed
};

/** Symbolic analysis of a fluid vector field. */
struct FluidODEs {
    std::vector<std::vector<std::string>> jacobian;         ///< d f_i / d x_j
    std::vector<std::string> latex;                         ///< LaTeX form of each right hand side
    std::vector<std::map<std::string, std::string>> equilibria;  ///< variable -> expression
    bool hasJacobian = false;
    bool hasLatex = false;
    bool hasEquilibria = false;
};

/** A computer algebra backend. */
class SymEngine {
public:
    virtual ~SymEngine() {}

    /** Name of the backing engine, e.g. "sage". */
    virtual std::string name() const = 0;

    /** True if the engine answers a health probe. */
    virtual bool isAvailable() const = 0;

    /**
     * Symbolic stationary distribution of a CTMC, pi Q = 0 with sum(pi) = 1.
     *
     * @param Q       generator entries as expression strings, row major and square
     * @param symbols the symbols appearing in Q, e.g. x1..xE
     * @return the solution, with pi also split over a common denominator
     */
    virtual CtmcSolution solveCTMC(const std::vector<std::vector<std::string>>& Q,
                                   const std::vector<std::string>& symbols) = 0;

    /**
     * Exact parametric sensitivity of a steady-state reward.
     *
     * @param Q       generator entries as expression strings, row major
     * @param symbols the symbols appearing in Q
     * @param theta   the symbol to differentiate with respect to
     * @param reward  reward rate per state, empty for the distribution alone
     * @return the sensitivity of the distribution and, if a reward is given, of its mean
     */
    virtual SymSensitivity ctmcSensitivity(const std::vector<std::vector<std::string>>& Q,
                                           const std::vector<std::string>& symbols,
                                           const std::string& theta,
                                           const std::vector<std::string>& reward) = 0;

    /**
     * Rewrites expressions into a normal form.
     *
     * @param exprs the expressions
     * @param form  one of simplify, factor, together, cancel, expand, latex
     * @return the rewritten expressions, in the input order
     */
    /**
     * Stationary distribution and a set of weighted sums of it.
     *
     * Every mean SolverCTMC reports is a NUMERIC linear functional of pi, so
     * the caller builds one weight vector per measure with no algebra at all
     * and this one round trip returns each `w . pi` as a rational function: a
     * throughput is `pi . depRates`, a queue length is `pi . stateSpaceAggr`, a
     * marginal probability is `pi . indicator` and a mean reward is `pi . r`.
     *
     * @param Q       generator entries as expression strings, row major
     * @param symbols the symbols appearing in Q
     * @param weights one row per measure, each of length Q.size(); entries are
     *                parsed in the same field as Q, so an expression is
     *                admissible and not only a number
     * @param names   a name per weight row, or empty for w0, w1, ...
     * @param ratios  ratios of two measures, or empty for none
     */
    virtual CtmcMeasures ctmcMeasures(const std::vector<std::vector<std::string>>& Q,
                                      const std::vector<std::string>& symbols,
                                      const std::vector<std::vector<std::string>>& weights,
                                      const std::vector<std::string>& names,
                                      const std::vector<RatioSpec>& ratios) = 0;

    /**
     * First passage time into a target set: the transform and the moments.
     *
     * With A the complement of the target, S = Q(A,A) the sub-generator,
     * s0 = -S*1 the exit vector and alpha the initial law restricted to A,
     * following Harrison and Knottenbelt (2002),
     *
     *     L(s) = alpha (sI - S)^-1 s0 + atom        Eqs. 1-2
     *     (-S) M(n) = n M(n-1),  M(0) = 1           Eq. 3
     *
     * THE MOMENTS COME FROM THE RECURSION, not from differentiating L: it
     * carries no transform symbol, stays in the same fraction field and is
     * `nmax` right solves against one matrix. Differentiating would leave the
     * fraction field, and the evaluation at s = 0 can hit 0/0 where a factor of
     * s failed to cancel.
     *
     * @param S       the sub-generator on the non-target states, row major
     * @param s0      the exit vector -S*1
     * @param alpha   the initial law restricted to the non-target states
     * @param atom    the initial mass already inside the target set
     * @param symbols the rate symbols appearing in S
     * @param svar    the transform symbol; adjoined to the field only when a
     *                transform is asked for, and required to differ from every
     *                rate symbol
     * @param want    any of lst, lstall, moments, momall
     * @param nmax    highest moment order, at least 1
     */
    virtual Passage ctmcPassage(const std::vector<std::vector<std::string>>& S,
                                const std::vector<std::string>& s0,
                                const std::vector<std::string>& alpha,
                                const std::string& atom,
                                const std::vector<std::string>& symbols,
                                const std::string& svar,
                                const std::vector<std::string>& want, int nmax) = 0;

    virtual std::vector<std::string> simplify(const std::vector<std::string>& exprs,
                                              const std::string& form) = 0;

    /**
     * Differentiates expressions.
     *
     * @param exprs    the expressions
     * @param variable the differentiation variable
     * @param order    the order of the derivative, at least 1
     * @return the derivatives, in the input order
     */
    virtual std::vector<std::string> diff(const std::vector<std::string>& exprs,
                                          const std::string& variable, int order) = 0;

    /**
     * Substitutes values for symbols and evaluates.
     *
     * @param exprs      the expressions
     * @param assignment value of each symbol
     * @return the numeric values, NaN where a free symbol remains
     */
    virtual std::vector<double> eval(const std::vector<std::string>& exprs,
                                     const std::map<std::string, double>& assignment) = 0;

    /**
     * Jacobian, LaTeX form and equilibria of a fluid vector field.
     *
     * @param rhs  the right hand side of dx/dt, one expression per state variable
     * @param vars the state variable names
     * @param want any of jacobian, latex, equilibria
     * @return the requested items
     */
    virtual FluidODEs fluidODEs(const std::vector<std::string>& rhs,
                                const std::vector<std::string>& vars,
                                const std::vector<std::string>& want) = 0;
};

}  // namespace sym
}  // namespace line

#endif  // LINE_API_SYM_SYM_ENGINE_H
