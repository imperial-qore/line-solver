/**
 * @file An exact symbolic scalar, as a `num_traits` numeric type.
 *
 * Copyright (c) 2012-2026, Imperial College London
 * All rights reserved.
 */
#ifndef LINE_NUM_SYM_NUMBER_H
#define LINE_NUM_SYM_NUMBER_H

#ifdef LINE_MP_USE_SYMENGINE

#include <cmath>
#include <limits>
#include <string>

#include <symengine/basic.h>
#include <symengine/eval_double.h>
#include <symengine/expression.h>
#include <symengine/symbol.h>
#include <symengine/visitor.h>

#include "line/num/number.h"

namespace line {

/**
 * A symbolic scalar: a SymEngine expression carrying `+ - * /` and powers.
 *
 * WHY THIS IS A `num_traits` TYPE AND NOT A NEW CODE PATH. Every routine that
 * is to accept symbolic input here is already `template <class T>` over
 * `Matrix<T>`, and `num_traits<T>` is the one place a numeric type declares
 * itself: `double`, `Rational`, `Real<D>` and `Complex` are already four such
 * declarations. Adding a fifth makes `pfqn_gld` and `pfqn_gldsingle` symbolic
 * WITHOUT EDITING EITHER OF THEM, which is also what keeps the symbolic and
 * numeric arms from drifting: there is only one arm.
 *
 * `is_exact` is true and `has_transcendental` is FALSE, matching `Rational`.
 * The second is what steers the gld family off its log-domain path, which is
 * correct here for the reason it is correct there: the log domain exists to stop
 * a double underflowing, and an exact expression cannot underflow. SymEngine
 * itself does carry `log`, so unlike the JAR's Rings backend this type CAN
 * represent `lG = log(G)`; `log_as_double` is nonetheless a double-valued
 * accessor by the trait's contract, so it evaluates and is defined only for a
 * fully substituted expression.
 *
 * COMPARISONS ARE STRUCTURAL, never decisions about the values a symbol may
 * take. `operator==` compares normal forms; there is deliberately no
 * `operator<`, because ordering a symbol has no meaning and the algorithms that
 * would want one (the lld threshold scan) are the very ones the reference
 * refuses on symbolic input.
 */
class Sym {
  public:
    Sym() : e_(SymEngine::integer(0)) {}
    explicit Sym(const SymEngine::Expression& e) : e_(e) {}
    Sym(int v) : e_(SymEngine::integer(v)) {}                       // NOLINT: implicit by design
    Sym(long v) : e_(SymEngine::integer(v)) {}                      // NOLINT
    Sym(const SymEngine::RCP<const SymEngine::Basic>& b) : e_(b) {}  // NOLINT

    /** A named symbol. */
    static Sym var(const std::string& name) {
        return Sym(SymEngine::Expression(SymEngine::symbol(name)));
    }

    /** A double carried in as the exact dyadic rational it already is. */
    static Sym from_double_exact(double v) {
        return Sym(SymEngine::Expression(SymEngine::real_double(v)));
    }

    const SymEngine::Expression& expr() const { return e_; }

    Sym operator+(const Sym& o) const { return Sym(e_ + o.e_); }
    Sym operator-(const Sym& o) const { return Sym(e_ - o.e_); }
    Sym operator*(const Sym& o) const { return Sym(e_ * o.e_); }
    Sym operator/(const Sym& o) const { return Sym(e_ / o.e_); }
    Sym operator-() const { return Sym(-e_); }

    Sym& operator+=(const Sym& o) { e_ = e_ + o.e_; return *this; }
    Sym& operator-=(const Sym& o) { e_ = e_ - o.e_; return *this; }
    Sym& operator*=(const Sym& o) { e_ = e_ * o.e_; return *this; }
    Sym& operator/=(const Sym& o) { e_ = e_ / o.e_; return *this; }

    /** Structural equality of the normal form, not a claim about values. */
    bool operator==(const Sym& o) const { return e_ == o.e_; }
    bool operator!=(const Sym& o) const { return !(*this == o); }

    std::string str() const { return SymEngine::str(e_); }

  private:
    SymEngine::Expression e_;
};

inline std::string to_string(const Sym& v) { return v.str(); }

template <>
struct num_traits<Sym> {
    using type = Sym;
    static constexpr bool is_exact = true;
    /**
     * FALSE, as for `Rational`, and for the same reason: it is what steers a
     * routine off its log-domain arm. SymEngine can represent `log`, but the
     * trait means "may this type be handed to std::log and give a value of its
     * own kind", and an expression cannot answer that without a substitution.
     */
    static constexpr bool has_transcendental = false;
    static const char* name() { return "symbolic"; }

    static Sym from_int(long v) { return Sym(v); }
    static Sym from_rational(long num, long den) {
        return Sym(SymEngine::Expression(SymEngine::integer(num)) /
                   SymEngine::Expression(SymEngine::integer(den)));
    }
    static Sym from_double(double v) { return Sym::from_double_exact(v); }

    /**
     * Defined only for a FULLY SUBSTITUTED expression: a symbol left in it has
     * no double, and SymEngine raises rather than guessing one. That raise is
     * the right behaviour for a caller who asked for a number.
     */
    static double to_double(const Sym& v) {
        return SymEngine::eval_double(*v.expr().get_basic());
    }

    /**
     * NaN ON A FREE SYMBOL, rather than a raise, and this is the one place the
     * two differ.
     *
     * `NcResult<T>` carries `lG` as a plain double and every return path of the
     * gld family fills it by calling this, unconditionally, on a value it has
     * just computed. Under `T = Sym` that value is usually symbolic, so a raise
     * here would make the SUCCESSFUL symbolic result unreachable: `G` is the
     * answer and `lG` is a convenience beside it. NaN says "there is no double
     * for this" in the field's own vocabulary and leaves `G` intact.
     *
     * A fully substituted expression still gets its real logarithm, so a caller
     * who substituted before asking loses nothing.
     */
    static double log_as_double(const Sym& v) {
        if (!SymEngine::free_symbols(*v.expr().get_basic()).empty()) {
            return std::numeric_limits<double>::quiet_NaN();
        }
        return std::log(to_double(v));
    }
    static std::string to_string(const Sym& v) { return v.str(); }
};

}  // namespace line

#endif  // LINE_MP_USE_SYMENGINE
#endif  // LINE_NUM_SYM_NUMBER_H
