/**
 * @file An exact symbolic scalar: a rational function over Q.
 *
 * @since LINE 3.0
 */
package jline.util.symbolic;

import java.util.List;
import java.util.Map;

import cc.redberry.rings.Rational;
import cc.redberry.rings.bigint.BigInteger;
import cc.redberry.rings.poly.multivar.MultivariatePolynomial;

/**
 * A rational function over Q in the variables of a {@link SymContext}.
 *
 * <p>This is the scalar the SYMBOLIC arm of the pfqn_gld family computes over.
 * It carries {@code +}, {@code -}, {@code *}, {@code /} and integer powers, and
 * NOTHING ELSE, which is exactly the operation set that arm uses: the recursion
 * is {@code g = g + L*g/mu} throughout, so its result is a polynomial over a
 * monomial denominator and every intermediate is exact.
 *
 * <p>WHAT IT DELIBERATELY DOES NOT CARRY IS A COMPARISON. {@code L > 0} and
 * {@code min(row) == 1} have no truth value on a symbol, which is the reason the
 * reference's symbolic arm skips its own load-independence scan and its
 * zero-demand guard rather than deciding them; the arms here do the same and
 * select on {@code N}, which is always concrete. An {@code equals} IS provided,
 * and is structural equality of the normal form Rings keeps, not a decision
 * about the values the variables might take.
 *
 * <p>NO TRANSCENDENTAL FUNCTION IS AVAILABLE, so there is no {@code log}. The
 * reference returns {@code lG = log(G)} as a symbolic expression and MATLAB,
 * sympy and the C++ SymEngine backend can all represent that; a field of
 * rational functions cannot. {@link jline.io.Ret.pfqnNcSym} therefore carries
 * the FORMAL expression as text and leaves the value to be taken after
 * substitution. See _kb/07-cross-language-parity.md.
 *
 * Copyright (c) 2012-2026, Imperial College London
 * All rights reserved.
 */
public final class SymExpr {

    private final SymContext ctx;
    private final Rational<MultivariatePolynomial<BigInteger>> value;

    SymExpr(SymContext ctx, Rational<MultivariatePolynomial<BigInteger>> value) {
        this.ctx = ctx;
        this.value = value;
    }

    /** @return the context this expression belongs to */
    public SymContext context() {
        return this.ctx;
    }

    /** @return the underlying Rings value */
    public Rational<MultivariatePolynomial<BigInteger>> value() {
        return this.value;
    }

    private void sameContext(SymExpr other) {
        if (other.ctx != this.ctx) {
            throw new RuntimeException("SymExpr: the two operands come from different symbolic contexts; "
                    + "Rings fixes the indeterminates when the ring is built, so they cannot be combined.");
        }
    }

    /**
     * @param other the addend
     * @return this + other
     */
    public SymExpr add(SymExpr other) {
        sameContext(other);
        return new SymExpr(this.ctx, this.ctx.field().add(this.value, other.value));
    }

    /**
     * @param other the subtrahend
     * @return this - other
     */
    public SymExpr subtract(SymExpr other) {
        sameContext(other);
        return new SymExpr(this.ctx, this.ctx.field().subtract(this.value, other.value));
    }

    /**
     * @param other the multiplicand
     * @return this * other
     */
    public SymExpr multiply(SymExpr other) {
        sameContext(other);
        return new SymExpr(this.ctx, this.ctx.field().multiply(this.value, other.value));
    }

    /**
     * @param other the divisor, non-zero
     * @return this / other
     */
    public SymExpr divide(SymExpr other) {
        sameContext(other);
        if (other.isZero()) {
            throw new RuntimeException("SymExpr: division by an expression that is identically zero.");
        }
        return new SymExpr(this.ctx, this.ctx.field().divideExact(this.value, other.value));
    }

    /**
     * @param k a non-negative exponent
     * @return this raised to k
     */
    public SymExpr pow(int k) {
        if (k < 0) {
            throw new RuntimeException("SymExpr: a negative exponent is a division; use divide instead.");
        }
        SymExpr acc = this.ctx.one();
        for (int i = 0; i < k; i++) {
            acc = acc.multiply(this);
        }
        return acc;
    }

    /** @return whether this is identically zero */
    public boolean isZero() {
        return this.value.isZero();
    }

    /** @return whether this is identically one */
    public boolean isOne() {
        return this.value.isOne();
    }

    /**
     * The value at a numeric assignment of every variable.
     *
     * <p>This is how a symbolic result is checked against the numeric routine:
     * substitute, then compare. Every variable of the context must be given a
     * value, since a partially substituted rational function is still symbolic
     * and has no double to return.
     *
     * @param at variable name to value
     * @return the value as a double
     */
    public double evaluate(Map<String, Double> at) {
        List<String> names = this.ctx.names();
        MultivariatePolynomial<BigInteger> num = this.value.numerator();
        MultivariatePolynomial<BigInteger> den = this.value.denominator();
        // Substitution is exact and done over a common power of ten, so that a
        // decimal assignment does not round before the division does.
        java.math.BigInteger scale = java.math.BigInteger.ONE;
        java.math.BigInteger[] vals = new java.math.BigInteger[names.size()];
        int[] scales = new int[names.size()];
        for (int i = 0; i < names.size(); i++) {
            Double v = at.get(names.get(i));
            if (v == null) {
                throw new RuntimeException("SymExpr: no value given for variable '" + names.get(i)
                        + "'; a partial substitution is still symbolic.");
            }
            java.math.BigDecimal d = new java.math.BigDecimal(v.doubleValue());
            vals[i] = d.unscaledValue();
            scales[i] = d.scale();
        }
        // Evaluate numerator and denominator as exact rationals, then divide once.
        java.math.BigInteger[] nd = evalExact(num, vals, scales);
        java.math.BigInteger[] dd = evalExact(den, vals, scales);
        java.math.BigDecimal numV = new java.math.BigDecimal(nd[0]).divide(
                new java.math.BigDecimal(nd[1]), java.math.MathContext.DECIMAL128);
        java.math.BigDecimal denV = new java.math.BigDecimal(dd[0]).divide(
                new java.math.BigDecimal(dd[1]), java.math.MathContext.DECIMAL128);
        if (denV.signum() == 0) {
            throw new RuntimeException("SymExpr: the denominator vanishes at this assignment.");
        }
        return numV.divide(denV, java.math.MathContext.DECIMAL128).doubleValue();
        }

    /** Exact value of a polynomial at the given scaled assignment, as num/den. */
    private java.math.BigInteger[] evalExact(MultivariatePolynomial<BigInteger> p,
                                             java.math.BigInteger[] vals, int[] scales) {
        java.math.BigInteger accNum = java.math.BigInteger.ZERO;
        java.math.BigInteger accDen = java.math.BigInteger.ONE;
        for (cc.redberry.rings.poly.multivar.Monomial<BigInteger> m : p) {
            java.math.BigInteger termNum = new java.math.BigInteger(m.coefficient.toString());
            java.math.BigInteger termDen = java.math.BigInteger.ONE;
            for (int i = 0; i < m.exponents.length; i++) {
                int e = m.exponents[i];
                if (e == 0) {
                    continue;
                }
                termNum = termNum.multiply(vals[i].pow(e));
                if (scales[i] > 0) {
                    termDen = termDen.multiply(java.math.BigInteger.TEN.pow(scales[i] * e));
                } else if (scales[i] < 0) {
                    termNum = termNum.multiply(java.math.BigInteger.TEN.pow(-scales[i] * e));
                }
            }
            // accNum/accDen + termNum/termDen
            accNum = accNum.multiply(termDen).add(termNum.multiply(accDen));
            accDen = accDen.multiply(termDen);
        }
        return new java.math.BigInteger[]{accNum, accDen};
    }

    /**
     * The expression printed with the context's own variable names.
     *
     * @return the printed form
     */
    @Override
    public String toString() {
        String[] vars = this.ctx.names().toArray(new String[0]);
        String num = this.value.numerator().toString(vars);
        if (this.value.denominator().isOne()) {
            return num;
        }
        return "(" + num + ")/(" + this.value.denominator().toString(vars) + ")";
    }

    @Override
    public boolean equals(Object o) {
        if (this == o) {
            return true;
        }
        if (!(o instanceof SymExpr)) {
            return false;
        }
        SymExpr other = (SymExpr) o;
        return this.ctx == other.ctx && this.value.equals(other.value);
    }

    @Override
    public int hashCode() {
        return this.value.hashCode();
    }
}
