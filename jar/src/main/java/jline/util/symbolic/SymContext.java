/**
 * @file The variable set a symbolic computation is carried out over.
 *
 * @since LINE 3.0
 */
package jline.util.symbolic;

import java.util.ArrayList;
import java.util.Arrays;
import java.util.Collections;
import java.util.HashMap;
import java.util.List;
import java.util.Map;

import cc.redberry.rings.Rational;
import cc.redberry.rings.Rationals;
import cc.redberry.rings.Rings;
import cc.redberry.rings.bigint.BigInteger;
import cc.redberry.rings.poly.MultivariateRing;
import cc.redberry.rings.poly.multivar.MultivariatePolynomial;

/**
 * The ordered variable set a family of {@link SymExpr} values shares.
 *
 * <p>Rings fixes the number of indeterminates when the ring is built, so every
 * expression that is to be added or multiplied has to come from ONE context.
 * That is what this class owns: the names, their positions, and the two Rings
 * objects derived from them.
 *
 * <p>THE FIELD IS THE FIELD OF RATIONAL FUNCTIONS over Q, not a polynomial ring,
 * because the gld recursion divides. Every divisor there is a bare {@code mu}
 * symbol, so a denominator is always a monomial and the normal form Rings keeps
 * is already the smallest one; nothing here has to call gcd for correctness.
 *
 * <p>Names are compared literally and are the strings the caller printed them
 * with, so a round trip through {@link SymExpr#toString()} names the same
 * symbols the caller supplied.
 *
 * Copyright (c) 2012-2026, Imperial College London
 * All rights reserved.
 */
public final class SymContext {

    private final List<String> names;
    private final Map<String, Integer> index;
    private final MultivariateRing<MultivariatePolynomial<BigInteger>> ring;
    private final Rationals<MultivariatePolynomial<BigInteger>> field;

    private SymContext(List<String> names) {
        if (names.isEmpty()) {
            throw new RuntimeException("SymContext: a symbolic context needs at least one variable.");
        }
        this.names = Collections.unmodifiableList(new ArrayList<String>(names));
        this.index = new HashMap<String, Integer>();
        for (int i = 0; i < this.names.size(); i++) {
            if (this.index.put(this.names.get(i), Integer.valueOf(i)) != null) {
                throw new RuntimeException("SymContext: duplicate variable '" + this.names.get(i) + "'.");
            }
        }
        this.ring = Rings.MultivariateRing(this.names.size(), Rings.Z);
        this.field = Rings.Frac(this.ring);
    }

    /**
     * A context over the given variable names, in the order given.
     *
     * @param names the variable names, each distinct
     * @return the context
     */
    public static SymContext of(String... names) {
        return new SymContext(Arrays.asList(names));
    }

    /**
     * A context over the given variable names, in the order given.
     *
     * @param names the variable names, each distinct
     * @return the context
     */
    public static SymContext of(List<String> names) {
        return new SymContext(names);
    }

    /**
     * The variable of the given name as an expression.
     *
     * @param name a name this context was built with
     * @return the variable
     */
    public SymExpr var(String name) {
        Integer i = this.index.get(name);
        if (i == null) {
            throw new RuntimeException("SymContext: '" + name + "' is not a variable of this context; it holds "
                    + this.names + ".");
        }
        return new SymExpr(this, new Rational<MultivariatePolynomial<BigInteger>>(
                this.ring, this.ring.variable(i.intValue())));
    }

    /**
     * An exact integer constant.
     *
     * @param v the value
     * @return the constant
     */
    public SymExpr constant(long v) {
        return new SymExpr(this, new Rational<MultivariatePolynomial<BigInteger>>(
                this.ring, this.ring.valueOf(v)));
    }

    /**
     * An exact rational constant.
     *
     * @param num numerator
     * @param den denominator, non-zero
     * @return the constant
     */
    public SymExpr rational(long num, long den) {
        if (den == 0) {
            throw new RuntimeException("SymContext: a rational constant cannot have a zero denominator.");
        }
        return new SymExpr(this, new Rational<MultivariatePolynomial<BigInteger>>(
                this.ring, this.ring.valueOf(num), this.ring.valueOf(den)));
    }

    /**
     * A double carried in EXACTLY, as the dyadic rational it already is.
     *
     * <p>No rounding happens here and none is wanted: a demand written 0.1 in
     * the caller's source is the double nearest 0.1, and turning it into the
     * rational 1/10 would answer for a model the caller did not state. The same
     * rule the C++ {@code num_traits<Rational>::from_double} follows.
     *
     * @param v the value, finite
     * @return the constant
     */
    public SymExpr constant(double v) {
        if (Double.isNaN(v) || Double.isInfinite(v)) {
            throw new RuntimeException("SymContext: cannot carry " + v + " into exact arithmetic.");
        }
        java.math.BigDecimal d = new java.math.BigDecimal(v);
        java.math.BigInteger unscaled = d.unscaledValue();
        int scale = d.scale();
        MultivariatePolynomial<BigInteger> num;
        MultivariatePolynomial<BigInteger> den;
        if (scale >= 0) {
            num = constantPoly(unscaled);
            den = constantPoly(java.math.BigInteger.TEN.pow(scale));
        } else {
            num = constantPoly(unscaled.multiply(java.math.BigInteger.TEN.pow(-scale)));
            den = this.ring.getOne();
        }
        return new SymExpr(this, new Rational<MultivariatePolynomial<BigInteger>>(this.ring, num, den));
    }

    /**
     * A constant polynomial holding an arbitrary integer.
     *
     * <p>{@code ring.valueOf} only takes a long, and the unscaled value of a
     * double routinely exceeds one: {@code new BigDecimal(0.1)} carries 55
     * digits, because the double nearest 0.1 is not 1/10.
     */
    private MultivariatePolynomial<BigInteger> constantPoly(java.math.BigInteger v) {
        return this.ring.getOne().createConstant(new BigInteger(v.toString()));
    }

    /** @return zero in this context */
    public SymExpr zero() {
        return new SymExpr(this, this.field.getZero());
    }

    /** @return one in this context */
    public SymExpr one() {
        return new SymExpr(this, this.field.getOne());
    }

    /** @return the variable names, in ring order */
    public List<String> names() {
        return this.names;
    }

    /** @return the position of a variable in ring order, or -1 */
    public int positionOf(String name) {
        Integer i = this.index.get(name);
        return i == null ? -1 : i.intValue();
    }

    /** @return the underlying polynomial ring */
    public MultivariateRing<MultivariatePolynomial<BigInteger>> ring() {
        return this.ring;
    }

    /** @return the underlying field of rational functions */
    public Rationals<MultivariatePolynomial<BigInteger>> field() {
        return this.field;
    }
}
