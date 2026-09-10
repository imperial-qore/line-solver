/*
 * Copyright (c) 2012-2026, QORE Lab, Imperial College London
 * All rights reserved.
 */

package jline.util;

import java.io.Serializable;

/**
 * A generic triple container that holds three objects of potentially different types.
 *
 * <p>This utility class provides a simple way to group three related values together.
 * It implements {@link Serializable} for persistence. It is the three-element
 * counterpart of {@link Pair} and exposes the same property-style accessors
 * ({@link #getFirst()}, {@link #getSecond()}, {@link #getThird()}) so that code
 * previously written against {@code Triple} works unchanged.</p>
 *
 * <p>Common use cases:
 * <ul>
 *   <li>Returning three values from methods</li>
 *   <li>Temporary grouping of related objects</li>
 * </ul>
 * </p>
 *
 * @param <A> the type of the first element
 * @param <B> the type of the second element
 * @param <C> the type of the third element
 * @since 1.0
 */
public class Triple<A, B, C> implements Serializable {
    private A first;
    private B second;
    private C third;

    /**
     * Constructs a new triple with the specified elements.
     *
     * @param first  the first element of the triple
     * @param second the second element of the triple
     * @param third  the third element of the triple
     */
    public Triple(A first, B second, C third) {
        this.first = first;
        this.second = second;
        this.third = third;
    }

    /**
     * Copy constructor that creates a new triple with the same elements as another triple.
     *
     * @param other the triple to copy from
     */
    public Triple(Triple<A, B, C> other) {
        this.first = other.first;
        this.second = other.second;
        this.third = other.third;
    }

    public A getFirst() {
        return this.first;
    }

    public void setFirst(A first) {
        this.first = first;
    }

    public B getSecond() {
        return this.second;
    }

    public void setSecond(B second) {
        this.second = second;
    }

    public C getThird() {
        return this.third;
    }

    public void setThird(C third) {
        this.third = third;
    }

    @Override
    public boolean equals(Object o) {
        if (this == o) {
            return true;
        }
        if (!(o instanceof Triple)) {
            return false;
        }
        Triple<?, ?, ?> other = (Triple<?, ?, ?>) o;
        return (this.first == null ? other.first == null : this.first.equals(other.first))
                && (this.second == null ? other.second == null : this.second.equals(other.second))
                && (this.third == null ? other.third == null : this.third.equals(other.third));
    }

    @Override
    public int hashCode() {
        int result = (this.first == null ? 0 : this.first.hashCode());
        result = 31 * result + (this.second == null ? 0 : this.second.hashCode());
        result = 31 * result + (this.third == null ? 0 : this.third.hashCode());
        return result;
    }

    @Override
    public String toString() {
        return "(" + this.first + ", " + this.second + ", " + this.third + ")";
    }
}
