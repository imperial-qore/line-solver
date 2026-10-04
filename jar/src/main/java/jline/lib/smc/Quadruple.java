/**
 * @file Quadruple data class for MG1 shift results
 *
 * @since LINE 3.1.0
 */
package jline.lib.smc;

import java.util.Objects;

/**
 * Helper class for returning four values.
 */
public final class Quadruple<A, B, C, D> {
    private final A first;
    private final B second;
    private final C third;
    private final D fourth;

    public Quadruple(A first, B second, C third, D fourth) {
        this.first = first;
        this.second = second;
        this.third = third;
        this.fourth = fourth;
    }

    public A getFirst() { return first; }
    public B getSecond() { return second; }
    public C getThird() { return third; }
    public D getFourth() { return fourth; }

    @Override
    public boolean equals(Object o) {
        if (this == o) return true;
        if (!(o instanceof Quadruple)) return false;
        Quadruple<?, ?, ?, ?> that = (Quadruple<?, ?, ?, ?>) o;
        return Objects.equals(first, that.first)
                && Objects.equals(second, that.second)
                && Objects.equals(third, that.third)
                && Objects.equals(fourth, that.fourth);
    }

    @Override
    public int hashCode() {
        return Objects.hash(first, second, third, fourth);
    }

    @Override
    public String toString() {
        return "Quadruple(" + first + ", " + second + ", " + third + ", " + fourth + ")";
    }
}
