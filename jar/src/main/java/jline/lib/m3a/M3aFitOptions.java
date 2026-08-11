/**
 * @file M3aFitOptions data class for M3A fitting options
 *
 * @since LINE 3.0
 */
package jline.lib.m3a;

import java.util.Objects;

/**
 * Options for M3A fitting algorithms.
 */
public final class M3aFitOptions {
    private final int method;
    private final int numStates;
    private final Double timescale;
    private final Double timescaleAsy;

    public M3aFitOptions(int method, int numStates, Double timescale, Double timescaleAsy) {
        this.method = method;
        this.numStates = numStates;
        this.timescale = timescale;
        this.timescaleAsy = timescaleAsy;
    }

    public M3aFitOptions(int numStates) {
        this(1, numStates, null, null);
    }

    public int getMethod() { return method; }
    public int getNumStates() { return numStates; }
    public Double getTimescale() { return timescale; }
    public Double getTimescaleAsy() { return timescaleAsy; }

    @Override
    public boolean equals(Object o) {
        if (this == o) return true;
        if (!(o instanceof M3aFitOptions)) return false;
        M3aFitOptions that = (M3aFitOptions) o;
        return method == that.method && numStates == that.numStates
                && Objects.equals(timescale, that.timescale)
                && Objects.equals(timescaleAsy, that.timescaleAsy);
    }

    @Override
    public int hashCode() {
        return Objects.hash(method, numStates, timescale, timescaleAsy);
    }

    @Override
    public String toString() {
        return "M3aFitOptions(method=" + method + ", numStates=" + numStates
                + ", timescale=" + timescale + ", timescaleAsy=" + timescaleAsy + ")";
    }
}
