/**
 * @file Options for M3A compression
 *
 * @since LINE 3.0
 */
package jline.lib.m3a;

import java.util.Objects;

/**
 * Options for M3A compression.
 */
public final class M3aCompressOptions {
    private final M3aCompressMethod method;
    private final int numStates;
    private final boolean verbose;

    public M3aCompressOptions(M3aCompressMethod method, int numStates, boolean verbose) {
        this.method = method;
        this.numStates = numStates;
        this.verbose = verbose;
    }

    public M3aCompressOptions() {
        this(M3aCompressMethod.AMAP_2STATE, 2, false);
    }

    public M3aCompressMethod getMethod() { return method; }
    public int getNumStates() { return numStates; }
    public boolean isVerbose() { return verbose; }

    @Override
    public boolean equals(Object o) {
        if (this == o) return true;
        if (!(o instanceof M3aCompressOptions)) return false;
        M3aCompressOptions that = (M3aCompressOptions) o;
        return numStates == that.numStates
                && verbose == that.verbose
                && method == that.method;
    }

    @Override
    public int hashCode() {
        return Objects.hash(method, numStates, verbose);
    }

    @Override
    public String toString() {
        return "M3aCompressOptions(method=" + method + ", numStates=" + numStates
                + ", verbose=" + verbose + ")";
    }
}
