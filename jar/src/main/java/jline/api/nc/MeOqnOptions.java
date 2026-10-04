/**
 * @file Maximum Entropy OQN Algorithm Options
 *
 * Configuration options for the ME OQN algorithm.
 *
 * @since LINE 3.0
 */
package jline.api.nc;

import java.util.Objects;

/**
 * Options for the ME OQN algorithm.
 */
public final class MeOqnOptions {
    private final double tol;
    private final int maxIter;
    private final boolean verbose;

    public MeOqnOptions() {
        this(1e-6, 1000, false);
    }

    public MeOqnOptions(double tol) {
        this(tol, 1000, false);
    }

    public MeOqnOptions(double tol, int maxIter) {
        this(tol, maxIter, false);
    }

    public MeOqnOptions(double tol, int maxIter, boolean verbose) {
        this.tol = tol;
        this.maxIter = maxIter;
        this.verbose = verbose;
    }

    public double getTol() { return tol; }
    public int getMaxIter() { return maxIter; }
    public boolean getVerbose() { return verbose; }

    @Override
    public boolean equals(Object o) {
        if (this == o) return true;
        if (!(o instanceof MeOqnOptions)) return false;
        MeOqnOptions that = (MeOqnOptions) o;
        return Double.compare(that.tol, tol) == 0
                && maxIter == that.maxIter
                && verbose == that.verbose;
    }

    @Override
    public int hashCode() {
        return Objects.hash(tol, maxIter, verbose);
    }

    @Override
    public String toString() {
        return "MeOqnOptions(tol=" + tol + ", maxIter=" + maxIter + ", verbose=" + verbose + ")";
    }
}
