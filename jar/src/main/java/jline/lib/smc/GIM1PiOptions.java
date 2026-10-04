/**
 * @file Options for GIM1_pi solver
 *
 * @since LINE 3.1.0
 */
package jline.lib.smc;

import java.util.Objects;

import jline.util.matrix.Matrix;

/**
 * Options for GIM1_pi solver.
 */
public final class GIM1PiOptions {
    public final Matrix boundary;
    public final int maxNumComp;
    public final int verbose;

    public GIM1PiOptions(Matrix boundary, int maxNumComp, int verbose) {
        this.boundary = boundary;
        this.maxNumComp = maxNumComp;
        this.verbose = verbose;
    }

    public GIM1PiOptions() {
        this(null, 500, 0);
    }

    public Matrix getBoundary() { return boundary; }
    public int getMaxNumComp() { return maxNumComp; }
    public int getVerbose() { return verbose; }

    public Matrix component1() { return boundary; }
    public int component2() { return maxNumComp; }
    public int component3() { return verbose; }

    @Override
    public boolean equals(Object o) {
        if (this == o) return true;
        if (!(o instanceof GIM1PiOptions)) return false;
        GIM1PiOptions that = (GIM1PiOptions) o;
        return maxNumComp == that.maxNumComp
                && verbose == that.verbose
                && Objects.equals(boundary, that.boundary);
    }

    @Override
    public int hashCode() {
        return Objects.hash(boundary, maxNumComp, verbose);
    }

    @Override
    public String toString() {
        return "GIM1PiOptions(boundary=" + boundary + ", maxNumComp=" + maxNumComp + ", verbose=" + verbose + ")";
    }
}
