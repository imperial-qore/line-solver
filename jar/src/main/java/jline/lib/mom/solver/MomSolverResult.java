package jline.lib.mom.solver;

import java.util.Arrays;

import org.apache.commons.math3.linear.RealMatrix;

/**
 * Result container for MOM solver computations.
 *
 * <ul>
 *   <li>X Throughput matrix (M&times;R) - throughput of each class at each station</li>
 *   <li>Q Queue length matrix (M&times;R) - average queue length of each class at each station</li>
 *   <li>G Normalizing constant vector - normalizing constants for the queueing network</li>
 * </ul>
 */
public final class MomSolverResult {
    public final RealMatrix X;
    public final RealMatrix Q;
    public final double[] G;

    public MomSolverResult(RealMatrix X, RealMatrix Q, double[] G) {
        this.X = X;
        this.Q = Q;
        this.G = G;
    }

    public RealMatrix getX() { return X; }
    public RealMatrix getQ() { return Q; }
    public double[] getG() { return G; }

    @Override
    public boolean equals(Object other) {
        if (this == other) return true;
        if (other == null || getClass() != other.getClass()) return false;

        MomSolverResult that = (MomSolverResult) other;

        if (X != null ? !X.equals(that.X) : that.X != null) return false;
        if (Q != null ? !Q.equals(that.Q) : that.Q != null) return false;
        return Arrays.equals(G, that.G);
    }

    @Override
    public int hashCode() {
        int result = (X != null) ? X.hashCode() : 0;
        result = 31 * result + ((Q != null) ? Q.hashCode() : 0);
        result = 31 * result + Arrays.hashCode(G);
        return result;
    }

    @Override
    public String toString() {
        return "MomSolverResult(X=" + X + ", Q=" + Q + ", G=" + Arrays.toString(G) + ")";
    }
}
