package jline.solvers.ctmc;

import jline.solvers.SolverResult;
import jline.util.matrix.Matrix;

/**
 * Result class for CTMC Marginal Aggregated solver.
 * Contains the marginal probabilities, steady-state distribution, runtime, and filename.
 */
public class ResultCTMCMargAggr extends SolverResult {
    
    private final Matrix Pnir;     // Marginal probabilities at stations
    private final Matrix pi;       // Steady-state distribution
    private final double runtime;  // Computation runtime in seconds
    private final String fname;    // Workspace filename (if saved)
    
    public ResultCTMCMargAggr(Matrix Pnir, Matrix pi, double runtime, String fname) {
        this.Pnir = Pnir;
        this.pi = pi;
        this.runtime = runtime;
        this.fname = fname;
    }
    
    /**
     * Get the marginal probabilities at stations
     * @return Matrix of marginal probabilities (1 x nstations)
     */
    public Matrix getPnir() {
        return Pnir;
    }
    
    /**
     * Get the steady-state distribution
     * @return Column vector of steady-state probabilities
     */
    public Matrix getPi() {
        return pi;
    }
    
    /**
     * Get the computation runtime
     * @return Runtime in seconds
     */
    public double getRuntime() {
        return runtime;
    }
    
    /**
     * Get the workspace filename
     * @return Filename where workspace was saved, or empty string if not saved
     */
    public String getFname() {
        return fname;
    }
}