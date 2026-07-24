package jline.lib.smc;

import jline.util.matrix.Matrix;

/**
 * Options for MG1_pi solver.
 */
public class MG1PiOptions {
    private final Matrix boundary;
    private final int maxNumComp;
    private final int precision;
    private final String solver;
    private final boolean verbose;
    private final String mode;

    public MG1PiOptions() {
        this(null, 500, 200, "FI", false, "ShiftPWCR");
    }

    public MG1PiOptions(int verboseFlag) {
        this(null, 500, 200, "FI", verboseFlag != 0, "ShiftPWCR");
    }

    public MG1PiOptions(int verboseFlag, boolean verbose) {
        this(null, 500, 200, "FI", verbose, "ShiftPWCR");
    }

    public boolean isVerbose() { return verbose; }

    public MG1PiOptions(Matrix boundary, int maxNumComp, int precision, String solver, boolean verbose, String mode) {
        this.boundary = boundary;
        this.maxNumComp = maxNumComp;
        this.precision = precision;
        this.solver = solver;
        this.verbose = verbose;
        this.mode = mode;
    }

    public Matrix getBoundary() { return boundary; }
    public int getMaxNumComp() { return maxNumComp; }
    public int getPrecision() { return precision; }
    public String getSolver() { return solver; }
    public boolean getVerbose() { return verbose; }
    public String getMode() { return mode; }
}
