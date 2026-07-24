package jline.lib.smc;

/**
 * Options for MG1_CR solver.
 */
public class MG1CROptions {
    private final String mode;
    private final int maxNumIt;
    private final int maxNumRoot;
    private final double epsilonValue;
    private final int verbose;
    private final String shiftType;

    public MG1CROptions() {
        this("ShiftPWCR", 50, 2048, 1e-16, 0, "one");
    }

    public MG1CROptions(String mode, int maxNumIt, int maxNumRoot, double epsilonValue, int verbose, String shiftType) {
        this.mode = mode;
        this.maxNumIt = maxNumIt;
        this.maxNumRoot = maxNumRoot;
        this.epsilonValue = epsilonValue;
        this.verbose = verbose;
        this.shiftType = shiftType;
    }

    public String getMode() { return mode; }
    public int getMaxNumIt() { return maxNumIt; }
    public int getMaxNumRoot() { return maxNumRoot; }
    public double getEpsilonValue() { return epsilonValue; }
    public int getVerbose() { return verbose; }
    public String getShiftType() { return shiftType; }
}
