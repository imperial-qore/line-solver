/**
 * @file Options for MG1_FI solver
 *
 * @since LINE 3.1.0
 */
package jline.lib.smc;

import java.util.Arrays;
import java.util.Objects;

import jline.util.matrix.Matrix;

public final class MG1FIOptions {
    private final String mode;
    private final int maxNumIt;
    private final int verbose;
    private final String shiftType;
    private final Matrix startValue;
    private final int[] nonZeroBlocks;

    public MG1FIOptions(String mode, int maxNumIt, int verbose, String shiftType,
                        Matrix startValue, int[] nonZeroBlocks) {
        this.mode = mode;
        this.maxNumIt = maxNumIt;
        this.verbose = verbose;
        this.shiftType = shiftType;
        this.startValue = startValue;
        this.nonZeroBlocks = nonZeroBlocks;
    }

    public MG1FIOptions() {
        this("U-Based", 10000, 0, "one", null, null);
    }

    public MG1FIOptions(int verbose) {
        this("U-Based", 10000, verbose, "one", null, null);
    }

    public MG1FIOptions(String mode, int maxNumIt, int verbose) {
        this(mode, maxNumIt, verbose, "one", null, null);
    }

    public String getMode() { return mode; }
    public int getMaxNumIt() { return maxNumIt; }
    public int getVerbose() { return verbose; }
    public String getShiftType() { return shiftType; }
    public Matrix getStartValue() { return startValue; }
    public int[] getNonZeroBlocks() { return nonZeroBlocks; }

    @Override
    public boolean equals(Object o) {
        if (this == o) return true;
        if (!(o instanceof MG1FIOptions)) return false;
        MG1FIOptions that = (MG1FIOptions) o;
        return maxNumIt == that.maxNumIt && verbose == that.verbose
                && Objects.equals(mode, that.mode)
                && Objects.equals(shiftType, that.shiftType)
                && Objects.equals(startValue, that.startValue)
                && Arrays.equals(nonZeroBlocks, that.nonZeroBlocks);
    }

    @Override
    public int hashCode() {
        int result = Objects.hash(mode, maxNumIt, verbose, shiftType, startValue);
        result = 31 * result + Arrays.hashCode(nonZeroBlocks);
        return result;
    }
}
