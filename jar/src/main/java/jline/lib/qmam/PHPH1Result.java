/**
 * @file Result of PH/PH/1 queue analysis
 *
 * @since LINE 3.1.0
 */
package jline.lib.qmam;

import java.util.Objects;

import jline.util.matrix.Matrix;

/**
 * Result of PH/PH/1 queue analysis.
 */
public final class PHPH1Result {
    public final Matrix queueLength; // Queue length distribution
    public final Matrix waitAlpha;   // Waiting time PH alpha vector
    public final Matrix waitT;       // Waiting time PH matrix

    public PHPH1Result(Matrix queueLength, Matrix waitAlpha, Matrix waitT) {
        this.queueLength = queueLength;
        this.waitAlpha = waitAlpha;
        this.waitT = waitT;
    }

    public Matrix getQueueLength() { return queueLength; }
    public Matrix getWaitAlpha() { return waitAlpha; }
    public Matrix getWaitT() { return waitT; }

    public Matrix component1() { return queueLength; }
    public Matrix component2() { return waitAlpha; }
    public Matrix component3() { return waitT; }

    @Override
    public boolean equals(Object o) {
        if (this == o) return true;
        if (!(o instanceof PHPH1Result)) return false;
        PHPH1Result that = (PHPH1Result) o;
        return Objects.equals(queueLength, that.queueLength)
                && Objects.equals(waitAlpha, that.waitAlpha)
                && Objects.equals(waitT, that.waitT);
    }

    @Override
    public int hashCode() {
        return Objects.hash(queueLength, waitAlpha, waitT);
    }

    @Override
    public String toString() {
        return "PHPH1Result(queueLength=" + queueLength + ", waitAlpha=" + waitAlpha + ", waitT=" + waitT + ")";
    }
}
