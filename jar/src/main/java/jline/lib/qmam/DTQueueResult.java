/**
 * @file Result of a discrete-time single-server queue analysis
 *
 * @since LINE 3.1.0
 */
package jline.lib.qmam;

import java.util.Objects;

import jline.util.matrix.Matrix;

/**
 * Result of a discrete-time (slotted) single-server queue analysis.
 *
 * <p>The queue length is read at slot boundaries under the late arrival system
 * with delayed access: within a slot the completion resolves first, the arrival
 * is appended at the end of the slot and cannot enter service before the next
 * one, and the level is taken after both. Entry {@code i} of
 * {@link #queueLength} is Prob[i customers in system], starting at i = 0.
 */
public final class DTQueueResult {
    public final Matrix queueLength;

    public DTQueueResult(Matrix queueLength) {
        this.queueLength = queueLength;
    }

    public Matrix getQueueLength() { return queueLength; }

    public Matrix component1() { return queueLength; }

    /** Mean number in system, sum_i i * P[i in system]. */
    public double getMeanQueueLength() {
        double q = 0;
        for (int i = 0; i < queueLength.getNumCols(); i++) {
            q += i * queueLength.get(0, i);
        }
        return q;
    }

    /** Utilization, i.e. the fraction of slots in which the server is busy. */
    public double getUtilization() {
        return 1.0 - queueLength.get(0, 0);
    }

    @Override
    public boolean equals(Object o) {
        if (this == o) return true;
        if (!(o instanceof DTQueueResult)) return false;
        DTQueueResult that = (DTQueueResult) o;
        return Objects.equals(queueLength, that.queueLength);
    }

    @Override
    public int hashCode() {
        return Objects.hash(queueLength);
    }

    @Override
    public String toString() {
        return "DTQueueResult(queueLength=" + queueLength + ")";
    }
}
