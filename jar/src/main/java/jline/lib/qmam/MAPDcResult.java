package jline.lib.qmam;

import jline.util.matrix.Matrix;

/** Result of MAP/D/c queue analysis. */
public final class MAPDcResult {
    public final Matrix queueLength;
    public final Matrix waitingTime;

    public MAPDcResult(Matrix queueLength, Matrix waitingTime) {
        this.queueLength = queueLength;
        this.waitingTime = waitingTime;
    }

    public Matrix getQueueLength() { return queueLength; }
    public Matrix getWaitingTime() { return waitingTime; }
}
