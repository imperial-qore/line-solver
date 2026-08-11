/**
 * @file Information about detected Fork-Join topology
 *
 * Copyright (c) 2012-2026, Imperial College London
 * All rights reserved.
 */
package jline.api.fj;

import java.util.Arrays;
import java.util.Objects;

/**
 * Information about detected Fork-Join topology.
 */
public final class FJInfo {
    public final int K;
    public final int forkIdx;
    public final int joinIdx;
    public final int[] queueIndices;
    public final int sourceIdx;
    public final int sinkIdx;
    public final boolean isValid;

    public FJInfo(int K, int forkIdx, int joinIdx, int[] queueIndices,
                  int sourceIdx, int sinkIdx, boolean isValid) {
        this.K = K;
        this.forkIdx = forkIdx;
        this.joinIdx = joinIdx;
        this.queueIndices = queueIndices;
        this.sourceIdx = sourceIdx;
        this.sinkIdx = sinkIdx;
        this.isValid = isValid;
    }

    public int getK() { return K; }
    public int getForkIdx() { return forkIdx; }
    public int getJoinIdx() { return joinIdx; }
    public int[] getQueueIndices() { return queueIndices; }
    public int getSourceIdx() { return sourceIdx; }
    public int getSinkIdx() { return sinkIdx; }
    public boolean isValid() { return isValid; }

    public int component1() { return K; }
    public int component2() { return forkIdx; }
    public int component3() { return joinIdx; }
    public int[] component4() { return queueIndices; }
    public int component5() { return sourceIdx; }
    public int component6() { return sinkIdx; }
    public boolean component7() { return isValid; }

    @Override
    public boolean equals(Object o) {
        if (this == o) return true;
        if (!(o instanceof FJInfo)) return false;
        FJInfo that = (FJInfo) o;
        return K == that.K
                && forkIdx == that.forkIdx
                && joinIdx == that.joinIdx
                && sourceIdx == that.sourceIdx
                && sinkIdx == that.sinkIdx
                && isValid == that.isValid
                && Arrays.equals(queueIndices, that.queueIndices);
    }

    @Override
    public int hashCode() {
        int result = Objects.hash(K, forkIdx, joinIdx, sourceIdx, sinkIdx, isValid);
        result = 31 * result + Arrays.hashCode(queueIndices);
        return result;
    }

    @Override
    public String toString() {
        return "FJInfo(K=" + K + ", forkIdx=" + forkIdx + ", joinIdx=" + joinIdx
                + ", queueIndices=" + Arrays.toString(queueIndices)
                + ", sourceIdx=" + sourceIdx + ", sinkIdx=" + sinkIdx
                + ", isValid=" + isValid + ")";
    }
}
