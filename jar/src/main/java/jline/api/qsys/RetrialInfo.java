/**
 * @file Information about a valid retrial queue topology
 *
 * @since LINE 3.1.0
 */
package jline.api.qsys;

import java.util.Objects;

/**
 * Information about a valid retrial queue topology.
 */
public final class RetrialInfo {
    private boolean isRetrial;
    private int stationIdx;
    private int nodeIdx;
    private int sourceIdx;
    private int classIdx;
    private String errorMsg;
    private int N;
    private double alpha;
    private double gamma;
    private double p;
    private int R;

    public RetrialInfo() {
        this(false, -1, -1, -1, -1, "", 0, 0.1, 0.0, 0.0, 0);
    }

    public RetrialInfo(boolean isRetrial, int stationIdx, int nodeIdx, int sourceIdx,
                       int classIdx, String errorMsg, int N, double alpha,
                       double gamma, double p, int R) {
        this.isRetrial = isRetrial;
        this.stationIdx = stationIdx;
        this.nodeIdx = nodeIdx;
        this.sourceIdx = sourceIdx;
        this.classIdx = classIdx;
        this.errorMsg = errorMsg;
        this.N = N;
        this.alpha = alpha;
        this.gamma = gamma;
        this.p = p;
        this.R = R;
    }

    public boolean isRetrial() { return isRetrial; }
    public void setRetrial(boolean isRetrial) { this.isRetrial = isRetrial; }

    public int getStationIdx() { return stationIdx; }
    public void setStationIdx(int stationIdx) { this.stationIdx = stationIdx; }

    public int getNodeIdx() { return nodeIdx; }
    public void setNodeIdx(int nodeIdx) { this.nodeIdx = nodeIdx; }

    public int getSourceIdx() { return sourceIdx; }
    public void setSourceIdx(int sourceIdx) { this.sourceIdx = sourceIdx; }

    public int getClassIdx() { return classIdx; }
    public void setClassIdx(int classIdx) { this.classIdx = classIdx; }

    public String getErrorMsg() { return errorMsg; }
    public void setErrorMsg(String errorMsg) { this.errorMsg = errorMsg; }

    public int getN() { return N; }
    public void setN(int N) { this.N = N; }

    public double getAlpha() { return alpha; }
    public void setAlpha(double alpha) { this.alpha = alpha; }

    public double getGamma() { return gamma; }
    public void setGamma(double gamma) { this.gamma = gamma; }

    public double getP() { return p; }
    public void setP(double p) { this.p = p; }

    public int getR() { return R; }
    public void setR(int R) { this.R = R; }

    @Override
    public boolean equals(Object o) {
        if (this == o) return true;
        if (!(o instanceof RetrialInfo)) return false;
        RetrialInfo that = (RetrialInfo) o;
        return isRetrial == that.isRetrial && stationIdx == that.stationIdx
                && nodeIdx == that.nodeIdx && sourceIdx == that.sourceIdx
                && classIdx == that.classIdx && N == that.N
                && Double.compare(that.alpha, alpha) == 0
                && Double.compare(that.gamma, gamma) == 0
                && Double.compare(that.p, p) == 0 && R == that.R
                && Objects.equals(errorMsg, that.errorMsg);
    }

    @Override
    public int hashCode() {
        return Objects.hash(isRetrial, stationIdx, nodeIdx, sourceIdx, classIdx,
                errorMsg, N, alpha, gamma, p, R);
    }

    @Override
    public String toString() {
        return "RetrialInfo(isRetrial=" + isRetrial + ", stationIdx=" + stationIdx
                + ", N=" + N + ", R=" + R + ", errorMsg=" + errorMsg + ")";
    }
}
