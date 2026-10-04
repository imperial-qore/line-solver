/**
 * @file AoI topology validation result type
 *
 * @since LINE 3.0
 */
package jline.api.aoi;

import java.util.Objects;

/**
 * Result of AoI topology validation.
 */
public final class AoiValidationResult {
    private final boolean isAoI;
    private final int sourceIdx;
    private final int queueIdx;
    private final int sinkIdx;
    private final int sourceStation;
    private final int queueStation;
    private final int capacity;
    private final String schedStrategy;
    private final String systemType;
    private final String errorMsg;

    public AoiValidationResult(boolean isAoI, int sourceIdx, int queueIdx, int sinkIdx,
                               int sourceStation, int queueStation, int capacity,
                               String schedStrategy, String systemType, String errorMsg) {
        this.isAoI = isAoI;
        this.sourceIdx = sourceIdx;
        this.queueIdx = queueIdx;
        this.sinkIdx = sinkIdx;
        this.sourceStation = sourceStation;
        this.queueStation = queueStation;
        this.capacity = capacity;
        this.schedStrategy = schedStrategy;
        this.systemType = systemType;
        this.errorMsg = errorMsg;
    }

    public boolean isAoI() { return isAoI; }
    public int getSourceIdx() { return sourceIdx; }
    public int getQueueIdx() { return queueIdx; }
    public int getSinkIdx() { return sinkIdx; }
    public int getSourceStation() { return sourceStation; }
    public int getQueueStation() { return queueStation; }
    public int getCapacity() { return capacity; }
    public String getSchedStrategy() { return schedStrategy; }
    public String getSystemType() { return systemType; }
    public String getErrorMsg() { return errorMsg; }

    @Override
    public boolean equals(Object o) {
        if (this == o) return true;
        if (!(o instanceof AoiValidationResult)) return false;
        AoiValidationResult that = (AoiValidationResult) o;
        return isAoI == that.isAoI
                && sourceIdx == that.sourceIdx
                && queueIdx == that.queueIdx
                && sinkIdx == that.sinkIdx
                && sourceStation == that.sourceStation
                && queueStation == that.queueStation
                && capacity == that.capacity
                && Objects.equals(schedStrategy, that.schedStrategy)
                && Objects.equals(systemType, that.systemType)
                && Objects.equals(errorMsg, that.errorMsg);
    }

    @Override
    public int hashCode() {
        return Objects.hash(isAoI, sourceIdx, queueIdx, sinkIdx, sourceStation, queueStation,
                capacity, schedStrategy, systemType, errorMsg);
    }

    @Override
    public String toString() {
        return "AoiValidationResult(isAoI=" + isAoI + ", sourceIdx=" + sourceIdx
                + ", queueIdx=" + queueIdx + ", sinkIdx=" + sinkIdx
                + ", sourceStation=" + sourceStation + ", queueStation=" + queueStation
                + ", capacity=" + capacity + ", schedStrategy=" + schedStrategy
                + ", systemType=" + systemType + ", errorMsg=" + errorMsg + ")";
    }
}
