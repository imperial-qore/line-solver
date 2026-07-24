/**
 * @file BMAP sample result
 *
 * @since LINE 3.0
 */
package jline.api.mam;

/**
 * Result of BMAP sampling containing inter-arrival time and batch size.
 */
public final class BmapSample {
    private final double interarrivalTime;
    private final int batchSize;

    public BmapSample(double interarrivalTime, int batchSize) {
        this.interarrivalTime = interarrivalTime;
        this.batchSize = batchSize;
    }

    public double getInterarrivalTime() { return interarrivalTime; }
    public int getBatchSize() { return batchSize; }

    @Override
    public boolean equals(Object o) {
        if (this == o) return true;
        if (!(o instanceof BmapSample)) return false;
        BmapSample that = (BmapSample) o;
        return Double.compare(that.interarrivalTime, interarrivalTime) == 0
                && batchSize == that.batchSize;
    }

    @Override
    public int hashCode() {
        int r = Double.hashCode(interarrivalTime);
        r = 31 * r + Integer.hashCode(batchSize);
        return r;
    }

    @Override
    public String toString() {
        return "BmapSample(interarrivalTime=" + interarrivalTime + ", batchSize=" + batchSize + ")";
    }
}
