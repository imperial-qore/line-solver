/**
 * @file CountStatistics data class for multi-class trace
 *
 * @since LINE 3.0
 */
package jline.api.trace;

import java.util.Arrays;
import java.util.Objects;

/**
 * Data class for count statistics
 */
public final class CountStatistics {
    private final double windowSize;
    private final int numWindows;
    private final SingleClassCountStats[] classStats;
    private final double totalTime;

    public CountStatistics(double windowSize, int numWindows,
                           SingleClassCountStats[] classStats, double totalTime) {
        this.windowSize = windowSize;
        this.numWindows = numWindows;
        this.classStats = classStats;
        this.totalTime = totalTime;
    }

    public double getWindowSize() { return windowSize; }
    public int getNumWindows() { return numWindows; }
    public SingleClassCountStats[] getClassStats() { return classStats; }
    public double getTotalTime() { return totalTime; }

    @Override
    public boolean equals(Object o) {
        if (this == o) return true;
        if (!(o instanceof CountStatistics)) return false;
        CountStatistics that = (CountStatistics) o;
        return Double.compare(that.windowSize, windowSize) == 0
                && numWindows == that.numWindows
                && Double.compare(that.totalTime, totalTime) == 0
                && Arrays.equals(classStats, that.classStats);
    }

    @Override
    public int hashCode() {
        int result = Objects.hash(windowSize, numWindows, totalTime);
        result = 31 * result + Arrays.hashCode(classStats);
        return result;
    }

    @Override
    public String toString() {
        return "CountStatistics(windowSize=" + windowSize + ", numWindows=" + numWindows
                + ", classStats=" + Arrays.toString(classStats) + ", totalTime=" + totalTime + ")";
    }
}
