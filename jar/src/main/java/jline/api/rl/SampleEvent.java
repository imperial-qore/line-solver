/**
 * @file Sampled event from RL environment
 *
 * @since LINE 3.0
 */
package jline.api.rl;

/**
 * Container for a sampled event from the RL environment.
 */
public final class SampleEvent {
    public final double t;
    public final int depNode;

    public SampleEvent(double t, int depNode) {
        this.t = t;
        this.depNode = depNode;
    }

    public double getT() { return t; }
    public int getDepNode() { return depNode; }

    public double component1() { return t; }
    public int component2() { return depNode; }

    @Override
    public boolean equals(Object o) {
        if (this == o) return true;
        if (!(o instanceof SampleEvent)) return false;
        SampleEvent that = (SampleEvent) o;
        return Double.compare(that.t, t) == 0 && depNode == that.depNode;
    }

    @Override
    public int hashCode() {
        int result = Double.hashCode(t);
        result = 31 * result + depNode;
        return result;
    }

    @Override
    public String toString() {
        return "SampleEvent(t=" + t + ", depNode=" + depNode + ")";
    }
}
