/**
 * @file Container for a sampled event from the general RL environment
 *
 * @since LINE 3.0
 */
package jline.api.rl;

import java.util.Objects;

import jline.io.Ret.SampleResult;

/**
 * Container for a sampled event from the general RL environment.
 */
public final class GeneralSampleEvent {
    public final double dt;
    public final int depNode;
    public final int arvNode;
    public final SampleResult sampleResult;

    public GeneralSampleEvent(double dt, int depNode, int arvNode, SampleResult sampleResult) {
        this.dt = dt;
        this.depNode = depNode;
        this.arvNode = arvNode;
        this.sampleResult = sampleResult;
    }

    public double getDt() { return dt; }
    public int getDepNode() { return depNode; }
    public int getArvNode() { return arvNode; }
    public SampleResult getSampleResult() { return sampleResult; }

    public double component1() { return dt; }
    public int component2() { return depNode; }
    public int component3() { return arvNode; }
    public SampleResult component4() { return sampleResult; }

    @Override
    public boolean equals(Object o) {
        if (this == o) return true;
        if (!(o instanceof GeneralSampleEvent)) return false;
        GeneralSampleEvent that = (GeneralSampleEvent) o;
        return Double.compare(that.dt, dt) == 0
                && depNode == that.depNode
                && arvNode == that.arvNode
                && Objects.equals(sampleResult, that.sampleResult);
    }

    @Override
    public int hashCode() {
        return Objects.hash(dt, depNode, arvNode, sampleResult);
    }

    @Override
    public String toString() {
        return "GeneralSampleEvent(dt=" + dt + ", depNode=" + depNode
                + ", arvNode=" + arvNode + ", sampleResult=" + sampleResult + ")";
    }
}
