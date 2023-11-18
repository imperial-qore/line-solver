package jline.solvers.ssa.metrics;

import jline.solvers.ssa.events.Event;
import jline.solvers.ssa.state.SSAStateMatrix;

public class DummyMetric extends Metric<Double, Double>{
    public DummyMetric(int nodeIdx, int classIdx, int nServers) {
        super("None", Double.NaN,Double.NaN, nodeIdx, classIdx, nServers);
    }

    protected void addSample(double currentTime, Double metric) {
    }

    public void fromStateMatrix(double t, SSAStateMatrix networkState) {

    }

    public Double getMetric() {
        return Double.NaN;
    }

    public boolean useMatrix() {
        return false;
    }

    public void fromEvent(double t, Event e) {

    }


    public void fromEvent(double t, Event e, int n) {

    }
}
