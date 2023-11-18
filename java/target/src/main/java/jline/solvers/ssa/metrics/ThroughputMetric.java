package jline.solvers.ssa.metrics;

import java.util.ArrayList;

import jline.solvers.ssa.events.ArrivalEvent;
import jline.solvers.ssa.events.DepartureEvent;
import jline.solvers.ssa.events.Event;
import jline.solvers.ssa.state.SSAStateMatrix;
import jline.util.Pair;

public class ThroughputMetric extends Metric<Double, Double> {
    public ThroughputMetric(int nodeIdx, int classIdx, int nServers, boolean record) {
        super("Throughput", 0.0, 0.0, record, nodeIdx, classIdx, nServers);
    }
    public ThroughputMetric(int nodeIdx, int classIdx, int nServers) {
        this(nodeIdx, classIdx, nServers, false);
    }

    protected void addSample(double currentTime, Double metric) {
    }

    public Double getMetric() {
        if (this.useMSER5) {
            return this.cutoffMSER5().getRight();
        } else if (this.useR5) {
            return this.cutoffR5().getRight();
        }
        if (Double.isNaN(this.metricValue)) {
            return 0.0;
        }
        return this.metricValue;
    }

    public void fromStateMatrix(double t, SSAStateMatrix networkState) {
        /*double stateValue = stateMatrix.getState(i,j);
        double throughput = Math.max(stateMatrix.getState(i,j) - this.prevState, 0);
        this.addSample(t, throughput);
        this.prevState = stateValue;*/
    }

    public void resetHistory() {
        this.metricHistory = new ArrayList<Pair<Double, Double>>(5000);
        this.metricValue = this.initialMetric;
        this.currentValue = this.initialValue;
        this.cutoffTime = this.time;
    }

    @SuppressWarnings("unchecked")
    public void fromEvent (double t, Event e) {
        if (e instanceof DepartureEvent) {
            DepartureEvent de = (DepartureEvent) e;
            if ((de.getNode().getStatefulIdx() == this.nodeIdx) &&
                    (de.getJobClass().getJobClassIdx() == this.classIdx)){
                double tDelta = t - this.time;
                double totalThroughput = this.metricValue * (this.time-this.cutoffTime);
                totalThroughput += 1.0;
                if (this.record) {
                    this.metricHistory.add(new Pair(this.time, 1.0/tDelta));
                } else {
                    if (this.useR5 && (this.crossCount < this.r5Value)) {
                        if (this.belowAverage && ((1.0/tDelta) >= this.metricValue)) {
                            this.crossCount++;
                            this.belowAverage = false;

                            if (this.crossCount >= this.r5Value) {
                                this.cutoffTime = t;
                                this.metricValue = 0.0;
                                return;
                            }
                        } else if (!this.belowAverage && ((1.0/tDelta) < this.metricValue)) {
                            this.crossCount++;
                            this.belowAverage = true;

                            if (this.crossCount >= this.r5Value) {
                                this.cutoffTime = t;
                                this.metricValue = 0.0;
                                return;
                            }
                        }
                    }
                }
                this.time = t;
                this.metricValue = totalThroughput/(this.time-this.cutoffTime);
            }
        }
    }
    @SuppressWarnings("unchecked")
    public void fromEvent(double t, Event e, int n) {
        if (e instanceof DepartureEvent) {
            DepartureEvent de = (DepartureEvent) e;
            if ((de.getNode().getStatefulIdx() == this.nodeIdx) &&
                    (de.getJobClass().getJobClassIdx() == this.classIdx)){
                double tDelta = t - this.time;
                double totalThroughput = this.metricValue * (this.time-this.cutoffTime);
                totalThroughput += n;

                if (this.record) {
                    this.metricHistory.add(new Pair(this.time, n/tDelta));
                } else {
                    if (this.useR5) {
                        if (this.belowAverage && ((n/tDelta) >= this.metricValue)) {
                            this.crossCount++;
                            this.belowAverage = false;

                            if (this.crossCount >= this.r5Value) {
                                this.cutoffTime = t;
                                this.metricValue = 0.0;
                            }
                        } else if (!this.belowAverage && ((n/tDelta) < this.metricValue)) {
                            this.crossCount++;
                            this.belowAverage = true;

                            if (this.crossCount >= this.r5Value) {
                                this.cutoffTime = t;
                                this.metricValue = 0.0;
                            }
                        }
                    }
                }
                this.time = t;
                this.metricValue = totalThroughput/(this.time-this.cutoffTime);
            }
        }
    }

    @Override
    @SuppressWarnings("unchecked")
    public void taper(double t) {
        double totalThroughput = this.metricValue * this.time;
        this.metricValue = totalThroughput/this.time;
        if (this.record) {
            this.metricHistory.add(new Pair(t, 0.0));
        }
        this.time = t;
    }

    @Override
    public boolean recognizeEvent(Event e) {
        if (e instanceof DepartureEvent) {
            DepartureEvent de = (DepartureEvent)e;
            return ((de.getNode().getStatefulIdx() == this.nodeIdx) &&
                    (de.getJobClass().getJobClassIdx() == this.classIdx));
        } else if (e instanceof ArrivalEvent) {
            ArrivalEvent ae = (ArrivalEvent)e;
            return ((ae.getNode().getStatefulIdx() == this.nodeIdx) &&
                    (ae.getJobClass().getJobClassIdx() == this.classIdx));
        }
        return false;
    }
}
