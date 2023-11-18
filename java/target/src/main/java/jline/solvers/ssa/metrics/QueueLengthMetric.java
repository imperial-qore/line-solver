package jline.solvers.ssa.metrics;

import jline.solvers.ssa.events.Event;
import jline.solvers.ssa.state.SSAStateMatrix;
import jline.util.Pair;

public class QueueLengthMetric extends Metric<Double, Double> {
    public QueueLengthMetric(int nodeIdx, int classIdx, int nServers, boolean record) {
        super("Queue Length", 0.0, 0.0, record, nodeIdx, classIdx, nServers);
        this.useMatrix = true;
    }

    public QueueLengthMetric (int nodeIdx, int classIdx, int nServers) {
        this(nodeIdx, classIdx, nServers, false);
    }

    @SuppressWarnings("unchecked")
    protected void addSample(double currentTime, Double metric) {
        double timeDelta = currentTime-this.time;
        this.metricValue = ((this.metricValue * this.time) +
                (this.currentValue*timeDelta))/(currentTime);
        this.currentValue = metric;
        this.time = currentTime;
        if (this.record) {
            this.metricHistory.add(new Pair(currentTime, metric));
        }
    }

    public Double getMetric() {
        if (this.useMSER5) {
            return this.cutoffMSER5().getRight();
        } else if (this.useR5) {
            return this.cutoffR5().getRight();
        }

        if (this.record) {
            double cumulativeQT = 0;
            double prevTime = 0;
            double timeIter;
            double prevValue = 0;
            if (this.metricHistory.size() == 0) {
                return 0.0;
            }
            for (Pair<Double, Double> timelineEntry : this.metricHistory) {
                timeIter = timelineEntry.getLeft();
                cumulativeQT += (timeIter-prevTime)*prevValue;
                prevTime = timeIter;
                prevValue = timelineEntry.getRight();
            }
            cumulativeQT += (this.time-prevTime)*prevValue;
            prevTime = this.time;
            return cumulativeQT/prevTime;
        } else {
            return this.metricValue;
        }
    }
    @SuppressWarnings("unchecked")
    public void fromStateMatrix(double t, SSAStateMatrix networkState) {
        double Q = Math.min(networkState.getState(this.nodeIdx,this.classIdx),
                networkState.getCapacity(this.nodeIdx,this.classIdx));

        if (Q == this.metricValue) {
            return;
        }
        if (this.record) {
            this.metricValue = Q;
            this.metricHistory.add(new Pair(t, Q));
        } else {
            if (this.useR5 && (this.crossCount < this.r5Value)) {
                if (this.belowAverage && (Q >= this.metricValue)) {
                    this.crossCount++;
                    this.belowAverage = false;

                    if (this.crossCount >= this.r5Value) {
                        this.cutoffTime = t;
                        this.metricValue = 0.0;
                        return;
                    }
                } else if (!this.belowAverage && (Q < this.metricValue)) {
                    this.crossCount++;
                    this.belowAverage = true;

                    if (this.crossCount >= this.r5Value) {
                        this.cutoffTime = t;
                        this.metricValue = 0.0;
                        return;
                    }
                }
            }
            double timeDelta = t-this.time;
            this.metricValue = ((this.metricValue * (this.time-this.cutoffTime)) +
                    (this.currentValue*timeDelta))/(t-this.cutoffTime);
            this.currentValue = Q;
            this.time = t;
        }
    }

    public void fromEvent(double t, Event e) {

    }
    public void fromEvent(double t, Event e, int n) {

    }
}
