package jline.solvers.ssa.metrics;

import jline.solvers.ssa.events.Event;
import jline.solvers.ssa.state.SSAStateMatrix;
import jline.util.Pair;

public class UtilizationMetric extends Metric<Double, Double> {
    protected boolean isDelay;

    public UtilizationMetric(int nodeIdx, int classIdx, int nServers, boolean record, boolean isDelay) {
        super("Utilization", 0.0, 0.0, record, nodeIdx, classIdx, nServers);
        this.useMatrix = true;
        this.isDelay = isDelay;
    }

    public UtilizationMetric(int nodeIdx, int classIdx, int nServers, boolean record) {
        super("Utilization", 0.0, 0.0, record, nodeIdx, classIdx, nServers);
        this.useMatrix = true;
        this.isDelay = false;
    }

    public UtilizationMetric(int nodeIdx, int classIdx, int nServers) {
        this(nodeIdx, classIdx, nServers, false);
    }

    public void addSample(double currentTime, Double metric) {
    }

    public Double getMetric() {
        if (this.useMSER5) {
            return this.cutoffMSER5().getRight();
        } else if (this.useR5) {
            return this.cutoffR5().getRight();
        }

        if (this.record) {
            double cumulativeUT = 0;
            double prevTime = 0;
            double timeIter;
            double prevValue = 0;
            if (this.metricHistory.size() == 0) {
                return 0.0;
            }
            for (Pair<Double, Double> timelineEntry : this.metricHistory) {
                timeIter = timelineEntry.getLeft();
                cumulativeUT += (timeIter-prevTime)*prevValue;
                prevTime = timeIter;
                prevValue = timelineEntry.getRight();
            }
            cumulativeUT += (this.time-prevTime)*prevValue;
            prevTime = this.time;
            if (prevTime == 0) {
                return 0.0;
            }
            return cumulativeUT/prevTime;
        } else {
            return this.metricValue;
        }
    }

    @SuppressWarnings("unchecked")
    public void fromStateMatrix(double t, SSAStateMatrix networkState) {
        double inProcess = networkState.inProcess(this.nodeIdx, this.classIdx);
        double utilization;

        if (this.isDelay) {
            utilization = inProcess;
        } else if (nServers == 0) {
            return;
        } else {
            utilization = inProcess/nServers;
        }

        if (utilization == this.metricValue) {
            return;
        }
        if (this.record) {
            this.metricValue = utilization;
            this.metricHistory.add(new Pair(t, utilization));
        } else {
            if (this.useR5 && (this.crossCount < this.r5Value)) {
                if (this.belowAverage && (utilization >= this.metricValue)) {
                    this.crossCount++;
                    this.belowAverage = false;

                    if (this.crossCount >= this.r5Value) {
                        this.cutoffTime = t;
                        this.metricValue = 0.0;
                        return;
                    }
                } else if (!this.belowAverage && (utilization < this.metricValue)) {
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
            this.currentValue = utilization;
            this.time = t;
        }
    }

    public void fromEvent(double t, Event e) {
    }
    public void fromEvent(double t, Event e, int n) {

    }
}
