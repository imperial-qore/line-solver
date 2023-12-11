package jline.solvers.ssa.metrics;

import jline.lang.constant.SchedStrategy;
import jline.solvers.ssa.events.ArrivalEvent;
import jline.solvers.ssa.events.DepartureEvent;
import jline.solvers.ssa.events.Event;
import jline.solvers.ssa.state.SSAStateMatrix;
import jline.util.Pair;

import java.util.ArrayList;
import java.util.Deque;
import java.util.LinkedList;

public class ResponseTimeMetric extends Metric<Double, Double> {
    int nDepartures;
    Deque<Double> buffer;
    SchedStrategy schedStrategy;
    public ResponseTimeMetric(int nodeIdx, int classIdx, int nServers, SchedStrategy schedStrategy, boolean record) {
        super("Response Time", 0.0, 0.0, record, nodeIdx, classIdx, nServers);
        this.shortName = "RpT";
        this.nDepartures = 0;
        this.buffer = new LinkedList<Double>();
        this.schedStrategy = schedStrategy;
    }

    public ResponseTimeMetric (int nodeIdx, int classIdx, int nServers, SchedStrategy schedStrategy) {
        this(nodeIdx, classIdx, nServers, schedStrategy, true);
    }
    protected void addSample(double currentTime, Double metric) {
    }

    public Double getMetric() {
        if (this.useMSER5) {
            return this.cutoffMSER5().getRight();
        } else if (this.useR5) {
            return this.cutoffR5().getRight();
        }

        return this.metricValue;
    }

    public void fromStateMatrix(double t, SSAStateMatrix networkState) {
    }

    @SuppressWarnings("unchecked")
    public void fromEvent(double t, Event e) {
        if (e instanceof DepartureEvent) {
            if (this.buffer.isEmpty()) {
                return;
            }
            DepartureEvent de = (DepartureEvent) e;
            if ((de.getNode().getStatefulIdx() == this.nodeIdx) &&
                    (de.getJobClass().getJobClassIdx() == this.classIdx)){
                double timePassed = t-this.buffer.pop();
                this.metricValue = ((this.metricValue * this.nDepartures) + timePassed)/(this.nDepartures+1);
                this.nDepartures++;

                if (this.record) {
                    this.metricHistory.add(new Pair(this.time, timePassed));
                } else {
                    if (this.useR5 && (this.crossCount < this.r5Value)) {
                        if (this.belowAverage && ((timePassed) >= this.metricValue)) {
                            this.crossCount++;
                            this.belowAverage = false;

                            if (this.crossCount >= this.r5Value) {
                                this.cutoffTime = t;
                                this.metricValue = 0.0;
                                this.nDepartures = 0;
                            }
                        } else if (!this.belowAverage && ((timePassed) < this.metricValue)) {
                            this.crossCount++;
                            this.belowAverage = true;

                            if (this.crossCount >= this.r5Value) {
                                this.cutoffTime = t;
                                this.metricValue = 0.0;
                                this.nDepartures = 0;
                            }
                        }
                    }
                }
                this.time = t;
            }
        } else if (e instanceof ArrivalEvent) {
            ArrivalEvent ae = (ArrivalEvent) e;
            if ((ae.getNode().getStatefulIdx() == this.nodeIdx) &&
                    (ae.getJobClass().getJobClassIdx() == this.classIdx)){
                if (this.schedStrategy == SchedStrategy.FCFS) {
                    this.buffer.addLast(t);
                } else if (this.schedStrategy == SchedStrategy.LCFS) {
                    this.buffer.addFirst(t);
                } else {
                    this.buffer.addLast(t);
                }
            }
        }
    }

    public void resetHistory() {
        this.metricHistory = new ArrayList<Pair<Double, Double>>(5000);
        this.metricValue = this.initialMetric;
        this.currentValue = this.initialValue;
        this.nDepartures = 0;
        this.cutoffTime = this.time;
    }
    @SuppressWarnings("unchecked")
    public void fromEvent(double t, Event e, int n) {
        if (e instanceof DepartureEvent) {
            if (this.buffer.isEmpty()) {
                return;
            }
            DepartureEvent de = (DepartureEvent) e;
            if ((de.getNode().getStatefulIdx() == this.nodeIdx) &&
                    (de.getJobClass().getJobClassIdx() == this.classIdx)){
                double timePassed = 0;
                for (int i = 0; i < n; i++) {
                    if (buffer.isEmpty()) {
                        n -= i;
                        break;
                    }
                    timePassed += t-buffer.pop();
                }
                this.metricValue = ((this.metricValue * this.nDepartures) + (timePassed*n))/(this.nDepartures+n);
                this.nDepartures += n;

                if (this.record) {
                    this.metricHistory.add(new Pair(this.time, timePassed));
                } else {
                    if (this.useR5 && (this.crossCount < this.r5Value)) {
                        if (this.belowAverage && ((timePassed) >= this.metricValue)) {
                            this.crossCount++;
                            this.belowAverage = false;

                            if (this.crossCount >= this.r5Value) {
                                this.cutoffTime = t;
                                this.metricValue = 0.0;
                                this.nDepartures = 0;
                            }
                        } else if (!this.belowAverage && ((timePassed) < this.metricValue)) {
                            this.crossCount++;
                            this.belowAverage = true;

                            if (this.crossCount >= this.r5Value) {
                                this.cutoffTime = t;
                                this.metricValue = 0.0;
                                this.nDepartures = 0;
                            }
                        }
                    }
                }
                this.time = t;
            }
        } else if (e instanceof ArrivalEvent) {
            ArrivalEvent ae = (ArrivalEvent) e;
            if ((ae.getNode().getStatefulIdx() == this.nodeIdx) &&
                    (ae.getJobClass().getJobClassIdx() == this.classIdx)){
                for (int i = 0; i < n; i++) {
                    if (this.schedStrategy == SchedStrategy.FCFS) {
                        this.buffer.addLast(t);
                    } else if (this.schedStrategy == SchedStrategy.LCFS) {
                        this.buffer.addFirst(t);
                    } else {
                        this.buffer.addLast(t);
                    }
                }
            }
        }
    }

    @Override
    public boolean recognizeEvent(Event e) {
        if (e instanceof DepartureEvent) {
            DepartureEvent de = (DepartureEvent)e;
            return ((de.getNode().getStatefulIdx() == this.nodeIdx) &&
                    (de.getJobClass().getJobClassIdx() == this.classIdx));
        } else if (e instanceof ArrivalEvent) {
            ArrivalEvent ae = (ArrivalEvent) e;
            return ((ae.getNode().getStatefulIdx() == this.nodeIdx) &&
                    (ae.getJobClass().getJobClassIdx() == this.classIdx));
        }
        return true;
    }
}
