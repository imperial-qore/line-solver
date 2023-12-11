package jline.solvers.ssa.metrics;

import java.util.*;

import jline.solvers.ssa.events.Event;
import jline.solvers.ssa.state.SSAStateMatrix;
import jline.util.Pair;

public abstract class Metric<T extends Number, U> {
    protected double time;
    protected T metricValue;
    protected T currentValue;
    protected String name;
    protected String shortName;
    protected List<Pair<Double, T>> metricHistory;
    protected boolean record;
    protected int nodeIdx;
    protected int classIdx;
    protected int nServers;
    public boolean useMatrix;

    protected int r5Value;
    protected boolean useR5;
    protected boolean useMSER5;
    protected boolean belowAverage;
    protected int crossCount;
    protected double cutoffTime;
    protected T initialMetric;
    protected T initialValue;

    public Metric(String name, T initialMetric, T initialValue, boolean record, int nodeIdx, int classIdx, int nServers) {
        this.name = name;
        this.shortName = Character.toString(this.name.toUpperCase().charAt(0));
        this.metricValue = initialMetric;
        this.initialMetric = initialMetric;
        this.currentValue = initialValue;
        this.initialValue = initialValue;
        this.record = record;
        this.useMatrix = false;

        if (this.record) {
            this.metricHistory = new ArrayList<Pair<Double, T>>(5000);
        }

        this.r5Value = 19;
        this.useR5 = false;
        this.useMSER5 = false;
        this.belowAverage = false;
        this.crossCount = 0;

        this.nodeIdx = nodeIdx;
        this.classIdx = classIdx;
        this.nServers = nServers;
        this.cutoffTime = 0;
    }

    public Metric(String name, T initialMetric, T initialValue, int nodeIdx, int classIdx, int nServers) {
        this(name, initialMetric, initialValue, true, nodeIdx, classIdx, nServers);
    }

    public Metric(String name) {
        this.name = name;
        this.shortName = Character.toString(this.name.toUpperCase().charAt(0));
    }

    public void configureR5(int k) {
        this.useR5 = true;
        this.useMSER5 = false;
        this.r5Value = k;
    }

    public void configureMSER5() {
        if (!this.record) {
            throw new RuntimeException("MSER5 requires a timeline!");
        }
        this.useMSER5 = true;
        this.useR5 = false;
    }

    protected abstract void addSample(double currentTime, T metric);
    public abstract void fromStateMatrix(double t, SSAStateMatrix networkState);
    public abstract void fromEvent(double t, Event e);
    public abstract void fromEvent(double t, Event e, int n);
    public abstract T getMetric();

    public String getName() {
        return this.name;
    }
    public String getShortName() {
        return this.shortName;
    }

    public void setShortName(String shortName) {
        this.shortName = shortName;
    }

    public void taper(double t) {
        this.time = t;
    }

    public void resetHistory() {
        this.metricHistory = new ArrayList<Pair<Double, T>>(5000);
        this.metricValue = this.initialMetric;
        this.currentValue = this.initialValue;
        this.cutoffTime = this.time;
    }

    public Pair<Integer, Double> cutoffR5() {
        if (!this.record) {
            this.useR5 = false;
            Pair<Integer,Double> outPair =  new Pair<Integer, Double>(0,this.getMetric().doubleValue());
            this.useR5 = true;
            return outPair;
        }
        Iterator<Pair<Double, T>> pairIterator = metricHistory.listIterator();
        double movingAverage = 0;
        double curTime = 0;
        boolean isTransient = false;
        double cumulativeMT = 0;
        double cutoffTime = 0;
        int cutoffPoint = 0;
        int nCrosses = 0;

        double prevValue = 0;
        double val = 0;
        while(pairIterator.hasNext()) {
            Pair<Double, T> curPair = pairIterator.next();

            if (isTransient) {
                cutoffPoint++;
                boolean below = val < movingAverage;

                val = curPair.getLeft();

                movingAverage *= curTime;
                movingAverage += val;
                movingAverage /= curPair.getLeft();

                if (below && (val >= movingAverage)) {
                    nCrosses++;
                } else if (!below && (val < movingAverage)) {
                    nCrosses++;
                }

                if (nCrosses >= this.r5Value) {
                    isTransient = false;
                    cutoffTime = curPair.getLeft();
                    prevValue = val;
                }
            } else {
                cumulativeMT += prevValue*(curPair.getLeft()-curTime);
                prevValue = curPair.getRight().doubleValue();
            }

            curTime = curPair.getLeft();
        }

        return new Pair<Integer, Double>(cutoffPoint, cumulativeMT/(curTime-cutoffTime));
    }

    public Pair<Integer, Double> cutoffMSER5() {
        if (!this.record) {
            return new Pair<Integer, Double>(0, 0.0);
        }

        Iterator<Pair<Double, T>> pairIterator = metricHistory.listIterator();

        this.useMSER5 = false; // prevent infinite regress..
        double metricMean = this.getMetric().doubleValue();
        this.useMSER5 = true;

        double currentD = 0;
        double currentBatchValue = 0;
        int batchCtr = 0;
        double batchT = 0;
        double minD = Double.POSITIVE_INFINITY;

        int n = metricHistory.size();
        double curT = 0;

        while (pairIterator.hasNext()) {
            Pair<Double, T> candidateCutoff = pairIterator.next();
            double prevT = curT;
            curT = candidateCutoff.getLeft();
            currentBatchValue += candidateCutoff.getRight().doubleValue()*(curT-prevT);
            batchCtr++;

            if ((batchCtr == 4) || (!pairIterator.hasNext())) {
                currentD += Math.pow((currentBatchValue/(curT-batchT))-metricMean,2);

                currentBatchValue = 0;
                batchCtr = 0;
                batchT = curT;
            }
        }

        currentD /= Math.pow(n,2);
        int iAtMinD = 0;

        int d = 0;
        batchCtr = 0;
        currentBatchValue = 0;
        pairIterator = metricHistory.listIterator();
        curT = 0;

        while (pairIterator.hasNext()) {
            Pair<Double, T> candidateCutoff = pairIterator.next();
            double prevT = curT;
            curT = candidateCutoff.getLeft();
            currentBatchValue += candidateCutoff.getRight().doubleValue() * (curT-prevT);
            batchCtr++;

            if ((batchCtr == 4) || (!pairIterator.hasNext())) {
                currentD *= Math.pow(n-d,2);
                d += 5;
                currentD -= Math.pow((currentBatchValue/(curT-batchT))-metricMean, 2);
                currentD /= Math.pow(n-d,2);
                if (currentD < minD) {
                    minD = currentD;
                    iAtMinD = d;
                }

                currentBatchValue = 0;
                batchCtr = 0;
                batchT = curT;
            }
        }

        if (iAtMinD > (n/2)) {
            iAtMinD = n/2;
        }

        double steadyStateMean = 0;
        double startTime = 0;
        double prevTime = 0;
        double curTime = 0;
        pairIterator = metricHistory.listIterator();
        for (int i = 0; i < iAtMinD; i++) {
            curTime = pairIterator.next().getLeft();
        }

        double prevValue;
        Pair<Double, T> curPair;

        if (pairIterator.hasNext()) {
            curPair = pairIterator.next();
            curTime = curPair.getLeft();
            prevValue = curPair.getRight().doubleValue();
        } else {
            return new Pair <Integer, Double>(0, 0.0);
        }

        startTime = curTime;
        double cumulativeMT = 0;
        prevTime = 0;
        double timeIter;
        prevValue = 0;

        while (pairIterator.hasNext()) {
            Pair<Double, T> timelineEntry = pairIterator.next();
            timeIter = timelineEntry.getLeft();
            cumulativeMT += (timeIter-prevTime)*prevValue;
            prevTime = timeIter;
            prevValue = timelineEntry.getRight().doubleValue();
        }

        return new Pair<Integer, Double>(iAtMinD, cumulativeMT/(prevTime-startTime));
    }

    public boolean recognizeEvent(Event e) {
        return true;
    }

    public void setRecord(boolean record) {
        this.record = record;
    }
}
