/*
 * Copyright (c) 2012-2026, QORE Lab, Imperial College London
 * All rights reserved.
 */

package jline.lang.workflow;

import jline.GlobalConstants;
import jline.lang.Element;
import jline.lang.processes.APH;
import jline.lang.processes.Distribution;
import jline.lang.processes.DistributionScaling;
import jline.lang.processes.Exp;
import jline.lang.processes.Immediate;
import jline.lang.processes.Markovian;
import jline.util.Pair;
import jline.util.matrix.Matrix;

import java.util.HashMap;
import java.util.Map;

/**
 * A computational activity in a Workflow.
 * <p>
 * Unlike Activity in LayeredNetwork, it carries no call list: an external call
 * is represented as an activity whose host demand is the law of the call
 * response time, so that a synchronous call and a local computation compose in
 * the same way. An asynchronous call blocks the caller for no time and is
 * simply left out of the workflow.
 * </p>
 */
public class WorkflowActivity extends Element {

    private Distribution hostDemand;
    private double hostDemandMean;
    private double hostDemandSCV;
    private Workflow workflow;
    private int index = -1;
    private Map<String, Object> metadata;

    public WorkflowActivity(Workflow workflow, String name, double meanServiceTime) {
        super(name);
        this.workflow = workflow;
        setHostDemand(meanServiceTime);
    }

    public WorkflowActivity(Workflow workflow, String name, Distribution hostDemand) {
        super(name);
        this.workflow = workflow;
        setHostDemand(hostDemand);
    }

    public void setHostDemand(double meanServiceTime) {
        if (meanServiceTime <= GlobalConstants.FineTol) {
            this.hostDemand = new Immediate();
            this.hostDemandMean = GlobalConstants.FineTol;
            this.hostDemandSCV = GlobalConstants.FineTol;
        } else {
            this.hostDemand = Exp.fitMean(meanServiceTime);
            this.hostDemandMean = meanServiceTime;
            this.hostDemandSCV = 1.0;
        }
        invalidateParent();
    }

    public void setHostDemand(Distribution hostDemand) {
        this.hostDemand = hostDemand;
        this.hostDemandMean = hostDemand.getMean();
        this.hostDemandSCV = hostDemand.getSCV();
        invalidateParent();
    }

    /**
     * Change the mean, preserving the shape.
     * <p>
     * Scales the current law in time rather than refitting it, so the SCV, the
     * skewness and the order are preserved and the cached series-parallel tree
     * keeps its shape.
     * </p>
     *
     * @param meanValue new mean, positive and finite
     */
    public void setHostDemandMean(double meanValue) {
        if (!(meanValue > 0) || Double.isInfinite(meanValue)) {
            throw new IllegalArgumentException("The activity mean must be a positive finite scalar.");
        }

        double oldMean = hostDemandMean;
        if (hostDemand == null || hostDemand instanceof Immediate || !(oldMean > 0)
                || Double.isInfinite(oldMean)) {
            setHostDemand(Exp.fitMean(meanValue));
            return;
        }

        double factor = oldMean / meanValue;
        this.hostDemand = DistributionScaling.scaleRate(hostDemand, factor);
        this.hostDemandMean = meanValue;
        // The SCV is invariant under a time scaling
        if (workflow != null && index >= 0) {
            workflow.rescaleActivityLeaf(index, factor);
        }
    }

    private void invalidateParent() {
        // The parent caches the composed law, so the leaf must be marked dirty
        // here as well as on a topology change
        if (workflow != null && index >= 0) {
            workflow.invalidateActivity(index);
        }
    }

    public Distribution getHostDemand() {
        return hostDemand;
    }

    public double getHostDemandMean() {
        return hostDemandMean;
    }

    public double getHostDemandSCV() {
        return hostDemandSCV;
    }

    public Workflow getWorkflow() {
        return workflow;
    }

    public int getIndex() {
        return index;
    }

    public void setIndex(int index) {
        this.index = index;
    }

    public Pair<Matrix, Matrix> getPHRepresentation() {
        if (hostDemand instanceof Immediate) {
            Matrix alpha = Matrix.singleton(1.0);
            Matrix T = Matrix.singleton(-GlobalConstants.Immediate);
            return new Pair<Matrix, Matrix>(alpha, T);
        }

        if (hostDemand instanceof Markovian) {
            Markovian markov = (Markovian) hostDemand;
            Matrix alpha = markov.getInitProb();
            Matrix T = markov.D(0);

            if (alpha.getNumRows() > 1 && alpha.getNumCols() == 1) {
                alpha = alpha.transpose();
            }

            return new Pair<Matrix, Matrix>(alpha, T);
        }

        double mean = hostDemandMean;
        double scv = hostDemandSCV;
        if (scv < GlobalConstants.FineTol) {
            scv = 1.0;
        }
        APH aph = APH.fitMeanAndSCV(mean, scv);
        Matrix alpha = aph.getInitProb();
        Matrix T = aph.D(0);

        if (alpha.getNumRows() > 1 && alpha.getNumCols() == 1) {
            alpha = alpha.transpose();
        }

        return new Pair<Matrix, Matrix>(alpha, T);
    }

    public int getNumberOfPhases() {
        if (hostDemand instanceof Immediate) {
            return 1;
        }
        if (hostDemand instanceof Markovian) {
            return (int) ((Markovian) hostDemand).getNumberOfPhases();
        }
        Pair<Matrix, Matrix> ph = getPHRepresentation();
        return ph.getRight().getNumRows();
    }

    /**
     * Get optional metadata (e.g., from WfCommons).
     * @return Metadata map or null if not set
     */
    public Map<String, Object> getMetadata() {
        return metadata;
    }

    /**
     * Set optional metadata.
     * @param metadata Metadata map
     */
    public void setMetadata(Map<String, Object> metadata) {
        this.metadata = metadata;
    }

    /**
     * Check if metadata is present.
     * @return true if metadata is set
     */
    public boolean hasMetadata() {
        return metadata != null && !metadata.isEmpty();
    }

    /**
     * Get a metadata value by key.
     * @param key Metadata key
     * @return Value or null if not present
     */
    public Object getMetadataValue(String key) {
        if (metadata == null) {
            return null;
        }
        return metadata.get(key);
    }

    /**
     * Set a metadata value.
     * @param key Metadata key
     * @param value Metadata value
     */
    public void setMetadataValue(String key, Object value) {
        if (metadata == null) {
            metadata = new HashMap<String, Object>();
        }
        metadata.put(key, value);
    }
}
