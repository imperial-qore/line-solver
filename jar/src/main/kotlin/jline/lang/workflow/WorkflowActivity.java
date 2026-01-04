/*
 * Copyright (c) 2012-2026, QORE Lab, Imperial College London
 * All rights reserved.
 */

package jline.lang.workflow;

import jline.GlobalConstants;
import jline.lang.Element;
import jline.lang.processes.APH;
import jline.lang.processes.Distribution;
import jline.lang.processes.Exp;
import jline.lang.processes.Immediate;
import jline.lang.processes.Markovian;
import jline.util.Pair;
import jline.util.matrix.Matrix;

import java.util.HashMap;
import java.util.Map;

/**
 * A computational activity in a Workflow.
 */
public class WorkflowActivity extends Element {

    private Distribution hostDemand;
    private double hostDemandMean;
    private double hostDemandSCV;
    private Workflow workflow;
    private int index;
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
    }

    public void setHostDemand(Distribution hostDemand) {
        this.hostDemand = hostDemand;
        this.hostDemandMean = hostDemand.getMean();
        this.hostDemandSCV = hostDemand.getSCV();
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
            Matrix T = Matrix.singleton(-1e10);
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
