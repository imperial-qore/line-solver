/*
 * Copyright (c) 2012-2026, QORE Lab, Imperial College London
 * All rights reserved.
 */

package jline.lang.nodes;

import jline.lang.ClassSwitchMatrix;
import jline.lang.JobClass;
import jline.lang.Network;
import jline.lang.constant.RoutingStrategy;
import jline.lang.constant.SchedStrategy;
import jline.lang.constant.SchedStrategyType;
import jline.lang.sections.Buffer;
import jline.lang.sections.Dispatcher;
import jline.lang.sections.StatelessClassSwitcher;
import jline.util.matrix.Matrix;

import java.io.Serializable;
import java.util.List;

/**
 * A node that switches the class of an incoming job based on a probability table
 */
public class ClassSwitch extends Node implements Serializable {
    public boolean autoAdded;
    protected SchedStrategyType schedPolicy;
    protected SchedStrategy schedStrategy;

    public ClassSwitch(Network model, String name) {
        super(name);

        List<JobClass> jobClasses = model.getClasses();
        int m = jobClasses.size();

        this.input = new Buffer(jobClasses);
        this.output = new Dispatcher(jobClasses);

        Matrix csFun = new Matrix(m, m, m * m);

        this.setModel(model);
        this.model.addNode(this);

        this.schedPolicy = SchedStrategyType.NP;
        this.schedStrategy = SchedStrategy.FCFS;
        this.server = new StatelessClassSwitcher(jobClasses, csFun);
        this.autoAdded = false;
    }

    public ClassSwitch(Network model, String name, Matrix csFun) {
        super(name);

        List<JobClass> jobClasses = model.getClasses();
        this.input = new Buffer(jobClasses);
        this.output = new Dispatcher(jobClasses);


        this.setModel(model);
        this.model.addNode(this);

        this.schedPolicy = SchedStrategyType.NP;
        this.schedStrategy = SchedStrategy.FCFS;
        this.server = new StatelessClassSwitcher(jobClasses, csFun);
    }

    public ClassSwitchMatrix initClassSwitchMatrix() {
        int m = this.model.getNumberOfClasses();
        // Find the maximum class index to ensure matrix can accommodate 1-based indexing
        int maxIndex = 0;
        for (JobClass jc : this.model.getClasses()) {
            maxIndex = Math.max(maxIndex, jc.getIndex());
        }
        // Matrix size needs to be at least maxIndex+1 to accommodate 0-based access with 1-based indices
        int matrixSize = Math.max(m, maxIndex + 1);
        return new ClassSwitchMatrix(matrixSize, matrixSize, matrixSize * matrixSize);
    }

    public void setClassSwitchingMatrix(ClassSwitchMatrix csMatrix) {
        ((StatelessClassSwitcher) this.server).updateClasses(this.model.getClasses());
        ((StatelessClassSwitcher) this.server).updateClassSwitch(csMatrix.getMatrix());
    }

    public void setProbRouting(JobClass jobClass, Node destination, double probability) {
        this.setRouting(jobClass, RoutingStrategy.PROB, destination, probability);
    }

}
