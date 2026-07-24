/*
 * Copyright (c) 2012-2026, QORE Lab, Imperial College London
 * All rights reserved.
 */

package jline.lang.sections;

import jline.lang.JobClass;
import jline.util.SerializableFunction;
import jline.util.matrix.Matrix;

import java.io.Serializable;
import java.util.List;

/**
 * A class switcher that does not have a local state
 */
public class StatelessClassSwitcher extends ClassSwitcher implements Serializable {
    public StatelessClassSwitcher(List<JobClass> jobClasses, Matrix csMatrix) {
        super(jobClasses, "StatelessClassSwitcher");

        // Use anonymous inner class instead of lambda to avoid BootstrapMethodError
        // with CSFunInput inner class reference in Java 8/MATLAB classloader
        this.csFun = new SerializableFunction<CSFunInput, Double>() {
            @Override
            public Double apply(CSFunInput input) {
                int row = input.r;
                int col = input.s;
                return csMatrix.get(row, col);
            }
        };
    }

    public List<JobClass> getJobClasses() {
        return this.jobClasses;
    }

    public void updateClassSwitch(Matrix csMatrix) {
        // Use anonymous inner class instead of lambda to avoid BootstrapMethodError
        this.csFun = new SerializableFunction<CSFunInput, Double>() {
            @Override
            public Double apply(CSFunInput input) {
                int row = input.r;
                int col = input.s;
                return csMatrix.get(row, col);
            }
        };
    }

    public void updateClasses(List<JobClass> jobClasses) {
        this.jobClasses = jobClasses;
    }
}
