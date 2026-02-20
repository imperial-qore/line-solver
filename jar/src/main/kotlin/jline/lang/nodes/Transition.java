/*
 * Copyright (c) 2012-2026, QORE Lab, Imperial College London
 * All rights reserved.
 */

package jline.lang.nodes;

import jline.lang.JobClass;
import jline.lang.Mode;
import jline.lang.Network;
import jline.lang.constant.TimingStrategy;
import jline.lang.processes.Distribution;
import jline.lang.processes.Exp;
import jline.lang.processes.Immediate;
import jline.lang.sections.Enabling;
import jline.lang.sections.Firing;
import jline.lang.sections.Timing;
import jline.util.matrix.Matrix;

import java.util.ArrayList;
import java.util.HashMap;
import java.util.List;
import java.util.Map;

/**
 * Transition as in a stochastic Petri net model
 */
public class Transition extends ServiceNode {
    public Map<Mode, Matrix> enablingConditions;
    public Map<Mode, Matrix> inhibitingConditions;
    public Map<Mode, Matrix> firingOutcomes;
    public Map<Mode, TimingStrategy> timingStrategies;
    public Map<Mode, Distribution> distributions;
    public Matrix firingPriorities;
    public Matrix firingWeights;
    List<Mode> modes;
    Map<Mode, String> modeNames;
    Map<Mode, Integer> numberOfServers;
    int cap;

    public Transition(Network model, String nodeName) {
        super(nodeName);

        this.setModel(model);
        this.model.addNode(this);

        this.input = new Enabling();
        this.server = new Timing();
        this.output = new Firing(model.getClasses());

        this.cap = Integer.MAX_VALUE;
        this.enablingConditions = new HashMap<>();
        this.inhibitingConditions = new HashMap<>();
        this.modes = new ArrayList<>();
        this.modeNames = new HashMap<>();
        this.distributions = new HashMap<>();
        this.firingOutcomes = new HashMap<>();
        this.numberOfServers = new HashMap<>();
        this.timingStrategies = new HashMap<>();
        this.firingPriorities = new Matrix(0, 0);
        this.firingWeights = new Matrix(0, 0);
    }

    public Mode addMode(String modename) {
        Mode newmode = new Mode(this, modename);
        addMode(newmode);
        return newmode;
    }

    public void addMode(Mode mode) {
        int nclasses = this.model.getNumberOfClasses();
        int nnodes = this.model.getNumberOfNodes();
        this.modes.add(mode);
        this.modeNames.put(mode, mode.getName());
        this.enablingConditions.put(mode, new Matrix(nnodes, nclasses));
        for (Mode m : modes) {
            Matrix E = this.enablingConditions.get(mode);
            if (E.getNumRows() != nnodes) {
                Matrix.concatRows(E, Matrix.zeros(nnodes - E.getNumRows(), E.getNumCols()), E);
            }
            if (E.getNumCols() != nclasses) {
                Matrix.concatColumns(E, Matrix.zeros(E.getNumRows(), nclasses - E.getNumCols()), E);
            }
        }
        this.inhibitingConditions.put(mode, new Matrix(nnodes, nclasses));
        for (Mode m : modes) {
            Matrix I = this.inhibitingConditions.get(mode);
            if (I.getNumRows() != nnodes) {
                Matrix.concatRows(I, Matrix.zeros(nnodes - I.getNumRows(), I.getNumCols()), I);
            }
            if (I.getNumCols() != nclasses) {
                Matrix.concatColumns(I, Matrix.zeros(I.getNumRows(), nclasses - I.getNumCols()), I);
            }
        }
        this.numberOfServers.put(mode, 1); // Initialise transition with 1 server
        this.timingStrategies.put(mode, TimingStrategy.TIMED);
        if (this.firingWeights.isEmpty()) {
            this.firingWeights = Matrix.singleton(1.0);
        } else {
            this.firingWeights = this.firingWeights.concatCols(Matrix.singleton(1.0));
        }
        if (this.firingPriorities.isEmpty()) {
            this.firingPriorities = Matrix.singleton(1.0);
        } else {
            this.firingPriorities = this.firingPriorities.concatCols(Matrix.singleton(1.0));
        }
        this.distributions.put(mode, new Exp(1));
        this.firingOutcomes.put(mode, new Matrix(nnodes, nclasses));
    }

    public Distribution getFiringDistribution(Mode m) {
        return distributions.get(m);
    }

    public List<String> getModeNames() {
        return new ArrayList<>(modeNames.values());
    }

    public List<Mode> getModes() {
        return modes;
    }

    public Matrix getNumberOfModeServers() {
        Matrix modeServers = new Matrix(1, getNumberOfModes(), getNumberOfModes());
        int mctr = 0;
        for (Mode m : modes) {
            modeServers.set(0, mctr, numberOfServers.get(m));
            mctr++;
        }
        return modeServers;
    }

    public int getNumberOfModeServers(Mode m) {
        return numberOfServers.get(m);
    }

    public int getNumberOfModes() {
        return modes.size();
    }

    public void setDistribution(Mode mode, Distribution distribution) {
        this.distributions.put(mode, distribution);
    }

    public void setEnablingConditions(Mode mode, JobClass jobclass, Place inputPlace, int enablingCondition) {
        enablingConditions.get(mode).set(inputPlace.getNodeIndex(), jobclass.getIndex() - 1, enablingCondition);
    }

    public void setFiringOutcome(Mode mode, JobClass jobclass, Node node, int firingOutcome) {
        this.firingOutcomes.get(mode).set(node.getNodeIndex(), jobclass.getIndex() - 1, firingOutcome);
    }

    public void setFiringPriorities(Mode mode, int firingPriority) {
        this.firingPriorities.set(0, modes.indexOf(mode), firingPriority);
    }

    public void setFiringWeights(Mode mode, double firingWeight) {
        this.firingWeights.set(modes.indexOf(mode), firingWeight);
    }

    public void setInhibitingConditions(Mode mode, JobClass jobclass, Place inputPlace, int inhibitingCondition) {
        inhibitingConditions.get(mode).set(inputPlace.getNodeIndex(), jobclass.getIndex() - 1, inhibitingCondition);
    }

    public void setModeNames(Mode mode, String modeName) {
        this.modeNames.put(mode, modeName);
    }

    public void setNumberOfServers(Mode mode, Integer numberOfServers) {
        this.numberOfServers.put(mode, numberOfServers);
    }

    public void setTimingStrategy(Mode mode, TimingStrategy timingStrategy) {
        this.timingStrategies.put(mode, timingStrategy);
        if (timingStrategy == TimingStrategy.IMMEDIATE) {
            this.distributions.put(mode, Immediate.getInstance());
        }
    }

}
