/*
 * Copyright (c) 2012-2026, QORE Lab, Imperial College London
 * All rights reserved.
 */

package jline.solvers.fluid.petri;

import jline.lang.constant.TimingStrategy;
import jline.util.matrix.Matrix;

import java.util.ArrayList;
import java.util.List;

/**
 * One transition mode as a reaction record: its input arcs and their weights,
 * its inhibitor arcs and their thresholds, its incidence column, and the firing
 * process that times it. Java twin of the {@code modes} struct array built by
 * the MATLAB {@code fluid_petri_terms}.
 */
public class PetriMode {

    /** Node index of the owning Transition. */
    public int node;
    /** Mode index within that Transition. */
    public int mode;
    public TimingStrategy timing;
    /** Marking coordinates this mode consumes from. */
    public List<Integer> arcSlot = new ArrayList<Integer>();
    /** Multiplicity of each input arc, in the same order. */
    public List<Double> arcW = new ArrayList<Double>();
    /** Marking coordinates that inhibit this mode. */
    public List<Integer> inhSlot = new ArrayList<Integer>();
    /** Threshold of each inhibitor arc, in the same order. */
    public List<Double> inhThr = new ArrayList<Double>();
    /** Server count; infinite for an unbounded mode. */
    public double c = 1.0;
    /** Number of phases of the firing process; 0 for an immediate mode. */
    public int nph = 1;
    public Matrix D0;
    public Matrix D1;
    /** Row sums of D1, i.e. the completion rate out of each phase. */
    public double[] d1 = new double[]{1.0};
    /** Entry distribution over the phases. */
    public double[] pie = new double[]{1.0};
    /** Marking-dependent firing multiplier, null when the mode declares none. */
    public jline.util.SerializableFunction<Matrix, Double> dep;
    public int prio = 1;
    public double weight = 1.0;
    /** Incidence column over the whole state, post minus pre. */
    public double[] cvec;
    /** State coordinates of the phase block, empty for a single-phase mode. */
    public int[] zblk = new int[0];
    /**
     * Whether the min() this mode takes is worth closing. A mode whose enabling
     * degree cannot reach its server count has min() exact on its whole
     * support, so closing it is an error rather than an improvement -- exactly
     * as a station that cannot fill its servers is held first order.
     */
    public boolean closable;
    public String label = "";

    public PetriMode(int node, int mode, TimingStrategy timing, String label, int nm) {
        this.node = node;
        this.mode = mode;
        this.timing = timing;
        this.label = label;
        this.cvec = new double[nm];
    }
}
