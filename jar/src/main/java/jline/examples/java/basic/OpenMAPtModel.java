/*
 * Copyright (c) 2012-2026, QORE Lab, Imperial College London
 * All rights reserved.
 */

package jline.examples.java.basic;

import java.util.ArrayList;
import java.util.List;

import jline.lang.Network;
import jline.lang.OpenClass;
import jline.lang.constant.SolverType;
import jline.lang.nodes.Delay;
import jline.lang.nodes.Sink;
import jline.lang.nodes.Source;
import jline.lang.processes.Exp;
import jline.lang.processes.MAPt;
import jline.solvers.SolverOptions;
import jline.solvers.fluid.SolverFluid;
import jline.util.matrix.Matrix;

/**
 * Open queueing network with a MAP_t arrival process and the Ko-Pender limits.
 *
 * <p>MAPt is a time-inhomogeneous Markovian arrival process: segment k covers
 * [breakpoints[k], breakpoints[k+1]) and carries the pair (D0[k], D1[k]), so the stream is both
 * non-renewal, through the modulating phase, and non-stationary, through the schedule. Setting
 * one phase recovers an NHPP; setting one segment recovers an ordinary MAP.
 *
 * <p>The fluid solver's "kp" method integrates the fluid and diffusion limits of Ko and Pender,
 * "Diffusion limits for the (MAP_t/Ph_t/inf)^N queueing network", Oper. Res. Lett. 45 (2017)
 * 248-253: the mean and the covariance of the queue length are integrated jointly. For
 * infinite-server stations the rate functions are affine in the state, so both are exact rather
 * than asymptotic.
 */
public class OpenMAPtModel {

    private static Matrix matrix(double[][] entries) {
        Matrix out = new Matrix(entries.length, entries[0].length);
        for (int i = 0; i < entries.length; i++) {
            for (int j = 0; j < entries[0].length; j++) {
                out.set(i, j, entries[i][j]);
            }
        }
        return out;
    }

    /**
     * Builds a Source -> Delay -> Sink model whose arrivals follow a two-segment 2-phase MAP_t.
     */
    public static Network example() {
        Network model = new Network("model");

        Source source = new Source(model, "Source");
        Delay delay = new Delay(model, "Delay");
        Sink sink = new Sink(model, "Sink");

        OpenClass jobclass = new OpenClass(model, "OpenClass", 0);

        // Two segments of a 2-phase MAP, held for 1 and 1.5 time units and repeating.
        // The second segment runs the same phase graph at roughly twice the rate.
        List<Matrix> d0 = new ArrayList<Matrix>();
        List<Matrix> d1 = new ArrayList<Matrix>();
        d0.add(matrix(new double[][]{{-5.0, 1.0}, {2.0, -4.0}}));
        d0.add(matrix(new double[][]{{-12.0, 3.0}, {5.0, -9.0}}));
        d1.add(matrix(new double[][]{{3.0, 1.0}, {1.0, 1.0}}));
        d1.add(matrix(new double[][]{{7.0, 2.0}, {2.0, 2.0}}));

        source.setArrival(jobclass, new MAPt(new double[]{0.0, 1.0, 2.5}, d0, d1, true));
        delay.setService(jobclass, new Exp(2.0));

        model.link(Network.serialRouting(source, delay, sink));
        return model;
    }

    public static void main(String[] args) {
        Network model = example();

        SolverOptions options = new SolverOptions(SolverType.FLUID);
        options.method = "kp";
        options.tol = 1e-9;

        SolverFluid solver = new SolverFluid(model, options);
        // Steady state of a cyclic schedule is the average over one period.
        solver.getAvgTable();
    }
}
