/*
 * Copyright (c) 2012-2026, QORE Lab, Imperial College London
 * All rights reserved.
 */

package jline.examples.java.basic;

import java.util.ArrayList;
import java.util.Arrays;
import java.util.List;

import jline.lang.JobClass;
import jline.lang.Network;
import jline.lang.OpenClass;
import jline.lang.constant.SchedStrategy;
import jline.lang.nodes.Queue;
import jline.lang.nodes.Sink;
import jline.lang.nodes.Source;
import jline.lang.processes.Exp;
import jline.lang.processes.MMAPt;
import jline.solvers.ldes.LDESOptions;
import jline.solvers.ldes.SolverLDES;
import jline.util.matrix.Matrix;

/**
 * Open queueing network whose arrival stream is a MARKED, time-inhomogeneous MAP.
 *
 * <p>MMAPt crosses the two axes of MAPt and MarkedMAP: arrivals are labelled with one of K
 * marks, AND the matrices that generate them are functions of the wall clock. Segment j carries
 * D0[j] together with the K blocks D1k[c][j], and D0[j] plus the sum over c of D1k[c][j] is a
 * generator in every segment. The aggregate is the D1 of the underlying MAPt, so hiding the
 * marks recovers exactly that process.
 *
 * <p>At a Source the mark SELECTS THE CLASS of the arriving job: one modulating chain drives
 * every marked class, and only the first one carries the stream. That is what makes this
 * different from declaring K independent arrival processes: the classes are correlated through
 * the shared phase, and the class MIX shifts with the schedule even when the total rate does
 * not.
 *
 * <p>The schedule below is built so that only the mix moves. Both segments carry an aggregate
 * rate of 4, so the total arrival stream is statistically identical throughout; what changes is
 * the split, 9:1 towards the first class in the morning segment and 1:9 towards the second in
 * the evening one. A model that read the marks off the time-averaged matrices would report a
 * flat 2:2 split.
 *
 * <p>Reductions worth knowing: with K = 1 an MMAPt IS the MAPt with the same matrices, and with
 * identical segments it IS the stationary MMAP. MPHt is the phase-type twin, stored lowered to
 * this same form.
 *
 * <p>References: Q.-M. He, "The versatility of MMAP[K] and the MMAP[K]/G[K]/1 queue", Queueing
 * Systems 38(4), 2001, for the marked structure; Y. M. Ko and J. Pender, "Diffusion limits for
 * the (MAP_t/Ph_t/inf)^N queueing network", Oper. Res. Lett. 45(3), 2017, for the
 * time-inhomogeneous one.
 */
public class OpenMMAPtModel {

    private static Matrix one(double v) {
        Matrix out = new Matrix(1, 1);
        out.set(0, 0, v);
        return out;
    }

    private static List<Matrix> segments(double a, double b) {
        List<Matrix> out = new ArrayList<Matrix>();
        out.add(one(a));
        out.add(one(b));
        return out;
    }

    /** Builds a Source -&gt; Queue -&gt; Sink model whose arrivals follow a two-segment MMAPt. */
    public static Network example() {
        Network model = new Network("model");

        Source source = new Source(model, "Source");
        Queue queue = new Queue(model, "Queue", SchedStrategy.FCFS);
        Sink sink = new Sink(model, "Sink");

        OpenClass morningClass = new OpenClass(model, "Morning", 0);
        OpenClass eveningClass = new OpenClass(model, "Evening", 0);

        // Two twelve-hour segments repeating on a daily cycle. One phase, so the
        // aggregate stream is Poisson at rate 4 throughout and ONLY the mix moves.
        double[] breakpoints = new double[]{0.0, 12.0, 24.0};
        List<Matrix> d0 = segments(-4.0, -4.0);

        // D1k is MARK-MAJOR, then segment: mark 1 -> Morning, 9:1 by day and 1:9
        // by night; mark 2 -> Evening, the reverse.
        List<List<Matrix>> d1k = new ArrayList<List<Matrix>>();
        d1k.add(segments(3.6, 0.4));
        d1k.add(segments(0.4, 3.6));

        MMAPt arrival = new MMAPt(breakpoints, d0, d1k, true);
        source.setMarkedArrival(arrival,
                Arrays.asList((JobClass) morningClass, (JobClass) eveningClass));
        queue.setService(morningClass, new Exp(8.0));
        queue.setService(eveningClass, new Exp(8.0));

        model.link(Network.serialRouting(source, queue, sink));
        return model;
    }

    /** The arrival process of {@link #example()}, for the rate identities it satisfies. */
    public static MMAPt arrival() {
        double[] breakpoints = new double[]{0.0, 12.0, 24.0};
        List<List<Matrix>> d1k = new ArrayList<List<Matrix>>();
        d1k.add(segments(3.6, 0.4));
        d1k.add(segments(0.4, 3.6));
        return new MMAPt(breakpoints, segments(-4.0, -4.0), d1k, true);
    }

    public static void main(String[] args) {
        // LDES simulates the marked schedule directly: the mark is decided by WHICH
        // block's transition fired, inside the same competing-transitions draw that
        // ends the interval, so no extra randomness is introduced by the labelling.
        LDESOptions options = new LDESOptions();
        options.samples = 400000;
        options.seed = 23000;
        new SolverLDES(example(), options).getAvgTable().print();
    }
}
