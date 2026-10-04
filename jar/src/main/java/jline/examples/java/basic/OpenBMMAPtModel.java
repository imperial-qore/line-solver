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
import jline.lang.processes.BMMAPt;
import jline.lang.processes.Exp;
import jline.solvers.ldes.LDESOptions;
import jline.solvers.ldes.SolverLDES;
import jline.util.matrix.Matrix;

/**
 * Open queueing network whose arrival stream is a BATCH, MARKED, time-inhomogeneous MAP.
 *
 * <p>BMMAPt crosses the three axes LINE models for an arrival stream. A block is indexed by
 * segment, mark and batch size: D1kb[c][b][j] holds the rates that, in segment j, release a
 * batch of b+1 jobs ALL carrying mark c+1, and D0[j] plus the sum of every block is a generator
 * in each segment. At a Source the mark selects the CLASS of the arriving jobs and the batch
 * size says HOW MANY of them arrive at that instant.
 *
 * <p>Two derived levels come for free and are what keep the family legible: summing over b gives
 * the MMAPt of the marks alone, and summing over c as well gives the MAPt of the epochs alone. So
 * a consumer that ignores batches sees exactly the marked schedule, and one that ignores marks
 * too sees exactly the unmarked one.
 *
 * <p>The schedule below moves BOTH the class mix and the batch size, which is the point: it is
 * the combination no existing family expresses. Every segment fires epochs at rate 4, so the
 * epoch stream is statistically identical throughout; what changes is which class the batch
 * carries and how big it is. In the morning segment the Premium class arrives in PAIRS and
 * Economy singly; in the evening segment that reverses. The JOB rate therefore differs from the
 * EPOCH rate, and differs per class within a segment even though the epoch rate does not.
 *
 * <p>Reductions worth knowing: with every batch size 1 a BMMAPt IS the MMAPt with the same
 * blocks, sample path for sample path on one seed; with K = 1 it is the unmarked batch schedule;
 * with both it is the MAPt; and with identical segments it is the stationary BMAP.
 *
 * <p>References: Q.-M. He, "The versatility of MMAP[K] and the MMAP[K]/G[K]/1 queue", Queueing
 * Systems 38(4), 2001, for the marked structure; D. M. Lucantoni, "New results on the single
 * server queue with a batch Markovian arrival process", Stochastic Models 7(1), 1991, for the
 * batch structure; Y. M. Ko and J. Pender, "Diffusion limits for the (MAP_t/Ph_t/inf)^N queueing
 * network", Oper. Res. Lett. 45(3), 2017, for the time-inhomogeneous one.
 */
public class OpenBMMAPtModel {

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

    /** Builds a Source -> Queue -> Sink model whose arrivals follow a two-segment BMMAPt. */
    public static Network example() {
        Network model = new Network("model");

        Source source = new Source(model, "Source");
        Queue queue = new Queue(model, "Queue", SchedStrategy.FCFS);
        Sink sink = new Sink(model, "Sink");

        OpenClass premium = new OpenClass(model, "Premium", 0);
        OpenClass economy = new OpenClass(model, "Economy", 0);

        // Two twelve-hour segments repeating on a daily cycle. One phase, so the epoch
        // stream is Poisson at rate 4 throughout and only the labels move.
        double[] breakpoints = new double[]{0.0, 12.0, 24.0};
        List<Matrix> d0 = segments(-4.0, -4.0);

        // D1kb is MARK-MAJOR, then batch, then segment. The batch axis is DENSE: a mark
        // that never releases a batch of that size still declares a zero block for it,
        // exactly as BMAP's {D0, D1, ..., Dk} does.
        //
        //   morning (segment 1): Premium in PAIRS at 3.6, Economy singly at 0.4
        //   evening (segment 2): Premium singly at 0.4,   Economy in PAIRS at 3.6
        List<List<List<Matrix>>> d1kb = new ArrayList<List<List<Matrix>>>();
        List<List<Matrix>> mark1 = new ArrayList<List<Matrix>>();
        mark1.add(segments(0.0, 0.4));   // batch 1
        mark1.add(segments(3.6, 0.0));   // batch 2
        List<List<Matrix>> mark2 = new ArrayList<List<Matrix>>();
        mark2.add(segments(0.4, 0.0));   // batch 1
        mark2.add(segments(0.0, 3.6));   // batch 2
        d1kb.add(mark1);
        d1kb.add(mark2);

        BMMAPt arrival = new BMMAPt(breakpoints, d0, d1kb, true);
        source.setMarkedArrival(arrival, Arrays.asList((JobClass) premium, (JobClass) economy));
        queue.setService(premium, new Exp(16.0));
        queue.setService(economy, new Exp(16.0));

        model.link(Network.serialRouting(source, queue, sink));
        return model;
    }

    public static void main(String[] args) {
        Network model = example();

        // LDES simulates the batch marked schedule directly: the mark and the batch size
        // are both decided by WHICH block's transition fired, inside the same
        // competing-transitions draw that ends the interval, so neither label costs an
        // extra random draw and the reduction to MMAPt at batch size 1 is exact.
        LDESOptions options = new LDESOptions();
        options.samples = 400000;
        options.seed = 23000;
        new SolverLDES(model, options).getAvgTable().print();
    }
}
