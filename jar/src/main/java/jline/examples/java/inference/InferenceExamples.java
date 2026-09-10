/*
 * Copyright (c) 2012-2026, QORE Lab, Imperial College London
 * All rights reserved.
 */
package jline.examples.java.inference;

import java.util.Arrays;

import jline.inference.lang.ParamEstimator;
import jline.inference.lang.SampledMetric;
import jline.lang.ClosedClass;
import jline.lang.Network;
import jline.lang.constant.MetricType;
import jline.lang.constant.SchedStrategy;
import jline.lang.nodes.Delay;
import jline.lang.nodes.Queue;
import jline.lang.nodes.Station;
import jline.lang.processes.Exp;
import jline.solvers.mva.SolverMVA;
import jline.util.matrix.Matrix;

/**
 * Java twin of the inference example scripts (matlab/examples/inference,
 * python/examples/inference, cpp/examples/inference).
 */
public class InferenceExamples {

    /**
     * Variational inference for Markovian queueing networks, on a closed loop.
     *
     * Method "vi" (I. Perez, G. Casale, Adv. Appl. Prob. 53(3), 2021) infers
     * service rates from NOISY QUEUE-LENGTH READINGS taken over time: each
     * reading is exact with probability 1-epsilon and uniform over the remaining
     * feasible values otherwise. Unlike the other estimators it returns a
     * conjugate Gamma POSTERIOR per rate, not only a point estimate.
     *
     * It reads the queue lengths of EVERY station, not only the estimated one:
     * the transition counts the method is written in are pinned by the whole
     * picture.
     */
    public static void est_vi_closed() {
        int N = 10;

        // define model, with the queue rate to be estimated
        Network model = new Network("model");
        Delay delay = new Delay(model, "Delay");
        Queue queue = new Queue(model, "Queue1", SchedStrategy.FCFS);
        ClosedClass jobclass = new ClosedClass(model, "Class1", N, delay, 0);
        delay.setService(jobclass, new Exp(0.5));
        queue.setService(jobclass, new Exp(2.0));   // starting point of the estimate
        model.link(model.serialRouting(delay, queue));

        // queue-length readings, one per unit time, 10% of them faulty
        double[] qlen = {1, 2, 3, 2, 4, 3, 5, 4, 3, 4};
        double[] ts = new double[qlen.length];
        double[] dlen = new double[qlen.length];
        for (int k = 0; k < qlen.length; k++) {
            ts[k] = k + 1;
            dlen[k] = N - qlen[k];
        }

        // estimate the service rate of the queue
        ParamEstimator se = new ParamEstimator(model);
        se.options.method = "vi";
        se.options.epsilon = 0.1;      // probability that a reading is faulty
        se.options.priorShape = 2.0;   // Gamma prior shape; the rate is set from the model
        se.options.variational.ngrid = 51;      // time grid of the backward and forward passes
        se.options.variational.nsamples = 32;   // lattice points per marginal
        se.options.variational.ymax = 60;       // transition-count truncation
        se.options.variational.iterMax = 5;
        se.addSamples(new SampledMetric(MetricType.QLen, ts, dlen, delay, jobclass));
        se.addSamples(new SampledMetric(MetricType.QLen, ts, qlen, queue, jobclass));
        Matrix estVal = se.estimateAt(Arrays.asList((Station) queue));

        System.out.printf("Estimated demand: %.8f%n", estVal.get(0, 0));
        System.out.printf("posterior service rate ~ Gamma(%.4f, %.4f), mean %.4f%n",
                se.options.posteriorAlpha[0], se.options.posteriorBeta[0],
                se.options.posteriorAlpha[0] / se.options.posteriorBeta[0]);
        System.out.print("evidence lower bound over the iterations:");
        for (int i = 0; i < se.options.bound.length; i++) {
            System.out.printf(" %.3f", se.options.bound[i]);
        }
        System.out.println();

        // solve the model the estimate has been written into
        System.out.println("SOLVER: MVA");
        new SolverMVA(model).getAvgTable().print();
    }

    public static void main(String[] args) {
        est_vi_closed();
    }
}
