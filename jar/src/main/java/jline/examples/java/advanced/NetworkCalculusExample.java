/*
 * Copyright (c) 2012-2026, QORE Lab, Imperial College London
 * All rights reserved.
 */
package jline.examples.java.advanced;

import jline.api.snc.Snc_bound_delay;
import jline.api.snc.Snc_conv;
import jline.api.snc.Snc_env_poisson;
import jline.api.snc.Snc_perc_delay;
import jline.api.snc.Snc_srv_exp;
import jline.api.snc.SncEnvelope;
import jline.api.snc.SncResult;
import jline.lang.Network;
import jline.lang.OpenClass;
import jline.lang.RoutingMatrix;
import jline.lang.constant.SchedStrategy;
import jline.lang.nodes.Queue;
import jline.lang.nodes.Sink;
import jline.lang.nodes.Source;
import jline.lang.processes.Exp;
import jline.solvers.NetworkAvgTable;
import jline.solvers.NetworkPercTable;
import jline.solvers.ba.SolverBA;
import jline.util.matrix.Matrix;

/**
 * Stochastic network calculus: a delay quantile with a certified violation
 * probability, and what envelope propagation costs across a network.
 *
 * <p>Every other solver in LINE answers with a MEAN. The 'snc' family of
 * {@link SolverBA} answers with a TAIL: given a violation probability eps it
 * returns a delay d for which P{D &gt; d} &lt;= eps holds, and the guarantee is
 * valid for any work-conserving scheduling policy at the station. That is the
 * quantity a service-level objective is written against.</p>
 *
 * <p>The models are an M/M/1 and a tandem of M/M/1 stations, whose exact
 * answers are known in closed form, so every number printed can be checked.
 * Twin of matlab/examples/advanced/networkCalculus/snc_delay_quantile.m and
 * snc_tandem_multiclass.m.</p>
 */
public class NetworkCalculusExample {

    /** The quantile, how it tightens deep in the tail, and the loose mean. */
    public static void delayQuantile() {
        final double lambda = 0.6;
        final double mu = 1.0;
        Network model = new Network("SncQuantile");
        Source source = new Source(model, "Source");
        Queue queue = new Queue(model, "Queue", SchedStrategy.FCFS);
        Sink sink = new Sink(model, "Sink");
        OpenClass jobclass = new OpenClass(model, "Class1");
        source.setArrival(jobclass, new Exp(lambda));
        queue.setService(jobclass, new Exp(mu));
        RoutingMatrix P = model.initRoutingMatrix();
        P.set(jobclass, jobclass, Network.serialRouting(source, queue, sink));
        model.link(P);

        SolverBA solver = new SolverBA(model, "snc.upper");

        // getPercTable is to the snc family what getAvgTable is to a mean
        // solver: one row per station and class, reporting the response-time
        // and queue-length quantiles at the requested violation probability.
        NetworkPercTable percTable = solver.getPercTable(1e-3);
        percTable.print();

        // The exact M/M/1 sojourn tail is exp(-(mu-lambda)*d) and the exact
        // queue-length tail is rho^(n+1). The bound reproduces both DECAY RATES
        // exactly and pays a constant prefactor, so the ratio of the bounded
        // quantile to the exact one falls towards 1 as eps is tightened: the
        // family is at its best exactly where simulation is at its worst.
        System.out.printf("%n%-8s %10s %10s %8s   %10s %10s %8s%n",
                "eps", "d bound", "d exact", "ratio", "n bound", "n exact", "ratio");
        double[] epsList = {1e-2, 1e-3, 1e-6, 1e-9, 1e-12};
        for (int i = 0; i < epsList.length; i++) {
            Matrix d = solver.getDelayPerc(epsList[i]);
            Matrix b = solver.getBacklogPerc(epsList[i]);
            double dexact = -Math.log(epsList[i]) / (mu - lambda);
            double nexact = Math.log(epsList[i]) / Math.log(lambda / mu) - 1.0;
            System.out.printf("%-8.0e %10.4f %10.4f %8.3f   %10.4f %10.4f %8.3f%n",
                    epsList[i], d.get(1, 0), dexact, d.get(1, 0) / dexact,
                    b.get(1, 0), nexact, b.get(1, 0) / nexact);
        }

        // getAvgTable still works: the response time is the integral of the
        // tail bound, hence an upper bound on the mean. It is loose, and
        // deliberately so -- integrating over the whole axis is dominated by
        // the prefactor rather than by the decay rate the family gets right.
        NetworkAvgTable avgTable = solver.getAvgTable();
        avgTable.print();
        double exactR = 1.0 / (mu - lambda);
        double exactQ = (lambda / mu) / (1.0 - lambda / mu);
        double r = avgTable.getRespT().get(avgTable.getRespT().size() - 1);
        double q = avgTable.getQLen().get(avgTable.getQLen().size() - 1);
        System.out.printf("%nexact M/M/1: R = %.4f, Q = %.4f%n", exactR, exactQ);
        System.out.printf("snc.upper  : R = %.4f (%.1fx), Q = %.4f (%.1fx)%n",
                r, r / exactR, q, q / exactQ);

        // The solver is a thin wrapper over jline.api.snc. Note Snc_srv_exp,
        // not Snc_srv_rate: the work unit here is the JOB, so the server is the
        // counting process of an Exp(mu) service, and a constant-rate element
        // would model an M/D/1 and understate the delay.
        SncEnvelope arv = Snc_env_poisson.of(lambda);
        SncEnvelope srv = Snc_srv_exp.of(mu);
        SncResult d = Snc_perc_delay.snc_perc_delay(arv, srv, 1e-3);
        System.out.printf("%napi: d(1e-3) = %.4f at theta = %.4f%n", d.value, d.theta);
        System.out.printf("     the optimal theta approaches log(mu/lambda) = %.4f, which is%n",
                Math.log(mu / lambda));
        System.out.printf("     what makes the backlog decay rate exact%n");
        SncResult back = Snc_bound_delay.snc_bound_delay(arv, srv, d.value);
        System.out.printf("api: P{D > %.4f} <= %.3e (theta = %.4f), exact tail %.3e%n",
                d.value, back.value, back.theta, Math.exp(-(mu - lambda) * d.value));
    }

    /** Envelope propagation, pay-bursts-only-once, and blind multiplexing. */
    public static void tandemAndSharing() {
        final double lambda = 0.6;
        final double[] rates = {1.5, 1.2, 1.0};
        Network model = new Network("SncTandem");
        Source source = new Source(model, "Source");
        Queue q1 = new Queue(model, "Q1", SchedStrategy.FCFS);
        Queue q2 = new Queue(model, "Q2", SchedStrategy.FCFS);
        Queue q3 = new Queue(model, "Q3", SchedStrategy.FCFS);
        Sink sink = new Sink(model, "Sink");
        OpenClass jobclass = new OpenClass(model, "Class1");
        source.setArrival(jobclass, new Exp(lambda));
        q1.setService(jobclass, new Exp(rates[0]));
        q2.setService(jobclass, new Exp(rates[1]));
        q3.setService(jobclass, new Exp(rates[2]));
        RoutingMatrix P = model.initRoutingMatrix();
        P.set(jobclass, jobclass, Network.serialRouting(source, q1, q2, q3, sink));
        model.link(P);

        NetworkAvgTable tandemTable = new SolverBA(model, "snc.upper").getAvgTable();
        tandemTable.print();

        // Each hop replaces the arrival envelope by the DEPARTURE envelope of
        // the station upstream, which carries the burst the server has added.
        // The exact answer does not degrade this way -- by Burke's theorem the
        // departure process of an M/M/1 is again Poisson -- so the ratio to the
        // exact response time grows hop by hop. The bound stays valid; it is the
        // price of assuming nothing about the departure process.
        System.out.printf("%n%-8s %10s %10s %8s%n", "station", "R bound", "R exact", "ratio");
        for (int i = 0; i < 3; i++) {
            double exact = 1.0 / (rates[i] - lambda);
            double bound = tandemTable.getRespT().get(i + 1);
            System.out.printf("%-8s %10.4f %10.4f %8.2f%n",
                    tandemTable.getStationNames().get(i + 1), bound, exact, bound / exact);
        }

        // Summing the per-station bounds pays the burst term at every hop.
        // Concatenating the three service envelopes with Snc_conv first and
        // bounding the composed element once pays it only once, which is the
        // classical result of the network calculus.
        final SncEnvelope arv = Snc_env_poisson.of(lambda);
        SncEnvelope endToEnd = new SncEnvelope() {
            public double[] eval(double theta) {
                double[] acc = Snc_srv_exp.snc_srv_exp(rates[0], theta);
                for (int i = 1; i < rates.length; i++) {
                    double[] next = Snc_srv_exp.snc_srv_exp(rates[i], theta);
                    acc = Snc_conv.snc_conv(acc[0], acc[1], next[0], next[1], theta);
                }
                return acc;
            }
        };
        SncResult concat = Snc_perc_delay.snc_perc_delay(arv, endToEnd, 1e-3);
        double hopByHop = 0;
        for (int i = 0; i < rates.length; i++) {
            hopByHop += Snc_perc_delay.snc_perc_delay(arv, Snc_srv_exp.of(rates[i]), 1e-3 / 3).value;
        }
        System.out.printf("%nend-to-end delay quantile at eps=1e-3%n");
        System.out.printf("  concatenated (Snc_conv) : %8.4f  at theta = %.4f%n",
                concat.value, concat.theta);
        System.out.printf("  summed per hop          : %8.4f%n", hopByHop);
        System.out.printf("  pay bursts once saves   : %7.1f%%%n",
                100.0 * (1.0 - concat.value / hopByHop));

        // A class sharing a station sees the server minus whatever the other
        // classes take from it: Snc_leftover subtracts the cross-flow arrival
        // envelope from the service envelope. The result holds for ANY
        // work-conserving discipline, which is why it is well above the FCFS
        // answer -- it also covers the policy that serves the other class first
        // whenever it can.
        Network shared = new Network("SncShared");
        Source src2 = new Source(shared, "Source");
        Queue qs = new Queue(shared, "Shared", SchedStrategy.FCFS);
        Sink snk2 = new Sink(shared, "Sink");
        OpenClass classA = new OpenClass(shared, "ClassA");
        OpenClass classB = new OpenClass(shared, "ClassB");
        src2.setArrival(classA, new Exp(0.3));
        src2.setArrival(classB, new Exp(0.3));
        qs.setService(classA, new Exp(1.0));
        qs.setService(classB, new Exp(1.0));
        RoutingMatrix Ps = shared.initRoutingMatrix();
        Ps.set(classA, classA, Network.serialRouting(src2, qs, snk2));
        Ps.set(classB, classB, Network.serialRouting(src2, qs, snk2));
        shared.link(Ps);

        SolverBA sharedSolver = new SolverBA(shared, "snc.upper");
        sharedSolver.getAvgTable().print();
        sharedSolver.getPercTable(1e-3).print();
        System.out.printf("exact per-class response time (aggregate M/M/1, lambda=0.6, mu=1): %.4f%n",
                1.0 / (1.0 - 0.6));
    }

    /**
     * @param args unused
     */
    public static void main(String[] args) {
        delayQuantile();
        tandemAndSharing();
    }
}
