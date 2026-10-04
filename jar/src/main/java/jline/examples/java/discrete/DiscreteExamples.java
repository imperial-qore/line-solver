/*
 * Copyright (c) 2012-2026, QORE Lab, Imperial College London
 * All rights reserved.
 */

package jline.examples.java.discrete;

import jline.api.dpfqn.Dpfqn_nc;
import jline.api.dpfqn.Dpfqn_ncld;
import jline.api.dpfqn.DpfqnNcLdResult;
import jline.api.dpfqn.DpfqnNcResult;
import jline.api.dqsys.Bernoulli1Result;
import jline.api.dqsys.Dqsys_bernoulli1;
import jline.lang.Network;
import jline.solvers.SolverOptions;
import jline.solvers.nc.SolverNC;

import java.util.Arrays;

/**
 * Discrete-time (slotted) examples, one per feature of the discrete-time
 * normalizing-constant route.
 *
 * <p>{@link SolverNC} enters this route only when {@code options.config.slotted}
 * is true, the same switch {@link jline.solvers.ldes.SolverLDES} uses; it is
 * never inferred from the presence of a {@code Geometric} distribution. When
 * the switch is on and the model falls outside the discrete-time product form,
 * the solver errors instead of approximating.</p>
 *
 * <p>Reference: H. Daduna, Queueing Networks with Discrete Time Scale, LNCS
 * 2046, Springer, 2001.</p>
 */
public class DiscreteExamples {

    private static SolverOptions slotted() {
        SolverOptions options = new SolverOptions(jline.lang.constant.SolverType.NC);
        options.config.slotted = true;
        return options;
    }

    /** Geo/Geo/1 with an unbounded buffer. */
    public static void ex1_geogeo1() throws IllegalAccessException {
        System.out.println("\n--- ex1_geogeo1: Geo/Geo/1 on the slot lattice ---");
        double a = 0.2;
        double s = 0.5;
        Network model = DiscreteModel.dt_geogeo1(a, s);
        new SolverNC(model, slotted()).getAvgTable().print();
        System.out.printf("closed form: E[N] = %g, E[T] = %g slots%n",
                a * (1 - a) / (s - a), (1 - a) / (s - a));
    }

    /** Geo/Geo/1/L loss system and its blocking probability. */
    public static void ex2_geogeo1_loss() throws IllegalAccessException {
        System.out.println("\n--- ex2_geogeo1_loss: finite buffer, corollary 2.8 ---");
        double a = 0.2;
        double s = 0.5;
        int L = 4;
        Network model = DiscreteModel.dt_geogeo1_loss(a, s, L);
        new SolverNC(model, slotted()).getAvgTable().print();
        Bernoulli1Result r = Dqsys_bernoulli1.dqsys_bernoulli1(new double[]{a}, new double[]{s}, L);
        System.out.println("queue length law   : " + Arrays.toString(r.pmf));
        System.out.printf("loss probability   : %g%n", r.lossProb);
        System.out.printf("carried throughput : %g of the %g offered per slot%n", r.throughput, a);
    }

    /**
     * Load-dependent Bernoulli server, and the discrete-time arrival theorem.
     *
     * <p>Example 2.10 notes that a discrete-time M/M/c queue has no exactly
     * equivalent state dependent single server, but that p(n) = p min(n,c)
     * reproduces its conditional service intensity. The arrival law of theorem
     * 2.11 is printed alongside the time-stationary one: discrete time has no
     * PASTA analogue and the two differ.</p>
     */
    public static void ex3_bernoulli_loaddep() throws IllegalAccessException {
        System.out.println("\n--- ex3_bernoulli_loaddep: example 2.10 and theorem 2.11 ---");
        double a = 0.6;
        double s = 0.3;
        int c = 3;
        int L = 20;
        Network model = DiscreteModel.dt_bernoulli_loaddep(a, s, c, L);
        new SolverNC(model, slotted()).getAvgTable().print();
        double[] p = new double[L];
        for (int n = 1; n <= L; n++) {
            p[n - 1] = s * Math.min(n, c);
        }
        Bernoulli1Result r = Dqsys_bernoulli1.dqsys_bernoulli1(new double[]{a}, p, L);
        System.out.println("time-stationary law pi(0..6)  : "
                + Arrays.toString(Arrays.copyOf(r.pmf, 7)));
        System.out.println("arrival law         pi_1(0..6): "
                + Arrays.toString(Arrays.copyOf(r.arrivalPmf, 7)));
        System.out.println("no PASTA in discrete time: the two rows above are different laws");
    }

    /** Closed cycle of Bernoulli servers. */
    public static void ex4_cycle() throws IllegalAccessException {
        System.out.println("\n--- ex4_cycle: closed cycle, corollary 3.4 ---");
        double[] p = DiscreteModel.defaultCycleRates();
        int N = 5;
        Network model = DiscreteModel.dt_cycle(p, N);
        new SolverNC(model, slotted()).getAvgTable().print();
        DpfqnNcResult nc = Dpfqn_nc.dpfqn_nc(p, N);
        System.out.printf("log G(N,J)   = %g%n", nc.lG);
        System.out.printf("throughput   = %g jobs per slot%n", nc.throughput());
        double[] util = new double[p.length];
        for (int j = 0; j < p.length; j++) {
            util[j] = nc.throughput() / p[j];
        }
        System.out.println("utilizations = " + Arrays.toString(util));
    }

    /** Closed cycle with state dependent service probabilities. */
    public static void ex5_cycle_loaddep() throws IllegalAccessException {
        System.out.println("\n--- ex5_cycle_loaddep: closed cycle, theorem 3.2 ---");
        double[] p = DiscreteModel.defaultCycleRates();
        int N = 5;
        Network model = DiscreteModel.dt_cycle_loaddep(p, N, 1, 2);
        new SolverNC(model, slotted()).getAvgTable().print();
        double[][] P = new double[p.length][N];
        for (int j = 0; j < p.length; j++) {
            for (int n = 1; n <= N; n++) {
                P[j][n - 1] = (j == 1) ? p[j] * Math.min(n, 2) : p[j];
            }
        }
        DpfqnNcLdResult nc = Dpfqn_ncld.dpfqn_ncld(P, N);
        System.out.println("P(X_2 = 0..N) = " + Arrays.toString(nc.marginal(1)));
    }

    /** Multichain closed cycle. */
    public static void ex6_cycle_multiclass() throws IllegalAccessException {
        System.out.println("\n--- ex6_cycle_multiclass: multichain cycle, section 3.2 ---");
        double[] p = DiscreteModel.defaultCycleRates();
        int[] pops = new int[]{3, 2};
        Network model = DiscreteModel.dt_cycle_multiclass(p, pops);
        new SolverNC(model, slotted()).getAvgTable().print();
        DpfqnNcResult nc = Dpfqn_nc.dpfqn_nc(p, pops[0] + pops[1]);
        double X = nc.throughput();
        System.out.printf("aggregate throughput  = %g jobs per slot%n", X);
        System.out.printf("per-chain throughput  = %g and %g%n",
                X * pops[0] / (pops[0] + pops[1]), X * pops[1] / (pops[0] + pops[1]));
    }

    public static void main(String[] args) throws IllegalAccessException {
        ex1_geogeo1();
        ex2_geogeo1_loss();
        ex3_bernoulli_loaddep();
        ex4_cycle();
        ex5_cycle_loaddep();
        ex6_cycle_multiclass();
    }
}
