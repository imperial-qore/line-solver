/*
 * Copyright (c) 2012-2026, QORE Lab, Imperial College London
 * All rights reserved.
 */

package jline.examples.java.basic;

import jline.api.spn.Spn_lpbnd;
import jline.api.spn.Spn_metrics;
import jline.api.spn.Spn_pf;
import jline.io.LineCitations;
import jline.lang.Network;
import jline.lang.NetworkStruct;
import jline.solvers.NetworkAvgTable;
import jline.solvers.ba.BA;
import jline.solvers.ctmc.CTMC;
import jline.solvers.fluid.FLD;
import jline.solvers.fluid.FluidResult;
import jline.solvers.fluid.petri.PetriSolver;
import jline.solvers.ldes.LDES;
import jline.solvers.mva.MVA;
import jline.solvers.nc.NC;
import jline.solvers.ssa.SSA;
import jline.solvers.wrappers.jmt.JMT;
import jline.util.matrix.Matrix;
import java.util.Scanner;

/**
 * Stochastic Petri net examples mirroring the example notebooks in stochPetriNet.
 * <p>
 * This class contains Java implementations that mirror the example notebooks
 * found in jar/src/main/java/jline/examples/java/basic/stochPetriNet/. Each method
 * demonstrates a specific Petri net concept using models from the basic package.
 * <p>
 * The examples cover:
 * - Basic open and closed Petri net structures
 * - Multiple firing modes and batch processing
 * - Inhibiting conditions and complex token flows
 * - Various stochastic distributions in transitions
 * - Multi-class token systems
 */
public class StochPetriNetExamples {

    private static final Scanner scanner = new Scanner(System.in);

    private static void pauseForUser() {
        // Skip pause if running in non-interactive mode (e.g., Maven exec)
        if (System.console() == null) {
            System.out.println("\n[Running in non-interactive mode, continuing...]");
            return;
        }
        System.out.println("\nPress Enter to continue to next example...");
        try {
            scanner.nextLine();
        } catch (Exception e) {
            // Ignore scanner errors in case of pipe or redirection
        }
    }

    /**
     * Basic closed stochastic Petri net (spn_basic_closed.ipynb).
     * <p>
     * Demonstrates a simple closed Petri net with cyclic token flow
     * between two places using multiple solvers.
     * <p>
     * Features:
     * - 3 tokens initially in Place1, cycling between Place1 and Place2
     * - Exponential firing rates for transitions T1 and T2
     * - Multiple solver comparison (CTMC and JMT)
     * - Token conservation analysis
     */
    /**
     * Pareto-distributed transition firing (spn_pareto_service.ipynb).
     *
     * <p>Simulation only: the golden holds the JMT row at the reference's seed
     * and run length.
     */
    public static void spn_pareto_service() throws Exception {
        Network model = StochPetriNetModel.spn_pareto_service();

        try {
            new jline.solvers.wrappers.jmt.JMT(model, "seed", 23000, "samples", 10000)
                    .getAvgTable().print();
        } catch (Exception e) {
            System.out.println("JMT failed: " + e.getMessage());
        }
    }

    public static void spn_basic_closed() throws Exception {
        Network model = StochPetriNetModel.spn_basic_closed();
        
        // CTMC solver (exact for small Petri nets)
//        try {
//            CTMC solverCtmc = new CTMC(model);
//            solverCtmc.getAvgTable().print();
//        } catch (Exception e) {
//            System.out.println("CTMC solver error: " + e.getMessage());
//        }

        // JMT solver (simulation)
        try {
            JMT solverJmt = new JMT(model, "seed", 23000);
            solverJmt.getAvgTable().print();
        } catch (Exception e) {
            System.out.println("JMT solver not available: " + e.getMessage());
        }
        pauseForUser();
    }
    
    /**
     * Basic open stochastic Petri net (spn_basic_open.ipynb).
     * <p>
     * Shows fundamental open Petri net structure with
     * source, place, transition, and sink.
     */
    public static void spn_basic_open() throws Exception {
        Network model = StochPetriNetModel.spn_basic_open();
        JMT solver = new JMT(model, "seed", 23000);
        
        solver.getAvgTable().print();
        pauseForUser();
    }
    
    /**
     * Batch processing Petri net (spn_twomodes.ipynb).
     * <p>
     * Demonstrates batch token processing where transitions
     * require and produce multiple tokens.
     */
    public static void spn_twomodes() throws Exception {
        Network model = StochPetriNetModel.spn_twomodes();
        JMT solver = new JMT(model, "seed", 23000);
        
        solver.getAvgTable().print();
        pauseForUser();
    }
    
    /**
     * Competing transitions Petri net (spn_fourmodes.ipynb).
     * <p>
     * Shows resource competition between transitions
     * requiring different numbers of tokens.
     */
    public static void spn_fourmodes() throws Exception {
        Network model = StochPetriNetModel.spn_fourmodes();
        JMT solver = new JMT(model, "seed", 23000);
        
        solver.getAvgTable().print();
        pauseForUser();
    }
    
    /**
     * Inhibiting transitions Petri net (spn_inhibiting.ipynb).
     * <p>
     * Demonstrates inhibitor arcs that prevent firing
     * when tokens are present in certain places.
     */
    public static void spn_inhibiting() throws Exception {
        Network model = StochPetriNetModel.spn_inhibiting();
        JMT solver = new JMT(model, "seed", 23000);
        
        solver.getAvgTable().print();
        pauseForUser();
    }
    
    /**
     * Closed Petri net with two places (spn_closed_twoplaces.ipynb).
     * <p>
     * Simple closed system demonstrating token circulation
     * between two places with different service rates.
     */
    public static void spn_closed_twoplaces() throws Exception {
        Network model = StochPetriNetModel.spn_closed_twoplaces();
        
        // CTMC solver for exact solution
        try {
            CTMC solverCtmc = new CTMC(model);
            solverCtmc.getAvgTable().print();
        } catch (Exception e) {
            System.out.println("CTMC solver error: " + e.getMessage());
        }
        
        // JMT solver for simulation
        try {
            JMT solverJmt = new JMT(model, "seed", 23000);
            solverJmt.getAvgTable().print();
        } catch (Exception e) {
            System.out.println("JMT solver not available: " + e.getMessage());
        }
        pauseForUser();
    }
    
    /**
     * Closed Petri net with four places (spn_closed_fourplaces.ipynb).
     * <p>
     * More complex closed system with tokens cycling through
     * four different places using multiple transitions.
     */
    public static void spn_closed_fourplaces() throws Exception {
        Network model = StochPetriNetModel.spn_closed_fourplaces();
        
        // CTMC solver for exact solution
        try {
            CTMC solverCtmc = new CTMC(model);
            solverCtmc.getAvgTable().print();
        } catch (Exception e) {
            System.out.println("CTMC solver error: " + e.getMessage());
        }
        
        // JMT solver for simulation
        try {
            JMT solverJmt = new JMT(model, "seed", 23000);
            solverJmt.getAvgTable().print();
        } catch (Exception e) {
            System.out.println("JMT solver not available: " + e.getMessage());
        }
        pauseForUser();
    }
    
    /**
     * Open Petri net with seven places (spn_open_sevenplaces.ipynb).
     * <p>
     * Complex open system demonstrating token flow through
     * multiple places with varied routing probabilities.
     */
    public static void spn_open_sevenplaces() throws Exception {
        Network model = StochPetriNetModel.spn_open_sevenplaces();
        JMT solver = new JMT(model, "seed", 23000);
        
        solver.getAvgTable().print();
        pauseForUser();
    }

    /**
     * `spn_productform_nc.m`: solve a stochastic Petri net analytically with NC.
     *
     * <p>NC's {@code rec} method is the first ANALYTICAL route LINE offers for a
     * Petri net: CTMC builds the explicit generator, SSA and LDES simulate, FLD
     * fluidises. It works in three steps, each with its own reference:</p>
     *
     * <ul>
     *   <li>{@code Spn_pf} decides whether the net has a product form and
     *       derives the per-place factors g_l, by complex balance
     *       (Coleman-Henderson-Taylor, Perform. Eval. 26(3), 1996);</li>
     *   <li>{@code Mdd_rec} evaluates G as ONE memoised walk of the decision
     *       diagram holding the reachable set (Balsamo-Marin-Stojic, FGCS 111
     *       (2020) 475-490);</li>
     *   <li>{@code Spn_metrics} reads the mean tokens, the utilisations and the
     *       throughputs off masked walks of the same diagram.</li>
     * </ul>
     *
     * <p>Two nets are solved. The first is a closed cycle, which a queueing
     * network could also express. The second FORKS: Tf consumes one token and
     * produces two, so the marking is not a conserved job population and there
     * is no queueing-network counterpart, which is the limitation the MDD-rec
     * paper opens with.</p>
     *
     * @throws Exception if a solver encounters an error
     */
    public static void spn_productform_nc() throws Exception {
        // A closed cycle, where the exact CTMC gives the reference.
        System.out.println("\n--- 3-place cyclic net, N = 4 ---");
        new NC(StochPetriNetModel.spn_productform_cyclic()).getAvgTable().print();

        System.out.println("The same net through the explicit generator, for comparison:");
        new CTMC(StochPetriNetModel.spn_productform_cyclic()).getAvgTable().print();

        // The certificate the derivation produced. Deficiency zero plus weak
        // reversibility is what Feinberg's theorem needs for a positive
        // complex-balanced point to exist at ANY choice of rates.
        Spn_pf.SpnPfOptions pfOpt = new Spn_pf.SpnPfOptions();
        pfOpt.verbose = true;
        Spn_pf.SpnPfResult pf =
                Spn_pf.spn_pf(StochPetriNetModel.spn_productform_cyclic(), pfOpt);
        System.out.printf("product form: %s, deficiency %d, %d linkage classes, rank %d%n",
                pf.kind, pf.deficiency, pf.linkage, pf.srank);

        // A fork-join net, which has no queueing-network form at all.
        Network fj = StochPetriNetModel.spn_productform_forkjoin();
        System.out.println("\n--- fork-join net, P0 -> P1+P2 -> P3 -> P0, 3 tokens at P0 ---");
        new NC(fj).getAvgTable().print();

        Spn_pf.SpnPfResult pfj = Spn_pf.spn_pf(StochPetriNetModel.spn_productform_forkjoin());
        Spn_metrics.SpnMetricsResult met =
                Spn_metrics.spn_metrics(pfj.spn.mdds, pfj.g, pfj.spn.info);
        // Every token that forks must later join and return, so the three modes
        // share one throughput: a flow-conservation law nothing in the
        // derivation was told.
        System.out.println("mode throughputs: " + java.util.Arrays.toString(met.modeTput)
                + " (they must all agree)");
        // The place invariant of this net is 2*m0 + m1 + m2 + 2*m3 = 6, not the
        // token count, which is what "not a conserved population" means.
        double[] w = {2, 1, 1, 2};
        double inv = 0;
        for (int i = 0; i < w.length && i < met.tokens.length; i++) {
            inv += w[i] * met.tokens[i];
        }
        System.out.printf("place invariant 2*m0 + m1 + m2 + 2*m3 = %.6f%n", inv);
        pauseForUser();
    }

    /**
     * `spn_lpbounds.m`: bound a stochastic Petri net by linear programming.
     *
     * <p>BA's {@code spnlp} family is the first BOUNDING route LINE offers for a
     * Petri net, and it needs neither a generator nor a product form. It relaxes
     * the stationary chain to a MOMENT POLYTOPE (the uniformized evolution
     * equation written for E[X_p], E[X_p^2] and E[X_p1 X_p2], plus behavioural
     * and probabilistic inequalities) and then minimises and maximises each
     * reported measure over it. Every stationary point of the true chain
     * satisfies every row, so the two optima bracket the exact value.</p>
     *
     * <p>Reference: Z. Liu, "Performance Analysis of Stochastic Timed Petri Nets
     * Using Linear Programming Approach", IEEE Trans. Software Engineering
     * 24(11), 1998, 1014-1030.</p>
     *
     * @throws Exception if a solver encounters an error
     */
    public static void spn_lpbounds() throws Exception {
        // The reference's own Fig. 2b, and its Table 2. The one column not
        // reproduced is its u.b.1, the upper side further tightened by the
        // subnet-throughput theorems (its Thms 1 and 2); those are not
        // implemented, so u.b.2 is the column to compare against.
        double[][] mus = {{1, 1.25, 2, 0.5}, {1, 1.25, 2, 2.5}, {1, 1.25, 1.25, 2.5},
                          {1, 1.25, 1.25, 1}, {1.111, 1.111, 1.111, 1.111}};
        double[][] pub = {{1.165, 1.951, 2.000, 0.930, 2.000},
                          {1.829, 2.978, 3.529, 1.481, 4.000},
                          {1.581, 2.873, 3.333, 1.333, 4.000},
                          {1.359, 2.757, 3.333, 1.111, 4.000},
                          {1.350, 2.667, 2.963, 1.111, 4.444}};

        System.out.println("\nLiu (1998) Table 2: total throughput of the production line");
        System.out.printf("%-5s %-32s %-19s %-19s%n",
                "case", "Markovian LP", "published", "operational LP");
        System.out.printf("%-5s %9s %9s %9s   %9s %9s   %9s %9s%n",
                "", "lower", "simul", "upper", "l.b.", "u.b.2", "o.l.b.", "o.u.b.");
        for (int c = 0; c < mus.length; c++) {
            NetworkStruct sn = StochPetriNetModel.spn_lpbounds_prodline(mus[c]).getStruct();
            // The liveness rows of the reference's Table 1 are OPT-IN, because
            // they hold only on a live net and spn_lpbnd cannot certify
            // liveness. This one is live: a strongly connected marked graph
            // with a token on every cycle. They are the whole of the lower
            // side, so the published l.b. needs them.
            Spn_lpbnd.SpnLpOptions oLo = new Spn_lpbnd.SpnLpOptions();
            oLo.markovian = true;
            oLo.assumelive = true;
            Spn_lpbnd.SpnLpOptions oUp = new Spn_lpbnd.SpnLpOptions();
            oUp.markovian = true;
            Spn_lpbnd.SpnLpOptions oOp = new Spn_lpbnd.SpnLpOptions();
            oOp.markovian = false;
            oOp.assumelive = true;
            Spn_lpbnd.SpnLpBounds bLo = Spn_lpbnd.spn_lpbnd(sn, oLo);
            Spn_lpbnd.SpnLpBounds bUp = Spn_lpbnd.spn_lpbnd(sn, oUp);
            Spn_lpbnd.SpnLpBounds bOp = Spn_lpbnd.spn_lpbnd(sn, oOp);
            System.out.printf("%-5d %9.4f %9.4f %9.4f   %9.3f %9.3f   %9.4f %9.4f%n", c + 1,
                    rowSum(bLo.modeTput, 0), pub[c][1], rowSum(bUp.modeTput, 1),
                    pub[c][0], pub[c][2], rowSum(bOp.modeTput, 0), rowSum(bOp.modeTput, 1));
        }

        // Through BA, on a net with an inhibitor arc. The four method names are
        // spnlp2.upper/lower and their spnlp1.upper/lower counterparts, which drop
        // the second-moment, covariance and Little's-law families and so need
        // only a mean firing time rather than an exponential one. They are the
        // only family BA offers on a Petri net, and the only one it withholds
        // off a Petri net: every other family is parameterized by demands and a
        // population, which a marking is not.
        Network model = StochPetriNetModel.spn_inhibiting();
        System.out.println("\nmethods offered on this net: "
                + String.join(", ", new BA(model).listValidMethods()));

        NetworkAvgTable lo = new BA(StochPetriNetModel.spn_inhibiting(), "spnlp2.lower").getAvgTable();
        NetworkAvgTable up = new BA(StochPetriNetModel.spn_inhibiting(), "spnlp2.upper").getAvgTable();
        NetworkAvgTable ex = new CTMC(StochPetriNetModel.spn_inhibiting()).getAvgTable();
        System.out.println("\nMean tokens per place, exact between the two sides:");
        System.out.printf("%-8s %10s %10s %10s%n", "place", "lower", "exact", "upper");
        for (int i = 0; i < ex.getStationNames().size(); i++) {
            System.out.printf("%-8s %10.5f %10.5f %10.5f%n", ex.getStationNames().get(i),
                    lo.getQLen().get(i), ex.getQLen().get(i), up.getQLen().get(i));
        }

        // U = Q at a Place, which LINE models as an INF station, the same
        // convention CTMC and NC report. The paper's place utilization
        // 1 - P(m = 0) is a different quantity and is not this column.
        System.out.println("\nsolver.citations():");
        for (LineCitations.Citation c
                : new BA(StochPetriNetModel.spn_inhibiting(), "spnlp2.upper").citations()) {
            System.out.println(c);
        }
        pauseForUser();
    }

    /** Sum of one row of a bound block, which is how a total throughput is read. */
    private static double rowSum(double[][] m, int row) {
        double s = 0;
        if (m != null && m.length > row && m[row] != null) {
            for (int j = 0; j < m[row].length; j++) {
                s += m[row][j];
            }
        }
        return s;
    }

    /**
     * `spn_queueing_place.m`: a closed queueing Petri net (QPN) with two
     * queueing places.
     *
     * <p>A queueing place embeds a scheduling station inside a Petri-net place:
     * an arriving token is served by the place's embedded queue and, on
     * completion, moves to a depository from which the output transitions
     * consume it. A CPU (single-server FCFS) and a think stage (infinite
     * server) here exchange a fixed population of N tokens through two
     * immediate transitions, which is the queueing-Petri-net rendering of a
     * machine-repairman model, so MVA on the equivalent Delay+Queue network
     * cross-validates it exactly.</p>
     *
     * @throws Exception if a solver encounters an error
     */
    public static void spn_queueing_place() throws Exception {
        int N = 4;
        // Queueing places are simulated by LDES.
        new LDES(StochPetriNetModel.spn_queueing_place(N), "seed", 23000, "samples", 2e5)
                .getAvgTable().print();
        // Exact cross-check: the equivalent finite-population Delay + M/M/1.
        new MVA(StochPetriNetModel.spn_queueing_place_ref(N)).getAvgTable().print();
        pauseForUser();
    }

    /**
     * `spn_colored_gspn.m`: a colored generalized stochastic Petri net.
     *
     * <p>A GSPN mixes timed transitions, which fire after a random delay, with
     * immediate ones, which fire as soon as they are enabled and consume no
     * time. It is colored when the tokens carry a type, which in LINE is a job
     * class: a place holds one marking per color and a transition declares one
     * mode per color. Two colors circulate between a buffer and a server; the
     * immediate transition admits a token only into an empty server, through
     * inhibitor arcs on BOTH colors, and the firing weights arbitrate the
     * conflict when both colors are waiting. CTMC eliminates the vanishing
     * markings by stochastic complementation and solves the 8 tangible states
     * exactly; SSA reproduces them by simulation.</p>
     *
     * @throws Exception if a solver encounters an error
     */
    public static void spn_colored_gspn() throws Exception {
        Network model = StochPetriNetModel.spn_colored_gspn();
        new CTMC(model, "cutoff", 4).getAvgTable().print();
        new SSA(model, "seed", 23000, "samples", 2e5).getAvgTable().print();
        pauseForUser();
    }

    /**
     * `spn_fluid_dae.m`: fluid (mean-field) analysis of a stochastic Petri net.
     *
     * <p>A GSPN is a density-dependent Markov population process: the marking is
     * the population, a transition mode is a reaction, and the firing rate
     * {@code lambda*min(enabling degree, servers)} is the same min()
     * non-linearity the min-normal closure of FLD exists to smooth. The
     * {@code dae} method is the one that can carry it, because a Petri net needs
     * three things stated as EQUATIONS rather than integrated: the P-invariants,
     * which hold to solver tolerance instead of integrator tolerance; the firing
     * FLOW of an immediate transition, an algebraic unknown pinned by the
     * constraint that its input place holds no mass; and a bounded place, a
     * linear inequality on the marking.</p>
     *
     * <p>FLD resolves to {@code dae} on any model holding a Transition node, so
     * no method has to be named. Unlike every other solver of a Petri net in
     * LINE it also returns a SECOND MOMENT: the marking covariance of the linear
     * noise approximation, on {@code FluidResult.petri}.</p>
     *
     * @throws Exception if a solver encounters an error
     */
    public static void spn_fluid_dae() throws Exception {
        // A closed net whose fluid answer is EXACT.
        Network exact = StochPetriNetModel.spn_fluid_exact();
        FLD fld = new FLD(exact);
        fld.getAvgTable().print();
        new CTMC(StochPetriNetModel.spn_fluid_exact(), "cutoff", 6).getAvgTable().print();

        // The exact marking is Binomial(4, 3/5), so the variance is
        // 4*0.6*0.4 = 0.96.
        PetriSolver.PetriReport pr = ((FluidResult) fld.result).petri;
        System.out.printf("marking variance: %.6f %.6f  (exact 0.96)%n",
                pr.markingVar.get(0, 0), pr.markingVar.get(1, 0));
        System.out.printf("invariant \"%s\" = %g, error %.2e%n",
                pr.invariantLabel.get(0), pr.invariantValue[0], pr.invariantError[0]);

        // An immediate transition, as an algebraic flow. The vanishing place P2
        // holds exactly zero mass and the net answers as the reduced two-place
        // net does, which is what the algebraic flow buys: an approximation of
        // the immediate transition by a large finite rate would only approach it.
        new FLD(StochPetriNetModel.spn_fluid_immediate()).getAvgTable().print();
        new FLD(StochPetriNetModel.spn_fluid_reduced()).getAvgTable().print();
        pauseForUser();
    }

    /**
     * `test_spn_nrm_open.m`: the SSA Next-Reaction-Method path on OPEN nets.
     *
     * <p>A Source feeds a Place whose tokens drain through a Transition to a
     * Sink. Before the Source-arrival reaction was added the fed Place stayed
     * empty and the run threw "Deadlock: no transition is enabled".</p>
     *
     * <p>Each net asserts that the solver actually ran method {@code nrm} (never
     * a silent serial fallback) and that the simulated marking mean and
     * throughput match the analytic M/M/1 result: a Source Exp(lambda) feeding a
     * single-server Transition Exp(mu) is an M/M/1 queue at the Place, with mean
     * tokens rho/(1-rho) and throughput lambda. The canonical net is
     * cross-checked against JMT.</p>
     *
     * @throws Exception if a solver encounters an error
     */
    public static void test_spn_nrm_open() throws Exception {
        final double RTOL = 0.04;
        final double SAMPLES = 3e5;
        final int SEED = 23000;

        // Net 1: M/M/1 SPN, Source Exp(0.5) -> P1 -> T1 Exp(1.0) -> Sink.
        double lambda = 0.5;
        double mu = 1.0;
        double rho = lambda / mu;
        double qExact = rho / (1 - rho);   // = 1.0
        SSA solver = new SSA(StochPetriNetModel.spn_nrm_mm1(lambda, mu),
                "method", "nrm", "samples", SAMPLES, "seed", SEED);
        solver.getAvg();
        Matrix Qn = solver.result.QN;
        Matrix Tn = solver.result.TN;
        require(solver.result.method != null && solver.result.method.contains("nrm"),
                "net1 did not run NRM");
        // stations: Source(0), P1(1)
        require(Math.abs(Qn.get(1, 0) - qExact) / qExact < RTOL,
                String.format("net1 P1 tokens %g vs exact %g", Qn.get(1, 0), qExact));
        require(Math.abs(Tn.get(0, 0) - lambda) / lambda < RTOL,
                String.format("net1 Source tput %g vs %g", Tn.get(0, 0), lambda));
        require(Math.abs(Tn.get(1, 0) - lambda) / lambda < RTOL,
                String.format("net1 P1 tput %g vs %g", Tn.get(1, 0), lambda));

        // Cross-check mean tokens against JMT's simulation of the same net.
        JMT jmt = new JMT(StochPetriNetModel.spn_nrm_mm1(lambda, mu),
                "samples", SAMPLES, "seed", SEED);
        jmt.getAvg();
        Matrix Qj = jmt.result.QN;
        require(Math.abs(Qn.get(1, 0) - Qj.get(1, 0)) / Math.max(Qj.get(1, 0), 1e-9) < RTOL,
                String.format("net1 P1 tokens NRM %g vs JMT %g", Qn.get(1, 0), Qj.get(1, 0)));

        // Net 2: open tandem, two places in series.
        double mu1 = 1.0;
        double mu2 = 2.0;
        double q1 = (lambda / mu1) / (1 - lambda / mu1);   // = 1.0
        double q2 = (lambda / mu2) / (1 - lambda / mu2);   // = 1/3
        SSA solver2 = new SSA(StochPetriNetModel.spn_nrm_tandem(lambda, mu1, mu2),
                "method", "nrm", "samples", SAMPLES, "seed", SEED);
        solver2.getAvg();
        Matrix Qn2 = solver2.result.QN;
        Matrix Tn2 = solver2.result.TN;
        require(solver2.result.method != null && solver2.result.method.contains("nrm"),
                "net2 did not run NRM");
        // stations: Source(0), P1(1), P2(2)
        require(Math.abs(Qn2.get(1, 0) - q1) / q1 < RTOL,
                String.format("net2 P1 tokens %g vs exact %g", Qn2.get(1, 0), q1));
        require(Math.abs(Qn2.get(2, 0) - q2) / q2 < RTOL,
                String.format("net2 P2 tokens %g vs exact %g", Qn2.get(2, 0), q2));
        require(Math.abs(Tn2.get(1, 0) - lambda) / lambda < RTOL,
                String.format("net2 P1 tput %g vs %g", Tn2.get(1, 0), lambda));
        require(Math.abs(Tn2.get(2, 0) - lambda) / lambda < RTOL,
                String.format("net2 P2 tput %g vs %g", Tn2.get(2, 0), lambda));

        System.out.println("test_spn_nrm_open passed");
        pauseForUser();
    }

    /** The reference's {@code assert}: a failed check is the example's result. */
    private static void require(boolean condition, String message) {
        if (!condition) {
            throw new RuntimeException(message);
        }
    }

    /**
     * Main method demonstrating selected Petri net examples.
     */
    public static void main(String[] args) throws Exception {
        System.out.println("\n=== Running example: spn_basic_closed ===");
        try {
            spn_basic_closed();
        } catch (Exception e) {
            System.err.println("spn_basic_closed failed: " + e.getMessage());
            e.printStackTrace();
        }
        
        System.out.println("\n=== Running example: spn_basic_open ===");
        try {
            spn_basic_open();
        } catch (Exception e) {
            System.err.println("spn_basic_open failed: " + e.getMessage());
            e.printStackTrace();
        }
        
        System.out.println("\n=== Running example: spn_twomodes ===");
        try {
            spn_twomodes();
        } catch (Exception e) {
            System.err.println("spn_twomodes failed: " + e.getMessage());
            e.printStackTrace();
        }
        
        System.out.println("\n=== Running example: spn_fourmodes ===");
        try {
            spn_fourmodes();
        } catch (Exception e) {
            System.err.println("spn_fourmodes failed: " + e.getMessage());
            e.printStackTrace();
        }
        
        System.out.println("\n=== Running example: spn_inhibiting ===");
        try {
            spn_inhibiting();
        } catch (Exception e) {
            System.err.println("spn_inhibiting failed: " + e.getMessage());
            e.printStackTrace();
        }
        
        System.out.println("\n=== Running example: spn_closed_twoplaces ===");
        try {
            spn_closed_twoplaces();
        } catch (Exception e) {
            System.err.println("spn_closed_twoplaces failed: " + e.getMessage());
            e.printStackTrace();
        }
        
        System.out.println("\n=== Running example: spn_closed_fourplaces ===");
        try {
            spn_closed_fourplaces();
        } catch (Exception e) {
            System.err.println("spn_closed_fourplaces failed: " + e.getMessage());
            e.printStackTrace();
        }
        
        System.out.println("\n=== Running example: spn_open_sevenplaces ===");
        try {
            spn_open_sevenplaces();
        } catch (Exception e) {
            System.err.println("spn_open_sevenplaces failed: " + e.getMessage());
            e.printStackTrace();
        }
        
        System.out.println("\n=== Running example: spn_productform_nc ===");
        try {
            spn_productform_nc();
        } catch (Exception e) {
            System.err.println("spn_productform_nc failed: " + e.getMessage());
            e.printStackTrace();
        }
        
        System.out.println("\n=== Running example: spn_lpbounds ===");
        try {
            spn_lpbounds();
        } catch (Exception e) {
            System.err.println("spn_lpbounds failed: " + e.getMessage());
            e.printStackTrace();
        }
        
        System.out.println("\n=== Running example: spn_queueing_place ===");
        try {
            spn_queueing_place();
        } catch (Exception e) {
            System.err.println("spn_queueing_place failed: " + e.getMessage());
            e.printStackTrace();
        }
        
        System.out.println("\n=== Running example: spn_colored_gspn ===");
        try {
            spn_colored_gspn();
        } catch (Exception e) {
            System.err.println("spn_colored_gspn failed: " + e.getMessage());
            e.printStackTrace();
        }
        
        System.out.println("\n=== Running example: spn_fluid_dae ===");
        try {
            spn_fluid_dae();
        } catch (Exception e) {
            System.err.println("spn_fluid_dae failed: " + e.getMessage());
            e.printStackTrace();
        }
        
        System.out.println("\n=== Running example: test_spn_nrm_open ===");
        try {
            test_spn_nrm_open();
        } catch (Exception e) {
            System.err.println("test_spn_nrm_open failed: " + e.getMessage());
            e.printStackTrace();
        }
        
        scanner.close();
    }
}
