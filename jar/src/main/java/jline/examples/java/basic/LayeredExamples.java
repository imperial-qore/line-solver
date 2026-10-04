/*
 * Copyright (c) 2012-2026, QORE Lab, Imperial College London
 * All rights reserved.
 */

package jline.examples.java.basic;

import jline.VerboseLevel;
import jline.lang.constant.SolverType;
import jline.lang.Network;
import jline.lang.layered.LayeredNetwork;
import jline.api.infer.InferLqn;
import jline.api.infer.InferLqnOptions;
import jline.api.infer.InferLqnResult;
import jline.api.infer.ObsSpec;
import jline.api.infer.ParamSpec;
import jline.io.Ret.DistributionResult;
import jline.solvers.AvgTable;
import jline.solvers.LayeredNetworkAvgTable;
import jline.solvers.fluid.FLD;
import jline.solvers.ln.LN;
import jline.solvers.ln.LNTranAvgResult;
import jline.solvers.wrappers.lqns.LQNS;
import jline.solvers.mva.MVA;
import jline.solvers.nc.NC;
import jline.solvers.ssa.SSA;
import jline.solvers.SolverOptions;
import jline.solvers.SolverResult;
import jline.util.matrix.Matrix;
import java.util.Arrays;
import java.util.HashMap;
import java.util.List;
import java.util.Map;
import java.util.Scanner;

/**
 * Layered network examples mirroring the example notebooks in layeredModel.
 * <p>
 * This class contains Java implementations that mirror the example notebooks
 * found in jar/src/main/java/jline/examples/java/basic/layeredModel/. Each method
 * demonstrates a specific layered network concept using models from the basic package.
 * <p>
 * The examples cover:
 * - Basic layered network structures with processors and tasks
 * - Activity precedence patterns and synchronous calls
 * - Multi-solver approaches for layered networks
 * - Enterprise application modeling patterns
 * - BPMN-style workflow representations
 */
public class LayeredExamples {

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
     * Options that will print their table.
     *
     * <p>A solver built with {@code VerboseLevel.SILENT} writes that level into
     * {@code GlobalConstants}, and {@code NetworkAvgTable.print()} honours the
     * options the table carries, so an unqualified print after a silent solve
     * emits NOTHING. The reference silences the solver's per-layer chatter and
     * then displays the table regardless, which is what this asks for by name.
     *
     * @return LN options at {@code VerboseLevel.STD}
     */
    private static SolverOptions printing() {
        SolverOptions o = new SolverOptions(SolverType.LN);
        o.verbose = VerboseLevel.STD;
        return o;
    }


    /**
     * Basic layered network (lqn_basic.ipynb).
     * <p>
     * Demonstrates fundamental layered network concepts with
     * processors, tasks, and activity precedence.
     */
    /**
     * Round-robin call dispatch over a set of target tasks (lqn_rrobin.m).
     *
     * <p>Only the squashed ("flat") layering can express this: under "srvn" each
     * server task lives in its own submodel and is replaced, in the client's
     * submodel, by a surrogate delay, so no node ever has arcs to more than one
     * of them. The layer solver must also implement state-dependent routing (SSA
     * here); MVA, NC and FLD are rejected rather than silently returning the
     * probabilistic split.</p>
     *
     * @throws Exception if the solver fails
     */
    public static void lqn_rrobin() throws Exception {
        LayeredNetwork model = LayeredModel.lqn_rrobin();

        SolverOptions options = LN.defaultOptions();
        options.config.layering = "flat";
        options.verbose = VerboseLevel.SILENT;
        AvgTable avgTable = new LN(model, m -> new SSA(m, "verbose", false), options,
                SolverType.SSA).getAvgTable();
        avgTable.setOptions(printing());
        avgTable.print();
        pauseForUser();
    }

    /**
     * Join-the-shortest-queue call dispatch over a set of target tasks (lqn_jsq.m).
     *
     * <p>The twin of {@link #lqn_rrobin()}, with the cyclic pointer replaced by
     * the least loaded target at dispatch time. Like round-robin it needs the
     * squashed layering and a layer solver carrying state-dependent routing.</p>
     *
     * @throws Exception if the solver fails
     */
    public static void lqn_jsq() throws Exception {
        LayeredNetwork model = LayeredModel.lqn_jsq();

        SolverOptions options = LN.defaultOptions();
        options.config.layering = "flat";
        options.verbose = VerboseLevel.SILENT;
        AvgTable avgTable = new LN(model, m -> new SSA(m, "verbose", false), options,
                SolverType.SSA).getAvgTable();
        avgTable.setOptions(printing());
        avgTable.print();
        pauseForUser();
    }

    /**
     * Iteration trace of the layered fixed point (lqn_bpmn_trace.m).
     *
     * <p>Drives the layered solver by hand (init, pre, analyze, post, converged)
     * instead of through the analyzer, printing the largest throughput and queue
     * length each layer carries at each iteration. That is what makes a
     * non-converging or a silently-zero layer visible, which an aggregate table
     * hides.</p>
     *
     * @throws Exception if the solver fails
     */
    public static void lqn_bpmn_trace() throws Exception {
        LayeredNetwork model = LayeredModel.lqn_bpmn_trace();

        SolverOptions lnoptions = LN.defaultOptions();
        lnoptions.verbose = VerboseLevel.SILENT;
        lnoptions.iter_max = 10;
        SolverOptions options = MVA.defaultOptions();
        options.verbose = VerboseLevel.SILENT;
        LN solverLN = new LN(model, m -> new MVA(m, options), lnoptions);

        System.out.println("=== Running LN solver with tracing ===");
        solverLN.init();

        Map<Integer, Map<Integer, SolverResult>> results =
                new HashMap<Integer, Map<Integer, SolverResult>>();
        for (int it = 1; it <= 5; it++) {
            System.out.println("=== Iteration " + it + " ===");
            solverLN.pre(it);
            Map<Integer, SolverResult> perLayer = new HashMap<Integer, SolverResult>();
            // the layer index is 0-based here and 1-based in the reference, and
            // post() reads the map back by that index, so only the printed
            // number is shifted
            for (int e = 0; e < solverLN.nlayers; e++) {
                SolverResult res = solverLN.analyze(it, e);
                if (res == null) {
                    continue;
                }
                perLayer.put(e, res.deepCopy());
                double maxTN = res.TN == null || res.TN.isEmpty() ? 0 : res.TN.elementMax();
                double maxQN = res.QN == null || res.QN.isEmpty() ? 0 : res.QN.elementMax();
                if (maxTN > 1e-10 || maxQN > 1e-10) {
                    System.out.printf("  Layer %d: max(TN)=%.4e, max(QN)=%.4e%n",
                            e + 1, maxTN, maxQN);
                }
            }
            results.put(it, perLayer);
            solverLN.setEnsembleResults(results);
            solverLN.post(it);

            if (solverLN.converged(it)) {
                System.out.println("Converged at iteration " + it);
                break;
            }
        }

        System.out.println("=== Final check: Getting AvgTable ===");
        LayeredNetworkAvgTable avgTable = (LayeredNetworkAvgTable) solverLN.getAvgTable();
        String[] taskNames = {"R1_Task", "R2_Task", "R3_Task", "R1A_Task", "R1B_Task",
                "R2A_Task", "R2B_Task"};
        List<String> nodeNames = avgTable.getNodeNames();
        List<Double> tput = avgTable.getTput();
        for (String taskName : taskNames) {
            int idx = nodeNames.indexOf(taskName);
            if (idx >= 0) {
                System.out.printf("%s: Tput=%.6f%n", taskName, tput.get(idx));
            }
        }
        pauseForUser();
    }

    /**
     * Identify hidden LQN parameters from measured performance (lqn_paramident.m).
     *
     * <p>Demonstrates {@code InferLqn.inferLqn}, an Extended Kalman Filter that
     * tracks hidden LQN parameters (host demands, think times) from measurable
     * performance data, following Zheng, Yang, Woodside, Litoiu, Iszlai,
     * "Tracking Time-Varying Parameters in Software Systems with Extended Kalman
     * Filters", CASCON 2005.</p>
     *
     * <p>Two parameters are hidden: the reference-task think time (the paper's Z)
     * and the P2 host demand of activity AS3 (the paper's service demand S_d).
     * The measurable vector is [R(E1), U(P1), U(P2)]. A measurement sequence with
     * a step change plus noise is synthesised, then the parameter trajectory is
     * recovered.</p>
     *
     * @throws Exception if the solver fails
     */
    public static void lqn_paramident() throws Exception {
        LayeredNetwork model = LayeredModel.lqn_paramident();

        // what to estimate (paramSpec) and what is observed (obsSpec)
        List<ParamSpec> paramSpec = Arrays.asList(
                ParamSpec.think("T1"), ParamSpec.hostDemand("AS3"));
        List<ObsSpec> obsSpec = Arrays.asList(
                ObsSpec.respT("E1"), ObsSpec.util("P1"), ObsSpec.util("P2"));

        // ground-truth parameter trajectory with a step change
        int nsteps = 24;
        Matrix aTrueSeq = new Matrix(2, nsteps);
        for (int k = 0; k < nsteps; k++) {
            aTrueSeq.set(0, k, k >= 12 ? 1.0 : 0.5);          // step: Z doubles at step 13
            aTrueSeq.set(1, k, (k >= 6 && k <= 17) ? 1.0 / 25 : 1.0 / 50); // S_d pulse
        }

        // synthesise the measurement matrix Z from the true model plus noise
        SolverOptions solveropts = LN.defaultOptions();
        solveropts.verbose = VerboseLevel.SILENT;
        int no = obsSpec.size();
        Matrix Z = new Matrix(no, nsteps);
        java.util.Random rng = new java.util.Random(12345);
        for (int k = 0; k < nsteps; k++) {
            Matrix a = new Matrix(2, 1);
            a.set(0, 0, aTrueSeq.get(0, k));
            a.set(1, 0, aTrueSeq.get(1, k));
            InferLqn.setParams(model, paramSpec, a);
            LayeredNetworkAvgTable table =
                    (LayeredNetworkAvgTable) new LN(model, solveropts).getEnsembleAvg();
            Matrix ztrue = InferLqn.getObs(table, obsSpec);
            for (int o = 0; o < no; o++) {
                // ~2% measurement noise
                Z.set(o, k, ztrue.get(o, 0) * (1 + 0.02 * rng.nextGaussian()));
            }
        }

        // run the EKF identification
        InferLqnOptions opt = new InferLqnOptions();
        opt.a0 = new Matrix(2, 1);
        opt.a0.set(0, 0, 0.7);          // deliberately wrong initial guess
        opt.a0.set(1, 0, 1.0 / 40);
        opt.QFac = 0.1;                 // eq 9a drift-noise factor
        opt.RFac = 0.2;                 // eq 9b measurement-noise factor
        opt.gammaT = 1.0;               // T = T* (measurement interval == system constant)
        opt.aTrue = new Matrix(2, 1);
        opt.aTrue.set(0, 0, aTrueSeq.get(0, nsteps - 1));
        opt.aTrue.set(1, 0, aTrueSeq.get(1, nsteps - 1));
        InferLqnResult info = InferLqn.inferLqn(model, paramSpec, obsSpec, Z, opt);

        System.out.printf("%nStep | Z_true  Z_hat  | Sd_true  Sd_hat  | ||e||%n");
        for (int k = 0; k < nsteps; k++) {
            double enorm = 0;
            for (int o = 0; o < no; o++) {
                enorm += info.e.get(o, k) * info.e.get(o, k);
            }
            System.out.printf("%4d | %6.3f  %6.3f | %7.4f  %7.4f | %.3g%n", k + 1,
                    aTrueSeq.get(0, k), info.ahat.get(0, k),
                    aTrueSeq.get(1, k), info.ahat.get(1, k), Math.sqrt(enorm));
        }
        System.out.printf(
                "%nFinal estimate:  Z = %.4f (true %.4f),  Sd = %.5f (true %.5f)%n",
                info.ahat.get(0, nsteps - 1), aTrueSeq.get(0, nsteps - 1),
                info.ahat.get(1, nsteps - 1), aTrueSeq.get(1, nsteps - 1));
        System.out.printf(
                "RMS parameter tracking error Ea = %.4g,  prediction error Er = %.4g%n",
                info.Ea, info.Er);
        pauseForUser();
    }

    /**
     * Transient (time-dependent) analysis of a layered network (lqn_transient.m).
     *
     * <p>{@code getTranAvg} returns the transient mean queue length, utilization
     * and throughput of every ensemble layer over time. The traces are assembled
     * block-diagonally: layer e occupies a disjoint block of rows (its stations)
     * and columns (its classes). Only transient-capable per-layer solvers
     * (Fluid, CTMC, SSA) produce them; here the layers are solved with the fluid
     * ODE solver. Do NOT set a timespan on the per-layer factory: the
     * steady-state fixed point rejects it, and getTranAvg auto-selects the
     * timespan per layer.</p>
     *
     * <p>The initial point is the layer's default state (all closed jobs at the
     * reference station), NOT the converged occupancy, so each curve relaxes
     * from all-at-reference to the layer steady state.</p>
     *
     * @throws Exception if the solver fails
     */
    public static void lqn_transient() throws Exception {
        LayeredNetwork model = LayeredModel.lqn_transient();

        SolverOptions options = LN.defaultOptions();
        options.verbose = VerboseLevel.SILENT;
        LN solver = new LN(model, m -> new FLD(m, "verbose", false), options);
        LNTranAvgResult tran = solver.getTranAvg();

        int E = solver.nlayers;
        System.out.printf("LN.getTranAvg returned transient traces for %d layers.%n", E);

        // The block-diagonal offsets (r0,c0) advance by each layer's station and
        // class counts, mirroring how getTranAvg stacks the per-layer blocks.
        int r0 = 0;
        int c0 = 0;
        for (int e = 0; e < E; e++) {
            Network layer = solver.getEnsemble().get(e);
            int M = layer.getNumberOfStations();
            int K = layer.getNumberOfClasses();
            System.out.printf("%n--- Layer %d (%d stations, %d classes)%n", e + 1, M, K);
            for (int i = 0; i < M; i++) {
                for (int r = 0; r < K; r++) {
                    Matrix trace = tran.QNt[r0 + i][c0 + r];
                    Matrix t = tran.t[r0 + i][c0 + r];
                    if (trace == null || trace.isEmpty() || trace.elementMax() <= 1e-6) {
                        continue;
                    }
                    int last = trace.getNumRows() - 1;
                    System.out.printf("  %-10s [%-10s] E[N](0) = %8.5f -> E[N](%.3f) = %8.5f%n",
                            layer.getStations().get(i).getName(),
                            layer.getClasses().get(r).getName(),
                            trace.get(0, 0), t.get(last, 0), trace.get(last, 0));
                }
            }
            r0 += M;
            c0 += K;
        }
        pauseForUser();
    }

    /**
     * Initializing a layered run from the LQNS solution (lqn_init.m).
     *
     * <p>Runs the external solver first, when it is installed, then the layered
     * solver on its own, and finally reads the response-time CDF off the layer
     * the calls terminate in. The layer solver must be transient-capable for
     * that last step, which is why it is the fluid one.</p>
     *
     * @throws Exception if the solver fails
     */
    public static void lqn_init() throws Exception {
        System.out.println(
                "This example illustrates the initialization of LN using the output of LQNS.");
        LayeredNetwork model = LayeredModel.lqn_serial();

        // LQNS, whose solution is what an LN run can be initialized from
        if (LQNS.isAvailable()) {
            SolverOptions options = LQNS.defaultOptions();
            options.keep = true; // keep the intermediate XML files of the translation
            AvgTable lqnsTable = new LQNS(model, options).getLNAvgTable();
            System.out.println("\nLQNS Results:");
            lqnsTable.setOptions(printing());
            lqnsTable.print();
        } else {
            System.out.println("\nLQNS solver not available - skipping.");
        }

        // LN without initialization, for the elapsed time the initialization saves
        System.out.println("\nSolve with LN without initialization:");
        long t0 = System.nanoTime();
        AvgTable lnTable = new LN(model, m -> new MVA(m)).getAvgTable();
        double elapsed = (System.nanoTime() - t0) / 1e9;
        lnTable.setOptions(printing());
        lnTable.print();
        System.out.printf("Time elapsed: %.3fs%n", elapsed);

        // CDF of response times, taken on the layer the calls terminate in
        System.out.println("\nWe now obtain the CDF of response times:");
        java.util.List<Network> ensemble = model.getEnsemble();
        if (ensemble.size() >= 3) {
            DistributionResult RD = new FLD(ensemble.get(2)).getCdfRespT();
            System.out.println("RD (CDF of response times):");
            for (int i = 0; i < RD.cdfData.size(); i++) {
                for (int r = 0; r < RD.cdfData.get(i).size(); r++) {
                    Matrix cdf = RD.cdfData.get(i).get(r);
                    if (cdf == null || cdf.isEmpty()) {
                        continue;
                    }
                    System.out.printf("  station %d, class %d: %d points, F(%.4f) = %.6f%n",
                            i + 1, r + 1, cdf.getNumRows(),
                            cdf.get(cdf.getNumRows() - 1, 1), cdf.get(cdf.getNumRows() - 1, 0));
                }
            }
        } else {
            System.out.println("Model ensemble has fewer than three layers.");
        }
        pauseForUser();
    }

    /**
     * Method 'flat.ph': the squashed layering with the composed encoding (lqn_flatph.m).
     *
     * <p>A method name of the layered solver carries TWO decisions. The part
     * before the dot is the LAYERING, which fixes what a submodel is; the part
     * after it is the ENCODING, which fixes how an activity graph is written
     * into that submodel. The four combinations are "srvn.cs", "srvn.ph",
     * "flat.cs" and "flat.ph".</p>
     *
     * <p>"flat" is the ALIAS of "flat.cs" and resolves unconditionally rather
     * than probing "flat.ph", because a model is squashed in order to express
     * the routed call groups that only the routing encoding dispatches; naming
     * "flat.ph" is therefore not the same as naming "flat".</p>
     *
     * @throws Exception if the solver fails
     */
    public static void lqn_flatph() throws Exception {
        LayeredNetwork model = LayeredModel.lqn_flatph();

        String[] methods = {"srvn.cs", "srvn.ph", "flat.cs", "flat.ph"};
        for (String method : methods) {
            SolverOptions options = LN.defaultOptions();
            options.method = method;
            options.verbose = VerboseLevel.SILENT;
            LN solver = new LN(model, options);
            AvgTable avgTable = solver.getAvgTable();
            System.out.printf("%n--- method=%s (built as %s, %d submodel(s))%n",
                    method, solver.lnmethod, solver.getEnsemble().size());
            avgTable.setOptions(printing());
            avgTable.print();
        }
        pauseForUser();
    }

    /**
     * Method 'srvn.ph': an entry's activity graph as a phase-type law (lqn_srvnph.m).
     *
     * <p>The default method turns every activity graph into routing: one class
     * per task, entry, activity and call, plus Fork, Join, Router and
     * ClassSwitch nodes. Method "srvn.ph" composes each entry graph into a
     * single phase-type law by the exact series-parallel reduction, so a layer
     * becomes a two-station cycle, Delay plus Queue, with one closed class per
     * caller task. The sequencing survives as a distribution rather than as
     * routing.</p>
     *
     * @throws Exception if the solver fails
     */
    public static void lqn_srvnph() throws Exception {
        LayeredNetwork model = LayeredModel.lqn_srvnph();

        SolverOptions lnoptions = LN.defaultOptions();
        lnoptions.verbose = VerboseLevel.SILENT;
        SolverOptions solveroptions = MVA.defaultOptions();
        solveroptions.verbose = VerboseLevel.SILENT;

        long t0 = System.nanoTime();
        LN defaultSolver = new LN(model, m -> new MVA(m, solveroptions), lnoptions);
        AvgTable defaultTable = defaultSolver.getAvgTable();
        double tDefault = (System.nanoTime() - t0) / 1e9;
        System.out.printf("%nLN(default) Results [%.3f s]:%n", tDefault);
        defaultTable.setOptions(printing());
        defaultTable.print();

        SolverOptions lnoptions2 = LN.defaultOptions();
        lnoptions2.verbose = VerboseLevel.SILENT;
        // the method name; bare "srvnph" is unrecognised and falls back to srvn.cs
        lnoptions2.method = "srvn.ph";
        t0 = System.nanoTime();
        LN phSolver = new LN(model, m -> new MVA(m, solveroptions), lnoptions2);
        AvgTable phTable = phSolver.getAvgTable();
        double tSrvnph = (System.nanoTime() - t0) / 1e9;
        System.out.printf("%nLN(srvn.ph) Results [%.3f s]:%n", tSrvnph);
        phTable.setOptions(printing());
        phTable.print();

        int nclassesDefault = 0;
        for (Network layer : defaultSolver.getEnsemble()) {
            nclassesDefault += layer.getNumberOfClasses();
        }
        int nclassesSrvnph = 0;
        for (Network layer : phSolver.getEnsemble()) {
            nclassesSrvnph += layer.getNumberOfClasses();
        }
        System.out.printf(
                "%nLayer classes: default %d, srvn.ph %d. Runtime: %.3fs vs %.3fs (%.2fx).%n",
                nclassesDefault, nclassesSrvnph, tDefault, tSrvnph, tDefault / tSrvnph);
        pauseForUser();
    }

    /**
     * A processor whose servers are not interchangeable (lqn_server_pools).
     * <p>
     * Prints the fully-compatible pool, which is the neutral declaration, and
     * then the compatibility graph, under which neither task reaches more than
     * two of the three servers.
     */
    public static void lqn_server_pools() throws Exception {
        SolverOptions opt = LN.defaultOptions();
        // a compatibility declaration is a station rate law, which only the
        // class-switching layerings can carry
        opt.method = "srvn.cs";

        System.out.println("--- homogeneous pool on P1 ---");
        new LN(LayeredModel.lqn_server_pools(false), opt).getAvgTable().print();

        System.out.println("--- compatibility pool on P1 ---");
        new LN(LayeredModel.lqn_server_pools(true), opt).getAvgTable().print();

        pauseForUser();
    }

    public static void lqn_basic() throws Exception {
        LayeredNetwork model = LayeredModel.lqn_basic();

        // The reference solves this with the layered solver on its default layers.
        new LN(model).getAvgTable().print();

        pauseForUser();
    }
    
    /**
     * Serial layered network (lqn_serial.ipynb).
     * <p>
     * Shows serial activity precedence in layered networks
     * with synchronous call patterns.
     */
    public static void lqn_serial() throws Exception {
        LayeredNetwork model = LayeredModel.lqn_serial();
        
        // Create LQNS solver as in MATLAB (lines 32-34)
        SolverOptions options = LQNS.defaultOptions();
        options.keep = true;
        LQNS solver = new LQNS(model, options);
        AvgTable avgTable = solver.getLNAvgTable();
        avgTable.print();
        
        // Get raw avg tables as in MATLAB (lines 37-39)
        AvgTable[] rawTables = solver.getRawAvgTables();
        rawTables[0].print();
        if (rawTables.length > 1) {
            rawTables[1].print();
        }
        
        // LN as in MATLAB (line 41)
        LN solverLN = new LN(model, SolverType.MVA);
        AvgTable avgTableLN = solverLN.getAvgTable();
        avgTableLN.print();
        
        pauseForUser();
    }
    
    /**
     * Multi-solver layered network (lqn_multi_solvers.ipynb).
     * <p>
     * Demonstrates different solver approaches for layered networks
     * with infinite server capacity.
     */
    public static void lqn_multi_solvers() throws Exception {
        LayeredNetwork model = LayeredModel.lqn_multi_solvers();
        
        // LQNS solver as in MATLAB (lines 27-29)
        SolverOptions lqnsOptions = LQNS.defaultOptions();
        lqnsOptions.keep = true;
        lqnsOptions.verbose = VerboseLevel.STD;
        LQNS lqnsSolver = new LQNS(model, lqnsOptions);
        AvgTable avgTableLQNS = lqnsSolver.getLNAvgTable();
        avgTableLQNS.print();
        
        // LN with MVA solver in each layer (lines 37-39)
        SolverOptions lnOptions = LN.defaultOptions();
        lnOptions.verbose = VerboseLevel.SILENT;
        lnOptions.seed = 2300;
        SolverOptions mvaOptions = MVA.defaultOptions();
        mvaOptions.verbose = VerboseLevel.SILENT;
        LN solverLN_MVA = new LN(model, (subModel) -> new MVA(subModel, mvaOptions), lnOptions);
        AvgTable avgTableLN_MVA = solverLN_MVA.getAvgTable();
        avgTableLN_MVA.print();
        
        // LN with NC solver in each layer (lines 47-49)
        SolverOptions ncOptions = NC.defaultOptions();
        ncOptions.verbose = VerboseLevel.SILENT;
        LN solverLN_NC = new LN(model, (subModel) -> new NC(subModel, ncOptions), lnOptions);
        AvgTable avgTableLN_NC = solverLN_NC.getAvgTable();
        avgTableLN_NC.print();
        
        pauseForUser();
    }
    
    /**
     * Two-task layered network (lqn_twotasks.ipynb).
     * <p>
     * Shows interaction between multiple tasks with
     * synchronous call patterns.
     */
    public static void lqn_twotasks() throws Exception {
        LayeredNetwork model = LayeredModel.lqn_twotasks();

        try {
            new LQNS(model).getLNAvgTable().print();
        } catch (Exception e) {
            System.out.println("LQNS failed: " + e.getMessage());
        }
        // NC LAYERS, NOT THE DEFAULT: MVA layers and NC layers are different fixed
        // points, so the layer solver the reference pins is part of the golden.
        new LN(model, (subModel) -> new NC(subModel, "verbose", false)).getAvgTable().print();

        pauseForUser();
    }
    
    /**
     * BPMN-style layered network (lqn_bpmn.ipynb).
     * <p>
     * Shows fork-join patterns with OrFork, AndFork,
     * OrJoin, and AndJoin activity precedence.
     */
    public static void lqn_bpmn() throws Exception {
        LayeredNetwork model = LayeredModel.lqn_bpmn();
        
        LQNS solver = new LQNS(model);
        AvgTable avgTable = solver.getLNAvgTable();
        avgTable.print();
        
        pauseForUser();
    }
    
    /**
     * Layered network with a setup / delay-off task (lqn_setup.ipynb).
     * <p>
     * The servers of task F2 switch off when idle and pay a setup
     * time when a request reactivates them.
     */
    public static void lqn_setup() throws Exception {
        LayeredNetwork model = LayeredModel.lqn_setup();

        // The reference drives MVA layers, which is what its golden holds.
        new LN(model, SolverType.MVA).getAvgTable().print();

        pauseForUser();
    }
    
    /**
     * Workflow layered network (lqn_workflows.ipynb).
     * <p>
     * Demonstrates loop and fork-join precedence patterns
     * with nested synchronous calls.
     */
    public static void lqn_workflows() throws Exception {
        LayeredNetwork model = LayeredModel.lqn_workflows();

        try {
            new LQNS(model).getLNAvgTable().print();
        } catch (Exception e) {
            System.out.println("LQNS failed: " + e.getMessage());
        }
        new LN(model).getAvgTable().print();

        pauseForUser();
    }
    
    /**
     * OFBiz-style layered network (lqn_ofbiz.ipynb).
     * <p>
     * Enterprise application model inspired by Apache OFBiz
     * with database and application layers.
     */
    public static void lqn_ofbiz() throws Exception {
        LayeredNetwork model = LayeredModel.lqn_ofbiz();

        try {
            new LQNS(model).getLNAvgTable().print();
        } catch (Exception e) {
            System.out.println("LQNS failed: " + e.getMessage());
        }
        // NC layers: the golden keys this row LN(NC).
        new LN(model, (subModel) -> new NC(subModel, "verbose", false)).getAvgTable().print();

        pauseForUser();
    }

    /**
     * Sock Shop microservice layered network (lqn_sockshop).
     * <p>
     * Demonstrates a multi-tier microservice architecture with
     * processor replication, fan-in/fan-out, and PS scheduling.
     */
    public static void lqn_sockshop() throws Exception {
        LayeredNetwork model = LayeredModel.lqn_sockshop();

        LN solverLN = new LN(model, SolverType.MVA);
        AvgTable avgTable = solverLN.getAvgTable();
        avgTable.print();

        pauseForUser();
    }

    /**
     * Main method demonstrating selected layered network examples.
     */
    public static void main(String[] args) throws Exception {
        try {
            lqn_basic();
        } catch (Exception e) {
            System.err.println("lqn_basic failed: " + e.getMessage());
            e.printStackTrace();
        }
        
        try {
            lqn_serial();
        } catch (Exception e) {
            System.err.println("lqn_serial failed: " + e.getMessage());
            e.printStackTrace();
        }
        
        try {
            lqn_multi_solvers();
        } catch (Exception e) {
            System.err.println("lqn_multi_solvers failed: " + e.getMessage());
            e.printStackTrace();
        }
        
        try {
            lqn_twotasks();
        } catch (Exception e) {
            System.err.println("lqn_twotasks failed: " + e.getMessage());
            e.printStackTrace();
        }
        
        try {
            lqn_bpmn();
        } catch (Exception e) {
            System.err.println("lqn_bpmn failed: " + e.getMessage());
            e.printStackTrace();
        }
        
        try {
            lqn_setup();
        } catch (Exception e) {
            System.err.println("lqn_setup failed: " + e.getMessage());
            e.printStackTrace();
        }
        
        try {
            lqn_workflows();
        } catch (Exception e) {
            System.err.println("lqn_workflows failed: " + e.getMessage());
            e.printStackTrace();
        }
        
        try {
            lqn_ofbiz();
        } catch (Exception e) {
            System.err.println("lqn_ofbiz failed: " + e.getMessage());
            e.printStackTrace();
        }

        try {
            lqn_sockshop();
        } catch (Exception e) {
            System.err.println("lqn_sockshop failed: " + e.getMessage());
            e.printStackTrace();
        }

        try {
            lqn_rrobin();
        } catch (Exception e) {
            System.err.println("lqn_rrobin failed: " + e.getMessage());
            e.printStackTrace();
        }

        try {
            lqn_jsq();
        } catch (Exception e) {
            System.err.println("lqn_jsq failed: " + e.getMessage());
            e.printStackTrace();
        }

        try {
            lqn_flatph();
        } catch (Exception e) {
            System.err.println("lqn_flatph failed: " + e.getMessage());
            e.printStackTrace();
        }

        try {
            lqn_srvnph();
        } catch (Exception e) {
            System.err.println("lqn_srvnph failed: " + e.getMessage());
            e.printStackTrace();
        }

        try {
            lqn_init();
        } catch (Exception e) {
            System.err.println("lqn_init failed: " + e.getMessage());
            e.printStackTrace();
        }

        try {
            lqn_transient();
        } catch (Exception e) {
            System.err.println("lqn_transient failed: " + e.getMessage());
            e.printStackTrace();
        }

        try {
            lqn_bpmn_trace();
        } catch (Exception e) {
            System.err.println("lqn_bpmn_trace failed: " + e.getMessage());
            e.printStackTrace();
        }

        try {
            lqn_paramident();
        } catch (Exception e) {
            System.err.println("lqn_paramident failed: " + e.getMessage());
            e.printStackTrace();
        }

        scanner.close();
    }
}
