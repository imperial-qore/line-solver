package jline.solvers.ssa.analyzers;

import java.util.HashMap;
import java.util.HashSet;
import java.util.List;
import java.util.Map;
import java.util.Set;

import jline.api.sn.SnNonmarkovToPh;
import jline.io.InputOutput;
import jline.lang.NetworkStruct;
import jline.lang.constant.BalkingStrategy;
import jline.lang.constant.DropStrategy;
import jline.lang.constant.ImpatienceType;
import jline.lang.constant.NodeType;
import jline.lang.constant.ProcessType;
import jline.lang.constant.SchedStrategy;
import jline.lang.nodeparam.TransitionNodeParam;
import jline.lang.nodes.StatefulNode;
import jline.solvers.SolverOptions;
import jline.util.SerializableFunction;
import jline.solvers.ssa.SSAResult;
import jline.solvers.ssa.SolverSSA;
import jline.util.matrix.Matrix;

public final class Solver_ssa_analyzer {
    private Solver_ssa_analyzer() {}

    /**
     * Rejects marking-dependent transition firing rates.
     *
     * <p>Neither SSA engine (NRM or serial) applies the g(marking) multiplier set by
     * {@code Transition.setFiringRateDependence}, unlike SolverCTMC and SolverLDES, so
     * the model is refused rather than silently simulated at the nominal rate. Mirrors
     * matlab/src/solvers/SSA/solver_ssa_nrm.m and the native-Python
     * api/solvers/ssa/{nrm,serial}.py guards.
     *
     * @param sn network structure of the analyzed model
     */
    private static void assertNoFiringDependence(NetworkStruct sn) {
        for (int ind = 0; ind < sn.nnodes; ind++) {
            if (sn.nodetype.get(ind) != NodeType.Transition) {
                continue;
            }
            Object np = sn.nodeparam.get(sn.nodes.get(ind));
            if (!(np instanceof TransitionNodeParam)) {
                continue;
            }
            List<SerializableFunction<Matrix, Double>> firingdep = ((TransitionNodeParam) np).firingdep;
            if (firingdep == null) {
                continue;
            }
            for (int m = 0; m < firingdep.size(); m++) {
                if (firingdep.get(m) != null) {
                    throw new RuntimeException(String.format(
                            "Transition %s mode %d uses a marking-dependent firing rate "
                            + "(setFiringRateDependence), which SolverSSA does not support; "
                            + "use SolverCTMC or SolverLDES.", sn.nodenames.get(ind), m + 1));
                }
            }
        }
    }

    public static SSAResult solver_ssa_analyzer(NetworkStruct snInput, SolverOptions options, SolverSSA solverSSA) {
        // see _kb/06-solver-catalog.md for rationale
        // SSA draws a sample path, so the surrogate must be a genuine phase-type
        String phfit0 = options.config.phfit;
        options.config.phfit = "ph";
        NetworkStruct sn = SnNonmarkovToPh.snNonmarkovToPh(snInput, options, false);
        options.config.phfit = phfit0;
        long Tstart = System.nanoTime();
        Map<StatefulNode, Matrix> init_state = new HashMap<StatefulNode, Matrix>();
        for (StatefulNode statefulNode : sn.state.keySet()) {
            if (sn.isfjaugmented) {
                // see _kb/06-solver-catalog.md for rationale
                init_state.put(statefulNode, sn.state.get(statefulNode));
                continue;
            }
            int node = statefulNode.getStatefulIndex();
            if (node == -1) {
                String nodeName = statefulNode.getName() != null ? statefulNode.getName() : "unknown";
                String nodeType = statefulNode.getClass().getSimpleName();
                Integer nodeByName = null;
                for (int idx = 0; idx < solverSSA.model.getStatefulNodes().size(); idx++) {
                    StatefulNode sf = solverSSA.model.getStatefulNodes().get(idx);
                    if (sf != null && nodeName.equals(sf.getName())) { nodeByName = idx; break; }
                }
                if (nodeByName != null) {
                    init_state.put(solverSSA.model.getStatefulNodes().get(nodeByName), sn.state.get(statefulNode));
                } else {
                    Integer nodeByType = null;
                    for (int idx = 0; idx < solverSSA.model.getStatefulNodes().size(); idx++) {
                        StatefulNode sf = solverSSA.model.getStatefulNodes().get(idx);
                        if (sf != null && nodeType.equals(sf.getClass().getSimpleName())) { nodeByType = idx; break; }
                    }
                    if (nodeByType != null) {
                        init_state.put(solverSSA.model.getStatefulNodes().get(nodeByType), sn.state.get(statefulNode));
                    } else {
                        if (nodeType.contains("Cache")) {
                            InputOutput.line_warning("solver_ssa_analyzer",
                                    "Skipping Cache node '%s' with invalid stateful index. This may be expected for cache nodes with class switching.",
                                    nodeName);
                        } else {
                            throw new RuntimeException("solver_ssa_analyzer: StatefulNode '" + nodeName + "' (type: " + nodeType + ") has invalid stateful index -1.");
                        }
                    }
                }
            } else {
                init_state.put(solverSSA.model.getStatefulNodes().get(node), sn.state.get(statefulNode));
            }
        }

        SSAResult res = new SSAResult();
        String method;
        SolverOptions actualOptions = options;
        boolean isDefaultMethod = false;

        // see _kb/06-solver-catalog.md for rationale
        boolean isSPN = false;
        for (NodeType nt : sn.nodetype) {
            if (nt == NodeType.Transition) { isSPN = true; break; }
        }
        if (isSPN) {
            assertNoFiringDependence(sn);
            boolean spnExplicitSerial = "serial".equals(actualOptions.method)
                    || "ssa".equals(actualOptions.method)
                    || !gdNrmOK(sn);
            if (spnExplicitSerial) {
                // see _kb/06-solver-catalog.md for rationale
                actualOptions.method = "serial";
            } else {
                boolean spnDefault = "default".equals(actualOptions.method);
                actualOptions.method = "nrm";
                res = Solver_ssa_analyzer_nrm.solver_ssa_analyzer_nrm(sn, init_state, actualOptions);
                res.method = spnDefault ? "default/nrm" : "nrm";
                res.runtime = (System.nanoTime() - Tstart) / 1.0e9;
                return res;
            }
        }

        if ("default".equals(actualOptions.method)) {
            boolean nrmSupported = nrmEligibleReason(sn).isEmpty();

            if (nrmSupported) {
                actualOptions.method = "nrm";
                res = Solver_ssa_analyzer_nrm.solver_ssa_analyzer_nrm(sn, init_state, actualOptions);
                res.method = "default/nrm";
                res.runtime = (System.nanoTime() - Tstart) / 1.0e9;
                return res;
            } else {
                actualOptions = options.copy();
                actualOptions.method = "serial";
                isDefaultMethod = true;
            }
        } else if ("nrm".equals(actualOptions.method)) {
            if (!renegeNrmOK(sn)) {
                // Phase-type patience would need the remaining-patience phase of
                // each waiting job, which the reaction network does not carry
                InputOutput.line_warning(InputOutput.mfilename(new Object() {}),
                        "NRM supports only exponential (memoryless) patience for reneging; falling back to the serial method.");
                actualOptions = options.copy();
                actualOptions.method = "serial";
            } else if (!gdNrmOK(sn)) {
                // A global (Whittle) dependence reads the whole population
                // matrix, which the per-station propensity closures do not get
                InputOutput.line_warning(InputOutput.mfilename(new Object() {}),
                        "NRM does not support a global dependence (setGlobalDependence); falling back to the serial method.");
                actualOptions = options.copy();
                actualOptions.method = "serial";
            } else if (!balkNrmOK(sn)) {
                // EXPECTED_WAIT / COMBINED balking depend on the mean waiting
                // time, which is not a function of the state vector
                InputOutput.line_warning(InputOutput.mfilename(new Object() {}),
                        "NRM only supports QUEUE_LENGTH balking; falling back to the serial method.");
                actualOptions = options.copy();
                actualOptions.method = "serial";
            } else {
                // Cache nodes (including the retrieval/delayed-hit system) are
                // fully supported by the NRM engine; see _kb/06-solver-catalog.md
                res = Solver_ssa_analyzer_nrm.solver_ssa_analyzer_nrm(sn, init_state, actualOptions);
                res.method = "nrm";
                res.runtime = (System.nanoTime() - Tstart) / 1.0e9;
                return res;
            }
        } else if ("ssa".equals(actualOptions.method)) {
            actualOptions = options.copy();
            actualOptions.method = "serial";
            isDefaultMethod = true;
        }

        // 'ssa' is the reference's alias for the serial engine, NOT for the NRM:
        // solver_ssa_analyzer.m:139 sets options.method='serial' under it. It was
        // advertised by SolverSSA.listValidMethods and reached the "Unknown
        // analysis method" arm below, so a listed name always threw.
        if ("ssa".equals(actualOptions.method)) {
            actualOptions.method = "serial";
        }
        if ("serial".equals(actualOptions.method)) {
            res = Solver_ssa_analyzer_serial.solver_ssa_analyzer_serial(sn, false, init_state, actualOptions, solverSSA);
            method = isDefaultMethod ? "default/serial" : "serial";
        } else if ("para".equals(actualOptions.method) || "parallel".equals(actualOptions.method)
                || "ssa.parallel".equals(actualOptions.method)) {
            try {
                res = Solver_ssa_analyzer_parallel.solver_ssa_analyzer_parallel(sn, init_state, actualOptions, solverSSA);
                method = isDefaultMethod ? "default/parallel" : "parallel";
            } catch (Exception e) {
                System.out.println("Parallel execution failed - falling back to serial SSA.");
                res = Solver_ssa_analyzer_serial.solver_ssa_analyzer_serial(sn, true, init_state, actualOptions, solverSSA);
                method = isDefaultMethod ? "default/serial" : "serial";
            }
        } else {
            throw new RuntimeException("solver_ssa_analyzer:UnknownMethod - Unknown analysis method: " + actualOptions.method);
        }

        res.method = method;
        res.runtime = (System.nanoTime() - Tstart) / 1.0e9;
        if (jline.solvers.Solver.timeExceeded(Tstart, options.timeout)) {
            res.timedOut = true;
            InputOutput.line_warning("Solver_ssa_analyzer",
                    "Solver stopped after the wall-clock time budget (options.timeout=" + options.timeout + "s) was exceeded; returning the interim solution.");
        }

        if (options.confint > 0 && res.tranSysState != null && res.tranSysState.size() > 1) {
            computeBatchMeansCI(res, sn, options.confint, options);
        }
        return res;
    }

    private static void computeBatchMeansCI(SSAResult result, NetworkStruct sn, double confintLevel, SolverOptions options) {
        int M = sn.nstations;
        int K = sn.nclasses;
        Matrix QNCI = new Matrix(M, K); QNCI.fill(0.0);
        Matrix UNCI = new Matrix(M, K); UNCI.fill(0.0);
        Matrix RNCI = new Matrix(M, K); RNCI.fill(0.0);
        Matrix TNCI = new Matrix(M, K); TNCI.fill(0.0);
        Matrix ANCI = new Matrix(M, K); ANCI.fill(0.0);
        Matrix WNCI = new Matrix(M, K); WNCI.fill(0.0);

        Matrix times = result.tranSysState.get(0);
        if (times == null) return;
        int nSamples = times.getNumRows();
        if (nSamples < 20) {
            result.QNCI = QNCI; result.UNCI = UNCI; result.RNCI = RNCI;
            result.TNCI = TNCI; result.ANCI = ANCI; result.WNCI = WNCI;
            return;
        }
        int numBatches = Math.min(20, nSamples / 10);
        if (numBatches < 2) {
            result.QNCI = QNCI; result.UNCI = UNCI; result.RNCI = RNCI;
            result.TNCI = TNCI; result.ANCI = ANCI; result.WNCI = WNCI;
            return;
        }
        int batchSize = nSamples / numBatches;
        // Discard initial transient before batch means: use
        // options.config.warmupfrac when set (> 0), else the legacy 10% discard
        double warmupFrac = (options != null && options.config != null
                && options.config.warmupfrac != null && options.config.warmupfrac > 0)
                ? options.config.warmupfrac : 0.1;
        int transientCutoff = Math.max(1, (int) (nSamples * warmupFrac));

        for (int ist = 0; ist < M; ist++) {
            int isf = (int) sn.stationToStateful.get(ist);
            if (isf >= 0 && isf < result.tranSysState.size() - 1) {
                Matrix stateData = result.tranSysState.get(1 + isf);
                if (stateData == null || stateData.isEmpty()) continue;
                for (int k = 0; k < K; k++) {
                    double[] batchMeans = new double[numBatches];
                    int validBatches = 0;
                    for (int b = 0; b < numBatches; b++) {
                        int startIdx = transientCutoff + b * batchSize;
                        int endIdx = Math.min(transientCutoff + (b + 1) * batchSize, nSamples);
                        if (startIdx >= nSamples || startIdx >= endIdx) continue;
                        double weightedSum = 0.0;
                        double totalTime = 0.0;
                        for (int i = startIdx; i < endIdx; i++) {
                            double dt = (i > 0) ? times.get(i, 0) - times.get(i - 1, 0) : times.get(i, 0);
                            double qLen;
                            if (stateData.getNumCols() > k) qLen = stateData.get(i, k);
                            else {
                                double sum = 0.0;
                                for (int c = 0; c < stateData.getNumCols(); c++) sum += stateData.get(i, c);
                                qLen = sum;
                            }
                            weightedSum += qLen * dt;
                            totalTime += dt;
                        }
                        if (totalTime > 0) {
                            batchMeans[validBatches] = weightedSum / totalTime;
                            validBatches++;
                        }
                    }
                    if (validBatches >= 2) {
                        double batchMean = 0.0;
                        for (int i = 0; i < validBatches; i++) batchMean += batchMeans[i];
                        batchMean /= validBatches;
                        double variance = 0.0;
                        for (int i = 0; i < validBatches; i++) {
                            double diff = batchMeans[i] - batchMean;
                            variance += diff * diff;
                        }
                        variance /= (validBatches - 1);
                        double stdErr = Math.sqrt(variance / validBatches);
                        double alpha = 1.0 - confintLevel;
                        double tCrit = getTCriticalValue(alpha, validBatches - 1);
                        QNCI.set(ist, k, tCrit * stdErr);
                    }
                }
            }
        }
        for (int i = 0; i < M; i++) {
            for (int k = 0; k < K; k++) {
                UNCI.set(i, k, QNCI.get(i, k));
                RNCI.set(i, k, QNCI.get(i, k));
                TNCI.set(i, k, QNCI.get(i, k));
            }
        }
        result.QNCI = QNCI; result.UNCI = UNCI; result.RNCI = RNCI;
        result.TNCI = TNCI; result.ANCI = ANCI; result.WNCI = WNCI;
    }

    private static double getTCriticalValue(double alpha, int df) {
        double halfAlpha = alpha / 2.0;
        if (df > 30) {
            if (halfAlpha <= 0.005) return 2.576;
            if (halfAlpha <= 0.025) return 1.96;
            if (halfAlpha <= 0.05) return 1.645;
            return 1.28;
        }
        Map<Integer, Double> t95 = new HashMap<Integer, Double>();
        t95.put(1, 12.71); t95.put(2, 4.30); t95.put(3, 3.18); t95.put(4, 2.78); t95.put(5, 2.57);
        t95.put(6, 2.45); t95.put(7, 2.36); t95.put(8, 2.31); t95.put(9, 2.26); t95.put(10, 2.23);
        t95.put(11, 2.20); t95.put(12, 2.18); t95.put(13, 2.16); t95.put(14, 2.14); t95.put(15, 2.13);
        t95.put(16, 2.12); t95.put(17, 2.11); t95.put(18, 2.10); t95.put(19, 2.09); t95.put(20, 2.09);
        t95.put(25, 2.06); t95.put(30, 2.04);
        if (halfAlpha >= 0.02 && halfAlpha <= 0.03) {
            Double v = t95.get(df);
            if (v != null) return v;
            v = t95.get(Math.min(df, 30));
            return v != null ? v : 2.0;
        }
        Double base = t95.get(df);
        if (base == null) base = t95.get(Math.min(df, 30));
        double baseT = (base != null) ? base : 2.0;
        if (halfAlpha <= 0.005) return baseT * 1.32;
        if (halfAlpha <= 0.05) return baseT * 0.84;
        return baseT * 0.65;
    }

    private static boolean isPopulationModel(NetworkStruct sn) {
        double totalJobs = sn.njobs.sumCols().get(0, 0);
        if (totalJobs <= 0) return false;
        for (int i = 0; i < sn.njobs.getNumCols(); i++) {
            if (sn.njobs.get(0, i) < 0) return false;
        }
        return true;
    }

    /**
     * True: finite capacity regions are supported under both rules. DROP is
     * reproduced by censoring the refused transition; WAITQ parks the refused job
     * in a per-region FIFO carried explicitly by the NRM engine
     * (Solver_ssa_nrm.fcrReleaseCascade) and admits it head-of-line as capacity
     * frees. No region rule forces a serial fallback.
     */
    /**
     * True when no station uses a balking strategy the NRM cannot evaluate.
     * QUEUE_LENGTH is a pure function of the state vector, so the NRM draws it at
     * firing time; EXPECTED_WAIT and COMBINED depend on the mean waiting time and
     * need the serial engine (State.afterEventStation rejects them likewise).
     */
    /**
     * Whether the NRM engine can run this model, as a SENTENCE.
     *
     * <p>ONE PREDICATE, TWO CALLERS. The dispatch above asks it to decide whether
     * to PREFER the NRM on the {@code default} path; {@code SolverSSA
     * .supportsModelMethod} asks it to decide whether {@code nrm} may be OFFERED
     * at all. Inlined in the dispatch it could answer only the first, so
     * {@code findSolver} reported {@code ssa.nrm} runnable on every model and an
     * explicit {@code SolverSSA(model,"nrm")} then raised from
     * {@code Solver_ssa_analyzer_nrm} on the very models this test excludes.</p>
     *
     * @param sn the network structure
     * @return empty string when the NRM can run it, else the reason it cannot
     */
    public static String nrmMethodRefusal(NetworkStruct sn) {
        // NOT nrmEligibleReason, and the difference is why there are two.
        // nrmEligibleReason asks whether the NRM should be PREFERRED on the
        // 'default' path and is deliberately wide. This asks what a GATE must
        // ask: will the name the caller typed produce an answer?
        //
        // The explicit "nrm" arm above FALLS BACK to the serial engine, with a
        // warning, for reneging patience, global dependence and balking, so a
        // model carrying any of those still gets a correct answer under 'nrm'
        // and must not be refused here. What raises instead of falling back is
        // exactly two things, both in Solver_ssa_analyzer_nrm: a scheduling
        // policy the reaction network has no form for, and phase-type service at
        // a station whose in-service multiset it does not record.
        //
        // MATLAB gates on the scheduling rule alone (its NRM analyzer has no
        // phase check) and C++ on all six (its NRM has no fallback arm). Each
        // gate states its own engine's reach; that is not a divergence.
        Set<SchedStrategy> allowedSched = nrmAllowedSched();
        for (SchedStrategy s : sn.sched.values()) {
            if (!allowedSched.contains(s)) {
                return "The 'nrm' engine has no reaction form for the " + s
                        + " scheduling policy. Use options.method='serial'.";
            }
        }
        if (!phaseNrmOK(sn)) {
            return "The 'nrm' engine expands phase-type service only at INF/PS and "
                    + "non-preemptive buffered stations. Use options.method='serial'.";
        }
        return "";
    }

    /** The scheduling policies the NRM has a reaction form for. */
    private static Set<SchedStrategy> nrmAllowedSched() {
        Set<SchedStrategy> allowedSched = new HashSet<SchedStrategy>();
        allowedSched.add(SchedStrategy.INF);
        allowedSched.add(SchedStrategy.EXT);
        allowedSched.add(SchedStrategy.PS);
        allowedSched.add(SchedStrategy.LPS);
        allowedSched.add(SchedStrategy.DPS);
        allowedSched.add(SchedStrategy.GPS);
        allowedSched.add(SchedStrategy.PSPRIO);
        allowedSched.add(SchedStrategy.DPSPRIO);
        allowedSched.add(SchedStrategy.GPSPRIO);
        allowedSched.add(SchedStrategy.SIRO);
        allowedSched.add(SchedStrategy.HOL);
        allowedSched.add(SchedStrategy.SEPT);
        allowedSched.add(SchedStrategy.LEPT);
        allowedSched.add(SchedStrategy.FCFS);
        allowedSched.add(SchedStrategy.LCFS);
        allowedSched.add(SchedStrategy.LCFSPR);
        allowedSched.add(SchedStrategy.PAS);
        allowedSched.add(SchedStrategy.POLLING);
        return allowedSched;
    }

    public static String nrmEligibleReason(NetworkStruct sn) {
        Set<SchedStrategy> allowedSched = nrmAllowedSched();
        for (SchedStrategy s : sn.sched.values()) {
            if (!allowedSched.contains(s)) {
                return "The 'nrm' engine has no reaction form for the " + s
                        + " discipline. Use options.method='serial'.";
            }
        }
        // see _kb/06-solver-catalog.md for rationale
        if (!cacheNrmOK(sn)) {
            return "The 'nrm' engine cannot reproduce this cache configuration. "
                    + "Use options.method='serial'.";
        }
        // A global (Whittle) dependence reads the whole population matrix, which
        // the per-station propensity closures do not receive.
        if (!gdNrmOK(sn)) {
            return "The 'nrm' engine builds each reaction propensity from one station's "
                    + "population slice, so it cannot evaluate a global (Whittle) scaling "
                    + "handle. Use options.method='serial'.";
        }
        // see _kb/06-solver-catalog.md for rationale
        for (NodeType nt : sn.nodetype) {
            if (nt == NodeType.Fork || nt == NodeType.Join) {
                return "The 'nrm' engine has no Fork/Join node handling. "
                        + "Use options.method='serial'.";
            }
        }
        // see _kb/06-solver-catalog.md for rationale
        if (!fcrNrmOK(sn)) {
            return "The 'nrm' engine cannot evaluate this finite capacity region rule. "
                    + "Use options.method='serial'.";
        }
        // Only the state-based QUEUE_LENGTH balking strategy is a function of the
        // state vector; reneging needs memoryless patience.
        if (!balkNrmOK(sn)) {
            return "The 'nrm' engine evaluates only QUEUE_LENGTH balking, which is a function "
                    + "of the state vector; EXPECTED_WAIT and COMBINED need the mean waiting "
                    + "time. Use options.method='serial'.";
        }
        if (!renegeNrmOK(sn)) {
            return "The 'nrm' engine abandons at the aggregate rate (waiting count)*mu, which "
                    + "is exact only for memoryless patience. Use options.method='serial'.";
        }
        // see _kb/06-solver-catalog.md for rationale
        if (!phaseNrmOK(sn)) {
            return "The 'nrm' engine expands non-exponential service exactly only where every "
                    + "job present is in service; this model buffers it under a discipline "
                    + "whose in-service phase multiset it does not record. "
                    + "Use options.method='serial'.";
        }
        return "";
    }

    private static boolean balkNrmOK(NetworkStruct sn) {
        if (sn.balkingStrategy == null || sn.balkingStrategy.isEmpty()) {
            return true;
        }
        for (Map<jline.lang.JobClass, BalkingStrategy> smap : sn.balkingStrategy.values()) {
            if (smap == null) continue;
            for (BalkingStrategy bs : smap.values()) {
                if (bs != null && bs != BalkingStrategy.QUEUE_LENGTH) {
                    return false;
                }
            }
        }
        return true;
    }

    /**
     * True when no station reneges with non-exponential patience. The NRM abandons
     * at the aggregate rate (waiting count)*mu, which is only correct when patience
     * is memoryless; phase-type patience would need each waiting job's remaining
     * phase. The serial engine rejects the same combination outright.
     */
    private static boolean renegeNrmOK(NetworkStruct sn) {
        if (sn.impatienceClass == null || sn.impatienceClass.isEmpty()) {
            return true;
        }
        for (Map.Entry<jline.lang.nodes.Station, Map<jline.lang.JobClass, ImpatienceType>> e
                : sn.impatienceClass.entrySet()) {
            Map<jline.lang.JobClass, ImpatienceType> cmap = e.getValue();
            if (cmap == null) continue;
            for (Map.Entry<jline.lang.JobClass, ImpatienceType> ce : cmap.entrySet()) {
                if (ce.getValue() != ImpatienceType.RENEGING) continue;
                Map<jline.lang.JobClass, ProcessType> tmap =
                        (sn.impatienceType != null) ? sn.impatienceType.get(e.getKey()) : null;
                ProcessType pt = (tmap != null) ? tmap.get(ce.getKey()) : null;
                if (pt != null && pt != ProcessType.EXP) {
                    return false;
                }
            }
        }
        return true;
    }

    /**
     * True when every non-exponential service sits at a station whose rate law the
     * NRM expands exactly. The INF/PS family (INF/EXT/PS/LPS/DPS/GPS and their PRIO
     * variants) needs only the per-phase populations, since every job present is in
     * service. The non-preemptive buffered family (FCFS/LCFS/SIRO/HOL/SEPT/LEPT)
     * tracks the phases of the jobs actually in service in the auxiliary multiset
     * svcph. A preemptive (LCFSPR) or polling station with phase-type service still
     * needs the serial engine.
     */
    static boolean phaseNrmOK(NetworkStruct sn) {
        Set<SchedStrategy> exact = new HashSet<SchedStrategy>();
        exact.add(SchedStrategy.INF);
        // see _kb/06-solver-catalog.md for rationale
        exact.add(SchedStrategy.PS);
        exact.add(SchedStrategy.LPS);
        exact.add(SchedStrategy.DPS);
        exact.add(SchedStrategy.GPS);
        exact.add(SchedStrategy.PSPRIO);
        exact.add(SchedStrategy.DPSPRIO);
        exact.add(SchedStrategy.GPSPRIO);
        exact.add(SchedStrategy.FCFS);
        exact.add(SchedStrategy.LCFS);
        exact.add(SchedStrategy.SIRO);
        exact.add(SchedStrategy.HOL);
        exact.add(SchedStrategy.SEPT);
        exact.add(SchedStrategy.LEPT);
        for (Map.Entry<jline.lang.nodes.Station, Map<jline.lang.JobClass, ProcessType>> e
                : sn.procid.entrySet()) {
            SchedStrategy sched = sn.sched.get(e.getKey());
            Map<jline.lang.JobClass, ProcessType> cmap = e.getValue();
            if (cmap == null) continue;
            for (ProcessType pt : cmap.values()) {
                if (pt != ProcessType.DISABLED && pt != ProcessType.EXP && !exact.contains(sched)) {
                    return false;
                }
            }
        }
        return true;
    }

    private static boolean fcrNrmOK(NetworkStruct sn) {
        // see _kb/06-solver-catalog.md for rationale
        return true;
    }

    /**
     * True when no global (Whittle) dependence is declared. The NRM builds one
     * propensity closure per reaction from the per-station population slice; a
     * global handle reads the whole population matrix, which that closure does
     * not receive. The serial engine carries the factor, so the model runs there.
     */
    private static boolean gdNrmOK(NetworkStruct sn) {
        return sn.gdscaling == null;
    }

    private static boolean cacheNrmOK(NetworkStruct sn) {
        // see _kb/06-solver-catalog.md for rationale
        return true;
    }
}
