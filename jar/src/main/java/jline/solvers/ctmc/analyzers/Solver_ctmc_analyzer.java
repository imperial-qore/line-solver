package jline.solvers.ctmc.analyzers;

import jline.api.mam.*;
import jline.api.mc.Ctmc_makeinfgen;
import jline.api.mc.Ctmc_solve;
import jline.api.pfqn.ld.CdPeakScaling;
import jline.api.sn.SnNonmarkovToPh;
import jline.io.InputOutput;
import jline.io.InputOutput;
import jline.lang.NetworkStruct;
import jline.lang.JobClass;
import jline.lang.constant.DropStrategy;
import jline.lang.constant.NodeType;
import jline.lang.constant.SchedStrategy;
import jline.lang.nodeparam.CacheNodeParam;
import jline.lang.state.ToMarginal;
import jline.solvers.SolverOptions;
import jline.solvers.ctmc.SolverCTMC;
import jline.solvers.ctmc.handlers.Solver_ctmc;
import jline.util.Maths;
import jline.util.MatFileUtils;
import jline.util.matrix.Matrix;
import jline.util.matrix.MatrixCell;

import java.util.ArrayList;
import java.util.List;
import java.util.Map;

public class Solver_ctmc_analyzer {

    private final SolverCTMC solverCTMC;

    public Solver_ctmc_analyzer(SolverCTMC solverCTMC) {
        this.solverCTMC = solverCTMC;
    }

    public static SolverCTMC.AnalyzerResult solver_ctmc_analyzer(NetworkStruct snInput, SolverOptions options) {
        // see _kb/06-solver-catalog.md for rationale
        NetworkStruct sn = SnNonmarkovToPh.snNonmarkovToPh(snInput, options, false);
        int M = sn.nstations;
        int K = sn.nclasses;
        Matrix S = sn.nservers;
        Matrix NK = sn.njobs.transpose();
        java.util.Map<jline.lang.nodes.Station, jline.lang.constant.SchedStrategy> schedid = sn.sched;
        long Tstart = System.nanoTime();
        Object PH = sn.proc;
        // sn.proc is stored as Object (erased); this cast to its concrete nested-map
        // type is inherently unchecked. Hoist it once and reuse.
        @SuppressWarnings("unchecked")
        java.util.Map<Object, java.util.Map<Object, MatrixCell>> PHmap =
                (java.util.Map<Object, java.util.Map<Object, MatrixCell>>) PH;

        jline.solvers.ctmc.ResultCTMC solverCTMCResult = Solver_ctmc.solver_ctmc(sn, options);
        Matrix InfGen = solverCTMCResult.getQ();
        Matrix StateSpace = solverCTMCResult.getStateSpace();
        Matrix StateSpaceAggr = solverCTMCResult.getStateSpaceAggr();
        jline.util.matrix.MatrixCell EventFiltration = solverCTMCResult.getDfilt();
        double[][][] arvRates = solverCTMCResult.getArvRates();
        double[][][] depRates = solverCTMCResult.getDepRates();
        sn = solverCTMCResult.getSn();

        for (int isf = 0; isf < sn.nstateful; isf++) {
            if (sn.state.get(sn.stateful.get(isf)).getNumCols() < sn.space.get(sn.stateful.get(isf)).getNumCols()) {
                Matrix state_matrix = new Matrix(1, sn.space.get(sn.stateful.get(isf)).getNumCols());
                state_matrix.zero();

                int startIdx = sn.space.get(sn.stateful.get(isf)).getNumCols() - sn.state.get(sn.stateful.get(isf)).getNumCols();
                int endIdx = state_matrix.getNumCols();
                for (int col = startIdx; col < endIdx; col++) {
                    state_matrix.set(0, col, sn.state.get(sn.stateful.get(isf)).get(col - startIdx));
                }
                sn.state.replace(sn.stateful.get(isf), state_matrix);
            }
        }

        NetworkStruct sncopy = sn;
        String fname = "";
        if (options.keep) {
            try {
                MatFileUtils.ensureWorkspaceDirectoryExists();
                fname = MatFileUtils.genFilename("workspace");
                MatFileUtils.saveCTMCWorkspace(StateSpace, InfGen, null, fname);
            } catch (Exception e) {
                InputOutput.line_warning("solver_ctmc_analyzer", "Could not save workspace to .mat file: %s", e.getMessage());
                fname = "";
            }
        } else {
            fname = "";
        }

        // wsetMap maps new indices to old indices; identity mapping unless reducible
        int[] wsetMap = null;
        Matrix wset = new Matrix(1, InfGen.length());
        for (int col = 0; col < wset.getNumCols(); col++) {
            wset.set(0, col, col);
        }

        // see _kb/06-solver-catalog.md for rationale
        Matrix probSysState = null;
        Matrix StateSpaceWork = StateSpace;
        Matrix StateSpaceAggrWork = StateSpaceAggr;
        // see _kb/06-solver-catalog.md for rationale
        Matrix InfGenWork = InfGen;

        // Detect connected components via symmetrized adjacency
        Matrix Bsym = InfGen.add(1.0, InfGen.transpose());
        Bsym.absEq();
        for (int i = 0; i < Bsym.getNumRows(); i++) {
            for (int j = 0; j < Bsym.getNumCols(); j++) {
                if (Bsym.get(i, j) > 0) Bsym.set(i, j, 1.0);
            }
        }
        java.util.Set<java.util.Set<Integer>> componentSets = Matrix.weaklyConnect(Bsym, null);
        List<List<Integer>> components = new java.util.ArrayList<List<Integer>>();
        for (java.util.Set<Integer> comp : componentSets) {
            components.add(new java.util.ArrayList<Integer>(comp));
        }
        int nConnComp = components.size();

        // see _kb/06-solver-catalog.md for rationale
        boolean pasModel = sn.sched != null && (sn.sched.containsValue(SchedStrategy.PAS) || sn.sched.containsValue(SchedStrategy.OI));
        if (pasModel) {
            List<Double> s0p = new ArrayList<Double>();
            boolean s0ok = true;
            for (int isf = 0; isf < sn.nstateful; isf++) {
                Matrix sm = sn.state.get(sn.stateful.get(isf));
                if (sm == null) { s0ok = false; break; }
                for (int ri = 0; ri < sm.getNumRows(); ri++) {
                    for (int ci = 0; ci < sm.getNumCols(); ci++) {
                        s0p.add(sm.get(ri, ci));
                    }
                }
            }
            int initPas = -1;
            if (s0ok) {
                Matrix s0m = new Matrix(1, s0p.size());
                for (int i = 0; i < s0p.size(); i++) s0m.set(0, i, s0p.get(i));
                initPas = Matrix.matchrow(StateSpace, s0m);
            }
            if (initPas >= 0) {
                List<Integer> reach = forwardReachable(InfGen, initPas);
                if (reach.size() < InfGen.length()) {
                    java.util.Set<Integer> reachSet = new java.util.HashSet<Integer>(reach);
                    List<List<Integer>> twoComp = new java.util.ArrayList<List<Integer>>();
                    twoComp.add(new java.util.ArrayList<Integer>(reach));
                    List<Integer> rest = new java.util.ArrayList<Integer>();
                    for (int i = 0; i < InfGen.length(); i++) {
                        if (!reachSet.contains(i)) rest.add(i);
                    }
                    if (!rest.isEmpty()) twoComp.add(rest);
                    components = twoComp;
                    nConnComp = components.size();
                }
            }
        }

        if (nConnComp > 1) {
            InputOutput.line_debug(options.verbose, String.format("CTMC is reducible: %d connected components", nConnComp));

            List<Double> s0parts = new ArrayList<Double>();
            for (int isf = 0; isf < sn.nstateful; isf++) {
                Matrix stateMatrix = sn.state.get(sn.stateful.get(isf));
                if (stateMatrix != null) {
                    for (int ri = 0; ri < stateMatrix.getNumRows(); ri++) {
                        for (int ci = 0; ci < stateMatrix.getNumCols(); ci++) {
                            s0parts.add(stateMatrix.get(ri, ci));
                        }
                    }
                }
            }
            Matrix s0 = new Matrix(1, s0parts.size());
            for (int i = 0; i < s0parts.size(); i++) {
                s0.set(0, i, s0parts.get(i));
            }
            int initStateIdx = Matrix.matchrow(StateSpace, s0);

            int[] connComp = new int[InfGen.length()];
            int compId = 0;
            for (List<Integer> comp : components) {
                for (int idx : comp) {
                    connComp[idx] = compId;
                }
                compId++;
            }

            List<Integer> wsetList;
            if (initStateIdx < 0) {
                int[] compSizes = new int[nConnComp];
                for (int idx = 0; idx < connComp.length; idx++) {
                    compSizes[connComp[idx]]++;
                }
                int largestComp = 0;
                int largestSize = 0;
                for (int ci = 0; ci < compSizes.length; ci++) {
                    if (compSizes[ci] > largestSize) {
                        largestSize = compSizes[ci];
                        largestComp = ci;
                    }
                }
                wsetList = new ArrayList<Integer>();
                for (int idx = 0; idx < connComp.length; idx++) {
                    if (connComp[idx] == largestComp) wsetList.add(idx);
                }
                InputOutput.line_debug(options.verbose, String.format("Using largest component with %d states (initial state removed by stochcomp)", wsetList.size()));
            } else {
                int targetComp = connComp[initStateIdx];
                wsetList = new ArrayList<Integer>();
                for (int idx = 0; idx < connComp.length; idx++) {
                    if (connComp[idx] == targetComp) wsetList.add(idx);
                }
                InputOutput.line_debug(options.verbose, String.format("Using component %d with %d states (from initial state)", targetComp, wsetList.size()));
            }

            int nw = wsetList.size();
            Matrix InfGenSub = new Matrix(nw, nw);
            for (int i = 0; i < nw; i++) {
                for (int j = 0; j < nw; j++) {
                    InfGenSub.set(i, j, InfGen.get(wsetList.get(i), wsetList.get(j)));
                }
            }
            probSysState = Ctmc_solve.ctmc_solve(InfGenSub);
            InfGenWork = InfGenSub;

            Matrix StateSpaceSub = new Matrix(nw, StateSpace.getNumCols());
            for (int i = 0; i < nw; i++) {
                for (int j = 0; j < StateSpace.getNumCols(); j++) {
                    StateSpaceSub.set(i, j, StateSpace.get(wsetList.get(i), j));
                }
            }
            StateSpaceWork = StateSpaceSub;

            Matrix StateSpaceAggrSub = new Matrix(nw, StateSpaceAggr.getNumCols());
            for (int i = 0; i < nw; i++) {
                for (int j = 0; j < StateSpaceAggr.getNumCols(); j++) {
                    StateSpaceAggrSub.set(i, j, StateSpaceAggr.get(wsetList.get(i), j));
                }
            }
            StateSpaceAggrWork = StateSpaceAggrSub;

            wsetMap = new int[nw];
            for (int i = 0; i < nw; i++) {
                wsetMap[i] = wsetList.get(i);
            }

            wset = new Matrix(1, nw);
            for (int i = 0; i < nw; i++) {
                wset.set(0, i, i);
            }
        } else {
            InputOutput.line_debug(options.verbose, "CTMC is irreducible, using full state space");
            probSysState = Ctmc_solve.ctmc_solve(InfGen);
        }

        if (probSysState.hasNaN() || probSysState.isEmpty()) {
            throw new RuntimeException("CTMC solver failed to compute steady-state probabilities for this cache model. " +
                    "This may indicate numerical instability or an invalid model configuration.");
        }

        for (int row = 0; row < probSysState.getNumRows(); row++) {
            for (int col = 0; col < probSysState.getNumCols(); col++) {
                if (probSysState.get(row, col) < 0) {
                    probSysState.set(row, col, 0);
                }
            }
        }

        double sum = probSysState.sumSubMatrix(0, probSysState.getNumRows(), 0, probSysState.getNumCols());
        if (sum > 0) {
            probSysState.divide(sum, probSysState, true);
        } else {
            throw new RuntimeException("CTMC solver computed zero total probability. This indicates an invalid model configuration.");
        }

        final int[] finalWsetMap = wsetMap;

        Matrix XN = new Matrix(1, K);
        XN.zero();
        Matrix UN = new Matrix(M, K);
        UN.zero();
        Matrix QN = new Matrix(M, K);
        QN.zero();
        Matrix RN = new Matrix(M, K);
        RN.zero();
        Matrix TN = new Matrix(M, K);
        TN.zero();
        Matrix CN = new Matrix(1, K);
        CN.zero();

        Matrix istSpaceShift = new Matrix(1, M);
        istSpaceShift.zero();

        for (int i = 0; i < M; i++) {
            if (i == 0) {
                istSpaceShift.set(0, i, 0);
            } else {
                double temp = istSpaceShift.get(0, i - 1) + sn.space.get(sn.stateful.get(i - 1)).getNumCols();
                istSpaceShift.set(0, i, temp);
            }
        }

        double refsf;
        for (int k = 0; k < K; k++) {
            refsf = sn.stationToStateful.get((int) sn.refstat.get(k));
            XN.set(0, k, refsf);
            double sumValue = 0.0;
            for (int i = 0; i < wset.getNumCols(); i++) {
                int index = (int) wset.get(i);
                int origIdx = (finalWsetMap != null) ? finalWsetMap[index] : index;
                sumValue += probSysState.get(index) * arvRates[origIdx][(int) refsf][k];
            }
            XN.set(0, k, sumValue);
        }

        // see _kb/06-solver-catalog.md for rationale
        boolean[] inDropRegion = new boolean[M];
        if (sn.nregions > 0 && sn.region != null && sn.regionrule != null) {
            for (int f = 0; f < sn.nregions; f++) {
                boolean isDropRegion = false;
                for (int r = 0; r < sn.regionrule.getNumCols(); r++) {
                    if (sn.regionrule.get(f, r) == DropStrategy.Drop.getID()) {
                        isDropRegion = true;
                        break;
                    }
                }
                if (!isDropRegion) {
                    continue;
                }
                Matrix regf = sn.region.get(f);
                for (int i = 0; i < M && i < regf.getNumRows(); i++) {
                    for (int c = 0; c < regf.getNumCols(); c++) {
                        if (regf.get(i, c) != -1) {
                            inDropRegion[i] = true;
                            break;
                        }
                    }
                }
            }
        }

        for (int i = 0; i < M; i++) {
            int isf = (int) sn.stationToStateful.get(i);
            int ind = (int) sn.stationToNode.get(i);
            for (int k = 0; k < K; k++) {
                double sumTN = 0.0;
                double sumQN = 0.0;
                for (int index = 0; index < wset.getNumCols(); index++) {
                    int wstIdx = (int) wset.get(index);
                    int origIdx = (finalWsetMap != null) ? finalWsetMap[wstIdx] : wstIdx;
                    double depRate = depRates[origIdx][isf][k];
                    double probState = probSysState.get(wstIdx);
                    sumTN += probState * depRate;
                    double ssaValue = StateSpaceAggrWork.get(wstIdx, i * K + k);
                    double product = probState * ssaValue;
                    sumQN += product;
                }
                TN.set(i, k, sumTN);
                QN.set(i, k, sumQN);
            }
            if (sn.nodetype.get(ind) != NodeType.Source) {
                // see _kb/06-solver-catalog.md for rationale
                boolean stationCapFinite = !Double.isInfinite(sn.cap.get(i)) && sn.cap.get(i) < Integer.MAX_VALUE;
                boolean[] canDropClass = new boolean[K];
                for (int r = 0; r < K; r++) {
                    boolean classCapFinite = !Double.isInfinite(sn.classcap.get(i, r)) && sn.classcap.get(i, r) < Integer.MAX_VALUE;
                    canDropClass[r] = Double.isInfinite(sn.njobs.get(r)) && (stationCapFinite || classCapFinite || inDropRegion[i]);
                }
                // see _kb/06-solver-catalog.md for rationale
                double[] arvAtStation = new double[K];
                for (int r = 0; r < K; r++) {
                    for (int idx = 0; idx < wset.length(); idx++) {
                        int wsetIdxSig = (int) wset.get(idx);
                        int origIdxSig = (finalWsetMap != null) ? finalWsetMap[wsetIdxSig] : wsetIdxSig;
                        arvAtStation[r] += probSysState.get(idx) * arvRates[origIdxSig][isf][r];
                    }
                }
                boolean[] signalLossy = jline.solvers.ctmc.handlers.CtmcSignalLossy.signalLossyClasses(sn, arvAtStation);
                for (int r = 0; r < K; r++) {
                    canDropClass[r] = canDropClass[r] || signalLossy[r];
                }
                SchedStrategy schedStrategy = schedid.get(sn.stations.get(i));
                if (schedStrategy == SchedStrategy.INF) {
                    int k = 0;
                    while (k < K) {
                        UN.set(i, k, QN.get(i, k));
                        k++;
                    }
                } else if (schedStrategy == SchedStrategy.PS || schedStrategy == SchedStrategy.DPS || schedStrategy == SchedStrategy.GPS) {
                    if (sn.lldscaling.isEmpty() && sn.cdscaling.isEmpty() && (sn.jdscaling == null || sn.jdscaling.isEmpty())) {
                        int k = 0;
                        while (k < K) {
                            if (!PHmap.get(sn.stations.get(i)).get(sn.jobclasses.get(k)).isEmpty()) {
                                MatrixCell value = PHmap.get(sn.stations.get(i)).get(sn.jobclasses.get(k));
                                double mean = Map_mean.map_mean(value) / S.get(i);
                                double UNarv_ik = 0.0;
                                int idx = 0;
                                while (idx < wset.length()) {
                                    int wsetIdx = (int) wset.get(idx);
                                    int origIdx = (finalWsetMap != null) ? finalWsetMap[wsetIdx] : wsetIdx;
                                    UNarv_ik += probSysState.get(idx) * arvRates[origIdx][isf][k];
                                    idx++;
                                }
                                UNarv_ik = UNarv_ik * mean;
                                double UNdep_ik = TN.get(i, k) * mean;
                                UN.set(i, k, canDropClass[k] ? UNdep_ik : Maths.max(UNarv_ik, UNdep_ik));
                            }
                            k++;
                        }
                    } else {
                        // see _kb/06-solver-catalog.md for rationale
                        ind = (int) sn.stationToNode.get(i);
                        double ceffPs = S.get(i);
                        if (sn.lldscaling != null && !sn.lldscaling.isEmpty() && i < sn.lldscaling.getNumRows()) {
                            for (int j = 0; j < sn.lldscaling.getNumCols(); j++) {
                                ceffPs = Math.max(ceffPs, sn.lldscaling.get(i, j));
                            }
                        }
                        int col = 0;
                        while (col < K) {
                            UN.set(i, col, 0);
                            col++;
                        }
                        int index = 0;
                        while (index < wset.getNumCols()) {
                            int st = (int) wset.get(index);
                            int StateSpaceColStart = (int) istSpaceShift.get(i);
                            int StateSpaceColEnd = (int) istSpaceShift.get(i) + sn.space.get(sn.stateful.get(i)).getNumCols();
                            jline.lang.state.State.StateMarginalStatistics toMarginalResult = ToMarginal.toMarginal(sn,
                                    ind,
                                    Matrix.extract(StateSpaceWork, st, st + 1, StateSpaceColStart, StateSpaceColEnd),
                                    null, null, null, null, null);
                            Matrix ni = toMarginalResult.ni;
                            Matrix nir = toMarginalResult.nir;
                            boolean checkNi = true;
                            double totJobsPs = 0.0;
                            int niIndex = 0;
                            while (niIndex < ni.length()) {
                                if (ni.get(niIndex) <= 0) {
                                    checkNi = false;
                                }
                                totJobsPs += ni.get(niIndex);
                                niIndex++;
                            }
                            if (checkNi) {
                                double lldNowPs = 1.0;
                                if (sn.lldscaling != null && !sn.lldscaling.isEmpty() && i < sn.lldscaling.getNumRows()) {
                                    int lldIdx = (int) Math.min(Math.max(totJobsPs, 1.0) - 1, sn.lldscaling.getNumCols() - 1);
                                    lldNowPs = sn.lldscaling.get(i, lldIdx);
                                }
                                int k = 0;
                                while (k < K) {
                                    double v = probSysState.get(st) * nir.get(k) * sn.schedparam.get(i, k);
                                    Matrix dividend = nir.mult(sn.schedparam.getRow(i).transpose());
                                    double addSum = 0.0;
                                    int divIndex = 0;
                                    while (divIndex < dividend.length()) {
                                        addSum += v / dividend.get(divIndex);
                                        divIndex++;
                                    }
                                    UN.set(i, k, addSum * lldNowPs / ceffPs + UN.get(i, k));
                                    k++;
                                }
                            }
                            index++;
                        }
                    }
                } else if (schedStrategy == SchedStrategy.PAS || schedStrategy == SchedStrategy.OI) {
                    // see _kb/06-solver-catalog.md for rationale
                    int col = 0;
                    while (col < K) {
                        UN.set(i, col, 0);
                        col++;
                    }
                    int index = 0;
                    while (index < wset.length()) {
                        int st = (int) wset.get(index);
                        int StateSpaceColStart = (int) istSpaceShift.get(i);
                        int StateSpaceColEnd = (int) istSpaceShift.get(i) + sn.space.get(sn.stateful.get(i)).getNumCols();
                        jline.lang.state.State.StateMarginalStatistics toMarginalResult = ToMarginal.toMarginal(sn,
                                ind,
                                Matrix.extract(StateSpaceWork, st, st + 1, StateSpaceColStart, StateSpaceColEnd),
                                null, null, null, null, null);
                        Matrix sir = toMarginalResult.sir;
                        int k = 0;
                        while (k < K) {
                            double v = UN.get(i, k) + probSysState.get(st) * sir.get(k) / S.get(i);
                            UN.set(i, k, v);
                            k++;
                        }
                        index++;
                    }
                } else {
                    if ((sn.lldscaling == null || sn.lldscaling.isEmpty()) && (sn.cdscaling == null || sn.cdscaling.isEmpty()) && (sn.jdscaling == null || sn.jdscaling.isEmpty())) {
                        int k = 0;
                        while (k < K) {
                            if (!PHmap.get(sn.stations.get(i)).get(sn.jobclasses.get(k)).isEmpty()) {
                                double UNarv_ik = 0.0;
                                MatrixCell value = PHmap.get(sn.stations.get(i)).get(sn.jobclasses.get(k));
                                double mean = Map_mean.map_mean(value);
                                int idx = 0;
                                while (idx < wset.length()) {
                                    int wsetIdx = (int) wset.get(idx);
                                    int origIdx = (finalWsetMap != null) ? finalWsetMap[wsetIdx] : wsetIdx;
                                    UNarv_ik += probSysState.get(idx) * arvRates[origIdx][isf][k];
                                    idx++;
                                }
                                UNarv_ik = UNarv_ik * mean / S.get(i);
                                double UNdep_ik = TN.get(i, k) * Map_mean.map_mean(value) / S.get(i);
                                UN.set(i, k, canDropClass[k] ? UNdep_ik : Maths.max(UNarv_ik, UNdep_ik));
                            }
                            k++;
                        }
                    } else {
                        ind = (int) sn.stationToNode.get(i);
                        int col = 0;
                        while (col < K) {
                            UN.set(i, col, 0);
                            col++;
                        }
                        int index = 0;
                        while (index < wset.length()) {
                            int st = (int) wset.get(index);
                            int StateSpaceColStart = (int) istSpaceShift.get(i);
                            int StateSpaceColEnd = (int) istSpaceShift.get(i) + sn.space.get(sn.stateful.get(i)).getNumCols();
                            jline.lang.state.State.StateMarginalStatistics toMarginalResult = ToMarginal.toMarginal(sn,
                                    ind,
                                    Matrix.extract(StateSpaceWork, st, st + 1, StateSpaceColStart, StateSpaceColEnd),
                                    null, null, null, null, null);
                            Matrix ni = toMarginalResult.ni;
                            Matrix sir = toMarginalResult.sir;
                            boolean checkNir = true;
                            int niIndex = 0;
                            while (niIndex < ni.length()) {
                                if (ni.get(niIndex) <= 0) {
                                    checkNir = false;
                                }
                                niIndex++;
                            }
                            if (checkNir) {
                                // see _kb/06-solver-catalog.md for rationale
                                double totJobs = 0.0;
                                for (int nidx = 0; nidx < ni.length(); nidx++) totJobs += ni.get(nidx);
                                double lldNow = 1.0;
                                double ceff = S.get(i);
                                if (sn.lldscaling != null && !sn.lldscaling.isEmpty() && i < sn.lldscaling.getNumRows()) {
                                    int colIdx = (int) Maths.min(Maths.max(totJobs, 1.0) - 1, sn.lldscaling.getNumCols() - 1);
                                    lldNow = sn.lldscaling.get(i, colIdx);
                                    for (int cidx = 0; cidx < sn.lldscaling.getNumCols(); cidx++) {
                                        ceff = Maths.max(ceff, sn.lldscaling.get(i, cidx));
                                    }
                                }
                                double sirTot = 0.0;
                                for (int k2 = 0; k2 < K; k2++) sirTot += sir.get(k2);
                                int k = 0;
                                while (k < K) {
                                    double share = sirTot > 0 ? sir.get(k) / sirTot : 0.0;
                                    double v = UN.get(i, k) + probSysState.get(st) * share * lldNow / ceff;
                                    UN.set(i, k, v);
                                    k++;
                                }
                            }
                            index++;
                        }
                    }
                }
                // see _kb/06-solver-catalog.md for rationale
                boolean anySignalLossy = false;
                for (int r = 0; r < K; r++) {
                    anySignalLossy = anySignalLossy || signalLossy[r];
                }
                if (anySignalLossy && schedStrategy != SchedStrategy.INF
                        && sn.lldscaling.isEmpty() && sn.cdscaling.isEmpty() && (sn.jdscaling == null || sn.jdscaling.isEmpty())) {
                    double[] UNb = jline.solvers.ctmc.handlers.CtmcSignalBusy.busyFraction(sn,
                            (int) sn.stationToNode.get(i), i, schedStrategy, S.get(i),
                            StateSpaceWork, istSpaceShift, wset, probSysState, K);
                    for (int k = 0; k < K; k++) {
                        if (signalLossy[k]) {
                            UN.set(i, k, UNb[k]);
                        }
                    }
                }
            }
        }

        // see _kb/06-solver-catalog.md for rationale
        if (sn.cdscaling != null && !sn.cdscaling.isEmpty()) {
            boolean allFinite = true;
            for (int k = 0; k < K; k++) {
                if (Double.isInfinite(sn.njobs.get(k))) { allFinite = false; break; }
            }
            if (allFinite) {
                for (int ist = 0; ist < M; ist++) {
                    jline.lang.nodes.Station stat = sn.stations.get(ist);
                    jline.util.SerializableFunction<Matrix, Matrix> beta = sn.cdscaling.get(stat);
                    if (beta == null) continue;
                    Matrix peakVec = sn.cdscalingpeak != null ? sn.cdscalingpeak.get(stat) : null;
                    for (int k = 0; k < K; k++) {
                        double rate = sn.rates.get(ist, k);
                        double bmax = (peakVec != null) ? peakVec.get(0, k) : 1.0;
                        if (Double.isFinite(rate) && rate > 0 && bmax > 0) {
                            UN.set(ist, k, TN.get(ist, k) / rate / bmax);
                        } else {
                            UN.set(ist, k, 0.0);
                        }
                    }
                }
            }
        }

        // joint-dependence utilization normalization (U = T*S/peak using the
        // declared sn.jdscalingpeak), mirroring the class-dependence block.
        if (sn.jdscaling != null && !sn.jdscaling.isEmpty()) {
            boolean allFinite = true;
            for (int k = 0; k < K; k++) {
                if (Double.isInfinite(sn.njobs.get(k))) { allFinite = false; break; }
            }
            if (allFinite) {
                for (int ist = 0; ist < M; ist++) {
                    jline.lang.nodes.Station stat = sn.stations.get(ist);
                    jline.util.SerializableFunction<Matrix, Matrix> eta = sn.jdscaling.get(stat);
                    if (eta == null) continue;
                    Matrix peakVec = sn.jdscalingpeak != null ? sn.jdscalingpeak.get(stat) : null;
                    for (int k = 0; k < K; k++) {
                        double rate = sn.rates.get(ist, k);
                        double bmax = (peakVec != null) ? peakVec.get(0, k) : 1.0;
                        if (Double.isFinite(rate) && rate > 0 && bmax > 0) {
                            UN.set(ist, k, TN.get(ist, k) / rate / bmax);
                        } else {
                            UN.set(ist, k, 0.0);
                        }
                    }
                }
            }
        }

        // see _kb/06-solver-catalog.md for rationale
        if (sn.isbasblocking != null && sn.connmatrix != null && !sn.connmatrix.isEmpty()) {
            for (int i = 0; i < M; i++) {
                int indB = (int) sn.stationToNode.get(i);
                // see _kb/06-solver-catalog.md for rationale
                if (indB >= sn.isbasblocking.length() || sn.isbasblocking.get(indB) != 1) {
                    continue; // no blocked marker at this station
                }
                int dest = -1;
                boolean ambiguous = false;
                for (int j = 0; j < sn.connmatrix.getNumCols(); j++) {
                    if (sn.connmatrix.get(indB, j) != 1 || sn.isstation.get(j) != 1) {
                        continue;
                    }
                    if (dest >= 0) {
                        ambiguous = true;
                        break;
                    }
                    dest = j;
                }
                if (ambiguous || dest < 0) {
                    continue; // ambiguous destination: leave the job where it sits
                }
                int jst = (int) sn.nodeToStation.get(dest);
                if (jst < 0) {
                    continue;
                }
                int bcol = (int) istSpaceShift.get(0, i)
                        + sn.space.get(sn.stateful.get(i)).getNumCols() - 1;
                for (int k = 0; k < K; k++) {
                    double shift = 0.0;
                    for (int index = 0; index < wset.getNumCols(); index++) {
                        int wstIdx = (int) wset.get(index);
                        // see _kb/06-solver-catalog.md for rationale
                        if (StateSpaceWork.get(wstIdx, bcol) != 1.0) {
                            continue;
                        }
                        // Only the held job itself moves, not the whole queue at i: a
                        // blocked state holds exactly one completed job, so cap at 1.
                        double nk = StateSpaceAggrWork.get(wstIdx, i * K + k);
                        shift += probSysState.get(wstIdx) * Math.min(nk, 1.0);
                    }
                    if (shift > 0) {
                        QN.set(i, k, QN.get(i, k) - shift);
                        QN.set(jst, k, QN.get(jst, k) + shift);
                    }
                }
            }
        }

        for (int k = 0; k < K; k++) {
            for (int i = 0; i < M; i++) {
                if (TN.get(i, k) > 0) {
                    RN.set(i, k, QN.get(i, k) / TN.get(i, k));
                } else {
                    RN.set(i, k, 0);
                }
            }
            CN.set(k, NK.get(k) / XN.get(k));
        }
        QN.setNaNToZero();
        CN.setNaNToZero();
        RN.setNaNToZero();
        UN.setNaNToZero();
        XN.setNaNToZero();
        TN.setNaNToZero();

        long Tstop = System.nanoTime();
        double runtime = ((double) (Tstop - Tstart)) / 1000000000.0;

        Matrix TNcache = new Matrix(sn.nstateful, K);
        Matrix XNcache = new Matrix(sn.nstateful, K);
        for (int k = 0; k < K; k++) {
            for (int isf = 0; isf < sn.nstateful; isf++) {
                int ind = (int) sncopy.statefulToNode.get(isf);
                if (sncopy.nodetype.get(ind) == NodeType.Cache) {
                    double TNcacheSum = 0.0;
                    double XNcacheSum = 0.0;
                    for (int i = 0; i < wset.getNumCols(); i++) {
                        int index = (int) wset.get(i);
                        int origIdx = (finalWsetMap != null) ? finalWsetMap[index] : index;
                        TNcacheSum += probSysState.get(i) * depRates[origIdx][isf][k];
                        XNcacheSum += probSysState.get(i) * arvRates[origIdx][isf][k];
                    }
                    TNcache.set(isf, k, TNcacheSum);
                    XNcache.set(isf, k, XNcacheSum);
                }
            }
        }

        boolean retrievalLatencyWarned = false;
        for (int k = 0; k < K; k++) {
            for (int isf = 0; isf < sncopy.nstateful; isf++) {
                int ind = (int) sncopy.statefulToNode.get(isf);
                if (sncopy.nodetype.get(ind) == NodeType.Cache) {
                    CacheNodeParam cacheParam = (CacheNodeParam) sncopy.nodeparam.get(sn.nodes.get(ind));
                    if (cacheParam.hitclass.length() > k) {
                        int h = (int) cacheParam.hitclass.get(k);
                        int m = (int) cacheParam.missclass.get(k);

                        if (cacheParam.actualhitprob == null) {
                            cacheParam.actualhitprob = new Matrix(1, K);
                            cacheParam.actualhitprob.fill(Double.NaN);
                        }
                        if (cacheParam.actualmissprob == null) {
                            cacheParam.actualmissprob = new Matrix(1, K);
                            cacheParam.actualmissprob.fill(Double.NaN);
                        }

                        double actualmissprobValue = Double.NaN;
                        double actualhitprobValue = Double.NaN;
                        if (h != -1 && m != -1) {
                            actualhitprobValue = TNcache.get(isf, h) / (TNcache.get(isf, h) + TNcache.get(isf, m));
                            actualmissprobValue = TNcache.get(isf, m) / (TNcache.get(isf, h) + TNcache.get(isf, m));
                        }
                        cacheParam.actualhitprob.set(k, actualhitprobValue);
                        cacheParam.actualmissprob.set(k, actualmissprobValue);

                        // see _kb/06-solver-catalog.md for rationale
                        double actualresidt = Double.NaN;
                        List<Integer> retrievalSystemQueueIndices =
                                (cacheParam.retrievalSystemQueueIndices != null)
                                        ? cacheParam.retrievalSystemQueueIndices.get(k) : null;
                        if (retrievalSystemQueueIndices != null && !retrievalSystemQueueIndices.isEmpty()) {
                            if (!retrievalLatencyWarned) {
                                InputOutput.line_warning("solver_ctmc_analyzer",
                                        "Retrieval-system expected latency is not currently implemented; "
                                        + "reporting NaN.");
                                retrievalLatencyWarned = true;
                            }
                        }

                        if (cacheParam.actualresidt == null) {
                            cacheParam.actualresidt = new Matrix(1, K);
                            cacheParam.actualresidt.fill(Double.NaN);
                        }
                        cacheParam.actualresidt.set(k, actualresidt);
                    }
                }
            }
        }
        return new SolverCTMC.AnalyzerResult(QN, UN, RN, TN, CN, XN,
                InfGen, StateSpace, StateSpaceAggr, EventFiltration,
                runtime, fname, sncopy,
                probSysState, StateSpaceWork, StateSpaceAggrWork, InfGenWork);
    }

    /**
     * Returns the set of state indices forward-reachable from {@code start} in the
     * directed graph of positive off-diagonal entries of the generator {@code Q}.
     * For a closed pass-and-swap network whose initial placement lies in a
     * recurrent class, this is exactly that recurrent class.
     */
    private static List<Integer> forwardReachable(Matrix Q, int start) {
        int n = Q.length();
        boolean[] visited = new boolean[n];
        java.util.ArrayDeque<Integer> queue = new java.util.ArrayDeque<Integer>();
        visited[start] = true;
        queue.add(start);
        while (!queue.isEmpty()) {
            int u = queue.poll();
            for (int v = 0; v < n; v++) {
                if (v != u && !visited[v] && Math.abs(Q.get(u, v)) > 1e-12) {
                    visited[v] = true;
                    queue.add(v);
                }
            }
        }
        List<Integer> result = new ArrayList<Integer>();
        for (int i = 0; i < n; i++) {
            if (visited[i]) result.add(i);
        }
        return result;
    }
}
