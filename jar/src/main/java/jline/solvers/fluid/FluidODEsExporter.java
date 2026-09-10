/*
 * Copyright (c) 2012-2026, QORE Lab, Imperial College London
 * All rights reserved.
 */
package jline.solvers.fluid;

import java.util.ArrayList;
import java.util.List;
import java.util.Locale;

import jline.GlobalConstants;
import jline.lang.JobClass;
import jline.lang.NetworkStruct;
import jline.lang.constant.SchedStrategy;
import jline.lang.nodes.Station;
import jline.solvers.SolverOptions;
import jline.api.mc.Dtmc_stochcomp;
import jline.util.matrix.Matrix;
import jline.util.matrix.MatrixCell;
import jline.util.matrix.MatrixEquation;

import static jline.api.mam.Map_pie.map_pie;

/**
 * Symbolic export of the mean-field ODE system integrated by SolverFluid.
 * Mirrors the MATLAB implementation in solver_fluid_symodes.m and
 * SolverFLD/exportODEs.m so that all codebases emit the same LaTeX document
 * for a given model. Two representations are produced:
 *
 * form = "W": dx/dt = W'*theta(x) + lambda (methods: default, matrix, pnorm)
 * form = "J": dx/dt = J*r(x) (methods: closing, statedep, softmin)
 */
public class FluidODEsExporter {

    /**
     * Structural description of the exported ODE system. State variable,
     * station and class indices are 0-based; they are rendered 1-based.
     */
    public static class SymODEs {
        public String form;
        public String method;
        public int nstates;
        public int[] stateStation;
        public int[] stateClass;
        public int[] statePhase; // 1-based phase, 0 marks removed placeholders
        public List<String> stationNames;
        public List<String> classNames;
        public List<String> schedNames;
        public SchedStrategy[] sched;
        public double[] S;
        /** S holds the population at an INF station, so isInfinite(S) cannot tell one apart. */
        public boolean[] isInfStation;
        // form W
        public Matrix W;
        public double[] Alambda;
        public boolean[] isSource;
        public String smoothing;
        public double[] pstar;
        public int[] keep;
        // form J
        public int nevents;
        public int[] eventFrom;
        public int[] eventTo;
        public double[] coeff;
        public int[] eventVar;
        public String[] factorType;
        public int[] factorStation;
        public int[] factorClass;
        public int[][] factorOthers;
        public double[][] dpsw;
        public double[] fcfsPhaseW;
        public double alpha;
        // initial condition
        public double[] x0;
    }

    /**
     * Build the structural description of the ODE system for the method set
     * in the solver options.
     */
    public static SymODEs build(NetworkStruct sn, SolverOptions options) {
        int M = sn.nstations;
        int K = sn.nclasses;

        String method = options.method.replace("fluid.", "");
        String form;
        if (method.equals("default") || method.equals("matrix") || method.equals("pnorm")) {
            form = "W";
            if (method.equals("default")) {
                method = "matrix";
            }
        } else if (method.equals("closing") || method.equals("statedep") || method.equals("softmin")) {
            form = "J";
        } else {
            throw new RuntimeException("Symbolic ODE export is unsupported for method '" + method
                    + "'. Supported methods: default, matrix, pnorm, closing, statedep, softmin.");
        }

        SymODEs sys = new SymODEs();
        sys.form = form;
        sys.method = method;
        sys.stationNames = new ArrayList<String>();
        sys.classNames = new ArrayList<String>();
        sys.schedNames = new ArrayList<String>();
        sys.sched = new SchedStrategy[M];
        for (int i = 0; i < M; i++) {
            sys.stationNames.add(sn.nodenames.get((int) sn.stationToNode.get(i)));
            sys.sched[i] = sn.sched.get(sn.stations.get(i));
            sys.schedNames.add(SchedStrategy.toText(sys.sched[i]));
        }
        for (int r = 0; r < K; r++) {
            sys.classNames.add(sn.classnames.get(r));
        }

        if (form.equals("W")) {
            buildWForm(sys, sn, options, M, K);
        } else {
            buildJForm(sys, sn, M, K, method);
        }
        buildX0(sys, sn, options, M, K);
        return sys;
    }

    private static boolean isDisabled(Matrix mu) {
        if (mu == null || mu.isEmpty()) {
            return true;
        }
        for (int k = 0; k < mu.length(); k++) {
            if (!Double.isNaN(mu.get(k))) {
                return false;
            }
        }
        return true;
    }

    private static void buildWForm(SymODEs sys, NetworkStruct sn, SolverOptions options, int M, int K) {
        double[] S = new double[M];
        double njobsSum = sn.njobs.elementSum();
        for (int i = 0; i < M; i++) {
            double si = sn.nservers.get(i, 0);
            S[i] = Double.isInfinite(si) ? njobsSum : si;
        }
        sys.isInfStation = new boolean[M];
        for (int i = 0; i < M; i++) {
            sys.isInfStation[i] = Double.isInfinite(sn.nservers.get(i, 0));
        }

        // station-to-station routing matrix via stochastic complementation
        List<Integer> stationIndices = new ArrayList<Integer>();
        for (int ist = 0; ist < M; ist++) {
            int isf = (int) sn.stationToStateful.get(ist);
            for (int r = 0; r < K; r++) {
                stationIndices.add(isf * K + r);
            }
        }
        Matrix P = Dtmc_stochcomp.dtmc_stochcomp(sn.rt, stationIndices);

        // remove Sink->Source feedback routing for open classes
        for (int srcIst = 0; srcIst < M; srcIst++) {
            if (sys.sched[srcIst] == SchedStrategy.EXT) {
                for (int r = 0; r < K; r++) {
                    if (!Double.isNaN(sn.rates.get(srcIst, r)) && sn.rates.get(srcIst, r) > 0) {
                        int srcCol = srcIst * K + r;
                        for (int fromIst = 0; fromIst < M; fromIst++) {
                            if (fromIst != srcIst) {
                                for (int fromR = 0; fromR < K; fromR++) {
                                    P.set(fromIst * K + fromR, srcCol, 0.0);
                                }
                            }
                        }
                    }
                }
            }
        }

        // W = Psi + B*P*A'
        Matrix psi = new Matrix(0, 0);
        Matrix A = new Matrix(0, 0);
        Matrix B = new Matrix(0, 0);
        for (int ist = 0; ist < M; ist++) {
            Station station = sn.stations.get(ist);
            for (int r = 0; r < K; r++) {
                JobClass jobClass = sn.jobclasses.get(r);
                if (sn.phases.get(ist, r) == 0.0) {
                    Matrix zeroMatrix = new Matrix(1, 1, 1);
                    Matrix nanMatrix = new Matrix(1, 1, 1);
                    nanMatrix.set(0, 0, Double.NaN);
                    psi = psi.createBlockDiagonal(zeroMatrix);
                    A = A.createBlockDiagonal(nanMatrix);
                    B = B.createBlockDiagonal(zeroMatrix);
                } else {
                    psi = psi.createBlockDiagonal(sn.proc.get(station).get(jobClass).get(0));
                    A = A.createBlockDiagonal(sn.pie.get(station).get(jobClass).transpose());
                    B = B.createBlockDiagonal(sn.proc.get(station).get(jobClass).get(1).sumRows());
                }
            }
        }
        MatrixEquation calc = new MatrixEquation();
        calc.alias(psi, "psi", A, "A", P, "P", B, "B");
        calc.process("W = psi + B*P*A'");
        Matrix W = new Matrix(calc.lookupSimple("W"));

        // exogenous arrival rates into queue phases
        Matrix sourceArrivals = new Matrix(M, K);
        for (int srcIst = 0; srcIst < M; srcIst++) {
            if (sys.sched[srcIst] == SchedStrategy.EXT) {
                for (int r = 0; r < K; r++) {
                    if (!Double.isNaN(sn.rates.get(srcIst, r)) && sn.rates.get(srcIst, r) > 0) {
                        sourceArrivals.set(srcIst, r, sn.rates.get(srcIst, r));
                    }
                }
            }
        }
        int totalStates = W.getNumRows();
        double[] AlambdaFull = new double[totalStates];
        int state = 0;
        for (int ist = 0; ist < M; ist++) {
            Station station = sn.stations.get(ist);
            for (int r = 0; r < K; r++) {
                JobClass jobClass = sn.jobclasses.get(r);
                int nPhases = (int) sn.phases.get(ist, r);
                if (nPhases > 0) {
                    if (sys.sched[ist] == SchedStrategy.EXT) {
                        state += nPhases;
                    } else {
                        double arrivalRateToQueue = 0.0;
                        for (int srcIst = 0; srcIst < M; srcIst++) {
                            if (sourceArrivals.get(srcIst, r) > 0) {
                                arrivalRateToQueue += sourceArrivals.get(srcIst, r) * P.get(srcIst * K + r, ist * K + r);
                            }
                        }
                        if (arrivalRateToQueue > 0) {
                            Matrix pieVec = sn.pie.get(station).get(jobClass);
                            for (int k = 0; k < nPhases; k++) {
                                AlambdaFull[state] = pieVec.get(k) * arrivalRateToQueue;
                                state++;
                            }
                        } else {
                            state += nPhases;
                        }
                    }
                } else {
                    state++;
                }
            }
        }

        // state metadata over the same enumeration, then keep-filter
        int[] stStation = new int[totalStates];
        int[] stClass = new int[totalStates];
        int[] stPhase = new int[totalStates];
        state = 0;
        for (int ist = 0; ist < M; ist++) {
            for (int r = 0; r < K; r++) {
                int nPhases = (int) sn.phases.get(ist, r);
                if (nPhases == 0) {
                    stStation[state] = ist;
                    stClass[state] = r;
                    stPhase[state] = 0; // placeholder, removed by keep filter
                    state++;
                } else {
                    for (int k = 0; k < nPhases; k++) {
                        stStation[state] = ist;
                        stClass[state] = r;
                        stPhase[state] = k + 1;
                        state++;
                    }
                }
            }
        }

        List<Integer> keepList = new ArrayList<Integer>();
        Matrix sumWCols = W.sumCols();
        for (int col = 0; col < W.getNumCols(); col++) {
            if (!Double.isNaN(sumWCols.get(0, col))) {
                keepList.add(col);
            }
        }
        int n = keepList.size();
        Matrix Wk = new Matrix(n, n);
        sys.keep = new int[n];
        sys.Alambda = new double[n];
        sys.stateStation = new int[n];
        sys.stateClass = new int[n];
        sys.statePhase = new int[n];
        sys.isSource = new boolean[n];
        for (int i = 0; i < n; i++) {
            int ki = keepList.get(i);
            sys.keep[i] = ki;
            for (int j = 0; j < n; j++) {
                Wk.set(i, j, W.get(ki, keepList.get(j)));
            }
            sys.Alambda[i] = AlambdaFull[ki];
            sys.stateStation[i] = stStation[ki];
            sys.stateClass[i] = stClass[ki];
            sys.statePhase[i] = stPhase[ki];
            sys.isSource[i] = (sys.sched[stStation[ki]] == SchedStrategy.EXT);
        }
        sys.W = Wk;
        sys.nstates = n;
        sys.S = S;

        // smoothing selection mirrors the pstar gate in the matrix analyzer
        if (options.config.pstar != null && options.config.pstar.size() > 0) {
            sys.smoothing = "pnorm";
            sys.pstar = new double[M];
            for (int i = 0; i < M; i++) {
                if (options.config.pstar.size() == 1) {
                    sys.pstar[i] = options.config.pstar.get(0);
                } else {
                    sys.pstar[i] = options.config.pstar.get(i);
                }
            }
        } else {
            sys.smoothing = "min";
            sys.pstar = null;
        }
    }

    private static void buildJForm(SymODEs sys, NetworkStruct sn, int M, int K, String method) {
        double[] S = new double[M];
        for (int i = 0; i < M; i++) {
            double si = sn.nservers.get(i, 0);
            S[i] = Double.isInfinite(si) ? sn.nclosedjobs : si;
        }
        sys.isInfStation = new boolean[M];
        for (int i = 0; i < M; i++) {
            sys.isInfStation[i] = Double.isInfinite(sn.nservers.get(i, 0));
        }

        boolean weighted = method.equals("statedep") || method.equals("softmin");
        if (weighted) {
            for (int i = 0; i < M; i++) {
                if (sys.sched[i] == SchedStrategy.EXT) {
                    throw new RuntimeException("The '" + method + "' ODE method does not support open models, "
                            + "so their ODE system cannot be exported. Use the 'matrix' or 'closing' method instead.");
                }
            }
        }

        // state indexing as in the closing/statedep/softmin ODE handlers
        Matrix[][] Mu = new Matrix[M][K];
        Matrix[][] Phi = new Matrix[M][K];
        int[][] qIndices = new int[M][K];
        int[][] Kic = new int[M][K];
        boolean[][] enabled = new boolean[M][K];
        int cs = 0;
        for (int i = 0; i < M; i++) {
            Station station = sn.stations.get(i);
            for (int c = 0; c < K; c++) {
                JobClass jobClass = sn.jobclasses.get(c);
                Matrix mu = sn.mu.get(station).get(jobClass);
                Matrix phi = sn.phi.get(station).get(jobClass);
                if (isDisabled(mu)) {
                    Mu[i][c] = null;
                    Phi[i][c] = null;
                    Kic[i][c] = 0;
                    enabled[i][c] = false;
                } else {
                    Mu[i][c] = mu;
                    Phi[i][c] = phi;
                    Kic[i][c] = mu.length();
                    enabled[i][c] = true;
                }
                qIndices[i][c] = cs;
                cs += Kic[i][c];
            }
        }
        int nstates = cs;

        sys.stateStation = new int[nstates];
        sys.stateClass = new int[nstates];
        sys.statePhase = new int[nstates];
        for (int i = 0; i < M; i++) {
            for (int c = 0; c < K; c++) {
                for (int k = 0; k < Kic[i][c]; k++) {
                    sys.stateStation[qIndices[i][c] + k] = i;
                    sys.stateClass[qIndices[i][c] + k] = c;
                    sys.statePhase[qIndices[i][c] + k] = k + 1;
                }
            }
        }

        // normalized DPS weights
        double[][] dpsw = new double[M][K];
        for (int i = 0; i < M; i++) {
            for (int c = 0; c < K; c++) {
                dpsw[i][c] = 1.0;
            }
            if (sys.sched[i] == SchedStrategy.DPS) {
                double sum = 0.0;
                for (int c = 0; c < K; c++) {
                    dpsw[i][c] = sn.schedparam.get(i, c);
                    sum += dpsw[i][c];
                }
                for (int c = 0; c < K; c++) {
                    dpsw[i][c] /= sum;
                }
            }
        }

        // FCFS phase weights used by statedep/softmin: w = -1/D0(k,k)
        double[] fcfsPhaseW = new double[nstates];
        if (weighted) {
            for (int i = 0; i < M; i++) {
                if (sys.sched[i] == SchedStrategy.FCFS) {
                    Station station = sn.stations.get(i);
                    for (int c = 0; c < K; c++) {
                        if (enabled[i][c]) {
                            MatrixCell ph = sn.proc.get(station).get(sn.jobclasses.get(c));
                            for (int k = 0; k < Kic[i][c]; k++) {
                                fcfsPhaseW[qIndices[i][c] + k] = -1.0 / ph.get(0).get(k, k);
                            }
                        }
                    }
                }
            }
        }

        // event enumeration (departures first, then phase changes); zero-rate events dropped
        List<Integer> evFrom = new ArrayList<Integer>();
        List<Integer> evTo = new ArrayList<Integer>();
        List<Double> evCoeff = new ArrayList<Double>();
        List<Integer> evVar = new ArrayList<Integer>();
        List<String> evType = new ArrayList<String>();
        List<Integer> evStation = new ArrayList<Integer>();
        List<Integer> evClass = new ArrayList<Integer>();
        List<int[]> evOthers = new ArrayList<int[]>();

        for (int i = 0; i < M; i++) {
            for (int c = 0; c < K; c++) {
                if (enabled[i][c]) {
                    int xic = qIndices[i][c];
                    for (int j = 0; j < M; j++) {
                        for (int l = 0; l < K; l++) {
                            if (sn.rt.get(i * K + c, j * K + l) > 0) {
                                MatrixCell phTarget = sn.proc.get(sn.stations.get(j)).get(sn.jobclasses.get(l));
                                Matrix pieVec;
                                if (phTarget == null || phTarget.isEmpty()) {
                                    pieVec = new Matrix(1, 1, 1);
                                    pieVec.set(0, 0, 1.0);
                                } else {
                                    pieVec = map_pie(phTarget.get(0), phTarget.get(1));
                                }
                                int xjl = qIndices[j][l];
                                for (int ki = 0; ki < Kic[i][c]; ki++) {
                                    for (int kj = 0; kj < Kic[j][l]; kj++) {
                                        if (weighted) {
                                            if (sys.sched[i] == SchedStrategy.INF && j == i) {
                                                continue; // self-loop departures skipped at INF stations
                                            }
                                            if (!isHandledSched(sys.sched[i])) {
                                                continue;
                                            }
                                        }
                                        double base = Phi[i][c].get(ki) * Mu[i][c].get(ki)
                                                * sn.rt.get(i * K + c, j * K + l) * pieVec.get(kj);
                                        if (base > 0) {
                                            addEvent(evFrom, evTo, evCoeff, evVar, evType, evStation, evClass, evOthers,
                                                    xic + ki, xjl + kj, base, method, sys.sched[i], i, c, ki,
                                                    qIndices, Kic, S, K, dpsw, fcfsPhaseW);
                                        }
                                    }
                                }
                            }
                        }
                    }
                }
            }
        }
        for (int i = 0; i < M; i++) {
            for (int c = 0; c < K; c++) {
                if (enabled[i][c]) {
                    if (weighted && !isHandledSched(sys.sched[i])) {
                        continue;
                    }
                    int xic = qIndices[i][c];
                    MatrixCell ph = sn.proc.get(sn.stations.get(i)).get(sn.jobclasses.get(c));
                    for (int ki = 0; ki < Kic[i][c] - 1; ki++) {
                        for (int kip = 0; kip < Kic[i][c]; kip++) {
                            if (ki != kip) {
                                double base = ph.get(0).get(ki, kip);
                                if (base > 0) {
                                    addEvent(evFrom, evTo, evCoeff, evVar, evType, evStation, evClass, evOthers,
                                            xic + ki, xic + kip, base, method, sys.sched[i], i, c, ki,
                                            qIndices, Kic, S, K, dpsw, fcfsPhaseW);
                                }
                            }
                        }
                    }
                }
            }
        }

        int ne = evCoeff.size();
        sys.nevents = ne;
        sys.eventFrom = new int[ne];
        sys.eventTo = new int[ne];
        sys.coeff = new double[ne];
        sys.eventVar = new int[ne];
        sys.factorType = new String[ne];
        sys.factorStation = new int[ne];
        sys.factorClass = new int[ne];
        sys.factorOthers = new int[ne][];
        for (int e = 0; e < ne; e++) {
            sys.eventFrom[e] = evFrom.get(e);
            sys.eventTo[e] = evTo.get(e);
            sys.coeff[e] = evCoeff.get(e);
            sys.eventVar[e] = evVar.get(e);
            sys.factorType[e] = evType.get(e);
            sys.factorStation[e] = evStation.get(e);
            sys.factorClass[e] = evClass.get(e);
            sys.factorOthers[e] = evOthers.get(e);
        }
        sys.nstates = nstates;
        sys.S = S;
        sys.dpsw = dpsw;
        sys.fcfsPhaseW = fcfsPhaseW;
        if (method.equals("softmin")) {
            sys.alpha = 20.0; // softmin parameter, as in the softmin ODE handler
        }
    }

    private static boolean isHandledSched(SchedStrategy s) {
        return s == SchedStrategy.INF || s == SchedStrategy.EXT || s == SchedStrategy.PS
                || s == SchedStrategy.FCFS || s == SchedStrategy.DPS;
    }

    private static void addEvent(List<Integer> evFrom, List<Integer> evTo, List<Double> evCoeff,
                                 List<Integer> evVar, List<String> evType, List<Integer> evStation,
                                 List<Integer> evClass, List<int[]> evOthers,
                                 int from, int to, double base, String method, SchedStrategy schedi,
                                 int i, int c, int ki, int[][] qIndices, int[][] Kic, double[] S,
                                 int K, double[][] dpsw, double[] fcfsPhaseW) {
        String ftype;
        int[] others = null;
        double coeff = base;
        if (method.equals("closing")) {
            if (schedi == SchedStrategy.INF) {
                ftype = "lin";
            } else if (schedi == SchedStrategy.EXT) {
                if (ki == 0) {
                    ftype = "ext1";
                    others = new int[Kic[i][c] - 1];
                    for (int u = 1; u < Kic[i][c]; u++) {
                        others[u - 1] = qIndices[i][c] + u;
                    }
                } else {
                    ftype = "lin";
                }
            } else if (schedi == SchedStrategy.PS || schedi == SchedStrategy.FCFS) {
                ftype = "min";
            } else if (schedi == SchedStrategy.DPS) {
                // the share w_ir*x/ntilde_i of the capacity min(n_i,S_i), as in
                // the closing rate factors: no additive seed, not the full S_i
                ftype = "dpsmin";
                coeff = coeff * dpsw[i][c];
            } else {
                // strategies without a case in the closing rates keep rates = x
                ftype = "lin";
            }
        } else { // statedep, softmin
            if (schedi == SchedStrategy.INF) {
                ftype = "lin";
            } else if (schedi == SchedStrategy.PS) {
                ftype = "min";
            } else if (schedi == SchedStrategy.FCFS) {
                ftype = method.equals("softmin") ? "fcfsws" : "fcfsw";
                coeff = coeff * fcfsPhaseW[qIndices[i][c] + ki];
            } else { // DPS
                ftype = "dpspw";
            }
        }
        evFrom.add(from);
        evTo.add(to);
        evCoeff.add(coeff);
        evVar.add(qIndices[i][c] + ki);
        evType.add(ftype);
        evStation.add(i);
        evClass.add(c);
        evOthers.add(others);
    }

    private static void buildX0(SymODEs sys, NetworkStruct sn, SolverOptions options, int M, int K) {
        Matrix initSol = options.init_sol;
        if (initSol == null || initSol.isEmpty()) {
            sys.x0 = null;
            return;
        }
        if (sys.form.equals("J")) {
            sys.x0 = new double[initSol.length()];
            for (int i = 0; i < initSol.length(); i++) {
                sys.x0[i] = initSol.get(i);
            }
        } else {
            int totalStates = 0;
            for (int ist = 0; ist < M; ist++) {
                for (int r = 0; r < K; r++) {
                    int np = (int) sn.phases.get(ist, r);
                    totalStates += (np == 0) ? 1 : np;
                }
            }
            double[] x0Build = new double[totalStates];
            int state = 0;
            int initIdx = 0;
            for (int ist = 0; ist < M; ist++) {
                for (int r = 0; r < K; r++) {
                    int nPhases = (int) sn.phases.get(ist, r);
                    if (nPhases == 0) {
                        state++;
                    } else {
                        for (int k = 0; k < nPhases; k++) {
                            if (Double.isNaN(sn.rates.get(ist, r))) {
                                x0Build[state] = 0.0;
                            } else {
                                x0Build[state] = initSol.get(initIdx);
                                initIdx++;
                            }
                            state++;
                        }
                    }
                }
            }
            sys.x0 = new double[sys.nstates];
            for (int i = 0; i < sys.nstates; i++) {
                sys.x0[i] = x0Build[sys.keep[i]];
                if (sys.isSource[i]) {
                    sys.x0[i] = 0.0;
                }
            }
        }
    }

    /**
     * Render the LaTeX document for the given system.
     *
     * @param sys       structural description built by build()
     * @param options   solver options (hide_immediate remark)
     * @param modelName model name shown in the document
     * @param notation  "scalar" or "matrix"
     * @return LaTeX source
     */
    public static String render(SymODEs sys, SolverOptions options, String modelName, String notation) {
        if (!notation.equals("scalar") && !notation.equals("matrix")) {
            throw new RuntimeException("Unknown notation '" + notation + "'. Valid notations: scalar, matrix.");
        }
        int n = sys.nstates;
        List<String> L = new ArrayList<String>();
        L.add("% Mean-field fluid ODE system exported by LINE SolverFLD");
        L.add("% model: " + modelName);
        L.add("% method: " + sys.method);
        if (sys.form.equals("W")) {
            L.add("% form: dx/dt = W^T*theta(x) + lambda");
        } else {
            L.add("% form: dx/dt = J*r(x)");
        }
        L.add("% notation: " + notation);
        L.add("% nstates: " + n);
        if (sys.form.equals("J")) {
            L.add("% nevents: " + sys.nevents);
        }
        for (int s = 0; s < n; s++) {
            L.add(String.format(Locale.US, "%% STATE %d station=%s class=%s phase=%d", s + 1,
                    sys.stationNames.get(sys.stateStation[s]), sys.classNames.get(sys.stateClass[s]), sys.statePhase[s]));
        }
        if (sys.form.equals("J")) {
            for (int e = 0; e < sys.nevents; e++) {
                L.add(String.format(Locale.US, "%% EVENT %d var=%d type=%s coeff=%s", e + 1,
                        sys.eventVar[e] + 1, sys.factorType[e], cformat(sys.coeff[e], 15)));
            }
        }

        L.add("\\documentclass{article}");
        L.add("\\usepackage{amsmath}");
        L.add("\\usepackage[margin=2.5cm]{geometry}");
        L.add("\\allowdisplaybreaks");
        L.add("\\setcounter{MaxMatrixCols}{500}");
        L.add("\\begin{document}");
        L.add("\\section*{Mean-field fluid ODE system}");
        L.add(String.format(Locale.US,
                "\\noindent Model: \\texttt{%s}. Solver: \\texttt{SolverFLD}, method \\texttt{%s}, %s notation.",
                texesc(modelName), texesc(sys.method), notation));
        if (sys.form.equals("W")) {
            L.add(String.format(Locale.US,
                    "The system has %d state variables and reads $\\frac{\\mathrm{d}\\mathbf{x}}{\\mathrm{d}t} = W^{\\top}\\theta(\\mathbf{x}) + \\boldsymbol{\\lambda}$.",
                    n));
        } else {
            L.add(String.format(Locale.US,
                    "The system has %d state variables and %d events and reads $\\frac{\\mathrm{d}\\mathbf{x}}{\\mathrm{d}t} = J\\,r(\\mathbf{x})$.",
                    n, sys.nevents));
        }

        L.add("\\subsection*{State variables}");
        L.add("Each state variable $x_{s}$ is the mean number of jobs of a class in a service phase at a station:");
        L.add("\\begin{center}");
        int chunk = 48;
        for (int s0 = 0; s0 < n; s0 += chunk) {
            int s1 = Math.min(n, s0 + chunk);
            L.add("\\begin{tabular}{rlll}");
            L.add("\\hline");
            L.add("$s$ & station & class & phase\\\\");
            L.add("\\hline");
            for (int s = s0; s < s1; s++) {
                L.add(String.format(Locale.US, "%d & \\texttt{%s} & \\texttt{%s} & %d\\\\", s + 1,
                        texesc(sys.stationNames.get(sys.stateStation[s])),
                        texesc(sys.classNames.get(sys.stateClass[s])), sys.statePhase[s]));
            }
            L.add("\\hline");
            L.add("\\end{tabular}");
            if (s1 < n) {
                L.add("\\par\\medskip");
            }
        }
        L.add("\\end{center}");

        List<Integer> usedStations = new ArrayList<Integer>();
        for (int s = 0; s < n; s++) {
            if (!usedStations.contains(sys.stateStation[s])) {
                usedStations.add(sys.stateStation[s]);
            }
        }
        java.util.Collections.sort(usedStations);
        L.add("\\begin{center}");
        L.add("\\begin{tabular}{rlll}");
        L.add("\\hline");
        L.add("$i$ & station & scheduling & $S_{i}$\\\\");
        L.add("\\hline");
        for (int u = 0; u < usedStations.size(); u++) {
            int i = usedStations.get(u);
            L.add(String.format(Locale.US, "%d & \\texttt{%s} & %s & $%s$\\\\", i + 1,
                    texesc(sys.stationNames.get(i)), texesc(sys.schedNames.get(i)), fmtnum(sys.S[i])));
        }
        L.add("\\hline");
        L.add("\\end{tabular}");
        L.add("\\end{center}");

        double[][] T = new double[n][n];
        String[] varType = new String[n];
        int[] varStation = new int[n];
        int[] varClass = new int[n];
        int[][] varOthers = new int[n][];
        double[] constTerm = new double[n];
        buildTerms(sys, T, varType, varStation, varClass, varOthers, constTerm);

        List<String> defs = buildDefs(sys, varType, varStation, varClass);
        if (!defs.isEmpty()) {
            L.add("\\subsection*{Definitions}");
            L.add("\\begin{align*}");
            for (int d = 0; d < defs.size(); d++) {
                L.add(defs.get(d));
            }
            L.add("\\end{align*}");
        }

        if (notation.equals("scalar")) {
            L.add("\\subsection*{ODE system (scalar notation)}");
            L.add("\\begin{align}");
            for (int s = 0; s < n; s++) {
                L.add(renderEquation(sys, s, T, varType, varStation, varClass, varOthers, constTerm, s == n - 1));
            }
            L.add("\\end{align}");
        } else {
            L.add("\\subsection*{ODE system (matrix notation)}");
            if (sys.form.equals("W")) {
                boolean haveLambda = false;
                for (int s = 0; s < n; s++) {
                    if (sys.Alambda[s] != 0) {
                        haveLambda = true;
                    }
                }
                L.add("\\begin{equation}");
                if (haveLambda) {
                    L.add("\\frac{\\mathrm{d}\\mathbf{x}}{\\mathrm{d}t} = W^{\\top}\\,\\theta(\\mathbf{x}) + \\boldsymbol{\\lambda}");
                } else {
                    L.add("\\frac{\\mathrm{d}\\mathbf{x}}{\\mathrm{d}t} = W^{\\top}\\,\\theta(\\mathbf{x})");
                }
                L.add("\\end{equation}");
                L.add("with $\\theta_{s}(\\mathbf{x})$ given componentwise by");
                L.add("\\begin{equation*}");
                StringBuilder theta = new StringBuilder("\\theta(\\mathbf{x}) = \\begin{bmatrix} ");
                for (int v = 0; v < n; v++) {
                    if (v > 0) {
                        theta.append(" \\\\ ");
                    }
                    if (sys.isSource[v]) {
                        theta.append("0");
                    } else {
                        theta.append(factorTex(sys, v, varType[v], varStation[v], varClass[v], varOthers[v]));
                    }
                }
                theta.append(" \\end{bmatrix}");
                L.add(theta.toString());
                L.add("\\end{equation*}");
                L.add("and");
                L.add("\\begin{equation*}");
                Matrix WT = sys.W.transpose();
                double[][] wt = new double[n][n];
                for (int i = 0; i < n; i++) {
                    for (int j = 0; j < n; j++) {
                        wt[i][j] = WT.get(i, j);
                    }
                }
                L.add("W^{\\top} = " + renderNumMatrix(wt));
                L.add("\\end{equation*}");
                if (haveLambda) {
                    L.add("\\begin{equation*}");
                    L.add("\\boldsymbol{\\lambda} = " + renderNumVector(sys.Alambda) + "^{\\top}");
                    L.add("\\end{equation*}");
                }
            } else {
                L.add("\\begin{equation}");
                L.add("\\frac{\\mathrm{d}\\mathbf{x}}{\\mathrm{d}t} = J\\,r(\\mathbf{x})");
                L.add("\\end{equation}");
                L.add("with stoichiometry matrix");
                L.add("\\begin{equation*}");
                double[][] J = new double[n][sys.nevents];
                for (int e = 0; e < sys.nevents; e++) {
                    J[sys.eventFrom[e]][e] -= 1;
                    J[sys.eventTo[e]][e] += 1;
                }
                L.add("J = " + renderNumMatrix(J));
                L.add("\\end{equation*}");
                L.add("and event rate functions");
                L.add("\\begin{align*}");
                for (int e = 0; e < sys.nevents; e++) {
                    String fstr = factorTex(sys, sys.eventVar[e], sys.factorType[e], sys.factorStation[e],
                            sys.factorClass[e], sys.factorOthers[e]);
                    L.add(String.format(Locale.US, "r_{%d}(\\mathbf{x}) &= %s%s", e + 1,
                            termTex(sys.coeff[e], fstr), (e < sys.nevents - 1) ? "\\\\" : ""));
                }
                L.add("\\end{align*}");
            }
        }

        if (sys.x0 != null) {
            L.add("\\subsection*{Initial condition}");
            L.add("\\begin{equation*}");
            L.add("\\mathbf{x}(0) = " + renderNumVector(sys.x0) + "^{\\top}");
            L.add("\\end{equation*}");
        }

        L.add("\\subsection*{Remarks}");
        L.add("\\begin{itemize}");
        L.add("\\item For each station $i$, $n_{i}(\\mathbf{x})$ denotes the total mass at the station and $S_{i}$ the number of servers (infinite-server stations use the closed job population, $\\infty$ denotes infinity).");
        L.add("\\item The numerical solver regularizes vanishing denominators with a small positive constant; these regularizations are omitted here.");
        if (sys.form.equals("J")) {
            boolean anyFcfsw = false;
            boolean anyDps = false;
            for (int e = 0; e < sys.nevents; e++) {
                if (sys.factorType[e].equals("fcfsw") || sys.factorType[e].equals("fcfsws")) {
                    anyFcfsw = true;
                }
                if (sys.factorType[e].equals("dpsmin")) {
                    anyDps = true;
                }
            }
            if (anyFcfsw) {
                L.add("\\item At FCFS stations, the mean phase residence times $w_{u} = -1/[D_{0}]_{kk}$ weight the backlog $\\hat{n}_{i}$; the factors $w_{u}$ of the departing phases are folded into the rate coefficients.");
            }
            if (anyDps) {
                L.add("\\item At DPS stations, weights are normalized to sum to one and the weight $w_{ir}$ of the departing class is folded into the rate coefficient; the class shares $w_{ir}x/\\tilde{n}_{i}$ divide the station capacity $\\min(n_{i},S_{i})$, so they sum to one whenever the station is busy.");
            }
        }
        boolean anyFcfsStation = false;
        for (int s = 0; s < n; s++) {
            if (sys.sched[sys.stateStation[s]] == SchedStrategy.FCFS) {
                anyFcfsStation = true;
            }
        }
        if (anyFcfsStation && (sys.method.equals("matrix") || sys.method.equals("closing"))) {
            L.add("\\item For FCFS stations with non-exponential service, the solver may iteratively re-fit the service distributions (non-exponential approximation); the exported system uses the nominal model parameters.");
        }
        if (options != null && options.config != null && options.config.hide_immediate) {
            L.add("\\item \\texttt{hide\\_immediate} is enabled in the solver options: the numerical integration may further eliminate immediate transitions by state-space reduction; the exported system is the unreduced one.");
        }
        L.add("\\end{itemize}");
        L.add("\\end{document}");

        StringBuilder tex = new StringBuilder();
        for (int i = 0; i < L.size(); i++) {
            tex.append(L.get(i)).append("\n");
        }
        return tex.toString();
    }

    private static void buildTerms(SymODEs sys, double[][] T, String[] varType, int[] varStation,
                                   int[] varClass, int[][] varOthers, double[] constTerm) {
        int n = sys.nstates;
        if (sys.form.equals("W")) {
            for (int s = 0; s < n; s++) {
                for (int v = 0; v < n; v++) {
                    T[s][v] = sys.isSource[v] ? 0.0 : sys.W.get(v, s);
                }
                constTerm[s] = sys.Alambda[s];
            }
            for (int v = 0; v < n; v++) {
                if (!sys.isSource[v]) {
                    int i = sys.stateStation[v];
                    if (sys.isInfStation != null && sys.isInfStation[i]) {
                        varType[v] = "lin";
                    } else {
                        varType[v] = sys.smoothing;
                    }
                    varStation[v] = i;
                    varClass[v] = sys.stateClass[v];
                }
            }
        } else {
            for (int e = 0; e < sys.nevents; e++) {
                int v = sys.eventVar[e];
                T[sys.eventFrom[e]][v] -= sys.coeff[e];
                T[sys.eventTo[e]][v] += sys.coeff[e];
                if (varType[v] == null) {
                    varType[v] = sys.factorType[e];
                    varStation[v] = sys.factorStation[e];
                    varClass[v] = sys.factorClass[e];
                    varOthers[v] = sys.factorOthers[e];
                }
            }
        }
    }

    private static List<String> buildDefs(SymODEs sys, String[] varType, int[] varStation,
                                          int[] varClass) {
        List<String> defs = new ArrayList<String>();
        int n = sys.nstates;
        int M = sys.stationNames.size();
        int K = sys.classNames.size();
        boolean[] needN = new boolean[M];
        boolean[] needNT = new boolean[M];
        boolean[] needNH = new boolean[M];
        String[] gdef = new String[M];
        for (int v = 0; v < n; v++) {
            String f = varType[v];
            if (f == null) {
                continue;
            }
            int i = varStation[v];
            if (f.equals("min")) {
                needN[i] = true;
                gdef[i] = String.format(Locale.US,
                        "g_{%d}(\\mathbf{x}) &= \\frac{\\min(n_{%d}(\\mathbf{x}),\\, %s)}{n_{%d}(\\mathbf{x})}",
                        i + 1, i + 1, fmtnum(sys.S[i]), i + 1);
            } else if (f.equals("pnorm")) {
                needN[i] = true;
                gdef[i] = String.format(Locale.US,
                        "g_{%d}(\\mathbf{x}) &= \\Bigl(1 + \\bigl(n_{%d}(\\mathbf{x})/%s\\bigr)^{%s}\\Bigr)^{-1/%s}",
                        i + 1, i + 1, fmtnum(sys.S[i]), fmtnum(sys.pstar[i]), fmtnum(sys.pstar[i]));
            } else if (f.equals("dpsmin")) {
                needN[i] = true;
                needNT[i] = true;
                gdef[i] = String.format(Locale.US,
                        "g_{%d}(\\mathbf{x}) &= \\frac{\\min(n_{%d}(\\mathbf{x}),\\, %s)}{\\tilde{n}_{%d}(\\mathbf{x})}",
                        i + 1, i + 1, fmtnum(sys.S[i]), i + 1);
            } else if (f.equals("dpspw")) {
                needN[i] = true;
                needNT[i] = true;
            } else if (f.equals("fcfsw") || f.equals("fcfsws")) {
                needN[i] = true;
                needNH[i] = true;
                if (f.equals("fcfsw")) {
                    gdef[i] = String.format(Locale.US,
                            "g_{%d}(\\mathbf{x}) &= \\frac{\\min(n_{%d}(\\mathbf{x}),\\, %s)}{\\hat{n}_{%d}(\\mathbf{x})}",
                            i + 1, i + 1, fmtnum(sys.S[i]), i + 1);
                } else {
                    gdef[i] = String.format(Locale.US,
                            "g_{%d}(\\mathbf{x}) &= \\frac{\\mathrm{softmin}\\bigl(n_{%d}(\\mathbf{x}),\\, %s\\bigr)}{\\hat{n}_{%d}(\\mathbf{x})}",
                            i + 1, i + 1, fmtnum(sys.S[i]), i + 1);
                }
            }
        }
        for (int i = 0; i < M; i++) {
            if (needN[i]) {
                StringBuilder sum = new StringBuilder();
                for (int v = 0; v < n; v++) {
                    if (sys.stateStation[v] == i) {
                        if (sum.length() > 0) {
                            sum.append(" + ");
                        }
                        sum.append("x_{").append(v + 1).append("}");
                    }
                }
                defs.add(String.format(Locale.US, "n_{%d}(\\mathbf{x}) &= %s\\\\", i + 1, sum.toString()));
            }
        }
        for (int i = 0; i < M; i++) {
            if (needNT[i]) {
                StringBuilder parts = new StringBuilder();
                for (int r = 0; r < K; r++) {
                    StringBuilder sum = new StringBuilder();
                    for (int v = 0; v < n; v++) {
                        if (sys.stateStation[v] == i && sys.stateClass[v] == r) {
                            if (sum.length() > 0) {
                                sum.append(" + ");
                            }
                            sum.append("x_{").append(v + 1).append("}");
                        }
                    }
                    if (sum.length() > 0) {
                        if (parts.length() > 0) {
                            parts.append(" + ");
                        }
                        parts.append(String.format(Locale.US, "%s\\,(%s)", fmtnum(sys.dpsw[i][r]), sum.toString()));
                    }
                }
                defs.add(String.format(Locale.US, "\\tilde{n}_{%d}(\\mathbf{x}) &= %s\\\\", i + 1, parts.toString()));
            }
        }
        for (int i = 0; i < M; i++) {
            if (needNH[i]) {
                StringBuilder parts = new StringBuilder();
                for (int v = 0; v < n; v++) {
                    if (sys.stateStation[v] == i) {
                        if (parts.length() > 0) {
                            parts.append(" + ");
                        }
                        parts.append(String.format(Locale.US, "%s\\,x_{%d}", fmtnum(sys.fcfsPhaseW[v]), v + 1));
                    }
                }
                defs.add(String.format(Locale.US, "\\hat{n}_{%d}(\\mathbf{x}) &= %s\\\\", i + 1, parts.toString()));
            }
        }
        for (int i = 0; i < M; i++) {
            if (gdef[i] != null) {
                defs.add(gdef[i] + "\\\\");
            }
        }
        for (int v = 0; v < n; v++) {
            if (varType[v] != null && varType[v].equals("dpspw")) {
                int i = varStation[v];
                int r = varClass[v];
                String d = String.format(Locale.US,
                        "g_{%d,%d}(\\mathbf{x}) &= \\begin{cases} 1 & n_{%d}(\\mathbf{x}) \\le %s\\\\ \\dfrac{%s}{\\tilde{n}_{%d}(\\mathbf{x})} & n_{%d}(\\mathbf{x}) > %s \\end{cases}\\\\",
                        i + 1, r + 1, i + 1, fmtnum(sys.S[i]), fmtnum(sys.S[i] * sys.dpsw[i][r]), i + 1, i + 1, fmtnum(sys.S[i]));
                if (!defs.contains(d)) {
                    defs.add(d);
                }
            }
        }
        boolean anySoftmin = false;
        for (int v = 0; v < n; v++) {
            if (varType[v] != null && varType[v].equals("fcfsws")) {
                anySoftmin = true;
            }
        }
        if (anySoftmin) {
            defs.add(String.format(Locale.US,
                    "\\mathrm{softmin}(a,b) &= \\frac{a\\,e^{-\\alpha a} + b\\,e^{-\\alpha b}}{e^{-\\alpha a} + e^{-\\alpha b}}, \\qquad \\alpha = %s\\\\",
                    fmtnum(sys.alpha)));
        }
        if (!defs.isEmpty()) {
            String last = defs.get(defs.size() - 1);
            if (last.endsWith("\\\\")) {
                defs.set(defs.size() - 1, last.substring(0, last.length() - 2));
            }
        }
        return defs;
    }

    private static String renderEquation(SymODEs sys, int sidx, double[][] T, String[] varType,
                                         int[] varStation, int[] varClass, int[][] varOthers,
                                         double[] constTerm, boolean isLast) {
        List<double[]> termC = new ArrayList<double[]>();
        List<String> termF = new ArrayList<String>();
        for (int v = 0; v < sys.nstates; v++) {
            double c = T[sidx][v];
            if (c != 0) {
                termC.add(new double[]{c});
                termF.add(factorTex(sys, v, varType[v], varStation[v], varClass[v], varOthers[v]));
            }
        }
        if (constTerm[sidx] != 0) {
            termC.add(new double[]{constTerm[sidx]});
            termF.add("");
        }
        String rhs;
        if (termC.isEmpty()) {
            rhs = "0";
        } else {
            StringBuilder parts = new StringBuilder();
            for (int k = 0; k < termC.size(); k++) {
                double c = termC.get(k)[0];
                String body = termTex(Math.abs(c), termF.get(k));
                if (k == 0) {
                    parts.append(c < 0 ? "-" : "").append(body);
                } else {
                    parts.append(c < 0 ? " - " : " + ").append(body);
                }
                if ((k + 1) % 4 == 0 && k < termC.size() - 1) {
                    parts.append("\\nonumber\\\\\n&\\quad ");
                }
            }
            rhs = parts.toString();
        }
        return String.format(Locale.US, "\\frac{\\mathrm{d}x_{%d}}{\\mathrm{d}t} &= %s%s", sidx + 1, rhs,
                isLast ? "" : "\\\\");
    }

    private static String factorTex(SymODEs sys, int v, String ftype, int station, int classIdx, int[] others) {
        if (ftype.equals("lin")) {
            return String.format(Locale.US, "x_{%d}", v + 1);
        } else if (ftype.equals("min") || ftype.equals("pnorm") || ftype.equals("dpsmin")
                || ftype.equals("fcfsw") || ftype.equals("fcfsws")) {
            return String.format(Locale.US, "x_{%d}\\,g_{%d}(\\mathbf{x})", v + 1, station + 1);
        } else if (ftype.equals("dpspw")) {
            return String.format(Locale.US, "x_{%d}\\,g_{%d,%d}(\\mathbf{x})", v + 1, station + 1, classIdx + 1);
        } else { // ext1
            if (others == null || others.length == 0) {
                return ""; // single-phase source class: constant unit mass
            }
            StringBuilder sum = new StringBuilder();
            for (int u = 0; u < others.length; u++) {
                if (u > 0) {
                    sum.append(" - ");
                }
                sum.append("x_{").append(others[u] + 1).append("}");
            }
            return String.format(Locale.US, "\\bigl(1 - %s\\bigr)", sum.toString());
        }
    }

    private static String termTex(double c, String fstr) {
        if (fstr.isEmpty()) {
            return fmtnum(c);
        } else if (c == 1) {
            return fstr;
        } else {
            return String.format(Locale.US, "%s\\,%s", fmtnum(c), fstr);
        }
    }

    private static String renderNumMatrix(double[][] A) {
        int m = A.length;
        int ncol = (m > 0) ? A[0].length : 0;
        StringBuilder body = new StringBuilder();
        for (int r = 0; r < m; r++) {
            if (r > 0) {
                body.append(" \\\\ ");
            }
            for (int c = 0; c < A[r].length; c++) {
                if (c > 0) {
                    body.append(" & ");
                }
                body.append(fmtnum(A[r][c]));
            }
        }
        if (Math.max(m, ncol) > 12) {
            return String.format(Locale.US, "{\\scriptsize\\begin{bmatrix} %s \\end{bmatrix}}", body.toString());
        }
        return String.format(Locale.US, "\\begin{bmatrix} %s \\end{bmatrix}", body.toString());
    }

    private static String renderNumVector(double[] v) {
        StringBuilder body = new StringBuilder();
        for (int i = 0; i < v.length; i++) {
            if (i > 0) {
                body.append(" & ");
            }
            body.append(fmtnum(v[i]));
        }
        return String.format(Locale.US, "\\begin{pmatrix} %s \\end{pmatrix}", body.toString());
    }

    /** Compact LaTeX-safe number formatting, matching MATLAB's %.8g style. */
    static String fmtnum(double v) {
        if (Double.isInfinite(v)) {
            return v > 0 ? "\\infty" : "-\\infty";
        }
        if (v == Math.rint(v) && Math.abs(v) < 1e15) {
            return String.format(Locale.US, "%d", (long) v);
        }
        return cformat(v, 8);
    }

    /** C-style %g formatting (trailing zeros stripped), used for parity with MATLAB sprintf. */
    static String cformat(double v, int prec) {
        if (v == Math.rint(v) && Math.abs(v) < 1e15) {
            return String.format(Locale.US, "%d", (long) v);
        }
        String s = String.format(Locale.US, "%." + prec + "g", v);
        if (s.contains("e") || s.contains("E")) {
            String[] parts = s.split("[eE]");
            String mant = parts[0];
            if (mant.contains(".")) {
                mant = mant.replaceAll("0+$", "").replaceAll("\\.$", "");
            }
            s = mant + "e" + parts[1];
        } else if (s.contains(".")) {
            s = s.replaceAll("0+$", "").replaceAll("\\.$", "");
        }
        return s;
    }

    private static String texesc(String s) {
        return s.replaceAll("([_%&#])", "\\\\$1");
    }

    // -----------------------------------------------------------------------
    // Symbolic drift, for the computer algebra backend
    // -----------------------------------------------------------------------

    /** Smooth factor types; every other type carries a min or a branch. */
    private static final List<String> SMOOTH_FACTORS =
            java.util.Arrays.asList("lin", "ext1", "fcfsws");

    /**
     * State variable names of the exported drift, x1 ... xn.
     *
     * @param sys the ODE system
     * @return the variable names
     */
    public static List<String> stateVariables(SymODEs sys) {
        List<String> vars = new ArrayList<String>();
        for (int s = 0; s < sys.nstates; s++) {
            vars.add("x" + (s + 1));
        }
        return vars;
    }

    /**
     * Right-hand side of the ODE system as expression strings, one per state
     * variable, in the format the symbolic backend parses.
     *
     * <p>ONLY SMOOTH DRIFTS ARE EXPORTED. The default, matrix, closing and
     * statedep methods scale rates by min(n_i, S_i), which is not
     * differentiable at n_i = S_i, so their Jacobian does not exist there;
     * emitting a one-sided derivative would be a silent lie exactly at the
     * regime switch that matters. Use the p-norm smoothing
     * ({@code options.config.pstar}) or the softmin method.</p>
     *
     * <p>FineTol is carried in exactly the places the integrated systems put
     * it, and nowhere else: the p-norm drift offsets the station total, the
     * softmin drift offsets the phase-weighted total but not the plain station
     * total that feeds the softmin, and the closing rates offset neither.
     * Mirrors {@code @SolverFLD/getSymbolicDrift.m}.</p>
     *
     * @param sys the ODE system
     * @return one expression per state variable
     * @throws RuntimeException if the drift is not differentiable
     */
    public static List<String> symbolicDrift(SymODEs sys) {
        List<String> vars = stateVariables(sys);
        String eps0 = num(GlobalConstants.FineTol);
        if ("W".equals(sys.form)) {
            if (!"pnorm".equals(sys.smoothing)) {
                throw new RuntimeException(
                        "The drift of this method scales rates by min(n_i, S_i), which is not "
                                + "differentiable at n_i = S_i, so it has no Jacobian there. Set "
                                + "options.config.pstar to use the p-norm smoothing, or use the "
                                + "'softmin' method.");
            }
            return wformDrift(sys, vars, eps0);
        }
        if ("J".equals(sys.form)) {
            return jformDrift(sys, vars, eps0);
        }
        throw new RuntimeException("unsupported ODE form '" + sys.form + "'");
    }

    private static List<String> wformDrift(SymODEs sys, List<String> vars, String eps0) {
        int n = sys.nstates;
        String[] theta = new String[n];
        for (int s = 0; s < n; s++) {
            if (sys.isSource[s]) {
                theta[s] = "0";
                continue;
            }
            int i = sys.stateStation[s];
            String ni = stationSum(sys, i, vars, eps0);
            double S = sys.S[i];
            double p = sys.pstar[i];
            if (S <= 0 || p <= 0) {
                theta[s] = vars.get(s);
            } else {
                theta[s] = vars.get(s) + "/(1 + (" + ni + "/" + num(S) + ")^" + num(p)
                        + ")^(1/" + num(p) + ")";
            }
        }

        List<String> rhs = new ArrayList<String>();
        for (int s = 0; s < n; s++) {
            List<String> terms = new ArrayList<String>();
            for (int t = 0; t < n; t++) {
                double w = sys.W.get(t, s);
                if (w == 0 || "0".equals(theta[t])) {
                    continue;
                }
                terms.add("(" + num(w) + ")*(" + theta[t] + ")");
            }
            if (sys.Alambda[s] != 0) {
                terms.add(num(sys.Alambda[s]));
            }
            rhs.add(terms.isEmpty() ? "0" : join(terms, " + "));
        }
        return rhs;
    }

    private static List<String> jformDrift(SymODEs sys, List<String> vars, String eps0) {
        String[] rate = new String[sys.nevents];
        for (int e = 0; e < sys.nevents; e++) {
            String ftype = sys.factorType[e];
            if (!SMOOTH_FACTORS.contains(ftype)) {
                throw new RuntimeException("Event " + (e + 1) + " scales its rate by the "
                        + "non-smooth factor '" + ftype + "', which has no derivative where the "
                        + "regime switches, so the system has no Jacobian. Use the 'softmin' "
                        + "method, or the p-norm smoothing of the 'matrix' method.");
            }
            String v = vars.get(sys.eventVar[e]);
            String factor;
            if ("lin".equals(ftype)) {
                factor = v;
            } else if ("ext1".equals(ftype)) {
                int[] others = sys.factorOthers[e];
                if (others == null || others.length == 0) {
                    factor = "1";
                } else {
                    List<String> parts = new ArrayList<String>();
                    for (int k = 0; k < others.length; k++) {
                        parts.add(vars.get(others[k]));
                    }
                    factor = "(1 - (" + join(parts, " + ") + "))";
                }
            } else {
                int i = sys.factorStation[e];
                String ni = stationSum(sys, i, vars, "0");
                String nhat = phaseWeightedStationSum(sys, i, vars, eps0);
                factor = v + "*(" + softmin(ni, num(sys.S[i]), sys.alpha) + ")/(" + nhat + ")";
            }
            rate[e] = "(" + num(sys.coeff[e]) + ")*(" + factor + ")";
        }

        List<String> rhs = new ArrayList<String>();
        for (int s = 0; s < sys.nstates; s++) {
            List<String> terms = new ArrayList<String>();
            for (int e = 0; e < sys.nevents; e++) {
                double j = 0;
                if (sys.eventFrom[e] == s) {
                    j -= 1;
                }
                if (sys.eventTo[e] == s) {
                    j += 1;
                }
                if (j == 0) {
                    continue;
                }
                terms.add("(" + num(j) + ")*(" + rate[e] + ")");
            }
            rhs.add(terms.isEmpty() ? "0" : join(terms, " + "));
        }
        return rhs;
    }

    /**
     * Smooth minimum in its weighted-average form. {@code softmin.m} computes
     * the algebraically identical lo + gap*w/(1+w) purely to keep the exponent
     * argument non-positive, an overflow guard that means nothing symbolically
     * and would reintroduce the min/max branch this export exists to avoid.
     */
    private static String softmin(String x, String y, double alpha) {
        String a = num(alpha);
        return "((" + x + ")*exp(-(" + a + ")*(" + x + ")) + (" + y + ")*exp(-(" + a + ")*("
                + y + ")))/(exp(-(" + a + ")*(" + x + ")) + exp(-(" + a + ")*(" + y + ")))";
    }

    private static String stationSum(SymODEs sys, int i, List<String> vars, String offset) {
        List<String> parts = new ArrayList<String>();
        for (int k = 0; k < sys.nstates; k++) {
            if (sys.stateStation[k] == i) {
                parts.add(vars.get(k));
            }
        }
        return joinSum(parts, offset);
    }

    private static String phaseWeightedStationSum(SymODEs sys, int i, List<String> vars,
                                                  String offset) {
        List<String> parts = new ArrayList<String>();
        for (int k = 0; k < sys.nstates; k++) {
            if (sys.stateStation[k] == i) {
                double wt = sys.fcfsPhaseW[k];
                if (wt != 0) {
                    parts.add("(" + num(wt) + ")*" + vars.get(k));
                }
            }
        }
        return joinSum(parts, offset);
    }

    private static String joinSum(List<String> parts, String offset) {
        if (parts.isEmpty()) {
            return "(" + offset + ")";
        }
        if ("0".equals(offset)) {
            return "(" + join(parts, " + ") + ")";
        }
        return "(" + offset + " + " + join(parts, " + ") + ")";
    }

    private static String join(List<String> parts, String sep) {
        StringBuilder sb = new StringBuilder();
        for (int i = 0; i < parts.size(); i++) {
            if (i > 0) {
                sb.append(sep);
            }
            sb.append(parts.get(i));
        }
        return sb.toString();
    }

    /** Decimal text the symbolic backend reads as an exact rational. */
    private static String num(double v) {
        if (v == Math.rint(v) && !Double.isInfinite(v) && Math.abs(v) < 1e15) {
            return Long.toString((long) v);
        }
        return String.format(Locale.US, "%.17g", v);
    }

    private FluidODEsExporter() {
        // static utility class
    }
}
