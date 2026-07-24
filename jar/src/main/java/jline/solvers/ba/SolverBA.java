/*
 * Copyright (c) 2012-2026, QORE Lab, Imperial College London
 * All rights reserved.
 */
package jline.solvers.ba;

import jline.GlobalConstants;
import jline.lang.Network;
import jline.lang.NetworkStruct;
import jline.lang.constant.SchedStrategy;
import jline.lang.constant.SolverType;
import jline.solvers.NetworkBoundsTable;
import jline.solvers.NetworkSolver;
import jline.solvers.Solver;
import jline.solvers.SolverOptions;
import jline.solvers.SolverResult;
import jline.solvers.mva.MVAResult;
import jline.solvers.mva.analyzers.Solver_mva_bound_analyzer;
import jline.solvers.ctmc.SolverCTMC;
import jline.solvers.ctmc.analyzers.Solver_ctmc_qrf_analyzer;
import jline.api.mapqn.Mapqn_bnd_lr_pf;
import jline.util.matrix.Matrix;

import java.util.ArrayList;
import java.util.Arrays;
import java.util.List;

import static jline.io.InputOutput.line_error;
import static jline.io.InputOutput.mfilename;

/**
 * SolverBA is the dedicated bound-analysis solver for closed queueing networks.
 *
 * <p>Unlike SolverMVA (point estimates), each method returns an optimistic
 * ({@code .upper}) or pessimistic ({@code .lower}) throughput/queue-length
 * bound; the two sides of a family bracket the exact solution. Bounds need only
 * demands and populations, so the feature set is narrow (closed
 * product-form-parameterized models). Method families:
 * <ul>
 *   <li>Noniterative: aba, bjb, pb, gb, sb, mwba (Majumdar-Woodside).</li>
 *   <li>Hierarchical: pbh (Eager-Sevcik), pbk/bjbk (iterative PB(k)/BJB(k)),
 *       cbh (Dowdy), ssd (Suri-Dallery multiserver), cub (Kerola composite
 *       upper) with mbjb (multiclass BJB lower), sib (Srinivasan), ldbcmp
 *       (Anselmi-Cremonesi LD-BCMP lower).</li>
 *   <li>QRF: qrf.* (Quadratic Reduction Framework, LP-based, PH service), with
 *       aliases qr -&gt; qrf.mmi and lr -&gt; qrf.mmi.linear.</li>
 * </ul>
 */
public class SolverBA extends NetworkSolver {

    public SolverBA(Network model, String method) {
        super(model, "SolverBA", SolverBA.defaultOptions().method(method));
        this.sn = model.getStruct(false);
        this.result = new SolverResult();
    }

    public SolverBA(Network model, SolverOptions options) {
        super(model, "SolverBA", options);
        this.sn = model.getStruct(false);
        this.result = new SolverResult();
    }

    public SolverBA(Network model, Object... varargin) {
        this(model, SolverBA.defaultOptions());
        this.options = Solver.parseOptions(this.options, varargin);
    }

    public SolverBA(Network model) {
        super(model, "SolverBA", SolverBA.defaultOptions());
        this.sn = model.getStruct(false);
        this.result = new SolverResult();
    }

    public static SolverOptions defaultOptions() {
        SolverOptions options = new SolverOptions(SolverType.BA);
        options.method = "default";
        options.level = 2;
        return options;
    }

    public NetworkStruct getStruct() {
        if (this.sn == null) {
            this.sn = this.model.getStruct(false);
        }
        return this.sn;
    }

    @Override
    public void runAnalyzer() throws IllegalAccessException {
        if (this.options == null) {
            this.options = SolverBA.defaultOptions();
        }
        GlobalConstants.Verbose = options.verbose;
        if (this.sn == null) {
            this.sn = this.model.getStruct(false);
        }
        if (this.result == null) {
            this.result = new SolverResult();
        }
        if (this.sn.nclosedjobs <= 0) {
            line_error(mfilename(new Object() {}), "SolverBA supports closed queueing networks only.");
        }

        // 'default' selects the tightest noniterative upper bound (geometric);
        // 'qr' is a friendly alias for the QRF quadratic reduction bound. 'lr'
        // is the LP-based Linear Reduction bound (Mapqn_bnd_lr_pf, simplex) and
        // is NOT an alias of qrf.mmi.linear: that method is "linear" only in its
        // explicit Aeq/beq constraint representation, its objective being the
        // nonlinear MEM entropy. Bare 'lr' means 'lr.upper'.
        String method = this.options.method;
        if ("default".equals(method)) {
            method = "gb.upper";
        } else if ("lr".equals(method)) {
            method = "lr.upper";
        } else if ("qr".equals(method)) {
            method = "qrf.mmi";
        }

        SolverOptions bopts = this.options.copy();
        bopts.method = method;

        if ("lr.upper".equals(method) || "lr.lower".equals(method)) {
            MVAResult ret = lrBounds(this.sn, method);
            this.setAvgResults(ret.QN, ret.UN, ret.RN, ret.TN, new Matrix(0, 0), new Matrix(0, 0),
                    ret.CN, ret.XN, ret.runtime, method, 1);
        } else if (method.startsWith("qrf")) {
            SolverCTMC.AnalyzerResult r = Solver_ctmc_qrf_analyzer.solver_ctmc_qrf_analyzer(this.sn, bopts);
            this.setAvgResults(r.QN, r.UN, r.RN, r.TN, new Matrix(0, 0), new Matrix(0, 0),
                    r.CN, r.XN, r.runtime, method, 1);
        } else if (Arrays.asList(this.listValidMethods()).contains(method)) {
            MVAResult ret = Solver_mva_bound_analyzer.solver_mva_bound_analyzer(this.sn, bopts);
            this.setAvgResults(ret.QN, ret.UN, ret.RN, ret.TN, new Matrix(0, 0), new Matrix(0, 0),
                    ret.CN, ret.XN, ret.runtime, method, ret.iter);
        } else {
            line_error(mfilename(new Object() {}), "Unknown bound method '" + method
                    + "'. Valid methods: " + String.join(", ", this.listValidMethods()));
        }
    }

    /**
     * LP-based Linear Reduction bound. Solves the product-form linear-reduction
     * LP once per station (Mapqn_bnd_lr_pf, Apache SimplexSolver), minimizing
     * for "lr.lower" and maximizing for "lr.upper". Both senses are valid
     * bounds because the LP relaxation contains the exact solution. Unlike
     * qrf.mmi.linear -- whose objective is the nonlinear MEM entropy -- this
     * method is a pure LP end to end.
     *
     * @param sn     network structure (single-class closed, single-server)
     * @param method "lr.upper" or "lr.lower"
     * @return per-station bound metrics
     */
    private static MVAResult lrBounds(NetworkStruct sn, String method) {
        long tstart = System.currentTimeMillis();
        if (sn.nclasses != 1 || sn.nclosedjobs <= 0) {
            line_error(mfilename(new Object() {}),
                    "Method lr supports single-class closed networks only.");
        }
        int M = sn.nstations;
        for (int i = 0; i < M; i++) {
            if (sn.sched.get(sn.stations.get(i)) == SchedStrategy.INF) {
                line_error(mfilename(new Object() {}),
                        "Method lr does not support delay (infinite-server) stations.");
            }
            if (sn.nservers.get(i, 0) > 1) {
                line_error(mfilename(new Object() {}),
                        "Unsupported method for a model with multi-server stations.");
            }
        }
        String sense = "lr.lower".equals(method) ? "min" : "max";
        int N = (int) Math.round(sn.njobs.get(0, 0));

        Matrix visits = sn.visits.get(0);
        double[] V = new double[M];
        double totV = 0.0;
        for (int i = 0; i < M; i++) {
            V[i] = visits.get(i, 0);
            totV += V[i];
        }
        double[] mu = new double[M];
        for (int i = 0; i < M; i++) {
            mu[i] = sn.rates.get(i, 0);
        }
        double[][] r = new double[M][M];
        for (int i = 0; i < M; i++) {
            for (int j = 0; j < M; j++) {
                r[i][j] = totV > 0 ? V[j] / totV : 0.0;
            }
        }

        Mapqn_bnd_lr_pf.PFParameters params = new Mapqn_bnd_lr_pf.PFParameters(M, N, mu, r);
        Matrix UN = new Matrix(M, 1);
        double X = Double.POSITIVE_INFINITY;
        for (int i = 0; i < M; i++) {
            double u = Mapqn_bnd_lr_pf.solve(params, i + 1, sense).getObjectiveValue();
            UN.set(i, 0, u);
            if (V[i] > 0) {
                double cand = u * mu[i] / V[i];
                if (cand < X) {
                    X = cand;
                }
            }
        }
        if (!Double.isFinite(X)) {
            X = 0.0;
        }

        Matrix QN = new Matrix(M, 1);
        Matrix RN = new Matrix(M, 1);
        Matrix TN = new Matrix(M, 1);
        double sumD = 0.0;
        for (int i = 0; i < M; i++) {
            double t = V[i] * X;
            TN.set(i, 0, t);
            double res = "max".equals(sense) ? (1.0 / mu[i]) * N : 1.0 / mu[i];
            RN.set(i, 0, res);
            QN.set(i, 0, t * res);
            sumD += V[i] / mu[i];
        }
        Matrix CN = new Matrix(1, 1);
        CN.set(0, 0, "max".equals(sense) ? N * sumD : sumD);
        Matrix XN = new Matrix(1, 1);
        XN.set(0, 0, X);

        MVAResult ret = new MVAResult();
        ret.QN = QN;
        ret.UN = UN;
        ret.RN = RN;
        ret.TN = TN;
        ret.CN = CN;
        ret.XN = XN;
        ret.iter = 1;
        ret.runtime = (System.currentTimeMillis() - tstart) / 1000.0;
        return ret;
    }

    /**
     * Container for the {lower,upper} throughput/queue-length bracket of a
     * bound family. One-sided families (cub upper-only, mbjb/ldbcmp lower-only)
     * leave the missing side null.
     */
    public static class Bounds {
        public Matrix Tlower, Tupper, Qlower, Qupper;
    }

    /**
     * Returns the {lower,upper} bracket for the current method's family.
     *
     * @return the throughput/queue-length bracket
     */
    public Bounds getBounds() {
        String fam = this.options.method;
        int dot = fam.indexOf('.');
        if (dot >= 0) {
            fam = fam.substring(0, dot);
        }
        List<String> valid = Arrays.asList(this.listValidMethods());
        Bounds b = new Bounds();
        if (valid.contains(fam + ".lower")) {
            SolverResult rl = solveSide(fam + ".lower");
            b.Qlower = rl.QN;
            b.Tlower = rl.TN;
        }
        if (valid.contains(fam + ".upper")) {
            SolverResult ru = solveSide(fam + ".upper");
            b.Qupper = ru.QN;
            b.Tupper = ru.TN;
        }
        return b;
    }

    /**
     * Re-runs the solver for one side of the bracket.
     *
     * <p>The re-run instance inherits the CALLER'S FULL OPTION SET (level,
     * verbose, tolerances, ...) and only overrides the method. Constructing it
     * with just the method string would silently reset options.level to its
     * default of 2, so a hierarchical family (pbh/cbh/pbk/bjbk/sib) reached
     * through getBounds would never tighten as level is raised.
     *
     * @param method the fully qualified method for this side of the bracket
     * @return the solver result for that side
     */
    private SolverResult solveSide(String method) {
        SolverOptions opts = this.options.copy();
        opts.method = method;
        return new SolverBA(this.model, opts).getAvg();
    }

    /**
     * Returns a table of the {lower,upper} bracket per station and class, in
     * the layout of getAvgTable. Columns: Qlower, Qupper, Tlower, Tupper.
     *
     * <p>One-sided families (cub upper-only, mbjb/ldbcmp lower-only) carry NaN
     * on the missing side; NaN is preserved, never replaced by zero.
     *
     * @return the bracket table
     */
    public NetworkBoundsTable getBoundsTable() {
        return getBoundsTable(false);
    }

    /**
     * Returns a table of the {lower,upper} bracket per station and class.
     *
     * @param keepDisabled whether to retain station-class pairs whose bounds are all zero
     * @return the bracket table
     */
    public NetworkBoundsTable getBoundsTable(boolean keepDisabled) {
        Bounds b = getBounds();
        NetworkStruct s = getStruct();
        int M = s.nstations;
        int K = s.nclasses;
        List<Double> Qlow = new ArrayList<Double>();
        List<Double> Qupp = new ArrayList<Double>();
        List<Double> Tlow = new ArrayList<Double>();
        List<Double> Tupp = new ArrayList<Double>();
        List<String> className = new ArrayList<String>();
        List<String> stationName = new ArrayList<String>();
        for (int i = 0; i < M; i++) {
            for (int k = 0; k < K; k++) {
                double ql = boundAt(b.Qlower, i, k);
                double qu = boundAt(b.Qupper, i, k);
                double tl = boundAt(b.Tlower, i, k);
                double tu = boundAt(b.Tupper, i, k);
                // Mirror getAvgTable's drop of disabled station-class pairs, but
                // NaN-safe: keep the row when any value that is present is
                // nonzero, so an all-NaN side never removes the row.
                double[] vals = new double[]{ql, qu, tl, tu};
                boolean anyPresent = false;
                boolean anyNonZero = false;
                for (int v = 0; v < vals.length; v++) {
                    if (!Double.isNaN(vals[v])) {
                        anyPresent = true;
                        if (vals[v] != 0.0) {
                            anyNonZero = true;
                        }
                    }
                }
                if (keepDisabled || !anyPresent || anyNonZero) {
                    Qlow.add(ql);
                    Qupp.add(qu);
                    Tlow.add(tl);
                    Tupp.add(tu);
                    className.add(this.model.getClasses().get(k).getName());
                    stationName.add(this.model.getStations().get(i).getName());
                }
            }
        }
        NetworkBoundsTable t = new NetworkBoundsTable(Qlow, Qupp, Tlow, Tupp);
        t.setOptions(this.options);
        t.setClassNames(className);
        t.setStationNames(stationName);
        return t;
    }

    /**
     * Reads one entry of a bracket side, mapping an absent side to NaN.
     *
     * @param m the bracket side, null when the family is one-sided
     * @param i station index
     * @param k class index
     * @return the entry, or NaN when the side is absent
     */
    private static double boundAt(Matrix m, int i, int k) {
        if (m == null || m.isEmpty() || i >= m.getNumRows() || k >= m.getNumCols()) {
            return Double.NaN;
        }
        return m.get(i, k);
    }

    public String[] listValidMethods() {
        List<String> m = new ArrayList<String>(Arrays.asList(
                "default",
                "aba.upper", "aba.lower",
                "bjb.upper", "bjb.lower",
                "pb.upper", "pb.lower",
                "gb.upper", "gb.lower",
                "sb.upper", "sb.lower",
                "mwba.upper", "mwba.lower",
                "pbh.upper", "pbh.lower",
                "pbk.upper", "pbk.lower",
                "bjbk.upper", "bjbk.lower",
                "cbh.upper", "cbh.lower",
                "ssd.upper", "ssd.lower",
                "cub.upper", "mbjb.lower",
                "sib.upper", "sib.lower",
                "ldbcmp.lower",
                "qr", "lr", "lr.upper", "lr.lower",
                "qrf.mmi", "qrf.mem", "qrf.mmi.ld", "qrf.mmi.linear",
                "qrf.bas.mmi", "qrf.bas.mem", "qrf.bas", "qrf.rsrd"
        ));
        return m.toArray(new String[0]);
    }
}
