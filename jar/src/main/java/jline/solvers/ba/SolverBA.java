/*
 * Copyright (c) 2012-2026, QORE Lab, Imperial College London
 * All rights reserved.
 */
package jline.solvers.ba;

import jline.GlobalConstants;
import jline.lang.FeatureSet;
import jline.lang.Network;
import jline.lang.NetworkStruct;
import jline.lang.constant.ProcessType;
import jline.lang.constant.SchedStrategy;
import jline.lang.constant.SolverType;
import jline.solvers.NetworkBoundsTable;
import jline.solvers.NetworkPercTable;
import jline.solvers.NetworkSolver;
import jline.solvers.Solver;
import jline.solvers.SolverOptions;
import jline.solvers.SolverResult;
import jline.solvers.mva.MVAResult;
import jline.lang.constant.NodeType;
import jline.solvers.ba.analyzers.Solver_ba_bgt_analyzer;
import jline.solvers.ba.analyzers.Solver_ba_spnlp_analyzer;
import jline.solvers.ba.analyzers.Solver_ba_bpt_analyzer;
import jline.solvers.ba.analyzers.Solver_ba_snc_analyzer;
import jline.api.sn.SnHasBlocking;
import jline.api.sn.SnToQrfAlpha;
import jline.api.sn.SnToQrfBlocking;
import jline.api.snc.Snc_perc_backlog;
import jline.api.snc.Snc_perc_delay;
import jline.solvers.mva.analyzers.Solver_mva_bound_analyzer;
import jline.solvers.ctmc.SolverCTMC;
import jline.solvers.ctmc.analyzers.Solver_ctmc_qrf_analyzer;
import jline.api.mam.Map_mean;
import jline.api.mapqn.Mapqn_bnd_lr_mva;
import jline.api.mapqn.Mapqn_bnd_lr_pf;
import jline.api.mapqn.MVAVersionParameters;
import jline.lang.JobClass;
import jline.lang.nodes.Station;
import jline.util.matrix.MatrixCell;
import jline.api.pfqn.mva.Pfqn_ldbcmp;
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
 * product-form-parameterized models).
 *
 * <p>FINITE-BUFFER BLOCKING IS REFUSED, NOT BOUNDED. Needing only demands and a
 * population is the BCMP parameterization, which presumes UNBOUNDED buffers; a
 * buffer that binds couples the station occupancies and the resulting numbers
 * do not bracket the blocked model. {@code runAnalyzer} therefore gates on
 * {@code SnHasBlocking} and {@code listValidMethods} drops every blocking-blind
 * method. The exceptions are {@code qrf.bas*}/{@code qrf.rsrd}, which carry the
 * blocking tables explicitly. Use SolverMVA method {@code 'sqd'} for a point
 * estimate.
 *
 * <p>Method families:
 * <ul>
 *   <li>Noniterative: aba, bjb, pb, gb, sb, mwba (Majumdar-Woodside).</li>
 *   <li>Hierarchical: pbh (Eager-Sevcik), pbk/bjbk (iterative PB(k)/BJB(k)),
 *       cbh (Dowdy), ssd (Suri-Dallery multiserver), cub (Kerola composite
 *       upper) with mbjb (multiclass BJB lower), sib (Srinivasan), ldbcmp
 *       (Anselmi-Cremonesi LD-BCMP lower), scb (Dowdy et al. 1992 single-class
 *       bounds, the only family that does NOT bracket this model's own solution
 *       but the multiclass system the single-class model aggregates; kept out of
 *       auto.* for that reason).</li>
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
        // Every bound family here is parameterized by demands and a closed
        // population, except the three OPEN-network families 'bpt', 'bgt' and
        // 'snc', which are refused on a closed model instead, and the 'spnlp'
        // family, which is parameterized by a MARKING: whether that marking is
        // bounded is a question about the P-invariants of the net and not about
        // nclosedjobs, so Spn_lpbnd decides it rather than this gate.
        if (this.sn.nclosedjobs <= 0 && !this.options.method.startsWith("bpt")
                && !this.options.method.startsWith("bgt")
                && !this.options.method.startsWith("snc")
                && !this.options.method.startsWith("spnlp")) {
            line_error(mfilename(new Object() {}), "SolverBA supports closed queueing networks only.");
        }

        // 'default' selects the tightest noniterative upper bound (geometric);
        // 'qr' is a friendly alias for the QRF quadratic reduction bound. 'lr'
        // is the LP-based Linear Reduction bound (Mapqn_bnd_lr_pf, simplex) and
        // is NOT an alias of qrf.mmi.linear: that method is "linear" only in its
        // explicit Aeq/beq constraint representation, its objective being the
        // nonlinear MMI mutual information. Bare 'lr' means 'lr.upper'.
        String method = this.options.method;
        if ("default".equals(method)) {
            method = "gb.upper";
        } else if ("auto".equals(method)) {
            method = "auto.upper";
        } else if ("lr".equals(method)) {
            method = "lr.upper";
        } else if ("qr".equals(method)) {
            method = "qrf.mmi";
        }

        // Finite-buffer BLOCKING is outside the premises of every family here
        // except the QRF blocking bounds: the rest are parameterized by demands
        // and a population alone, which presumes unbounded buffers and a
        // product form that the truncation destroys. Refusing is not
        // conservatism -- on cqn_bas_blocking, gb.upper reports QLen 1.28 at a
        // station capped at 1 job. Gated after the aliases so 'default' is
        // judged as the gb.upper it resolves to.
        String blockingWhy = "";
        if (ignoresBlocking(method) && SnHasBlocking.snHasBlocking(this.sn)) {
            // A blocked model whose shape admits the QRF BAS bound gets it as
            // the DEFAULT rather than a refusal: 'qrf.bas' models the finite
            // buffer, and since SnToQrfBlocking derives its tables from the
            // model there is nothing left for the caller to supply.
            String req = this.options.method;
            if ("default".equals(req) || "auto".equals(req) || "auto.upper".equals(req)) {
                String[] routed = blockingDefault(this.sn);
                blockingWhy = routed[1];
                if (!routed[0].isEmpty()) {
                    method = routed[0];
                }
            }
        }

        if (ignoresBlocking(method) && SnHasBlocking.snHasBlocking(this.sn)) {
            String detail = blockingWhy.isEmpty() ? ""
                    : " The QRF blocking bounds do not apply here either: " + blockingWhy;
            line_error(mfilename(new Object() {}), "Method '" + this.options.method
                    + "' does not support finite-buffer blocking: every SolverBA bound family"
                    + " but the QRF blocking ones is parameterized by demands and a population"
                    + " alone, so it bounds the model as if its buffers were unbounded. Use"
                    + " SolverMVA with method 'sqd', an exact solver (CTMC, SSA, JMT, LDES), or"
                    + " the QRF blocking bounds 'qrf.bas'/'qrf.rsrd', which model the finite"
                    + " buffer." + detail);
        }

        // The STRUCTURAL premises of every bound family -- single-class closed,
        // fully closed, single-server -- are asked here and nowhere else, so
        // this run and the report gate supportsModelMethod give one answer
        // whichever a caller meets first. The premises a feature name CAN
        // express (a delay station, a closed class in an open family) live in
        // getMethodFeatureSet instead. Asked after the aliases so 'default' is
        // judged as the gb.upper it runs as.
        String structuralRefusal = methodRefusal(this.sn, method);
        if (!structuralRefusal.isEmpty()) {
            line_error(mfilename(new Object() {}), structuralRefusal);
        }

        SolverOptions bopts = this.options.copy();
        bopts.method = method;

        if ("lr.upper".equals(method) || "lr.lower".equals(method)) {
            MVAResult ret = lrBounds(this.sn, method);
            this.setAvgResults(ret.QN, ret.UN, ret.RN, ret.TN, new Matrix(0, 0), new Matrix(0, 0),
                    ret.CN, ret.XN, ret.runtime, method, 1);
        } else if (method.startsWith("mapamva")) {
            MVAResult ret = mapamvaBounds(this.sn, method);
            this.setAvgResults(ret.QN, ret.UN, ret.RN, ret.TN, new Matrix(0, 0), new Matrix(0, 0),
                    ret.CN, ret.XN, ret.runtime, method, 1);
        } else if (method.startsWith("bgt")) {
            // Piecewise-linear Lyapunov upper bound (Bertsimas-Gamarnik-
            // Tsitsiklis 2001) for multitype OPEN networks; the analyzer gates
            // itself, including the deterministic non-merging routing premise.
            MVAResult r = Solver_ba_bgt_analyzer.solver_ba_bgt_analyzer(this.sn, bopts);
            this.setAvgResults(r.QN, r.UN, r.RN, r.TN, new Matrix(0, 0), new Matrix(0, 0),
                    r.CN, r.XN, r.runtime, method, 1);
        } else if (method.startsWith("spnlp")) {
            // Moment-relaxation LP bounds (Liu 1998) for a stochastic Petri
            // net; the analyzer gates itself.
            MVAResult r = Solver_ba_spnlp_analyzer.solver_ba_spnlp_analyzer(this.sn, bopts);
            this.setAvgResults(r.QN, r.UN, r.RN, r.TN, new Matrix(0, 0), new Matrix(0, 0),
                    r.CN, r.XN, r.runtime, method, 1);
        } else if (method.startsWith("snc")) {
            // Stochastic network calculus upper bound (Fidler-Rizk 2015) for
            // feed-forward OPEN networks; the analyzer gates itself, including
            // the feed-forward and deterministic-routing premises.
            MVAResult r = Solver_ba_snc_analyzer.solver_ba_snc_analyzer(this.sn, bopts);
            this.setAvgResults(r.QN, r.UN, r.RN, r.TN, new Matrix(0, 0), new Matrix(0, 0),
                    r.CN, r.XN, r.runtime, method, 1);
        } else if (method.startsWith("bpt")) {
            // Achievable-region LP relaxation (Bertsimas-Paschalidis-Tsitsiklis
            // 1994) for multiclass OPEN networks; the analyzer gates itself.
            MVAResult r = Solver_ba_bpt_analyzer.solver_ba_bpt_analyzer(this.sn, bopts);
            this.setAvgResults(r.QN, r.UN, r.RN, r.TN, new Matrix(0, 0), new Matrix(0, 0),
                    r.CN, r.XN, r.runtime, method, 1);
        } else if (method.startsWith("qrf")) {
            SolverCTMC.AnalyzerResult r = Solver_ctmc_qrf_analyzer.solver_ctmc_qrf_analyzer(this.sn, bopts);
            this.setAvgResults(r.QN, r.UN, r.RN, r.TN, new Matrix(0, 0), new Matrix(0, 0),
                    r.CN, r.XN, r.runtime, method, 1);
            // listAllMethods, NOT listValidMethods: the latter drops what this
            // model cannot run, and dispatching on it would answer a direct
            // request for a gated method with "unknown method" instead of the
            // analyzer's own reason.
        } else if (Arrays.asList(listAllMethods()).contains(method)) {
            MVAResult ret = Solver_mva_bound_analyzer.solver_mva_bound_analyzer(this.sn, bopts);
            this.setAvgResults(ret.QN, ret.UN, ret.RN, ret.TN, new Matrix(0, 0), new Matrix(0, 0),
                    ret.CN, ret.XN, ret.runtime, method, ret.iter);
        } else {
            line_error(mfilename(new Object() {}), "Unknown bound method '" + method
                    + "'. Valid methods: " + String.join(", ", this.listValidMethods()));
        }
    }

    /**
     * MAP-AMVA LP bounds (G. Casale, E. Smirni, IEEE/IFIP DSN 2009).
     *
     * <p>The linear program over the EXACT mean-value balance equations of a
     * closed MAP queueing network, solved through {@link Mapqn_bnd_lr_mva}. It
     * is the only family in SolverBA that consumes the CORRELATION between
     * successive services rather than the service mean alone: its variables are
     * the per-phase queue lengths QN(i,k) and utilizations UN(i,k), so a
     * workload whose burstiness moves the bottleneck between stations is
     * bounded rather than averaged into a renewal process. That is why "MAP"
     * and "MMPP2" reach this family's feature set and no other.
     *
     * <p>The LP carries phases at ONE queue and requires it to be the LAST --
     * q(i,j,k,h) reads the scalar muM(i) for i &lt; M and the (D0,D1) pair
     * muMAP/v for i == M -- so a model whose phase-carrying station sits
     * elsewhere is PERMUTED rather than refused, and the results are permuted
     * back before they are returned.
     *
     * @param sn     network structure (single-class closed, single-server, no delay)
     * @param method "mapamva.upper" or "mapamva.lower"
     * @return per-station bound metrics
     */
    private static MVAResult mapamvaBounds(NetworkStruct sn, String method) {
        long tstart = System.currentTimeMillis();
        int M = sn.nstations;
        for (int i = 0; i < M; i++) {
            if (sn.sched.get(sn.stations.get(i)) == SchedStrategy.INF) {
                line_error(mfilename(new Object() {}),
                        "Method 'mapamva' does not support delay (infinite-server) stations: the "
                        + "MAP-AMVA program of Casale-Smirni (DSN 2009) is written for a network of "
                        + "queues and the paper names the delay extension as open work. Use a QRF "
                        + "method, which carries the load-dependent rate law.");
            }
        }
        boolean upper = method.endsWith(".upper");
        String sense = upper ? "max" : "min";
        int N = (int) Math.round(sn.njobs.get(0, 0));

        // Phase order per station. One phase is an exponential server, which
        // enters the LP as the scalar rate muM(i); more than one is the (D0,D1)
        // pair, which only queue M can hold.
        int[] kph = new int[M];
        Matrix[][] maps = new Matrix[M][];
        JobClass jobClass = sn.jobclasses.get(0);
        for (int i = 0; i < M; i++) {
            Station station = sn.stations.get(i);
            MatrixCell procCell = sn.proc.get(station) != null ? sn.proc.get(station).get(jobClass) : null;
            if (procCell != null && procCell.size() > 0) {
                maps[i] = new Matrix[]{procCell.get(0), procCell.get(1)};
                kph[i] = procCell.get(0).getNumRows();
            } else {
                maps[i] = null;
                kph[i] = 1;
            }
        }
        int mapIdx = -1;
        List<Integer> phased = new ArrayList<Integer>();
        for (int i = 0; i < M; i++) {
            if (kph[i] > 1) {
                phased.add(i + 1);
                mapIdx = i;
            }
        }
        if (phased.size() > 1) {
            line_error(mfilename(new Object() {}),
                    "Method 'mapamva' carries phases at ONE station: the LP gives queue M the "
                    + "(D0,D1) pair and every other queue a scalar rate. Stations " + phased
                    + " are all non-exponential. Use a QRF method, whose q carries a phase at "
                    + "every station.");
        }
        if (mapIdx < 0) {
            // Every station exponential: the program is still the right one, it
            // just degenerates to K = 1, where the per-phase variables collapse
            // and the balances become the product-form ones of Mapqn_bnd_lr_pf.
            mapIdx = M - 1;
        }
        int[] perm = new int[M];
        int pos = 0;
        for (int i = 0; i < M; i++) {
            if (i != mapIdx) {
                perm[pos++] = i;
            }
        }
        perm[M - 1] = mapIdx;
        int K = kph[mapIdx];

        Matrix muMAP = new Matrix(K, K);
        Matrix v = new Matrix(K, K);
        if (maps[mapIdx] == null) {
            muMAP.set(0, 0, 1.0);
            v.set(0, 0, 0.0);
        } else {
            Matrix d0 = maps[mapIdx][0];
            Matrix d1 = maps[mapIdx][1];
            for (int a = 0; a < K; a++) {
                for (int b = 0; b < K; b++) {
                    // muMAP(k,h) is the completion rate out of phase k landing
                    // in phase h, i.e. D1(k,h); v(k,h) is the background phase
                    // change that completes no job, i.e. D0 off the diagonal.
                    // Same (from,to) convention as Solver_ctmc_qrf_analyzer --
                    // writing either as its transpose is invisible for a
                    // reversible D0 and silently reverses the phase order of an
                    // Erlang.
                    muMAP.set(a, b, d1.get(a, b));
                    v.set(a, b, (a == b) ? 0.0 : d0.get(a, b));
                }
            }
        }

        double[] muM = new double[Math.max(M - 1, 0)];
        for (int a = 0; a < M - 1; a++) {
            muM[a] = sn.rates.get(perm[a], 0);
        }
        Matrix r = new Matrix(M, M);
        for (int a = 0; a < M; a++) {
            for (int b = 0; b < M; b++) {
                r.set(a, b, sn.rt.get(perm[a], perm[b]));
            }
        }

        Matrix visits = sn.visits.get(0);
        double[] Vp = new double[M];
        double[] Sp = new double[M];
        for (int a = 0; a < M; a++) {
            int i = perm[a];
            Vp[a] = visits.get(i, 0);
            MatrixCell procI = sn.proc.get(sn.stations.get(i)) != null
                    ? sn.proc.get(sn.stations.get(i)).get(jobClass) : null;
            Sp[a] = (procI != null && procI.size() > 0)
                    ? Map_mean.map_mean(procI.get(0), procI.get(1)) : 0.0;
        }

        MVAVersionParameters params = new MVAVersionParameters(M, N, K, muM, muMAP, r, v);

        // THREE SWEEPS OF THE SAME LP, and the utilization one runs in BOTH
        // senses on purpose. R_i = Q_i/(V_i*X) rises with Q_i and FALLS with X,
        // so an upper bound on the response time pairs Q_i^max with X^min;
        // dividing by X^max on both sides is what would report an upper R below
        // the exact value and break the bracket.
        double[] umax = new double[M];
        double[] umin = new double[M];
        double[] qbnd = new double[M];
        for (int a = 0; a < M; a++) {
            umax[a] = Mapqn_bnd_lr_mva.solve(params, a + 1, 0, "max", "UN").getObjectiveValue();
            umin[a] = Mapqn_bnd_lr_mva.solve(params, a + 1, 0, "min", "UN").getObjectiveValue();
            qbnd[a] = Mapqn_bnd_lr_mva.solve(params, a + 1, 0, sense, "QN").getObjectiveValue();
        }

        // Utilization law U_i = X*V_i*S_i, exact at a single server under ANY
        // service law, so each station turns its own utilization bound into a
        // throughput bound and the tightest of the M survives. A station with no
        // visits or no service time carries no information and is skipped.
        double xUp = Double.POSITIVE_INFINITY;
        double xLo = 0.0;
        boolean any = false;
        for (int a = 0; a < M; a++) {
            double load = Vp[a] * Sp[a];
            if (!Double.isFinite(load) || load <= 0) {
                continue;
            }
            any = true;
            xUp = Math.min(xUp, umax[a] / load);
            xLo = Math.max(xLo, umin[a] / load);
        }
        if (!any) {
            line_error(mfilename(new Object() {}), "Method '" + method
                    + "' found no station with both a positive visit ratio and a positive mean "
                    + "service time.");
        }
        double xb = upper ? xUp : xLo;
        double xopp = upper ? xLo : xUp;

        // Unpermute: the LP orders the stations with the phase-carrying one last.
        Matrix QN = new Matrix(M, 1);
        Matrix UN = new Matrix(M, 1);
        Matrix RN = new Matrix(M, 1);
        Matrix TN = new Matrix(M, 1);
        for (int a = 0; a < M; a++) {
            int i = perm[a];
            UN.set(i, 0, upper ? umax[a] : umin[a]);
            QN.set(i, 0, qbnd[a]);
            TN.set(i, 0, Vp[a] * xb);
            if (xopp > 0 && Vp[a] > 0) {
                RN.set(i, 0, qbnd[a] / (Vp[a] * xopp));
            } else {
                RN.set(i, 0, Double.POSITIVE_INFINITY);
            }
        }
        Matrix CN = new Matrix(1, 1);
        // Delay stations are refused above, so the closed-network response time
        // is N/X exactly and the throughput bracket transfers to it directly.
        // Summing the per-station R bounds instead would add M separately
        // attained maxima and report a looser number.
        CN.set(0, 0, xopp > 0 ? N / xopp : Double.POSITIVE_INFINITY);
        Matrix XN = new Matrix(1, 1);
        XN.set(0, 0, xb);

        MVAResult ret = new MVAResult();
        ret.QN = QN;
        ret.UN = UN;
        ret.RN = RN;
        ret.TN = TN;
        ret.CN = CN;
        ret.XN = XN;
        ret.runtime = (System.currentTimeMillis() - tstart) / 1000.0;
        return ret;
    }

    /**
     * LP-based Linear Reduction bound. Solves the product-form linear-reduction
     * LP once per station (Mapqn_bnd_lr_pf, Apache SimplexSolver), minimizing
     * for "lr.lower" and maximizing for "lr.upper". Both senses are valid
     * bounds because the LP relaxation contains the exact solution. Unlike
     * qrf.mmi.linear -- whose objective is the nonlinear MMI mutual
     * information -- this method is a pure LP end to end.
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
        // A blocking-blind family on a blocked model is dropped from
        // listValidMethods, so both sides below would silently come back null.
        // Refuse by name instead, with the reason runAnalyzer gives.
        if (ignoresBlocking(resolveMethod(this.options.method))
                && SnHasBlocking.snHasBlocking(this.getStruct())) {
            // runAnalyzer routes 'default'/'auto' to 'qrf.bas' on a blocked
            // model, and getBounds has to say so rather than contradict it --
            // but it still cannot BRACKET, because the analyzer solves qrf.bas
            // in the 'max' direction alone.
            String req = this.options.method;
            if ("default".equals(req) || "auto".equals(req) || "auto.upper".equals(req)) {
                String[] routed = blockingDefault(this.getStruct());
                if (!routed[0].isEmpty()) {
                    line_error(mfilename(new Object() {}), "'" + req + "' resolves to '"
                            + routed[0] + "' on this model, which has a binding finite buffer,"
                            + " and that bound is UPPER-only: there is no bracket to return."
                            + " Call getAvgTable/getAvg for the upper bound, or SolverMVA with"
                            + " method 'sqd' for a point estimate.");
                }
            }
            line_error(mfilename(new Object() {}), "Family '" + fam + "' does not support"
                    + " finite-buffer blocking; see the SolverBA method gate. Use"
                    + " 'qrf.bas'/'qrf.rsrd', or SolverMVA with method 'sqd' for a point"
                    + " estimate.");
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
     * Returns the per-station response-time QUANTILE at violation probability
     * eps: the smallest d for which {@code P{D_ir > d} <= eps} is certified by
     * the stochastic network calculus bound.
     *
     * <p>This is the native output of the 'snc' family, so the accessor runs
     * the SNC envelope propagation directly whatever {@code options.method}
     * says; a station-class pair carrying no traffic stays NaN. Every other
     * family bounds means only and has no counterpart.</p>
     *
     * @param eps violation probability, 0 &lt; eps &lt; 1
     * @return an (nstations x nclasses) matrix of quantiles, NaN where no traffic
     */
    public Matrix getDelayPerc(double eps) {
        return sncPerc(eps)[0];
    }

    /**
     * Returns the per-station queue-length QUANTILE, in jobs, at violation
     * probability eps. The counterpart of {@link #getDelayPerc}; see it for the
     * conventions.
     *
     * @param eps violation probability, 0 &lt; eps &lt; 1
     * @return an (nstations x nclasses) matrix of quantiles, NaN where no traffic
     */
    public Matrix getBacklogPerc(double eps) {
        return sncPerc(eps)[1];
    }

    /**
     * Returns the response-time and queue-length quantiles at violation
     * probability eps, in the layout of getAvgTable. Rows are the
     * station-class pairs that carry traffic.
     *
     * @param eps violation probability, 0 &lt; eps &lt; 1
     * @return the quantile table
     */
    public NetworkPercTable getPercTable(double eps) {
        Matrix[] db = sncPerc(eps);
        NetworkStruct s = getStruct();
        List<Double> dcol = new ArrayList<Double>();
        List<Double> bcol = new ArrayList<Double>();
        List<String> className = new ArrayList<String>();
        List<String> stationName = new ArrayList<String>();
        for (int i = 0; i < s.nstations; i++) {
            for (int k = 0; k < s.nclasses; k++) {
                double d = db[0].get(i, k);
                double b = db[1].get(i, k);
                if (Double.isNaN(d) && Double.isNaN(b)) {
                    continue;
                }
                dcol.add(d);
                bcol.add(b);
                className.add(this.model.getClasses().get(k).getName());
                stationName.add(this.model.getStations().get(i).getName());
            }
        }
        NetworkPercTable t = new NetworkPercTable(dcol, bcol);
        t.setOptions(this.options);
        t.setClassNames(className);
        t.setStationNames(stationName);
        return t;
    }

    /**
     * Both quantile matrices from one envelope propagation.
     *
     * <p>The analyzer returns the (arrival, service) envelopes it built, so the
     * two quantiles cost one pass over the network and one Chernoff search per
     * pair and metric.</p>
     *
     * @param eps violation probability
     * @return {delay quantiles, backlog quantiles}
     */
    private Matrix[] sncPerc(double eps) {
        NetworkStruct s = getStruct();
        Solver_ba_snc_analyzer.Envelopes env = Solver_ba_snc_analyzer.envelopes(s);
        Matrix D = new Matrix(s.nstations, s.nclasses);
        Matrix B = new Matrix(s.nstations, s.nclasses);
        for (int i = 0; i < s.nstations; i++) {
            for (int k = 0; k < s.nclasses; k++) {
                long key = env.key(i, k);
                if (env.arv.containsKey(key)) {
                    D.set(i, k, Snc_perc_delay.snc_perc_delay(
                            env.arv.get(key), env.srv.get(key), eps).value);
                    B.set(i, k, Snc_perc_backlog.snc_perc_backlog(
                            env.arv.get(key), env.srv.get(key), eps).value);
                } else {
                    D.set(i, k, Double.NaN);
                    B.set(i, k, Double.NaN);
                }
            }
        }
        return new Matrix[]{D, B};
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

    /**
     * Bound methods this MODEL can run: a narrowing of {@link #listAllMethods}.
     *
     * <p>A name that would always be refused on this model is not offered, so a
     * caller enumerating the list never asks for one. Ask for it by name anyway
     * and {@code runAnalyzer} still dispatches, so the analyzer's own reason is
     * what comes back.</p>
     */
    public String[] listValidMethods() {
        List<String> m = new ArrayList<String>(Arrays.asList(listAllMethods()));
        // The STRUCTURAL premises -- single-class closed, fully closed,
        // single-server -- come from methodRefusal, the same predicate
        // runAnalyzer raises on and supportsModelMethod reports. Asked first so
        // every later narrowing works on names this model could actually run:
        // before it, a two-class closed network was offered all 36
        // demand-parameterized bounds and 30 of them raised on contact.
        NetworkStruct snb = this.getStruct();
        List<String> structural = new ArrayList<String>();
        for (int i = 0; i < m.size(); i++) {
            // ... and methodDegenerate withholds the second kind of name: one
            // whose premises this model MEETS but whose formula says nothing
            // here. 'ldbcmp.lower' at N == Qhat is the only such case.
            if (methodRefusal(snb, m.get(i)).isEmpty()
                    && methodDegenerate(snb, m.get(i)).isEmpty()) {
                structural.add(m.get(i));
            }
        }
        m = structural;
        // 'bpt', 'bgt' and 'snc' are the mirror image of the reduction bounds:
        // all three are derived for an OPEN network of single-server
        // exponential stations, so every closed model, every delay station and
        // every multiserver station rules them out. Every other family rules
        // OUT the open model, so on an open network the list narrows to those
        // three.
        // 'spnlp.*' is the only family indexed by a MARKING rather than by
        // demands and a population, and the split is total in both directions:
        // on a Petri net nothing else has a representation of the model, and off
        // one spnlp has nothing to read. Two of the gates below already
        // half-cover this by accident -- a Place is an INF station, so
        // isReducible and isBptFeasible are both false on any Petri net -- but
        // the demand-parameterized families survive them and must be dropped by
        // name. Applied first so the open/closed narrowing cannot reinstate one.
        if (isPetri(snb)) {
            List<String> onlyPetri = new ArrayList<String>();
            for (int i = 0; i < m.size(); i++) {
                if (m.get(i).startsWith("spnlp")) {
                    onlyPetri.add(m.get(i));
                }
            }
            return onlyPetri.toArray(new String[0]);
        }
        List<String> noPetri = new ArrayList<String>();
        for (int i = 0; i < m.size(); i++) {
            if (!m.get(i).startsWith("spnlp")) {
                noPetri.add(m.get(i));
            }
        }
        m = noPetri;
        if (!isBptFeasible(snb)) {
            m.remove("bpt.lower");
            m.remove("bgt.upper");
            m.remove("snc.upper");
        }
        if (isFullyOpen(snb)) {
            List<String> onlyOpen = new ArrayList<String>();
            if (m.contains("bpt.lower")) {
                onlyOpen.add("bpt.lower");
            }
            if (m.contains("bgt.upper")) {
                onlyOpen.add("bgt.upper");
            }
            if (m.contains("snc.upper")) {
                onlyOpen.add("snc.upper");
            }
            m = onlyOpen;
        } else if (!isReducible(snb)) {
            // The QR/LR/QRF reduction bounds share one premise: a single-class
            // closed network of single-server stations. Naming them on a model
            // they cannot run turns a rejection into a method a caller is
            // invited to ask for, the same reason 'sqni' is gated in SolverMVA.
            // Mirrors SolverBA.m and the python/C++ twins.
            //
            // The LOAD-DEPENDENT arms survive where the rest cannot run:
            // alpha(i,n) is the rate law of a delay (alpha = n), of a c-server
            // station (alpha = min(n,c)) and of limited load dependence alike,
            // so 'qrf.mmi.ld' and 'qrf.mmi.linear' answer those models on the
            // model's own chain. SnToQrfAlpha owns the one restriction that
            // survives, exponential service wherever a station serves several
            // jobs at once. Dropping them with the rest would hide from a
            // caller the only two bound methods such a model has.
            boolean ldReducible = isLdReducible(snb);
            List<String> keep = new ArrayList<String>();
            for (int i = 0; i < m.size(); i++) {
                String name = m.get(i);
                boolean ldArm = name.equals("qrf.mmi.ld") || name.equals("qrf.mmi.linear");
                if (ldReducible && ldArm) {
                    keep.add(name);
                    continue;
                }
                // 'mapamva' shares that premise exactly -- single-class
                // closed, single-server, no delay -- so it narrows with them
                // rather than being offered on a model it would refuse on
                // contact. It is NOT one of the load-dependent arms: its q
                // carries no population index, so a delay or a c-server station
                // has nowhere to go.
                if (name.equals("qr") || name.equals("lr")
                        || name.startsWith("lr.") || name.startsWith("qrf.")
                        || name.startsWith("mapamva")) {
                    continue;
                }
                keep.add(name);
            }
            m = keep;
        }
        // A binding finite buffer rules out everything but the QRF blocking
        // bounds: the other families presume unbounded buffers, and runAnalyzer
        // refuses them by name on such a model. The list can legitimately come
        // back EMPTY -- a blocked model that is not single-class closed
        // single-server has no bound method at all, and offering one would be
        // the mis-selection this gate exists to prevent.
        if (SnHasBlocking.snHasBlocking(snb)) {
            List<String> keepBlk = new ArrayList<String>();
            for (int i = 0; i < m.size(); i++) {
                if (!ignoresBlocking(resolveMethod(m.get(i)))) {
                    keepBlk.add(m.get(i));
                }
            }
            // 'default' is offered back when it now MEANS one of the survivors:
            // runAnalyzer routes it to 'qrf.bas' on a blocked model of the right
            // shape, so a caller enumerating the list would otherwise be told
            // the model's own default is invalid.
            if (!blockingDefault(snb)[0].isEmpty()) {
                keepBlk.add(0, "default");
            }
            m = keepBlk;
        }
        return m.toArray(new String[0]);
    }

    /**
     * The families whose bound is a function of the single-chain demand vector
     * D = V./rates, the think time Z and the population N. A multiclass or open
     * model simply does not have those, which is why the rule is total.
     */
    private static final String[] SINGLE_CLASS_FAMS = {
            "auto", "aba", "bjb", "pb", "sb", "gb", "harel", "lr",
            "pbh", "cbh", "pbk", "bjbk", "ssd", "sib", "scb", "ldbcmp",
            // 'mapamva' is single-class for a different reason from the rest --
            // its LP variables QN(i,k)/UN(i,k) are indexed by station and MAP
            // phase, with no class index at all -- but the premise it fails on
            // is the same one.
            "mapamva"};

    /**
     * The multiclass families: they take a per-chain demand MATRIX and a
     * population VECTOR, so several classes are fine and an infinite population
     * is not.
     */
    private static final String[] FULLY_CLOSED_FAMS = {"mwba", "cub", "mbjb", "looping"};

    /**
     * The single-server families whose alternative on a multiserver model IS
     * 'ssd', which is why the reason names it. 'ssd' is the multiserver bound
     * itself, 'ldbcmp' is parameterized by the limiting demand of a
     * load-dependent station and 'auto' composes whichever candidates survive,
     * so all three are absent.
     */
    private static final String[] SINGLE_SERVER_FAMS = {
            "aba", "bjb", "pb", "sb", "gb", "harel", "lr",
            "pbh", "cbh", "pbk", "bjbk", "sib", "scb", "mapamva"};

    private static boolean isIn(String[] set, String fam) {
        for (int i = 0; i < set.length; i++) {
            if (set[i].equals(fam)) {
                return true;
            }
        }
        return false;
    }

    /** The family prefix of a method name: everything before the first dot. */
    private static String familyOf(String method) {
        int dot = method.indexOf('.');
        return (dot < 0) ? method : method.substring(0, dot);
    }

    /**
     * The STRUCTURAL premises of the SolverBA bound families, in one place: the
     * reason METHOD cannot bound the model SN, or "" when it can.
     *
     * <p>ONE PREDICATE, TWO CALLERS. {@code runAnalyzer} asks it before
     * dispatching and raises what it returns; {@code supportsModelMethod} asks
     * it after the feature gate and reports the same sentence, which is what
     * findSolver, listValidMethods and SolverAUTO's ranked choice all read. A
     * second copy of any rule below is how the report and the run drift apart:
     * the report offers a pair that raises the moment it is run, which is the
     * defect this method exists to remove.
     *
     * <p>WHAT BELONGS HERE AND WHAT DOES NOT. Only the rules the feature
     * registry cannot name. {@link FeatureSet} has no entry for "one class",
     * for a server count or for a station count, so those are structural and
     * live here. Rules of the form "this family does not accept a delay
     * station" ARE nameable and belong in {@link #getMethodFeatureSet}, which
     * drops SchedStrategy_INF from the offending method's set instead: a
     * feature set can refuse a model for HAVING a construct, never for lacking
     * one.
     *
     * <p>METHOD is taken as the caller spells it and resolved through
     * {@code resolveMethod}, so 'default' is judged as the gb.upper it runs as
     * and the reason names that. The marking-parameterized spnlp and the QRF
     * reduction bounds carry no rule here: the QRF premise is the reducibility
     * test listValidMethods already applies. Of the three OPEN families, 'bpt'
     * and 'bgt' carry none either -- a closed model is refused by their feature
     * set and their analyzers walk the routing matrix for the rest -- while
     * 'snc' carries one, the SERVICE law.
     *
     * <p>WHY THE SNC SERVICE LAW IS HERE AND THE bpt/bgt ONE IS NOT. All three
     * analyzers refuse a non-exponential law at a queueing station. For bpt and
     * bgt that rule extends to the SOURCE and is registry-expressible, so it
     * rides in {@link #getMethodFeatureSet} as a dropped law: both are
     * invariant to the arrival law beyond its mean, so a non-exponential source
     * is not something they refuse, it is something they silently bound as if
     * it were Poisson. snc is the opposite: it CONSUMES the arrival law and its
     * analyzer branches on a non-exponential source deliberately. Its rule is
     * about the SERVICE only, and no feature name can say "Erlang at a Queue
     * but not at a Source", so it is structural.
     *
     * <p>Mirrors matlab/src/solvers/BA/ba_method_refusal.m and its native python
     * and C++ twins.
     *
     * @param sn     the network struct to judge
     * @param method the method name as the caller spells it
     * @return the reason the pair cannot run, or "" when it can
     */
    public static String methodRefusal(NetworkStruct sn, String method) {
        if (sn == null) {
            return "";
        }
        String resolved = resolveMethod(method);
        String fam = familyOf(resolved);
        boolean anyOpen = false;
        for (int r = 0; r < sn.njobs.getNumCols(); r++) {
            if (Double.isInfinite(sn.njobs.get(0, r))) {
                anyOpen = true;
                break;
            }
        }
        if (isIn(SINGLE_CLASS_FAMS, fam)) {
            if (sn.nclasses != 1 || sn.nclosedjobs <= 0) {
                return "Method '" + resolved + "' supports single-class closed networks only.";
            }
        } else if (isIn(FULLY_CLOSED_FAMS, fam)) {
            if (sn.nclosedjobs <= 0 || anyOpen) {
                return "Method '" + resolved + "' supports fully closed networks only.";
            }
        }
        boolean multiserver = false;
        for (int i = 0; i < sn.nstations; i++) {
            if (sn.sched.get(sn.stations.get(i)) == SchedStrategy.INF) {
                continue;
            }
            if (sn.nservers.get(i, 0) > 1) {
                multiserver = true;
                break;
            }
        }
        if (multiserver) {
            if (isIn(SINGLE_SERVER_FAMS, fam)) {
                return "Method '" + resolved
                        + "' does not support multi-server stations (use 'ssd').";
            }
            if (isIn(FULLY_CLOSED_FAMS, fam)) {
                return "Method '" + resolved + "' does not support multi-server stations.";
            }
        }

        // The SNC service law. Judged over the pairs a station COULD serve
        // rather than over the ones that carry traffic: the analyzer restricts
        // to the latter, which needs the traffic equations solved, and this is
        // their conservative outer approximation -- it never admits a pair the
        // analyzer refuses, and can only differ on a pair given a service time
        // at a station its class never visits. A Source is skipped, which is
        // the whole reason this is not a feature-set delta.
        if ("snc".equals(fam)) {
            for (int i = 0; i < sn.nstations; i++) {
                if (sn.nodetype.get((int) sn.stationToNode.get(i)) == NodeType.Source) {
                    continue;
                }
                for (int r = 0; r < sn.nclasses; r++) {
                    double mu = sn.rates.get(i, r);
                    if (Double.isNaN(mu) || Double.isInfinite(mu) || mu <= 0) {
                        continue;
                    }
                    ProcessType pid = sn.procid.get(sn.stations.get(i)).get(sn.jobclasses.get(r));
                    if (pid != ProcessType.EXP) {
                        return "Method '" + resolved + "' requires exponential service: station "
                                + (i + 1) + " class " + (r + 1) + " is " + pid + ".";
                    }
                }
            }
        }
        return "";
    }

    /**
     * Whether METHOD APPLIES to SN but its bound carries no information there,
     * and why. Empty when the bound is informative, and empty for every method
     * that has no such regime.
     *
     * <p>THIS IS A DIFFERENT QUESTION FROM {@link #methodRefusal}, which is why
     * it is a different method. That one answers "is this model outside the
     * method's domain", and its answer is what runAnalyzer raises. This one
     * answers "inside the domain, does the formula still say anything", and its
     * answer is NOT raised: a degenerate bound is a VALID bound, just a vacuous
     * one, so an analyzer asked for it by name is entitled to publish it. What
     * must not happen is OFFERING it: findSolver and listValidMethods exist to
     * name the pairs a caller can act on, and a table of zeros over a network
     * with jobs circulating in it is not something anyone can act on.
     *
     * <p>THE ONE METHOD WITH SUCH A REGIME IS 'ldbcmp.lower'. The
     * Anselmi-Cremonesi bound is built from the population SURPLUS a = N -
     * Qhat, where Qhat is the occupancy the non-bottleneck stations and the
     * think time would hold in the open network fed at the bottleneck's
     * saturation rate. Pfqn_ldbcmp returns NaN below the regime (a &lt; 0),
     * which the analyzer already refuses by name; AT the boundary a = 0 it
     * returns Xlo = 0, which is formally the trivial bound X &gt;= 0 and
     * propagates into a table whose queue lengths, utilizations and throughputs
     * are all zero. Every entry of that table is a true lower bound and none of
     * them is usable, and a caller cannot tell it from a real answer of zero.
     *
     * @param sn     the network struct to judge
     * @param method the method name as the caller spells it
     * @return the reason the bound says nothing here, or "" when it does
     */
    public static String methodDegenerate(NetworkStruct sn, String method) {
        if (sn == null || !"ldbcmp.lower".equals(resolveMethod(method))) {
            return "";
        }
        // The applicability rules come first and are not restated: a model this
        // method is outside the domain of has no bound to be degenerate about.
        if (!methodRefusal(sn, method).isEmpty()) {
            return "";
        }
        // TWO INDEX SPACES MEET HERE, and conflating them is what made this
        // predicate throw on a Petri net. sn.visits is STATEFUL-indexed --
        // SnRefreshVisits builds it (nstateful x nclasses) -- while sn.sched,
        // sn.rates and sn.stations are STATION-indexed. Every Station is a
        // StatefulNode, but not every StatefulNode is a Station: a Transition
        // extends ServiceNode and so is stateful WITHOUT being a station, and a
        // Place extends Station and is both. On the fork-join SPN of
        // SpnLpbndTest that is 4 stations against 7 stateful nodes, so walking
        // the VISIT rows and indexing sn.stations with the row threw
        // IndexOutOfBounds at row 4. Walk the STATIONS and convert with
        // stationToStateful, which is what SnRefreshVisits itself does.
        Matrix V = sn.visits.get(0);
        double Zt = 0.0;
        int nq = 0;
        for (int i = 0; i < sn.nstations; i++) {
            double v = V.get((int) sn.stationToStateful.get(i), 0);
            if (sn.sched.get(sn.stations.get(i)) == SchedStrategy.INF) {
                Zt += v / sn.rates.get(i, 0);
            } else {
                nq++;
            }
        }
        Matrix D = new Matrix(nq, 1);
        int idx = 0;
        for (int i = 0; i < sn.nstations; i++) {
            if (sn.sched.get(sn.stations.get(i)) != SchedStrategy.INF) {
                D.set(idx++, 0,
                        V.get((int) sn.stationToStateful.get(i), 0) / sn.rates.get(i, 0));
            }
        }
        // NO QUEUEING STATION, NO BOUND. Every station is a delay (or a Place,
        // on a Petri net, which is an INF station too), so there is no
        // bottleneck to build Qhat on and Pfqn_ldbcmp has no demand vector to
        // read. A PREDICATE MUST NOT THROW -- this one is asked once per name
        // by listValidMethods, before the Petri sieve has had a chance to drop
        // anything -- so the case is answered here. It is reached only once the
        // loops above are index-safe, which is the other half of the same fix.
        if (nq == 0) {
            return "Method 'ldbcmp.lower' has no queueing station to bound here: every station"
                    + " is an infinite server, so the bottleneck the open-network occupancy is"
                    + " built on does not exist.";
        }
        double N = (double) sn.nclosedjobs;
        double[] xb = Pfqn_ldbcmp.pfqn_ldbcmp(D, N, Zt, new Matrix(nq, 1));
        if (Double.isNaN(xb[2]) || Double.isInfinite(xb[2])) {
            return "Method 'ldbcmp.lower' does not apply here: the open-network occupancy Qhat"
                    + " the bound is built from does not exist, because a non-bottleneck station"
                    + " saturates at the bottleneck's arrival rate.";
        }
        if (Double.isNaN(xb[0]) || Double.isInfinite(xb[0]) || xb[0] <= 0) {
            return String.format("Method 'ldbcmp.lower' needs a population strictly above the"
                    + " open-network occupancy the bound is built from (Qhat=%.4f, N=%d): with no"
                    + " surplus it degenerates to the trivial bound X >= 0 and reports a table of"
                    + " zeros.", xb[2], (int) N);
        }
        return "";
    }

    /**
     * The base envelope with the per-method deltas the registry CAN name.
     *
     * <p>A feature set says "I accept this construct", so it can refuse a model
     * for HAVING one and never for lacking one; that is exactly the shape of
     * the delay-station and closed-class premises below, and exactly not the
     * shape of "one class" or "one server", which have no feature name and live
     * in {@link #methodRefusal} instead. Judged on the RESOLVED name so that
     * 'default' carries the envelope of the gb.upper it runs as.
     *
     * <p>DELAY STATIONS. 'sb' and 'lr' reject an infinite-server station
     * outright, and 'harel', 'sib' and 'scb' reject a nonzero think time, which
     * on these models is the same station: harel extrapolates the exact
     * normalizing constant of a delay-free network, SIB Section 3.2 is the
     * extension that would carry Z and is not implemented, and SCB Theorem 3
     * rests on the delay-free balanced-network throughput. The three OPEN
     * families reject one too, each being derived for one server per station.
     *
     * <p>CLASS TYPES. The three OPEN families drop ClosedClass, which is the
     * whole of their class premise. The MIRROR delta -- dropping OpenClass from
     * every demand-parameterized family -- is deliberately NOT applied:
     * "supports single-class closed networks only" is one rule, its
     * single-class half has no feature name, and splitting it across the two
     * mechanisms would report the closed half here and the single-class half in
     * methodRefusal for the same model. It is stated once, structurally.
     * 'spnlp' takes no delta at all: it is indexed by the marking, and whether
     * that marking is bounded is a question about the P-invariants of the net,
     * which Spn_lpbnd answers.
     *
     * @param method the concrete method name
     * @return the per-method feature envelope
     */
    @Override
    public FeatureSet getMethodFeatureSet(String method) {
        FeatureSet featSupported = SolverBA.getFeatureSet();
        // Model-aware: on a blocked model "default"/"auto"/"auto.upper" carry the
        // envelope of the 'qrf.bas' they are routed to, which is the same
        // resolution runAnalyzer dispatches on. Mirrors MATLAB
        // ba_resolve_model_method.
        String resolved = resolveModelMethod(method);
        String fam = familyOf(resolved);
        if ("bpt".equals(fam) || "bgt".equals(fam) || "snc".equals(fam)) {
            featSupported.setFalse(new String[]{
                    "ClosedClass", "Delay", "DelayStation", "SchedStrategy_INF"});
            if (!"snc".equals(fam)) {
                // SERVICE AND ARRIVAL LAWS. 'bpt' and 'bgt' are derived for a
                // MARKOVIAN open network and read the mean alone, so a
                // non-exponential law anywhere is not something they refuse at
                // run time -- it is something they silently bound as if it were
                // Poisson. Measured on the M/M/1 shape: replacing the Exp(1)
                // source by an Erlang of the same mean leaves bgt.upper at QLen
                // 32.6667 and bpt.lower at 1.0, digit for digit. That is a
                // bound on a DIFFERENT system, so the laws are dropped here
                // rather than left to a run-time check the analyzers do not
                // make: their own procid test covers the queueing stations only
                // and would miss exactly the source case.
                //
                // 'snc' is excluded: it CONSUMES the arrival law (the same
                // substitution moves it from 3.8244 to 3.0092) and its analyzer
                // branches on a non-exponential source deliberately. Its rule
                // is about the SERVICE only, which no feature name can say, so
                // it lives in methodRefusal instead.
                featSupported.setFalse(new String[]{
                        "APH", "Coxian", "Cox2", "Erlang", "HyperExp", "PH",
                        "Det", "Lognormal", "Pareto", "Uniform", "Weibull"});
            }
        } else if ("sb".equals(fam) || "harel".equals(fam) || "sib".equals(fam)
                || "scb".equals(fam) || "lr".equals(fam)) {
            featSupported.setFalse(new String[]{
                    "Delay", "DelayStation", "SchedStrategy_INF"});
        }
        // MODULATED SERVICE, THE ONE DELTA THAT RUNS THE OTHER WAY. "MAP" and
        // "MMPP2" sit in the base envelope so that 'mapamva' can accept them,
        // which means every OTHER family has to give them back: each of the rest
        // reads the service MEAN alone and would bound a correlated model as
        // though its services were independent, quietly returning a bracket for
        // a different system. Written as a delta on the complement rather than
        // as a grant because a feature set can refuse a model for HAVING a
        // construct and never for lacking one.
        //
        // 'mapamva' takes the delay delta instead: its LP is a network of
        // queues, and Casale-Smirni name the delay extension as open work.
        if ("mapamva".equals(fam)) {
            featSupported.setFalse(new String[]{
                    "Delay", "DelayStation", "SchedStrategy_INF"});
        } else {
            featSupported.setFalse(new String[]{"MAP", "MMPP2"});
        }
        // MULTISERVER (registry name since 2026-09-05) is out of the base
        // envelope: every demand-parameterized family reads one server per
        // station (methodRefusal names 'ssd' as the alternative), the alpha-free
        // QRF arms refuse a c-server station through SnToQrfAlpha and the open
        // families through their own refusal. What carries the count is granted
        // here: 'ssd' (the multiserver bound), 'ldbcmp' (its fixed-rate form runs
        // on the c-server rate law), 'auto' (which picks among them) and the two
        // load-dependent QRF arms, whose alpha(i,n) IS min(n,c).
        if ("ssd".equals(fam) || "ldbcmp".equals(fam) || "auto".equals(fam)
                || "qrf.mmi.ld".equals(resolved) || "qrf.mmi.linear".equals(resolved)) {
            featSupported.setTrue(new String[]{"MultiServer"});
        }
        // FINITECAPACITY (registry name since 2026-09-05): only the QRF blocking
        // bounds carry the buffer (the MM, MM1, ZZ, ZM, BB, F tables), which is
        // the same split ignoresBlocking makes; the structural refusal keeps
        // naming them. 'default'/'auto'/'auto.upper' resolve to 'qrf.bas' on a
        // blocked model of the right shape, so they are granted it through
        // RESOLVED. 'spnlp' is NOT granted: its polytope reads no Place capacity
        // (Spn_lpbnd reads sn.nodeparam only), so a capped place would be
        // relaxed away.
        if (resolved.startsWith("qrf.bas") || resolved.startsWith("qrf.rsrd")) {
            featSupported.setTrue(new String[]{"FiniteCapacity"});
        }
        return featSupported;
    }

    /**
     * The concrete bound method NAME runs as on THIS model: the model-free
     * aliases of {@link #resolveMethod}, then the finite-buffer routing of
     * {@link #blockingDefault}, under which "default", "auto" and "auto.upper"
     * run as 'qrf.bas' on a model whose buffer BINDS and whose shape that bound
     * needs. Mirrors MATLAB {@code ba_resolve_model_method}.
     *
     * <p>Only the UPPER side routes: the analyzer solves 'qrf.bas' in the 'max'
     * direction alone, so 'auto.lower' keeps refusing rather than being answered
     * with the wrong side.
     *
     * @param method the method name as given
     * @return the method the name runs as on this model
     */
    private String resolveModelMethod(String method) {
        boolean routable = "default".equals(method) || "auto".equals(method)
                || "auto.upper".equals(method);
        String resolved = resolveMethod(method);
        if (routable && ignoresBlocking(resolved) && this.model instanceof Network) {
            NetworkStruct sn = this.getStruct();
            if (SnHasBlocking.snHasBlocking(sn)) {
                String alt = blockingDefault(sn)[0];
                if (!alt.isEmpty()) {
                    resolved = alt;
                }
            }
        }
        return resolved;
    }

    /**
     * The feature gate above, plus the structural premises no feature name can
     * express.
     *
     * <p>{@link #methodRefusal} is the same predicate {@code runAnalyzer}
     * raises on, so a caller gets one answer whichever of the two it meets
     * first -- which is the point: this gate is what findSolver,
     * listValidMethods and SolverAUTO's ranked choice read, and while it was
     * silent about them a two-class closed network was reported as able to run
     * 36 bound methods of which 30 raised on contact.
     *
     * @param method the concrete method name
     * @return empty string if supported, else the offending reason
     */
    @Override
    public String supportsModelMethod(String method) {
        // STRUCTURAL PREDICATE FIRST, as SolverBA.m and the python twin do, and
        // for the reason this audit exists: runAnalyzer raises methodRefusal's
        // sentence, so the report has to lead with the same one or the gate and
        // the run disagree about WHY. Once MultiServer became a registry name in
        // 2026-09-05, the feature envelope started refusing the single-server
        // families first and the report answered "(feature: MultiServer)" where
        // the run said "Method 'gb.lower' does not support multi-server stations
        // (use 'ssd')", which is the more useful sentence and names the way out.
        if (this.model instanceof Network) {
            String structuralFirst = methodRefusal(this.getStruct(), method);
            if (!structuralFirst.isEmpty()) {
                return structuralFirst;
            }
        }
        String reason = super.supportsModelMethod(method);
        if (!reason.isEmpty()) {
            return reason;
        }
        if (!(this.model instanceof Network)) {
            return "";
        }
        // A bound that APPLIES but says nothing is not offered either; see
        // methodDegenerate for why that is a separate question and why the
        // analyzer is still allowed to answer it when asked by name.
        return methodDegenerate(this.getStruct(), method);
    }

    /**
     * The 'default'/'auto'/'qr'/'lr' aliases, as {@code runAnalyzer} resolves
     * them. Bare 'lr' means 'lr.upper' and 'qr' is the QRF quadratic reduction.
     *
     * @param method method name as given
     * @return the resolved name
     */
    private static String resolveMethod(String method) {
        if ("default".equals(method)) {
            return "gb.upper";
        } else if ("auto".equals(method)) {
            return "auto.upper";
        } else if ("lr".equals(method)) {
            return "lr.upper";
        } else if ("qr".equals(method)) {
            return "qrf.mmi";
        }
        return method;
    }

    /**
     * Whether METHOD bounds a model as if its buffers were unbounded.
     *
     * <p>Every family but the QRF BLOCKING bounds is parameterized by demands
     * (visits x service time) and a population alone, which is the BCMP
     * parameterization: unbounded buffers, and an equilibrium distribution that
     * factorizes. A finite buffer that BINDS breaks both premises, so the
     * numbers do not bracket the blocked model -- on cqn_bas_blocking (Queue2
     * capped at 1, N = 2) gb.upper reports QLen 1.28 at a station that can
     * never hold more than one job. 'qrf.bas*' and 'qrf.rsrd' carry the
     * blocking tables (MM, MM1, ZZ, ZM, BB, F) explicitly and are the
     * exceptions.</p>
     *
     * @param method the RESOLVED method name
     * @return true if the method assumes unbounded buffers
     */
    private static boolean ignoresBlocking(String method) {
        // 'spnlp.*' is an exception alongside the QRF blocking bounds: its
        // polytope is indexed by the marking itself, so a bounded place enters
        // it as a variable upper bound and as the P-invariant equality that
        // produced the bound. The buffer is modelled, not assumed away.
        return !(method.startsWith("qrf.bas") || method.startsWith("qrf.rsrd")
                || method.startsWith("spnlp"));
    }

    /** Whether the model holds Transition nodes, i.e. is a Petri net. */
    private static boolean isPetri(NetworkStruct sn) {
        for (int i = 0; i < sn.nodetype.size(); i++) {
            if (sn.nodetype.get(i) == NodeType.Transition) {
                return true;
            }
        }
        return false;
    }

    /** Whether every class of the model is open. */
    private static boolean isFullyOpen(NetworkStruct sn) {
        for (int r = 0; r < sn.njobs.length(); r++) {
            if (Double.isFinite(sn.njobs.get(r))) {
                return false;
            }
        }
        return true;
    }

    /**
     * Whether 'bpt' and 'bgt' apply: fully open, no INF station, one server
     * each. 'bgt.upper' additionally needs deterministic non-merging routes,
     * which its analyzer checks by walking the routing matrix -- too expensive
     * to repeat here, so it stays listed and refuses by name.
     */
    private static boolean isBptFeasible(NetworkStruct sn) {
        if (!isFullyOpen(sn)) {
            return false;
        }
        for (int i = 0; i < sn.nstations; i++) {
            if (sn.sched.get(sn.stations.get(i)) == SchedStrategy.INF) {
                return false;
            }
            if (sn.nservers.get(i, 0) > 1) {
                return false;
            }
        }
        return true;
    }

    /** Whether the QR/LR/QRF reduction applies: one closed class, no INF, c=1. */
    /**
     * What {@code 'default'}/{@code 'auto'} must mean on a model with finite
     * buffers: the QRF BAS bound, or nothing.
     *
     * <p>{@code 'default'} resolves to the geometric upper bound, which is
     * parameterized by demands and a population alone and therefore bounds a
     * blocked model as if its buffers were unbounded. {@code ignoresBlocking}
     * refuses that, which is right; but refusing is not the whole answer,
     * because {@code 'qrf.bas'} DOES model the finite buffer and
     * {@link SnToQrfBlocking} derives its tables from the model, so there is
     * nothing left for the caller to supply. A blocked model of the right shape
     * therefore gets {@code 'qrf.bas'} as its default, exactly as SolverMVA
     * routes a BAS model to {@code 'sqd'}.
     *
     * <p>The shape is the one {@link #isReducible} tests and the QRF analyzer
     * gates on. On top of it the tables must actually derive, which is asked of
     * {@code SnToQrfBlocking} rather than re-tested here -- it owns the
     * single-finite-buffer rule and the size guard, and a second copy of either
     * is how the two drift apart.
     *
     * <p>Only the UPPER side is routed: the analyzer solves qrf.bas in the
     * 'max' direction alone, so {@code 'auto.lower'} has no blocking
     * counterpart and keeps refusing rather than being answered with the wrong
     * side.
     *
     * @param sn network structure
     * @return {@code {method, why}}; method is empty when the routing does not
     *         apply, and why then carries the reason (empty for an unblocked model)
     */
    private static String[] blockingDefault(NetworkStruct sn) {
        if (!SnHasBlocking.snHasBlocking(sn)) {
            return new String[]{"", ""};
        }
        if (sn.nclasses != 1) {
            return new String[]{"", "the QRF blocking bounds are derived for a single-class"
                    + " closed network, which this model is not."};
        }
        for (int r = 0; r < sn.njobs.getNumCols(); r++) {
            if (Double.isInfinite(sn.njobs.get(0, r))) {
                return new String[]{"", "the QRF blocking bounds are derived for a single-class"
                        + " closed network, which this model is not."};
            }
        }
        for (int i = 0; i < sn.nstations; i++) {
            if (sn.sched.get(sn.stations.get(i)) == SchedStrategy.INF) {
                return new String[]{"", "the QRF blocking bounds model every station as a single"
                        + " server and have no infinite-server notion, so a delay station rules"
                        + " them out."};
            }
            if (sn.nservers.get(i, 0) > 1) {
                return new String[]{"", "the QRF blocking bounds model every station as a single"
                        + " server, so a multiserver station rules them out."};
            }
        }
        SnToQrfBlocking.Result derived =
                SnToQrfBlocking.snToQrfBlocking(sn, SnToQrfBlocking.DEFAULT_MAXVARS);
        if (!derived.msg.isEmpty()) {
            return new String[]{"", derived.msg};
        }
        return new String[]{"qrf.bas", ""};
    }

    /**
     * Whether the LOAD-DEPENDENT reduction applies where {@link #isReducible}
     * says no: one closed class, and a scaling alpha(i,n) that
     * {@link SnToQrfAlpha} can derive for every station.
     */
    private static boolean isLdReducible(NetworkStruct sn) {
        if (sn.nclasses != 1) {
            return false;
        }
        for (int r = 0; r < sn.njobs.getNumCols(); r++) {
            if (Double.isInfinite(sn.njobs.get(0, r))) {
                return false;
            }
        }
        return SnToQrfAlpha.snToQrfAlpha(sn).msg.isEmpty();
    }

    private static boolean isReducible(NetworkStruct sn) {
        if (sn.nclasses != 1) {
            return false;
        }
        for (int r = 0; r < sn.njobs.getNumCols(); r++) {
            if (Double.isInfinite(sn.njobs.get(0, r))) {
                return false;
            }
        }
        for (int i = 0; i < sn.nstations; i++) {
            if (sn.sched.get(sn.stations.get(i)) == SchedStrategy.INF) {
                return false;
            }
            if (sn.nservers.get(i, 0) > 1) {
                return false;
            }
        }
        return true;
    }

    /**
     * The feature envelope of the bounds.
     *
     * <p>SolverBA DECLARED NO FEATURE SET AT ALL, so it inherited
     * {@code Solver.supports}, which returns true: every model in the language
     * was accepted, including the ones whose constructs the analyzer has no
     * representation of. MATLAB had the mirror-image defect -- a set that named
     * no service distribution and therefore refused everything -- and both are
     * fixed to the same list.
     *
     * <p>WHAT MAKES A DISTRIBUTION ADMISSIBLE HERE IS ITS MEAN. The analyzer
     * reads sn.rates and sn.visits and nothing else: every bound in the
     * ABA/BJB/PB/GB/SB/Harel/MWBA families is a function of the demands
     * D = V./rates and the think time, so any renewal law with a finite mean is
     * admissible whatever its higher moments. The QRF reduction is the one that
     * needs more, and what it needs is a PH representation (the {D0,D1} pair out
     * of sn.proc), which the phase-type families below carry.
     *
     * <p>THE MODULATED LAWS "MAP" AND "MMPP2" ARE IN THE BASE ENVELOPE FOR ONE
     * FAMILY, 'mapamva', AND {@link #getMethodFeatureSet} STRIPS THEM FROM EVERY
     * OTHER. They were out entirely until mapamva landed, on the correct ground
     * that a renewal bound derived for a product-form network says nothing about
     * a correlated one: its mean rate exists, so the utilization law still holds
     * and the formula still returns a number, but that number brackets a
     * DIFFERENT system. MAP-AMVA (Casale-Smirni, DSN 2009) is derived FOR the
     * correlated model -- its variables are the per-phase QN(i,k) and UN(i,k) --
     * so the same reasoning that refuses the others admits it. The direction is
     * forced: a feature set refuses a model for HAVING a construct and never for
     * lacking one, so the only way to grant a law to one family is to put it in
     * the base envelope and take it away from the rest. MMAP and BMAP stay OUT
     * everywhere, and Cache and Fork/Join stay out for the older reason, that no
     * analyzer here has a representation of either.
     *
     * <p>"JobSink" IS DECLARED HERE AND NOT IN MATLAB, and the difference is in
     * the EMITTER, not in the envelope. {@code Network.getUsedLangFeatures}
     * marks a Sink with both "Sink" and "JobSink" (the node and its input
     * section); {@code MNetwork.getUsedLangFeatures} marks only "Sink". Naming
     * one and not the other therefore refused every OPEN model here while
     * MATLAB accepted it -- and the three open-network bounds 'bpt.lower',
     * 'bgt.upper' and 'snc.upper' are exactly the ones an open model is for, so
     * SolverAUTO.listValidMethods offered no 'ba.*' method name on a model whose
     * bounds this solver computes.
     *
     * @return the features the bounds can consume
     */
    public static FeatureSet getFeatureSet() {
        FeatureSet featSupported = new FeatureSet();
        featSupported.setTrue(new String[]{
                "ClassSwitch", "Delay", "DelayStation", "Queue",
                "Sink", "JobSink", "Source", "Router",
                "StatelessClassSwitcher",
                "ClosedClass", "OpenClass",
                // renewal service laws: the bounds need the mean, the QRF
                // reduction needs the PH form, and sn.proc carries both
                "APH", "Coxian", "Cox2", "Erlang", "Exp", "HyperExp", "PH",
                "Det", "Lognormal", "Pareto", "Uniform", "Weibull",
                // modulated service, for 'mapamva' alone; see the note above
                "MAP", "MMPP2",
                "SchedStrategy_INF", "SchedStrategy_PS",
                "SchedStrategy_FCFS", "SchedStrategy_LCFSPR",
                "RoutingStrategy_PROB", "RoutingStrategy_RAND",
                // Petri-net constructs, for the spnlp family. QueueingPlace is
                // deliberately absent: a place with an embedded queue has local
                // state the relaxation has no variable for, and the analyzer
                // refuses it by name. Same division SolverNC draws.
                "Place", "Transition", "Linkage", "Enabling", "Inhibiting",
                "Timing", "Firing", "Storage"});
        return featSupported;
    }

    /**
     * Whether SolverBA can bound this model.
     *
     * <p>Finite-buffer blocking is NOT decided here: there is no registry
     * feature name for a capacity, so a blocked model passes this gate and is
     * handled where it can be told apart -- {@code listValidMethods} drops the
     * blocking-blind bounds and keeps {@code qrf.bas*}/{@code qrf.rsrd}, which
     * carry the blocking tables explicitly.
     *
     * @param model the network model to check
     * @return true when every feature the model uses is one the bounds consume
     */
    @Override
    public boolean supports(Network model) {
        if (model == null) {
            return false;
        }
        return FeatureSet.supports(SolverBA.getFeatureSet(), model.getUsedLangFeatures());
    }

    /** Every bound method the solver implements, independently of the model. */
    public static String[] listAllMethods() {
        List<String> m = new ArrayList<String>(Arrays.asList(
                "default",
                "auto.upper", "auto.lower",
                "aba.upper", "aba.lower",
                "bjb.upper", "bjb.lower",
                "pb.upper", "pb.lower",
                "gb.upper", "gb.lower",
                "sb.upper", "sb.lower",
                "harel.upper", "harel.lower",
                "mwba.upper", "mwba.lower",
                "pbh.upper", "pbh.lower",
                "pbk.upper", "pbk.lower",
                "bjbk.upper", "bjbk.lower",
                "cbh.upper", "cbh.lower",
                "ssd.upper", "ssd.lower",
                "cub.upper", "mbjb.lower",
                "looping.upper", "looping.lower",
                "sib.upper", "sib.lower",
                "scb.upper", "scb.lower",
                "ldbcmp.lower",
                "bpt.lower", "bgt.upper", "snc.upper",
                // The three BARE ALIASES resolveMethod maps -- "auto" ->
                // auto.upper, "qr" -> qrf.mmi, "lr" -> lr.upper -- have to be
                // declared here, because checkDeclaredMethod validates the RAW
                // name against this list BEFORE resolveMethod ever runs. C++
                // and python resolve first and then validate, so they need no
                // such entry; this list is the JAR's compensation for that
                // asymmetry. "qr" and "lr" were already carried for exactly
                // this reason and "auto" was simply missed, which made the
                // class doc's "bare 'auto' means 'auto.upper'" untrue and the
                // "auto" branch of resolveMethod unreachable.
                "auto", "qr", "lr", "lr.upper", "lr.lower",
                "mapamva.upper", "mapamva.lower",
                "qrf.mmi", "qrf.mem", "qrf.bethe", "qrf.mmi.ld", "qrf.mmi.linear",
                "qrf.bas.mmi", "qrf.bas.mem", "qrf.bas.bethe", "qrf.bas", "qrf.rsrd",
                "spnlp.upper", "spnlp.lower", "spnlp.op.upper", "spnlp.op.lower"
        ));
        return m.toArray(new String[0]);
    }
}
