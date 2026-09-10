package jline.solvers.ctmc.analyzers;

import java.util.List;

import jline.api.mdd.MDD;
import jline.api.mdd.MddDescriptor;
import jline.api.mdd.MddMcdOptions;
import jline.api.mdd.MddMcdResult;
import jline.api.mdd.MddServiceLaw;
import jline.api.mdd.Mdd_descriptor;
import jline.api.mdd.Mdd_mcd;
import jline.api.mdd.Mdd_ps;
import jline.api.mdd.Mdd_reachset;
import jline.lang.JobClass;
import jline.lang.NetworkStruct;
import jline.lang.constant.NodeType;
import jline.lang.constant.SchedStrategy;
import jline.lang.nodes.Station;
import jline.solvers.SolverOptions;
import jline.util.matrix.Matrix;
import jline.util.matrix.MatrixCell;

/**
 * Stationary analysis by decision-diagram level aggregation.
 *
 * <p>Analyses a closed single-class network whose CTMC state space is held in a
 * decision diagram and solved by level aggregation, after A.S. Miner,
 * G. Ciardo, S. Donatelli, "Using the exact state space of a Markov model to
 * compute approximate stationary measures", SIGMETRICS 2000.</p>
 *
 * <p>This is the 'mdd' method of SolverCTMC. It never forms the |S|-state
 * generator: the reachable set is stored in an MDD and K coupled level-CTMCs are
 * iterated to a fixed point, so the memory cost is O(sum_k |M_k|) rather than
 * O(|S|). The saving grows with the number of stations, and is negative at K=3,
 * where the diagram compresses nothing.</p>
 *
 * <p><b>Exactness.</b> The single approximation is Pr{i_k | alpha} =
 * Pr{i_k | p}. It is EXACT on product-form networks (paper Sec. 5), which covers
 * exponential service under any work-conserving discipline and general service
 * at PS or IS stations (BCMP types 2 and 3). It is an approximation otherwise,
 * notably phase-type service at FCFS or LCFS, where errors of a fraction of a
 * percent on the mean queue lengths have been observed.</p>
 *
 * <p>MATLAB twin: {@code solver_ctmc_mdd_analyzer.m}. Python twin:
 * {@code api/solvers/ctmc/solver_ctmc_mdd_analyzer.py}.</p>
 */
public class Solver_ctmc_mdd_analyzer {

    private Solver_ctmc_mdd_analyzer() {}

    /** Metrics and diagram description returned by the aggregation analyzer. */
    public static class MddResult {
        public Matrix QN;
        public Matrix UN;
        public Matrix RN;
        public Matrix TN;
        public Matrix CN;
        public Matrix XN;
        /** The reachable set as stored. */
        public MDD mdd;
        /** The Kronecker descriptor that was built. */
        public MddDescriptor desc;
        /** |M_k| per paper level. */
        public int[] levelSizes;
        /** Fixed-point sweeps performed. */
        public int iters;
        /** |S|. */
        public long numStates;
        /** Local-state encoding chosen, "np" or "ps". */
        public String encoding;
        /** The SPN metadata, when the net route was taken; null otherwise. */
        public jline.api.spn.Spn_mdd.SpnInfo spn;
        /** Per-place-level marginal law, when the net route was taken. */
        public double[][] marginal;
    }

    /** Shared-server disciplines, where every job present holds its own phase. */
    private static boolean isShared(SchedStrategy s) {
        return s == SchedStrategy.PS || s == SchedStrategy.DPS || s == SchedStrategy.GPS
                || s == SchedStrategy.INF;
    }

    /** Non-preemptive disciplines the count-plus-phase encoding represents. */
    private static boolean isNonPreemptive(SchedStrategy s) {
        return s == SchedStrategy.FCFS || s == SchedStrategy.LCFS || s == SchedStrategy.SIRO
                || s == SchedStrategy.HOL;
    }

    /**
     * Analyse a closed single-class network by MDD level aggregation.
     *
     * @param sn the network structure
     * @param options solver options; the level knobs are read from
     *        options.config.mdd_tol / mdd_maxiter, never from iter_tol
     */
    public static MddResult solver_ctmc_mdd(NetworkStruct sn, SolverOptions options) {
        return solver_ctmc_mdd(sn, options, null);
    }

    /**
     * Can the mdd decision-diagram method be asked for this model?
     *
     * <p>The model-shape gate asked as a predicate rather than thrown. The
     * analyzer refuses with it before it builds anything, and
     * {@code SolverCTMC.supportsModelMethod} asks the very same call so that a
     * caller (model.help, findSolver, SolverAUTO) sees the verdict without
     * paying for a run. One predicate with two callers is what stops the report
     * and the analyzer from disagreeing about which models the method serves.</p>
     *
     * <p>A STOCHASTIC PETRI NET IS EXEMPT: a Place model is read through
     * {@link jline.api.spn.Spn_mdd}, which builds the reachable set and the
     * Kronecker descriptor from the marking rather than from the
     * (station,class) encoding, so neither the single-class rule nor the
     * closed-population rule applies to it.</p>
     *
     * <p>The deeper refusals the analyzer still raises -- a station-to-station
     * chain that is not stochastic, and a phase-type law at a discipline
     * neither local encoding represents -- are not restated here: they are
     * decided from quantities the analyzer computes on its way through, not
     * from the model shape, so a caller cannot be told about them without
     * doing the work.</p>
     *
     * @param sn the network structure
     * @return empty string when the method may run, else the refusal
     */
    public static String supportsReason(NetworkStruct sn) {
        if (sn.nodetype.contains(NodeType.Place)) {
            return "";
        }
        int R = sn.nclasses;
        if (R != 1) {
            return "the mdd method analyses single-class networks; this model has " + R
                    + " classes. The Kronecker descriptor would need one level per (station,class).";
        }
        for (NodeType nt : sn.nodetype) {
            if (nt == NodeType.Source || nt == NodeType.Sink) {
                return "the mdd method analyses CLOSED networks; an open stream makes the marking "
                        + "unbounded, so the reachable set has no finite decision diagram";
            }
        }
        double Nd = sn.njobs.get(0);
        if (Double.isInfinite(Nd) || Double.isNaN(Nd) || Nd <= 0) {
            return "the mdd method needs a finite positive closed population";
        }
        return "";
    }

    /**
     * Analyse a closed single-class network, or a stochastic Petri net, by MDD
     * level aggregation.
     *
     * <p>Passing {@code model} routes a net holding Places and Transitions through
     * {@link jline.api.spn.Spn_mdd} instead of {@link Mdd_descriptor}: the levels
     * are then (place, class) pairs plus one phase level per phase-type mode, and
     * the measures come back per place. The approximation is the same Eq. 5 as
     * for a queueing network, and it is exact on a product-form net, which
     * SolverNC's "rec" method solves exactly and far more cheaply -- the
     * aggregation earns its place on the nets that have NO product form.</p>
     *
     * <p>A net carries no per-station service rate, so Mdd_mcd returns only the
     * level marginals. The token throughput is then assembled here from the mode
     * rates and those marginals, under the SAME independence across levels that
     * the aggregation already assumes: it is the method's own approximation
     * applied once more, not a second one layered on top.</p>
     *
     * @param sn the network structure
     * @param options solver options
     * @param model the model object, needed only for the Petri-net route
     */
    public static MddResult solver_ctmc_mdd(NetworkStruct sn, SolverOptions options,
                                            jline.lang.Network model) {
        // The model-shape rules live in supportsReason, which
        // SolverCTMC.supportsModelMethod also asks: the analyzer must refuse
        // exactly what the report refuses, and one predicate with two callers
        // is what keeps the two from drifting apart.
        String shapeReason = supportsReason(sn);
        if (!shapeReason.isEmpty()) {
            throw new RuntimeException("solver_ctmc_mdd_analyzer: " + shapeReason);
        }
        if (sn.nodetype.contains(NodeType.Place)) {
            if (model == null) {
                throw new RuntimeException("solver_ctmc_mdd_analyzer: a stochastic Petri net is "
                        + "read from the model object, not from the network structure; call "
                        + "solver_ctmc_mdd(sn, options, model)");
            }
            return spn(model, sn, options);
        }
        int M = sn.nstations;
        int R = sn.nclasses;
        int N = (int) Math.round(sn.njobs.get(0));

        // ---- service laws and the station-to-station routing chain
        double[] mu = new double[M];
        double[] servers = new double[M];
        MddServiceLaw[] proc = new MddServiceLaw[M];
        SchedStrategy[] sched = new SchedStrategy[M];
        boolean anyPH = false;
        JobClass jobClass = sn.jobclasses.get(0);
        for (int i = 0; i < M; i++) {
            Station station = sn.stations.get(i);
            mu[i] = sn.rates.get(i, 0);                    // 1/E[S] by LINE convention
            servers[i] = sn.nservers.get(i, 0);
            sched[i] = sn.sched.get(station);
            int nph = 1;
            if (sn.phases != null && sn.phases.getNumRows() > i) {
                nph = (int) sn.phases.get(i, 0);
            }
            if (nph > 1) {
                MatrixCell pc = sn.proc.get(station).get(jobClass);
                if (pc == null || pc.size() < 2) {
                    throw new RuntimeException("solver_ctmc_mdd_analyzer: station " + (i + 1)
                            + " declares " + nph + " phases but carries no representable "
                            + "service law");
                }
                proc[i] = new MddServiceLaw(toArray(pc.get(0)), toArray(pc.get(1)));
                anyPH = true;
            }
        }
        double[][] P = routing(sn, M, R);
        String kind = encoding(sched, servers, proc, anyPH);

        MddDescriptor desc;
        if (kind.equals("ps")) {
            desc = Mdd_ps.mdd_ps(mu, P, servers, N, proc);
        } else {
            String[] names = new String[M];
            for (int i = 0; i < M; i++) {
                names[i] = sched[i] == null ? null : sched[i].toString().toUpperCase();
            }
            desc = Mdd_descriptor.mdd_descriptor(mu, P, servers, N, proc, names);
        }

        MDD mdd = Mdd_reachset.mdd_reachset(desc.domain, desc.init, desc.nextfun);

        // The level iteration is an INNER numerical solve and needs a far tighter
        // tolerance than the reported means: Mdd_mcd verifies the population
        // invariant at 1e-6, so a loose tolerance converges short of the fixed
        // point and trips that guard. options.iter_tol is the solver-level
        // fixed-point tolerance (sized for AMVA outer loops) and must NOT be
        // reused here; the level knobs come from options.config instead.
        MddMcdOptions mcdopt = new MddMcdOptions();
        if (options != null && options.config != null) {
            if (options.config.mdd_tol != null) {
                mcdopt.tol = options.config.mdd_tol.doubleValue();
            }
            if (options.config.mdd_maxiter != null) {
                mcdopt.maxiter = options.config.mdd_maxiter.intValue();
            }
        }
        MddMcdResult out = Mdd_mcd.mdd_mcd(mdd.toStruct(), desc, mcdopt);

        // ---- pack the analyzer contract
        MddResult res = new MddResult();
        res.QN = new Matrix(M, R);
        res.UN = new Matrix(M, R);
        res.RN = new Matrix(M, R);
        res.TN = new Matrix(M, R);
        for (int i = 0; i < M; i++) {
            res.QN.set(i, 0, out.QLen[i]);
            res.UN.set(i, 0, out.U[i]);
            res.TN.set(i, 0, out.X[i]);
            if (out.X[i] > 0) {
                res.RN.set(i, 0, out.QLen[i] / out.X[i]);   // Little's law at the station
            }
        }

        // system throughput at the reference station, per unit visit
        int ref = (int) sn.refstat.get(0);
        double vis = 1.0;
        if (sn.visits != null && sn.visits.containsKey(Integer.valueOf(0))) {
            Matrix v = sn.visits.get(Integer.valueOf(0));
            // visits is indexed by STATEFUL node, refstat by station
            int refsf = (int) sn.stationToStateful.get(ref);
            if (v != null && v.getNumRows() > refsf && v.get(refsf, 0) > 0) {
                vis = v.get(refsf, 0);
            }
        }
        res.XN = new Matrix(1, R);
        res.CN = new Matrix(1, R);
        res.XN.set(0, 0, res.TN.get(ref, 0) / vis);
        if (res.XN.get(0, 0) > 0) {
            res.CN.set(0, 0, N / res.XN.get(0, 0));
        }

        res.mdd = mdd;
        res.desc = desc;
        res.levelSizes = out.levelSizes;
        res.iters = out.iters;
        res.numStates = mdd.cardinality();
        res.encoding = kind;

        if (options != null && options.verbose != null
                && options.verbose != jline.VerboseLevel.SILENT) {
            int total = 0;
            for (int k = 0; k < out.levelSizes.length; k++) {
                total += out.levelSizes[k];
            }
            System.out.format("%nCTMC-mdd: %d levels, |S| = %d held as %d level states (%.1fx), "
                            + "%d fixed-point sweeps%n", desc.K, res.numStates, total,
                    (double) res.numStates / Math.max(total, 1), out.iters);
        }
        return res;
    }

    // -----------------------------------------------------------------------
    private static double[][] toArray(Matrix m) {
        int r = m.getNumRows();
        int c = m.getNumCols();
        double[][] a = new double[r][c];
        for (int i = 0; i < r; i++) {
            for (int j = 0; j < c; j++) {
                a[i][j] = m.get(i, j);
            }
        }
        return a;
    }

    /**
     * Station-to-station routing probabilities of the single class.
     *
     * <p>Read from sn.rt, which is indexed over stateful nodes, class-major.</p>
     */
    private static double[][] routing(NetworkStruct sn, int M, int R) {
        // stationToStateful is stored as a ROW vector in the JAR, so it is read
        // by linear index; get(i,0) is out of bounds past the first station
        double[][] P = new double[M][M];
        for (int i = 0; i < M; i++) {
            int ni = (int) sn.stationToStateful.get(i);
            double rs = 0;
            for (int j = 0; j < M; j++) {
                int nj = (int) sn.stationToStateful.get(j);
                P[i][j] = sn.rt.get(ni * R, nj * R);
                rs += P[i][j];
            }
            if (Math.abs(rs - 1.0) > 1e-8) {
                throw new RuntimeException("solver_ctmc_mdd_analyzer: the station-to-station "
                        + "routing chain is not stochastic; the mdd method needs every completion "
                        + "to move the job to another station");
            }
        }
        return P;
    }

    /**
     * Which local-state encoding represents these disciplines exactly.
     *
     * <p>Exponential service is discipline-insensitive for the queue-length law,
     * so the compact count encoding serves any work-conserving station.
     * Phase-type service is not: the count-plus-one-phase encoding is
     * non-preemptive, while processor sharing needs the per-phase counts of
     * every job present.</p>
     */
    private static String encoding(SchedStrategy[] sched, double[] servers,
                                   MddServiceLaw[] proc, boolean anyPH) {
        int M = sched.length;
        if (!anyPH) {
            return "np";
        }
        boolean allShared = true;
        for (int i = 0; i < M; i++) {
            if (!isShared(sched[i])) {
                allShared = false;
                break;
            }
        }
        if (allShared) {
            return "ps";
        }
        boolean allNP = true;
        int firstBad = -1;
        for (int i = 0; i < M; i++) {
            if (proc[i] == null) {
                continue;
            }
            boolean np = isNonPreemptive(sched[i]) && servers[i] == 1;
            if (!np) {
                allNP = false;
                if (firstBad < 0 && !isShared(sched[i])) {
                    firstBad = i;
                }
                if (firstBad < 0) {
                    firstBad = i;
                }
            }
        }
        if (allNP) {
            return "np";
        }
        throw new RuntimeException("solver_ctmc_mdd_analyzer: station " + (firstBad + 1)
                + " combines a phase-type service law with a discipline that neither local "
                + "encoding represents: the count-plus-phase encoding is non-preemptive, and the "
                + "per-phase-count encoding covers only shared servers (PS/DPS/GPS/INF). Mixing a "
                + "shared and a non-preemptive phase-type station in one model is likewise "
                + "unsupported.");
    }

    /**
     * The Petri-net route: Spn_mdd supplies the reachable set and the Kronecker
     * descriptor, Mdd_mcd aggregates, and the measures are read back per place.
     */
    private static MddResult spn(jline.lang.Network model, NetworkStruct sn,
                                 SolverOptions options) {
        jline.api.spn.Spn_mdd.SpnOptions mddopt = new jline.api.spn.Spn_mdd.SpnOptions();
        mddopt.verbose = options != null && options.verbose != null
                && options.verbose.ordinal() > 1;
        jline.api.spn.Spn_mdd.SpnResult spn = jline.api.spn.Spn_mdd.spn_mdd(model, mddopt);

        MddMcdOptions mcdopt = new MddMcdOptions();
        if (options != null && options.config != null) {
            if (options.config.mdd_maxiter != null) {
                mcdopt.maxiter = options.config.mdd_maxiter.intValue();
            }
            if (options.config.mdd_tol != null) {
                mcdopt.tol = options.config.mdd_tol.doubleValue();
            }
        }
        MddMcdResult out = Mdd_mcd.mdd_mcd(spn.mdds, spn.desc, mcdopt);

        int M = sn.nstations;
        int R = sn.nclasses;
        int L = spn.info.nplacelevels;
        Matrix QN = new Matrix(M, R);
        Matrix UN = new Matrix(M, R);
        Matrix RN = new Matrix(M, R);
        Matrix TN = new Matrix(M, R);
        QN.fill(0.0);
        UN.fill(0.0);
        RN.fill(0.0);
        TN.fill(0.0);

        int[] places = spn.info.places;
        for (int pp = 0; pp < places.length; pp++) {
            int ist = (int) sn.nodeToStation.get(places[pp]);
            if (ist < 0) {
                continue;
            }
            for (int k = 0; k < R; k++) {
                QN.set(ist, k, out.QLen[pp * R + k]);
                UN.set(ist, k, out.QLen[pp * R + k]);      // a Place is an INF station: U = Q
            }
        }

        // Mode throughputs from the level marginals. P(m_l = v) is read off the
        // level chain; the enabling degree of a mode is then treated as
        // independent across its input levels, which is Eq. 5 of the paper applied
        // once more rather than a fresh approximation.
        double[][] pl = levelMarginals(out, spn.mdds, L);
        for (int e = 0; e < spn.info.modes.size(); e++) {
            jline.api.spn.Spn_mdd.SpnMode mde = spn.info.modes.get(e);
            if (mde.nph > 1) {
                continue;                                  // no single rate; read the phase level
            }
            double x = mde.D1[0][0] * meanServers(pl, mde, L);
            // FIRING EVENTS, not tokens: snPnAvgRates converts the Place rows to a
            // token rate afterwards, exactly as it does for the explicit CTMC path,
            // and weighting here as well would count a weighted arc twice.
            for (int l = 0; l < L; l++) {
                if (mde.enab[l] > 0) {
                    int pp = l / R;
                    int k = l % R;
                    int ist = (int) sn.nodeToStation.get(places[pp]);
                    if (ist >= 0) {
                        TN.set(ist, k, TN.get(ist, k) + x);
                    }
                }
            }
        }
        for (int i = 0; i < M; i++) {
            for (int k = 0; k < R; k++) {
                if (TN.get(i, k) > 0) {
                    RN.set(i, k, QN.get(i, k) / TN.get(i, k));
                }
            }
        }

        Matrix XN = new Matrix(1, R);
        Matrix CN = new Matrix(1, R);
        XN.fill(0.0);
        CN.fill(0.0);
        for (int k = 0; k < R; k++) {
            int ref = (int) sn.refstat.get(k);
            if (ref >= 0 && ref < M) {
                XN.set(0, k, TN.get(ref, k));
            }
            double Nk = 0;
            for (int i = 0; i < M; i++) {
                Nk += QN.get(i, k);
            }
            if (XN.get(0, k) > 0 && Nk > 0) {
                CN.set(0, k, Nk / XN.get(0, k));
            }
        }

        MddResult res = new MddResult();
        res.QN = QN;
        res.UN = UN;
        res.RN = RN;
        res.TN = TN;
        res.CN = CN;
        res.XN = XN;
        res.mdd = spn.info.mdd;
        res.desc = spn.desc;
        res.levelSizes = out.levelSizes;
        res.iters = out.iters;
        res.numStates = spn.info.mdd.cardinality();
        res.encoding = "spn";
        res.spn = spn.info;
        res.marginal = pl;
        return res;
    }

    /**
     * P(level l = v) for each place level, from the converged level chains.
     * Mdd_mcd works in the paper's orientation, paper level k = level K+1-l.
     */
    private static double[][] levelMarginals(MddMcdResult out,
                                             jline.api.mdd.MddStruct mdds, int L) {
        int K = mdds.K;
        double[][] pl = new double[L][];
        for (int l = 0; l < L; l++) {
            int k = K - 1 - l;
            double[] p = new double[mdds.domain[l]];
            int[][] rows = out.Mrows[k];
            double[] pk = out.pik[k];
            for (int r = 0; r < rows.length; r++) {
                p[rows[r][1]] += pk[r];
            }
            double tot = 0;
            for (int v = 0; v < p.length; v++) {
                tot += p[v];
            }
            if (tot > 0) {
                for (int v = 0; v < p.length; v++) {
                    p[v] /= tot;
                }
            }
            pl[l] = p;
        }
        return pl;
    }

    /** E[min(enabling degree, servers)] under independence across input levels. */
    private static double meanServers(double[][] pl, jline.api.spn.Spn_mdd.SpnMode mde, int L) {
        java.util.List<Integer> lv = new java.util.ArrayList<Integer>();
        for (int l = 0; l < L; l++) {
            if (mde.enab[l] > 0) {
                lv.add(Integer.valueOf(l));
            }
        }
        if (lv.isEmpty()) {
            return 1.0;
        }
        double kmax = Double.POSITIVE_INFINITY;
        for (int i = 0; i < lv.size(); i++) {
            int l = lv.get(i).intValue();
            kmax = Math.min(kmax, Math.floor((pl[l].length - 1) / mde.enab[l]));
        }
        if (!Double.isInfinite(mde.srv)) {
            kmax = Math.min(kmax, mde.srv);
        }
        double n = 0;
        for (int k = 1; k <= (int) kmax; k++) {
            double ge = 1;                                 // P(deg >= k) = prod_l P(m_l >= k I_l)
            for (int i = 0; i < lv.size() && ge > 0; i++) {
                int l = lv.get(i).intValue();
                int thr = (int) (k * mde.enab[l]);
                if (thr >= pl[l].length) {
                    ge = 0;
                    break;
                }
                double s = 0;
                for (int v = thr; v < pl[l].length; v++) {
                    s += pl[l][v];
                }
                ge *= s;
            }
            n += ge;                                       // E[min] = sum_k P(min >= k)
        }
        return n;
    }
}
