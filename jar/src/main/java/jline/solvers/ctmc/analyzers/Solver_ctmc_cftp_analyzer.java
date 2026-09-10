package jline.solvers.ctmc.analyzers;

import java.util.Arrays;

import jline.api.pfqn.Pfqn_cftp;
import jline.api.sn.SnGetDemandsChain;
import jline.io.Ret;
import jline.lang.NetworkStruct;
import jline.lang.constant.NodeType;
import jline.lang.constant.RoutingStrategy;
import jline.lang.constant.SchedStrategy;
import jline.lang.nodes.Station;
import jline.solvers.SolverOptions;
import jline.util.matrix.Matrix;

/**
 * Perfect-sampling steady-state analysis of closed single-class product-form
 * networks, the sampling alternative to CTMC state-space enumeration.
 *
 * <p>States are drawn iid from the exact stationary distribution by monotone
 * Coupling From The Past ({@link Pfqn_cftp}), so the estimator carries Monte
 * Carlo error O(samples^(-1/2)) but never enumerates the state space.</p>
 *
 * <p>Reference: S. Kijima and T. Matsui, "Approximate/Perfect Samplers for
 * Closed Jackson Networks", Proc. Winter Simulation Conference, 2005.</p>
 */
public class Solver_ctmc_cftp_analyzer {

    private Solver_ctmc_cftp_analyzer() {}

    /** Metrics and sampled states returned by the perfect-sampling analyzer. */
    public static class CftpResult {
        public Matrix QN;
        public Matrix UN;
        public Matrix RN;
        public Matrix TN;
        public Matrix CN;
        public Matrix XN;
        /** Sampled states, one per row (samples x stations). */
        public Matrix samples;
        /** Per-sample coalescence horizon or mixing steps. */
        public Matrix horizon;
        /** Distinct sampled states. */
        public Matrix spaceAggr;
        /** Empirical probability of each row of spaceAggr. */
        public Matrix pi;
    }

    /** Maps the solver method string onto the Pfqn_cftp sampler name. */
    public static String sampler(String method) {
        String m = method == null ? "" : method.toLowerCase();
        if (m.equals("cftp") || m.equals("cftp.exact")) {
            return "cftp";
        }
        if (m.equals("cftp.approx")) {
            return "approx";
        }
        throw new RuntimeException("Unknown cftp variant '" + method + "'. Use 'cftp' or 'cftp.approx'.");
    }

    /**
     * Steady-state metrics of a closed single-class product-form network from
     * iid perfect samples of its stationary distribution.
     *
     * @param sn      network structure
     * @param options solver options carrying method, samples and timespan
     * @return the metrics together with the sampled states
     */
    public static CftpResult solver_ctmc_cftp(NetworkStruct sn, SolverOptions options) {
        String samplerName = sampler(options.method);
        assertSupported(sn, options);

        int M = sn.nstations;
        int K = sn.nclasses;

        Ret.snGetDemands dem = SnGetDemandsChain.snGetDemandsChain(sn);
        Matrix Lchain = dem.Dchain;
        Matrix STchain = dem.STchain;
        Matrix Vchain = dem.Vchain;
        int N = (int) Math.round(dem.Nchain.get(0));
        int iref = (int) Math.round(dem.refstatchain.get(0));

        Matrix L = new Matrix(1, M);
        Matrix S = new Matrix(1, M);
        for (int i = 0; i < M; i++) {
            L.set(0, i, Lchain.get(i, 0));
            Station sti = sn.stations.get(i);
            if (sn.sched.get(sti) == SchedStrategy.INF) {
                S.set(0, i, Double.POSITIVE_INFINITY);
            } else {
                S.set(0, i, sn.nservers.get(i));
            }
        }

        int nsamples = options.samples;
        if (nsamples < 1) {
            throw new RuntimeException("The cftp method requires a positive options.samples.");
        }

        Pfqn_cftp.PfqnCftpReturn draw = Pfqn_cftp.pfqn_cftp(L, N, S, nsamples, samplerName);

        double[] busy = new double[M];
        for (int i = 0; i < M; i++) {
            double acc = 0.0;
            for (int s = 0; s < nsamples; s++) {
                acc += Math.min(draw.X.get(s, i), S.get(0, i));
            }
            busy[i] = acc / nsamples;
        }

        CftpResult res = new CftpResult();
        res.QN = new Matrix(M, K);
        res.UN = new Matrix(M, K);
        res.RN = new Matrix(M, K);
        res.TN = new Matrix(M, K);
        res.CN = new Matrix(1, K);
        res.XN = new Matrix(1, K);

        // The utilization law X = mu_i*E[min(n_i,c_i)]/V_i holds at every station, but
        // each station estimates it with its own Monte Carlo error. Taking the estimate
        // at the reference station and propagating it through the visit ratios matches
        // the CTMC convention (XN is the arrival rate at the reference station) and
        // keeps flow balance, Little's law and C = N/X exact in the reported table.
        double X = 0.0;
        if (STchain.get(iref, 0) > 0 && Vchain.get(iref, 0) > 0) {
            X = busy[iref] / STchain.get(iref, 0) / Vchain.get(iref, 0);
        }
        res.XN.set(0, 0, X);
        if (X > 0) {
            res.CN.set(0, 0, N / X);
        }

        for (int i = 0; i < M; i++) {
            res.QN.set(i, 0, draw.Q.get(0, i));
            res.TN.set(i, 0, Vchain.get(i, 0) * X);
            // Utilization keeps its own estimator E[min(n_i,c_i)]/c_i: it is unbiased
            // and confined to [0,1] by construction, whereas deriving it from the
            // reference-station throughput lets Monte Carlo error push a saturated
            // station above 1.
            Station sti = sn.stations.get(i);
            if (sn.sched.get(sti) == SchedStrategy.INF) {
                res.UN.set(i, 0, res.QN.get(i, 0));
            } else {
                res.UN.set(i, 0, busy[i] / S.get(0, i));
            }
            if (res.TN.get(i, 0) > 0) {
                res.RN.set(i, 0, res.QN.get(i, 0) / res.TN.get(i, 0));
            }
        }

        res.samples = draw.X;
        res.horizon = draw.T;
        distinctStates(draw.X, nsamples, M, res);
        return res;
    }

    /** Collapses the sample matrix into distinct states and their frequencies. */
    private static void distinctStates(Matrix Xs, int nsamples, int M, CftpResult res) {
        double[][] rows = new double[nsamples][M];
        for (int s = 0; s < nsamples; s++) {
            for (int i = 0; i < M; i++) {
                rows[s][i] = Xs.get(s, i);
            }
        }
        Arrays.sort(rows, new java.util.Comparator<double[]>() {
            public int compare(double[] a, double[] b) {
                for (int i = 0; i < a.length; i++) {
                    int c = Double.compare(a[i], b[i]);
                    if (c != 0) {
                        return c;
                    }
                }
                return 0;
            }
        });
        int distinct = 0;
        for (int s = 0; s < nsamples; s++) {
            if (s == 0 || !Arrays.equals(rows[s], rows[s - 1])) {
                distinct++;
            }
        }
        res.spaceAggr = new Matrix(distinct, M);
        res.pi = new Matrix(distinct, 1);
        int r = -1;
        for (int s = 0; s < nsamples; s++) {
            if (s == 0 || !Arrays.equals(rows[s], rows[s - 1])) {
                r++;
                for (int i = 0; i < M; i++) {
                    res.spaceAggr.set(r, i, rows[s][i]);
                }
            }
            res.pi.set(r, 0, res.pi.get(r, 0) + 1.0 / nsamples);
        }
    }

    /**
     * Can the cftp perfect sampler be asked for this model?
     *
     * <p>The model-class gate asked as a predicate rather than thrown.
     * {@link #assertSupported} refuses with it, and
     * {@code SolverCTMC.supportsModelMethod} asks the very same call so that a
     * caller (model.help, findSolver, SolverAUTO) sees the verdict before
     * paying for a run. A second copy of the rules is how the report and the
     * run drift into two different answers.</p>
     *
     * <p>The sampler is exact only on the closed single-class product form its
     * balance function encodes; anything else must be refused, not
     * approximated. What the feature registry CAN name is also declared in
     * {@code SolverCTMC.getMethodFeatureSet}; this predicate is what carries
     * the structural rules the registry has no name for -- the class count, the
     * station count, the phase count and the steady-state restriction.</p>
     *
     * @param sn the network structure
     * @param options solver options, read for the timespan
     * @return empty string when the sampler may run, else the refusal
     */
    public static String supportsReason(NetworkStruct sn, SolverOptions options) {
        if (!Double.isInfinite(options.timespan[0])) {
            return ("The cftp method supports steady-state analysis only, not transient analysis.");
        }
        if (sn.nclasses != 1) {
            return ("The cftp method supports single-class models only, this model has "
                    + sn.nclasses + " classes.");
        }
        if (Double.isInfinite(sn.njobs.get(0)) || sn.njobs.get(0) < 1) {
            return ("The cftp method supports closed models only, with a finite positive population.");
        }
        if (sn.nstations < 2) {
            return ("The cftp method requires at least two stations.");
        }
        for (int ind = 0; ind < sn.nnodes; ind++) {
            NodeType nt = sn.nodetype.get(ind);
            if (nt != NodeType.Queue && nt != NodeType.Delay && nt != NodeType.Router) {
                return ("The cftp method supports Queue, Delay and Router nodes only, node "
                        + ind + " is of a different type.");
            }
        }
        for (int i = 0; i < sn.nstations; i++) {
            Station sti = sn.stations.get(i);
            SchedStrategy sched = sn.sched.get(sti);
            if (sched != SchedStrategy.INF && sched != SchedStrategy.PS && sched != SchedStrategy.FCFS
                    && sched != SchedStrategy.SIRO && sched != SchedStrategy.LCFSPR) {
                return ("The cftp method requires a product-form scheduling strategy "
                        + "(INF, PS, FCFS, SIRO, LCFSPR) at station " + i + ".");
            }
            if (sn.phases != null && sn.phases.get(i, 0) > 1) {
                return ("The cftp method requires exponential service times, station "
                        + i + " has " + (int) sn.phases.get(i, 0) + " phases.");
            }
            if (Double.isFinite(sn.cap.get(i)) && sn.cap.get(i) < sn.njobs.get(0)) {
                return ("The cftp method requires infinite buffers, station " + i
                        + " has capacity " + (int) sn.cap.get(i) + ".");
            }
            if (!Double.isFinite(sn.rates.get(i, 0)) || sn.rates.get(i, 0) <= 0) {
                return ("The cftp method requires a finite positive service rate at station "
                        + i + ".");
            }
        }
        if (sn.lldscaling != null && sn.lldscaling.getNumElements() > 0) {
            return ("The cftp method does not support load-dependent service rates.");
        }
        if (sn.cdscaling != null && !sn.cdscaling.isEmpty()) {
            return ("The cftp method does not support class-dependent service rates.");
        }
        if (sn.jdscaling != null && !sn.jdscaling.isEmpty()) {
            return ("The cftp method does not support joint-dependent service rates.");
        }
        if (sn.nregions > 0) {
            return ("The cftp method does not support finite capacity regions.");
        }
        for (int ind = 0; ind < sn.nnodes; ind++) {
            jline.lang.nodes.Node node = sn.nodes.get(ind);
            if (sn.routing == null || !sn.routing.containsKey(node)) {
                continue;
            }
            for (RoutingStrategy rs : sn.routing.get(node).values()) {
                if (rs != RoutingStrategy.PROB && rs != RoutingStrategy.RAND && rs != RoutingStrategy.DISABLED) {
                    return ("The cftp method requires Markovian routing (PROB, RAND), node "
                            + ind + " uses a state-dependent strategy.");
                }
            }
        }
        return "";
    }

    /** Refuses every model outside the closed single-class product form. */
    private static void assertSupported(NetworkStruct sn, SolverOptions options) {
        String reason = supportsReason(sn, options);
        if (!reason.isEmpty()) {
            throw new RuntimeException(reason);
        }
    }
}
