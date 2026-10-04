package jline.solvers.mam.handlers;

import java.util.HashMap;
import java.util.HashSet;
import java.util.Map;
import java.util.Set;

import jline.GlobalConstants;
import jline.api.mam.Map_cdf;
import jline.api.mam.Map_lambda;
import jline.api.mam.Map_mean;
import jline.api.mam.Map_pie;
import jline.api.mam.Map_scale;
import jline.api.mam.Map_var;
import jline.api.mam.Mmap_super;
import jline.api.map.MAPM1PSCdfRespT;
import jline.io.InputOutput;
import jline.io.InputOutput;
import jline.io.InputOutput;
import jline.lang.NetworkStruct;
import jline.lang.constant.SchedStrategy;
import jline.lib.butools.MMAPPH1FCFS;
import jline.lib.butools.MMAPPH1NPPR;
import jline.lib.butools.MMAPPH1PRPR;
import jline.solvers.SolverOptions;
import jline.util.matrix.Matrix;
import jline.util.matrix.MatrixCell;

public final class Solver_mam_passage_time {
    private Solver_mam_passage_time() {}

    public static PriorityAnalysis analyzePriorities(NetworkStruct sn, String context) {
        Matrix priorities = sn.classprio;
        double firstPriority = priorities.get(0);
        boolean identical = true;
        for (int i = 1; i < priorities.length(); i++) {
            if (priorities.get(i) != firstPriority) { identical = false; break; }
        }
        if (identical) {
            return new PriorityAnalysis(true, false, true, "Identical priorities - fully supported");
        }
        Set<Double> uniquePriorities = new HashSet<Double>();
        for (int i = 0; i < priorities.length(); i++) uniquePriorities.add(priorities.get(i));
        boolean allDistinct = uniquePriorities.size() == priorities.length();
        if (allDistinct) {
            return new PriorityAnalysis(false, true, true, "All distinct priorities - using NPPR analysis");
        }
        return new PriorityAnalysis(false, false, false, "Mixed priority configuration not supported - requires identical or all distinct priorities");
    }

    public static PriorityAnalysis analyzePriorities(NetworkStruct sn) {
        return analyzePriorities(sn, "basic");
    }

    public static Map<Integer, MatrixCell> solver_mam_passage_time(NetworkStruct sn,
                                                                    Map<jline.lang.nodes.Station, Map<jline.lang.JobClass, MatrixCell>> PH,
                                                                    SolverOptions options) {
        Map<Integer, MatrixCell> RD = new HashMap<Integer, MatrixCell>();
        int M = sn.nstations;
        int K = sn.nclasses;
        Matrix N = sn.njobs.transpose();

        boolean open = true;
        for (int i = 0; i < N.length(); i++) {
            if (Double.isFinite(N.get(i))) { open = false; break; }
        }

        MatrixCell A = new MatrixCell();
        int idx_arv = 0;
        int idx_q = 0;
        boolean is_ps = false;

        if (M == 2 && open) {
            Map<Integer, Matrix> pie = new HashMap<Integer, Matrix>();
            Map<Integer, Matrix> S = new HashMap<Integer, Matrix>();
            for (int i = 0; i < M; i++) {
                jline.lang.nodes.Station station = sn.stations.get(i);
                SchedStrategy schedI = sn.sched.get(station);
                if (schedI == SchedStrategy.EXT) {
                    Map<jline.lang.JobClass, MatrixCell> stationProc = PH.get(station);
                    if (stationProc != null) {
                        jline.lang.JobClass jobClass0 = sn.jobclasses.get(0);
                        MatrixCell proc0 = stationProc.get(jobClass0);
                        if (proc0 != null) {
                            A.set(0, proc0.get(0));
                            A.set(1, proc0.get(1));
                            A.set(2, proc0.get(1));
                            for (int k = 1; k < K; k++) {
                                jline.lang.JobClass jobClassK = sn.jobclasses.get(k);
                                MatrixCell procK = stationProc.get(jobClassK);
                                if (procK != null) {
                                    MatrixCell B = new MatrixCell();
                                    B.set(0, procK.get(0));
                                    B.set(1, procK.get(1));
                                    B.set(2, procK.get(1));
                                    A = Mmap_super.mmap_super(A, B);
                                }
                            }
                        }
                    }
                    idx_arv = i;
                } else if (schedI == SchedStrategy.FCFS || schedI == SchedStrategy.HOL
                        || schedI == SchedStrategy.FCFSPRIO || schedI == SchedStrategy.FCFSPRPRIO) {
                    Map<jline.lang.JobClass, MatrixCell> stationProc = PH.get(station);
                    if (stationProc != null) {
                        for (int k = 0; k < K; k++) {
                            jline.lang.JobClass jobClass = sn.jobclasses.get(k);
                            MatrixCell procK = stationProc.get(jobClass);
                            if (procK != null) {
                                double scaledMean = Map_mean.map_mean(procK.get(0), procK.get(1)) / sn.nservers.get(i);
                                MatrixCell scaledProc = Map_scale.map_scale(procK.get(0), procK.get(1), scaledMean);
                                pie.put(k, Map_pie.map_pie(scaledProc.get(0), scaledProc.get(1)));
                                S.put(k, scaledProc.get(0));
                            }
                        }
                    }
                    idx_q = i;
                    is_ps = false;
                } else if (schedI == SchedStrategy.PS) {
                    Map<jline.lang.JobClass, MatrixCell> stationProc = PH.get(station);
                    if (stationProc != null) {
                        for (int k = 0; k < K; k++) {
                            jline.lang.JobClass jobClass = sn.jobclasses.get(k);
                            MatrixCell procK = stationProc.get(jobClass);
                            if (procK != null) {
                                S.put(k, procK.get(0));
                                pie.put(k, Map_pie.map_pie(procK.get(0), procK.get(1)));
                            }
                        }
                    }
                    idx_q = i;
                    is_ps = true;
                } else {
                    throw new RuntimeException("Unsupported scheduling strategy");
                }
            }

            PriorityAnalysis priorityAnalysis = analyzePriorities(sn, "passage_time");
            SchedStrategy schedQ = sn.sched.get(sn.stations.get(idx_q));
            // Priorities select the law only under a priority DISCIPLINE: a
            // plain FCFS or PS queue serves in arrival or processor order
            // whatever the classprio column says
            boolean prioDiscipline = schedQ == SchedStrategy.HOL || schedQ == SchedStrategy.FCFSPRIO
                    || schedQ == SchedStrategy.FCFSPRPRIO;
            if (!priorityAnalysis.isSupported && !priorityAnalysis.isAllDistinct) {
                throw new RuntimeException(priorityAnalysis.message);
            } else if (priorityAnalysis.isAllDistinct && prioDiscipline) {
                // BUTools convention: D1=lowest .. DK=highest priority; LINE's
                // is lower value = higher priority, so the class order flips on
                // the way in and the outputs map back. Neither analyzer exports
                // the PH form, so the law is TABULATED: one solve for two
                // sojourn moments sizes a shared grid (mean + 10 sigma over the
                // classes), a second tabulates the CDF on it.
                final Matrix prios = sn.classprio;
                Integer[] order = new Integer[K];
                for (int k = 0; k < K; k++) order[k] = k;
                java.util.Arrays.sort(order, new java.util.Comparator<Integer>() {
                    @Override
                    public int compare(Integer a, Integer b) {
                        return Double.compare(prios.get(b), prios.get(a));
                    }
                });
                MatrixCell Db = new MatrixCell();
                Db.set(0, A.get(0));
                MatrixCell sigmaB = new MatrixCell();
                MatrixCell Sb = new MatrixCell();
                for (int b = 0; b < K; b++) {
                    Db.set(1 + b, A.get(2 + order[b]));
                    sigmaB.set(b, pie.get(order[b]));
                    Sb.set(b, S.get(order[b]));
                }
                boolean preemptive = schedQ == SchedStrategy.FCFSPRPRIO;
                Map<String, Map<Integer, Matrix>> moms = preemptive
                        ? MMAPPH1PRPR.MMAPPH1PRPR(Db, sigmaB, Sb, null, null, 2, null, null, null, null)
                        : MMAPPH1NPPR.MMAPPH1NPPR(Db, sigmaB, Sb, null, null, 2, null, null, null, null);
                double xmax = 0.0;
                for (int b = 0; b < K; b++) {
                    Matrix mk = moms.get("stMoms").get(b);
                    double m1 = mk.get(0);
                    double sig = Math.sqrt(Math.max(mk.get(1) - m1 * m1, 0.0));
                    xmax = Math.max(xmax, m1 + 10.0 * sig);
                }
                int n_pts = options.config.num_cdf_pts;
                // t = 0 is not evaluable (the Erlangization divides by t); a
                // sojourn time is strictly positive, so F(0) = 0 is prepended
                Matrix Xtail = new Matrix(1, n_pts - 1);
                for (int i = 1; i < n_pts; i++) Xtail.set(0, i - 1, xmax * i / (n_pts - 1.0));
                Map<String, Map<Integer, Matrix>> dist = preemptive
                        ? MMAPPH1PRPR.MMAPPH1PRPR(Db, sigmaB, Sb, null, null, null, Xtail, null, null, null)
                        : MMAPPH1NPPR.MMAPPH1NPPR(Db, sigmaB, Sb, null, null, null, Xtail, null, null, null);
                for (int b = 0; b < K; b++) {
                    int korig = order[b];
                    Matrix Fk = dist.get("stDistr").get(b);
                    Matrix F = new Matrix(n_pts, 1);
                    Matrix X = new Matrix(n_pts, 1);
                    for (int i = 1; i < n_pts; i++) {
                        F.set(i, 0, Fk.get(i - 1));
                        X.set(i, 0, xmax * i / (n_pts - 1.0));
                    }
                    if (!RD.containsKey(idx_arv)) RD.put(idx_arv, new MatrixCell());
                    RD.get(idx_arv).set(korig, new Matrix(0, 0));
                    if (!RD.containsKey(idx_q)) RD.put(idx_q, new MatrixCell());
                    RD.get(idx_q).set(korig, Matrix.concatColumns(F, X, null));
                }
            } else if (is_ps) {
                for (int k = 0; k < K; k++) {
                    if (S.get(k).getNumRows() != 1) {
                        InputOutput.line_error(InputOutput.mfilename(new Object()), "PS queue requires exponential (Markovian) service times");
                    }
                }
                if (K == 1) {
                    Matrix C_map = A.get(0);
                    Matrix D_map = A.get(1);
                    double mu = -S.get(0).get(0, 0);
                    double lambda = Map_lambda.map_lambda(C_map, D_map);
                    double rho = lambda / mu;
                    double approx_mean = 1.0 / (mu * (1.0 - rho));
                    int n_pts = options.config.num_cdf_pts;
                    double x_max = approx_mean * 10.0;
                    double[] x_vals = new double[n_pts];
                    for (int i = 0; i < n_pts; i++) x_vals[i] = (x_max * i) / (n_pts - 1.0);
                    double[] W_bar = MAPM1PSCdfRespT.computeCdf(C_map, D_map, mu, x_vals);
                    Matrix F = new Matrix(n_pts, 1);
                    Matrix X = new Matrix(n_pts, 1);
                    for (int i = 0; i < n_pts; i++) {
                        F.set(i, 0, 1.0 - W_bar[i]);
                        X.set(i, 0, x_vals[i]);
                    }
                    if (!RD.containsKey(idx_arv)) RD.put(idx_arv, new MatrixCell());
                    RD.get(idx_arv).set(0, new Matrix(0, 0));
                    if (!RD.containsKey(idx_q)) RD.put(idx_q, new MatrixCell());
                    RD.get(idx_q).set(0, Matrix.concatColumns(F, X, null));
                } else {
                    double[] mu_vec = new double[K];
                    for (int k = 0; k < K; k++) mu_vec[k] = -S.get(k).get(0, 0);
                    boolean allEqual = true;
                    for (int k = 1; k < K; k++) {
                        if (Math.abs(mu_vec[k] - mu_vec[0]) > GlobalConstants.FineTol) { allEqual = false; break; }
                    }
                    if (!allEqual) InputOutput.line_error(InputOutput.mfilename(new Object()), "Multi-class PS currently requires identical service rates");
                    double mu = mu_vec[0];
                    Matrix C_map = A.get(0);
                    // Sum the CLASS arrival matrices, A(2..): A(1) is already
                    // the aggregate the classes sum to, so starting the sum
                    // there counted every arrival twice and doubled the load
                    Matrix D_map_sum = A.get(2).copy();
                    for (int i = 3; i < A.size(); i++) D_map_sum = D_map_sum.add(1.0, A.get(i));
                    double lambda = Map_lambda.map_lambda(C_map, D_map_sum);
                    double rho = lambda / mu;
                    double approx_mean = 1.0 / (mu * (1.0 - rho));
                    int n_pts = options.config.num_cdf_pts;
                    double x_max = approx_mean * 10.0;
                    double[] x_vals = new double[n_pts];
                    for (int i = 0; i < n_pts; i++) x_vals[i] = (x_max * i) / (n_pts - 1.0);
                    double[] W_bar = MAPM1PSCdfRespT.computeCdf(C_map, D_map_sum, mu, x_vals);
                    Matrix F = new Matrix(n_pts, 1);
                    Matrix X = new Matrix(n_pts, 1);
                    for (int i = 0; i < n_pts; i++) {
                        F.set(i, 0, 1.0 - W_bar[i]);
                        X.set(i, 0, x_vals[i]);
                    }
                    for (int k = 0; k < K; k++) {
                        if (!RD.containsKey(idx_arv)) RD.put(idx_arv, new MatrixCell());
                        RD.get(idx_arv).set(k, new Matrix(0, 0));
                        if (!RD.containsKey(idx_q)) RD.put(idx_q, new MatrixCell());
                        RD.get(idx_q).set(k, Matrix.concatColumns(F, X, null));
                    }
                }
            } else {
                if (A.size() > 2) {
                    MatrixCell newA = new MatrixCell();
                    newA.set(0, A.get(0));
                    for (int i = 2; i < A.size(); i++) newA.set(i - 1, A.get(i));
                    A = newA;
                }
                // MMAPPH1FCFS.solve (butools) returns a raw Map; the assignment to
                // the parameterized type is inherently unchecked.
                @SuppressWarnings("unchecked")
                Map<String, Map<Integer, Matrix>> mmapResult = MMAPPH1FCFS.solve(A,
                        toMap(pie), toMap(S), null, null, null, null, false, true, null, null);
                Map<Integer, Matrix> alpha = mmapResult.get("stDistrPH_alpha");
                Map<Integer, Matrix> D0 = mmapResult.get("stDistrPH_A");
                for (int k = 0; k < K; k++) {
                    if (alpha.containsKey(k) && D0.containsKey(k)) {
                        Matrix alphaK = alpha.get(k);
                        Matrix D0K = D0.get(k);
                        Matrix negD0K = D0K.copy();
                        negD0K.scaleEq(-1.0);
                        Matrix D1K = negD0K.mult(Matrix.ones(alphaK.length(), 1)).mult(alphaK.transpose());
                        double variance = Map_var.map_var(D0K, D1K);
                        double meanResp = Map_mean.map_mean(D0K, D1K);
                        double sigma = Math.sqrt(variance);
                        int n = 5;
                        double maxTime = meanResp + n * sigma;
                        while (Map_cdf.map_cdf(D0K, D1K, Matrix.singleton(maxTime)).get(0) < 1 - GlobalConstants.FineTol) {
                            n++;
                            maxTime = meanResp + n * sigma;
                        }
                        int n_pts = options.config.num_cdf_pts;
                        Matrix F = new Matrix(n_pts, 1);
                        Matrix X = new Matrix(n_pts, 1);
                        for (int i = 0; i < n_pts; i++) {
                            double t = (maxTime * i) / (n_pts - 1.0);
                            X.set(i, 0, t);
                            F.set(i, 0, Map_cdf.map_cdf(D0K, D1K, Matrix.singleton(t)).get(0));
                        }
                        if (!RD.containsKey(idx_arv)) RD.put(idx_arv, new MatrixCell());
                        RD.get(idx_arv).set(k, new Matrix(0, 0));
                        if (!RD.containsKey(idx_q)) RD.put(idx_q, new MatrixCell());
                        RD.get(idx_q).set(k, Matrix.concatColumns(F, X, null));
                    }
                }
            }
        } else {
            InputOutput.line_warning(InputOutput.mfilename(new Object()),
                    "This model is not supported by SolverMAM yet. Returning with no result.");
        }
        return RD;
    }

    @SuppressWarnings({"unchecked", "rawtypes"})
    private static Map toMap(Map<?, ?> m) {
        return new HashMap(m);
    }

    @SuppressWarnings("unchecked")
    private static Map<Integer, Matrix> castMap(Object o) {
        Map<Integer, Matrix> result = new HashMap<Integer, Matrix>();
        if (o instanceof Map) {
            Map<?, ?> raw = (Map<?, ?>) o;
            for (Map.Entry<?, ?> e : raw.entrySet()) {
                if (e.getKey() instanceof Integer && e.getValue() instanceof Matrix) {
                    result.put((Integer) e.getKey(), (Matrix) e.getValue());
                }
            }
        }
        return result;
    }
}
