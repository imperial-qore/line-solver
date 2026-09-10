package jline.solvers.mva;

import jline.api.pfqn.mva.Pfqn_mvaoi;
import jline.lang.NetworkStruct;
import jline.lang.constant.SchedStrategy;
import jline.lang.nodeparam.QueueNodeParam;
import jline.lang.nodes.Station;
import jline.solvers.SolverOptions;
import jline.util.SerializableFunction;
import jline.util.matrix.Matrix;

import java.util.ArrayList;
import java.util.List;
import java.util.Map;
import java.util.function.ToDoubleFunction;

/**
 * Exact mean-value MVA for order-independent (OI) queueing networks.
 *
 * An OI station is a class-dependent load-dependent server whose total service
 * rate mu(n) is a permutation-invariant function of the per-class count vector n.
 * A closed network of infinite-server (delay) and load-independent (single-server,
 * product-form) stations plus ANY number of OI stations is product-form. This
 * analyzer aggregates the delay stations into a single think-time vector Z,
 * collects the load-independent (LI) queue demands and the OI-station rate
 * handles, and calls {@link Pfqn_mvaoi}, the mean-value Conditional-MVA (CMVA)
 * that carries one rate-shift vector per OI station and returns exact per-class
 * throughput and queue-lengths WITHOUT any normalizing constant or joint
 * marginal. The marginal-distribution counterpart is
 * {@link jline.api.pfqn.mva.Pfqn_mvaoi_marg}.
 *
 * Reference: Reiser, Lavenberg (1980), JACM 27(2); load-dependent extension Bruell,
 * Balbo, Afshari (1984); OI stations / CMVA Casale (2009); Casale, Comte, Dorsman (2026).
 */
public class SolverMVAOIAnalyzer {

    private final NetworkStruct sn;
    private final SolverOptions options;
    private final int R;
    private final int M;
    private final int[] njobs;
    private final double[][] demand;   // D_ir
    private final double[][] visits;   // V_ir
    private final boolean[] isDelay;
    private final boolean[] isOI;
    private final List<Integer> oiList;

    public SolverMVAOIAnalyzer(NetworkStruct sn, SolverOptions options) {
        this.sn = sn;
        this.options = options;
        this.R = sn.nclasses;
        this.M = sn.nstations;
        this.njobs = new int[R];
        for (int r = 0; r < R; r++) {
            this.njobs[r] = (int) Math.round(sn.njobs.get(0, r));
        }
        this.oiList = findOIStations();
        if (this.oiList.isEmpty()) {
            throw new RuntimeException("OI solver requires at least one order-independent station");
        }
        // ---- reject class switching (OI rank rates are per raw class) ------
        // The recursion is driven by the per-class population vector sn.njobs,
        // which class switching makes meaningless: a class that only ever
        // appears mid-chain carries njobs = 0, so the OI station would be
        // analyzed as if empty. Refuse it the way Solver_nc_oi does rather than
        // return that silently.
        for (int c = 0; c < sn.nchains; c++) {
            if (sn.inchain.get(c).getNumElements() > 1) {
                throw new RuntimeException("solver_mva_oi requires one class per chain (no class switching).");
            }
        }
        this.visits = new double[M][R];
        for (Map.Entry<Integer, Matrix> e : sn.visits.entrySet()) {
            Matrix vc = e.getValue();
            for (int i = 0; i < M; i++) {
                for (int r = 0; r < R; r++) {
                    this.visits[i][r] += vc.get(i, r);
                }
            }
        }
        this.demand = new double[M][R];
        for (int i = 0; i < M; i++) {
            for (int r = 0; r < R; r++) {
                double rate = sn.rates.get(i, r);
                if (!Double.isInfinite(rate) && !Double.isNaN(rate) && rate > 0) {
                    this.demand[i][r] = this.visits[i][r] / rate;
                }
            }
        }
        this.isDelay = new boolean[M];
        this.isOI = new boolean[M];
        for (int i = 0; i < M; i++) {
            this.isDelay[i] = sn.sched.get(sn.stations.get(i)) == SchedStrategy.INF;
        }
        for (int o : this.oiList) this.isOI[o] = true;
    }

    private List<Integer> findOIStations() {
        List<Integer> out = new ArrayList<Integer>();
        for (int i = 0; i < M; i++) {
            Station st = sn.stations.get(i);
            SchedStrategy s = sn.sched.get(st);
            if (s != SchedStrategy.PAS && s != SchedStrategy.OI) {
                continue;
            }
            Object p = sn.nodeparam.get(st);
            if (!(p instanceof QueueNodeParam)) {
                continue;
            }
            QueueNodeParam qp = (QueueNodeParam) p;
            if (qp.svcRateFun == null || qp.swapGraph == null) {
                continue;
            }
            if (qp.swapGraph.isEmpty() || qp.swapGraph.elementMaxAbs() == 0) {
                out.add(i);
            }
        }
        return out;
    }

    /** OI total service rate at count vector n via the (0-based) canonical microstate. */
    private static double muM(SerializableFunction<Matrix, Double> svcRateFun, int[] n, int R) {
        int tot = 0;
        for (int v : n) tot += v;
        if (tot == 0) {
            return 0.0;
        }
        Matrix micro = new Matrix(1, tot);
        int pos = 0;
        for (int r = 0; r < R; r++) {
            for (int c = 0; c < n[r]; c++) {
                micro.set(0, pos++, r);
            }
        }
        return svcRateFun.apply(micro);
    }

    /**
     * OI rate function reproducing the c-server BCMP station with per-class demands
     * Dq: mu(n) = (min(|n|,c)/|n|) * sum_{r: n_r&gt;0} n_r/Dq_r. For c = 1 this is the
     * familiar total completion rate of a multiclass single-server queue, and for a
     * single class it reduces to min(n,c)/Dq (M/M/c).
     */
    private static double msOIRate(int[] n, double[] Dq, double c) {
        int tot = 0;
        for (int v : n) tot += v;
        if (tot == 0) {
            return 0.0;
        }
        double acc = 0.0;
        for (int r = 0; r < n.length; r++) {
            if (n[r] > 0 && Dq[r] > 0) {
                acc += n[r] / Dq[r];
            }
        }
        return (Math.min(tot, c) / tot) * acc;
    }

    public AnalysisResults analyze() {
        long tstart = System.currentTimeMillis();

        // see _kb/06-solver-catalog.md for rationale
        double[] Z = new double[R];
        List<Integer> liList = new ArrayList<Integer>();
        List<Integer> msList = new ArrayList<Integer>();
        for (int i = 0; i < M; i++) {
            if (isOI[i]) {
                continue;
            } else if (isDelay[i]) {
                for (int r = 0; r < R; r++) Z[r] += demand[i][r];
            } else {
                double sv = sn.nservers.get(i, 0);
                if (!Double.isInfinite(sv) && !Double.isNaN(sv) && sv > 1) {
                    msList.add(i);
                } else {
                    liList.add(i);
                }
            }
        }
        double[][] Dli = new double[liList.size()][R];
        for (int j = 0; j < liList.size(); j++) {
            Dli[j] = demand[liList.get(j)].clone();
        }

        List<ToDoubleFunction<int[]>> muCell = new ArrayList<ToDoubleFunction<int[]>>();
        for (int o : oiList) {
            Station oiStation = sn.stations.get(o);
            QueueNodeParam qp = (QueueNodeParam) sn.nodeparam.get(oiStation);
            final SerializableFunction<Matrix, Double> f = qp.svcRateFun;
            muCell.add(new ToDoubleFunction<int[]>() {
                public double applyAsDouble(int[] n) {
                    return muM(f, n, R);
                }
            });
        }
        for (int j = 0; j < msList.size(); j++) {
            final double[] Dq = demand[msList.get(j)].clone();
            double svTmp = sn.nservers.get(msList.get(j), 0);
            if (Double.isInfinite(svTmp) || Double.isNaN(svTmp) || svTmp <= 0) {
                svTmp = 1;
            }
            final double cSv = svTmp;
            muCell.add(new ToDoubleFunction<int[]>() {
                public double applyAsDouble(int[] n) {
                    return msOIRate(n, Dq, cSv);
                }
            });
        }

        // Per-muCell-station visit vectors: genuine OI stations carry their
        // class visits (rate handle has none); ms-promoted stations pass unit
        // visits (already folded into demand by msOIRate).
        double[][] oivis = new double[oiList.size() + msList.size()][R];
        for (int o = 0; o < oiList.size(); o++) {
            for (int r = 0; r < R; r++) oivis[o][r] = visits[oiList.get(o)][r];
        }
        for (int j = 0; j < msList.size(); j++) {
            for (int r = 0; r < R; r++) oivis[oiList.size() + j][r] = 1.0;
        }
        Pfqn_mvaoi.Result res = Pfqn_mvaoi.pfqn_mvaoi(Z, njobs, muCell, Dli, oivis);

        double[] XN = res.X;

        // Row of res.Soi/res.Qoi holding each OI station (muCell order: oiList, msList).
        int[] oiRow = new int[M];
        for (int o = 0; o < oiList.size(); o++) oiRow[oiList.get(o)] = o;
        double[][] QN = new double[M][R];
        for (int o = 0; o < oiList.size(); o++) {
            QN[oiList.get(o)] = res.Qoi[o].clone();
        }
        for (int j = 0; j < msList.size(); j++) {
            QN[msList.get(j)] = res.Qoi[oiList.size() + j].clone();
        }
        for (int j = 0; j < liList.size(); j++) {
            QN[liList.get(j)] = res.Qli[j].clone();
        }
        for (int i = 0; i < M; i++) {
            if (isDelay[i]) {
                for (int r = 0; r < R; r++) QN[i][r] = XN[r] * demand[i][r];
            }
        }

        double[][] TN = new double[M][R];
        double[][] RN = new double[M][R];
        double[][] UN = new double[M][R];
        double[][] CN = new double[M][R];
        for (int i = 0; i < M; i++) {
            for (int r = 0; r < R; r++) {
                TN[i][r] = XN[r] * visits[i][r];
                if (XN[r] > 0) {
                    RN[i][r] = QN[i][r] / XN[r];
                }
                if (isOI[i]) {
                    // see _kb/06-solver-catalog.md for rationale
                    double sv = sn.nservers.get(i, 0);
                    if (Double.isInfinite(sv) || sv <= 0) sv = 1;
                    UN[i][r] = res.getSoi()[oiRow[i]][r] / sv;
                } else if (isDelay[i]) {
                    UN[i][r] = QN[i][r];
                } else {
                    double sv = sn.nservers.get(i, 0);
                    if (Double.isInfinite(sv) || sv <= 0) sv = 1;
                    UN[i][r] = XN[r] * demand[i][r] / sv;
                }
                CN[i][r] = RN[i][r];
            }
        }

        int Ntot = 0;
        for (int v : njobs) Ntot += v;
        long runtime = System.currentTimeMillis() - tstart;
        return new AnalysisResults(QN, UN, RN, TN, CN, XN, 0.0, runtime, Ntot, "oi");
    }

    public static class AnalysisResults {
        public double[][] QN, UN, RN, TN, CN;
        public double[] XN;
        public double lG;
        public long runtime;
        public int iter;
        public String method;

        public AnalysisResults(double[][] QN, double[][] UN, double[][] RN,
                double[][] TN, double[][] CN, double[] XN,
                double lG, long runtime, int iter, String method) {
            this.QN = QN;
            this.UN = UN;
            this.RN = RN;
            this.TN = TN;
            this.CN = CN;
            this.XN = XN;
            this.lG = lG;
            this.runtime = runtime;
            this.iter = iter;
            this.method = method;
        }
    }
}
