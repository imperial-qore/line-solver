package jline.examples.java.advanced;

import java.util.ArrayList;
import java.util.Arrays;
import java.util.List;

import jline.lang.ClosedClass;
import jline.lang.Network;
import jline.lang.RoutingMatrix;
import jline.lang.constant.SchedStrategy;
import jline.lang.nodes.Queue;
import jline.solvers.NetworkAvgTable;
import jline.solvers.SolverOptions;
import jline.solvers.ctmc.SolverCTMC;
import jline.solvers.ldes.SolverLDES;
import jline.util.matrix.Matrix;

/**
 * Validate the CTMC (exact) and LDES (simulation) solvers on a CLOSED CYCLIC
 * network of two pass-and-swap (PAS / order-independent) queues against the
 * exact product-form brute-force normalizing constant.
 *
 *   Topology:   --> Queue1 --> Queue2 -->   (cyclic, R closed classes)
 *
 * PAS queues are order-independent, so the network is product-form. The brute
 * force enumerates the whole ORDERED state space to obtain the normalizing
 * constant G, hence exact per-class throughputs and mean queue lengths, which
 * are compared against the solvers.
 *
 * Both stations must be valid OI queues: the TOTAL service rate mu_i(c) must
 * depend only on the customer multiset.
 *   station 1: M/M/k OI queue, class-independent  -> mu1(c) = min(n,k1)*s1
 *   station 2: infinite-server, class-dependent   -> mu2(c) = sum_j beta2[c_j]
 * For station 2 the prefix sums depend on order, so the ordered-state
 * enumeration is genuinely needed (it does not collapse to class counts).
 *
 * Class indices passed to mu(.) are 0-based (c(0) the oldest job).
 */
public class PassAndSwapCyclicExample {

    static final int[] K = {2, 2};        // per-class closed populations
    static final double S1 = 1.0;
    static final int K1 = 2;              // station 1: 2-server OI (M/M/2)
    static final double[] BETA2 = {1.5, 1.0}; // station 2: per-class inf-server rates
    static final int R = K.length;

    /** Total OI service rate as a function of an ordered 0-based prefix. */
    interface PrefixRate { double mu(int[] prefix, int len); }

    static final PrefixRate MU1 = (c, n) -> Math.min(n, K1) * S1;
    static final PrefixRate MU2 = (c, n) -> {
        double t = 0;
        for (int j = 0; j < n; j++) t += BETA2[c[j]];
        return t;
    };

    // ---- brute-force product form ------------------------------------
    static double oiDfs(int[] counts, int[] prefix, int depth, double w,
                        double[] e, PrefixRate rate) {
        boolean empty = true;
        for (int x : counts) if (x != 0) { empty = false; break; }
        if (empty) return w;
        double s = 0;
        for (int r = 0; r < counts.length; r++) {
            if (counts[r] > 0) {
                counts[r]--;
                prefix[depth] = r;
                s += oiDfs(counts, prefix, depth + 1,
                        w * e[r] / rate.mu(prefix, depth + 1), e, rate);
                counts[r]++;
            }
        }
        return s;
    }

    static double oiSum(int[] counts, double[] e, PrefixRate rate) {
        int n = 0;
        for (int x : counts) n += x;
        return oiDfs(counts.clone(), new int[n], 0, 1.0, e, rate);
    }

    static List<int[]> enumSplits(int[] kvec) {
        List<int[]> out = new ArrayList<>();
        int[] idx = new int[kvec.length];
        int total = 1;
        for (int k : kvec) total *= (k + 1);
        for (int s = 0; s < total; s++) {
            out.add(idx.clone());
            for (int r = 0; r < kvec.length; r++) {
                if (idx[r] < kvec[r]) { idx[r]++; break; }
                idx[r] = 0;
            }
        }
        return out;
    }

    static double normConst(int[] kvec, double[] e1, double[] e2,
                            PrefixRate r1, PrefixRate r2) {
        double G = 0;
        for (int[] m : enumSplits(kvec)) {
            int[] comp = new int[kvec.length];
            for (int r = 0; r < kvec.length; r++) comp[r] = kvec[r] - m[r];
            G += oiSum(m, e1, r1) * oiSum(comp, e2, r2);
        }
        return G;
    }

    /** Returns {G}, X[R], Q[2][R] packed as Object[]{Double, double[], double[][]}. */
    static Object[] bruteforce() {
        double[] e1 = new double[R];
        double[] e2 = new double[R];
        Arrays.fill(e1, 1.0);
        Arrays.fill(e2, 1.0);
        double G = normConst(K, e1, e2, MU1, MU2);
        double[] X = new double[R];
        for (int r = 0; r < R; r++) {
            int[] Km = K.clone();
            Km[r]--;
            X[r] = normConst(Km, e1, e2, MU1, MU2) / G;   // X_r = G(K-1_r)/G(K)
        }
        double[][] Q = new double[2][R];
        for (int[] m : enumSplits(K)) {
            int[] comp = new int[R];
            for (int r = 0; r < R; r++) comp[r] = K[r] - m[r];
            double w = oiSum(m, e1, MU1) * oiSum(comp, e2, MU2);
            for (int r = 0; r < R; r++) Q[0][r] += m[r] * w;
        }
        for (int r = 0; r < R; r++) {
            Q[0][r] /= G;
            Q[1][r] = K[r] - Q[0][r];
        }
        return new Object[]{G, X, Q};
    }

    // ---- LINE closed PAS network -------------------------------------
    static Network buildModel() {
        Network model = new Network("PAScyclic");
        Queue q1 = new Queue(model, "PASQueue1", SchedStrategy.PAS);
        Queue q2 = new Queue(model, "PASQueue2", SchedStrategy.PAS);
        ClosedClass[] cls = new ClosedClass[R];
        for (int r = 0; r < R; r++) {
            cls[r] = new ClosedClass(model, "Class" + (r + 1), K[r], q1);
        }
        q1.setService((Matrix c) -> Math.min(c.getNumCols(), K1) * S1);
        q2.setService((Matrix c) -> {
            double t = 0;
            for (int i = 0; i < c.getNumCols(); i++) t += BETA2[(int) c.get(0, i)];
            return t;
        });
        int total = 0;
        for (int k : K) total += k;
        q1.setSwapGraph(new Matrix(R, R));
        q1.setNumberOfServers(K1);
        q1.setCap(total);
        q2.setSwapGraph(new Matrix(R, R));
        q2.setNumberOfServers(total);
        q2.setCap(total);
        RoutingMatrix P = model.initRoutingMatrix();
        for (int r = 0; r < R; r++) {
            P.set(cls[r], cls[r], q1, q2, 1.0);
            P.set(cls[r], cls[r], q2, q1, 1.0);
        }
        model.link(P);
        return model;
    }

    /** Extract X[R] (chain throughput from q1) and Q[2][R] from an AvgTable. */
    static double[][] extractQ(NetworkAvgTable t) {
        List<Double> q = t.getQLen();
        double[][] Q = new double[2][R];
        for (int i = 0; i < 2; i++)
            for (int r = 0; r < R; r++)
                Q[i][r] = q.get(i * R + r);
        return Q;
    }

    static double[] extractX(NetworkAvgTable t) {
        List<Double> tp = t.getTput();
        double[] X = new double[R];
        for (int r = 0; r < R; r++) X[r] = tp.get(r);   // q1 rows
        return X;
    }

    /**
     * Brute-force normalizing constant of the closed cyclic PAS network
     * (pas_cyclic_normconst_bruteforce.m).
     *
     * <p>A pass-and-swap queue is an order-independent station: its total
     * service rate is a function mu_i(c) of the ordered class sequence c at the
     * station, and the stationary law is product-form and INVARIANT to the swap
     * graph, so the swap graph is not needed here. G is obtained by splitting
     * the population vector K into the per-class counts m at station 1 and K-m
     * at station 2, then summing the OI balance function over every ordered
     * arrangement of each station's count vector.</p>
     *
     * @return the normalizing constant G at the full population
     */
    public static double pas_cyclic_normconst_bruteforce() {
        double[] e1 = new double[R];
        double[] e2 = new double[R];
        Arrays.fill(e1, 1.0);
        Arrays.fill(e2, 1.0);
        double G = normConst(K, e1, e2, MU1, MU2);

        StringBuilder kstr = new StringBuilder();
        for (int r = 0; r < R; r++) {
            if (r > 0) {
                kstr.append(' ');
            }
            kstr.append(K[r]);
        }
        System.out.printf("Closed cyclic PAS network: R=%d classes, population K=[%s]%n",
                R, kstr);
        System.out.printf("Brute-force normalizing constant G = %.15g%n", G);

        // X_r = e_r * G(K - 1_r) / G(K); with e_r = 1 this is a ratio of
        // normalizing constants at populations K-1_r and K.
        StringBuilder xstr = new StringBuilder();
        for (int r = 0; r < R; r++) {
            double X = 0;
            if (K[r] > 0) {
                int[] Km = K.clone();
                Km[r]--;
                X = normConst(Km, e1, e2, MU1, MU2) / G;
            }
            if (r > 0) {
                xstr.append(' ');
            }
            xstr.append(String.format("%.6g", X));
        }
        System.out.printf("Per-class throughput X = [%s]%n", xstr);

        selftestSingleClass();
        return G;
    }

    /**
     * Single class, two single-server PAS stations reduces to a closed two-queue
     * Gordon-Newell network whose normalizing constant is the geometric sum
     * G = sum_n (1/mu1)^n (1/mu2)^(N-n). Verifies the enumeration reproduces it.
     */
    private static void selftestSingleClass() {
        final int N = 5;
        final double mu1 = 3.0;
        final double mu2 = 2.0;
        double[] e1 = {1.0};
        double[] e2 = {1.0};
        PrefixRate r1 = (c, n) -> pasRate(n, mu1, 1);
        PrefixRate r2 = (c, n) -> pasRate(n, mu2, 1);
        double Gbf = normConst(new int[]{N}, e1, e2, r1, r2);
        double Gcf = 0;
        for (int n = 0; n <= N; n++) {
            Gcf += Math.pow(1 / mu1, n) * Math.pow(1 / mu2, N - n);
        }
        if (Math.abs(Gbf - Gcf) > 1e-12 * Gcf) {
            throw new RuntimeException(String.format(
                    "self-test failed: brute-force %.15g vs closed-form %.15g", Gbf, Gcf));
        }
        System.out.printf("Self-test (single-class Gordon-Newell): PASS (G=%.12g)%n", Gbf);
    }

    /**
     * Total OI service rate of an ordered prefix of length n at a k-server PAS
     * station whose every class carries the same base rate.
     */
    private static double pasRate(int n, double beta, int k) {
        return Math.min(n, k) * beta;
    }

    /**
     * Validate CTMC and LDES against the exact brute-force product form
     * (pas_cyclic_ctmc_vs_bruteforce.m).
     */
    public static void pas_cyclic_ctmc_vs_bruteforce() {
        Object[] bf = bruteforce();
        double G = (Double) bf[0];
        double[] Xbf = (double[]) bf[1];
        double[][] Qbf = (double[][]) bf[2];
        System.out.printf("Brute-force normalizing constant G = %.12g%n%n", G);

        int total = 0;
        for (int k : K) total += k;

        // (3) CTMC -- exact
        SolverOptions opt = SolverCTMC.defaultOptions();
        opt.method = "exact";   // pin the state-space path
        opt.cutoff = Matrix.singleton(total);
        NetworkAvgTable tc = new SolverCTMC(buildModel(), opt).getAvgTable();
        tc.print();
        double[] Xc = extractX(tc);
        double[][] Qc = extractQ(tc);
        double errX = 0, errQ = 0;
        for (int r = 0; r < R; r++) errX = Math.max(errX, Math.abs(Xbf[r] - Xc[r]));
        for (int i = 0; i < 2; i++)
            for (int r = 0; r < R; r++) errQ = Math.max(errQ, Math.abs(Qbf[i][r] - Qc[i][r]));
        System.out.printf("%nCTMC:  max|dX| = %.3e   max|dQ| = %.3e%n", errX, errQ);
        if (errX > 1e-9 || errQ > 1e-9)
            throw new RuntimeException("CTMC vs brute-force mismatch");
        System.out.println("PASS: CTMC matches brute-force product form within 1.0e-9.");

        // (4) LDES -- simulation
        SolverOptions lopt = SolverLDES.defaultOptions();
        lopt.seed = 23000;
        lopt.samples = 200000;
        NetworkAvgTable tl = new SolverLDES(buildModel(), lopt).getAvgTable();
        double[] Xl = extractX(tl);
        double[][] Ql = extractQ(tl);
        double relX = 0, relQ = 0;
        for (int r = 0; r < R; r++) relX = Math.max(relX, Math.abs(Xbf[r] - Xl[r]) / Xbf[r]);
        for (int i = 0; i < 2; i++)
            for (int r = 0; r < R; r++) relQ = Math.max(relQ, Math.abs(Qbf[i][r] - Ql[i][r]) / Qbf[i][r]);
        System.out.printf("LDES:  max rel|dX| = %.2f%%   max rel|dQ| = %.2f%%%n",
                100 * relX, 100 * relQ);
        if (relX > 0.02 || relQ > 0.02)
            throw new RuntimeException("LDES vs brute-force exceeds 2%");
        System.out.println("PASS: LDES matches brute-force within 2% (simulation noise).");
    }

    public static void main(String[] args) {
        pas_cyclic_normconst_bruteforce();
        System.out.println();
        pas_cyclic_ctmc_vs_bruteforce();
    }
}
