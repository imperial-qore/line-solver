package jline.api.mdd;

import java.util.ArrayList;
import java.util.HashMap;
import java.util.List;
import java.util.Map;

/**
 * Kronecker rate descriptor for shared-server stations with phase-type service.
 *
 * <p>Under processor sharing every job at a station is in service at once, each
 * holding its own phase, so naming a single in-service phase (what
 * {@link Mdd_descriptor} does, which is non-preemptive semantics) cannot
 * represent the state. The local state here is instead the PER-PHASE COUNT
 * vector v = (v_1,...,v_h), v_a jobs in phase a, with n = sum(v) jobs present.
 * That is still a per-station quantity, so every event stays a product of
 * per-level terms and the Kronecker form of Eq. 1 survives.</p>
 *
 * <p>With one server shared by n jobs each job advances at rate 1/n, so from
 * local state v with n = sum(v):</p>
 * <pre>
 *   internal   v -&gt; v - e_a + e_b   at v_a * D0[a][b] / n     (a != b)
 *   departure  v -&gt; v - e_a         at v_a * t[a] / n * P[i][j]
 *   arrival    v -&gt; v + e_b         at pie[b]
 * </pre>
 * <p>An infinite-server (delay) station is the same without the 1/n scaling. For
 * h = 1 the departure rate collapses to n*mu/n = mu at PS and to n*mu at IS,
 * reproducing the usual single-server and delay rate laws.</p>
 *
 * <p>The local domain is the number of compositions of 0..N over h phases,
 * C(N+h,h), against 1+N*h for the non-preemptive encoding: the price of
 * tracking every job's phase rather than one.</p>
 */
public class Mdd_ps {

    private Mdd_ps() {}

    /**
     * Build the descriptor.
     *
     * @param mu station service rates, ignored where proc gives a law
     * @param P station-to-station routing matrix, row-stochastic
     * @param servers servers per station, 1 (PS) or infinite (IS); no other
     *        value has a per-phase-count encoding here
     * @param N closed population
     * @param proc per-station service law, or null for an exponential station
     */
    public static MddDescriptor mdd_ps(double[] mu, double[][] P, double[] servers, int N,
                                       MddServiceLaw[] proc) {
        int K = mu.length;
        double[][][] D0 = new double[K][][];
        double[][][] D1 = new double[K][][];
        double[][] entry = new double[K][];
        int[] h = new int[K];

        for (int i = 0; i < K; i++) {
            if (!(servers[i] == 1 || Double.isInfinite(servers[i]))) {
                throw new RuntimeException("mdd_ps: station " + (i + 1) + " has " + servers[i]
                        + " servers; only processor sharing (1) and infinite server have a "
                        + "per-phase-count encoding here");
            }
            MddServiceLaw law = (proc != null && i < proc.length) ? proc[i] : null;
            if (law == null) {
                D0[i] = new double[][]{{-mu[i]}};
                D1[i] = new double[][]{{mu[i]}};
                entry[i] = new double[]{1.0};
                h[i] = 1;
                continue;
            }
            D0[i] = law.D0;
            D1[i] = law.D1;
            h[i] = law.phases();
            entry[i] = Mdd_descriptor.entryLaw(law.pie, D1[i], h[i], i, "mdd_ps");
        }

        // ---- per-station composition state space
        int[][][] comp = new int[K][][];
        List<Map<CompKey, Integer>> lut = new ArrayList<Map<CompKey, Integer>>(K);
        int[] d = new int[K];
        for (int i = 0; i < K; i++) {
            comp[i] = compositions(h[i], N);
            Map<CompKey, Integer> table = new HashMap<CompKey, Integer>();
            for (int r = 0; r < comp[i].length; r++) {
                table.put(new CompKey(comp[i][r]), Integer.valueOf(r));
            }
            lut.add(table);
            d[i] = comp[i].length;
        }

        MddDescriptor desc = new MddDescriptor();
        desc.K = K;
        desc.N = N;
        desc.domain = d;
        desc.mu = mu;
        desc.servers = servers;
        desc.P = P;
        desc.nphases = h;
        desc.valuemap = new double[K][];
        for (int i = 0; i < K; i++) {
            double[] vm = new double[d[i]];
            for (int r = 0; r < d[i]; r++) {
                int s = 0;
                for (int a = 0; a < h[i]; a++) {
                    s += comp[i][r][a];
                }
                vm[r] = s;
            }
            desc.valuemap[i] = vm;
        }

        // ---- initial state: all jobs at station 1, entered in its entry phase
        int[] init = new int[K];
        for (int i = 0; i < K; i++) {
            int[] v = new int[h[i]];
            if (i == 0) {
                int first = 0;
                for (int a = 0; a < h[0]; a++) {
                    if (entry[0][a] > 0) {
                        first = a;
                        break;
                    }
                }
                v[first] = N;
            }
            init[i] = lut.get(i).get(new CompKey(v)).intValue();
        }
        desc.init = init;

        final double[][][] fD0 = D0;
        final double[][][] fD1 = D1;
        final double[][] fPie = entry;
        final int[] fh = h;
        final int[][][] fComp = comp;
        final List<Map<CompKey, Integer>> fLut = lut;
        final int fN = N;
        final double[][] fP = P;
        desc.nextfun = new MddNextState() {
            @Override
            public int[][] next(int[] state) {
                return successors(state, fD0, fD1, fPie, fh, fComp, fLut, fN, fP);
            }
        };

        // ---- events
        List<MddEvent> events = new ArrayList<MddEvent>();
        for (int i = 0; i < K; i++) {
            if (h[i] == 1) {
                continue;
            }
            MddLocalMatrix Wi = internal(D0[i], comp[i], lut.get(i), h[i], d[i], servers[i]);
            if (Wi.nnz == 0) {
                continue;
            }
            events.add(new MddEvent(i, i, new int[]{i}, new MddLocalMatrix[]{Wi}));
        }
        for (int a = 0; a < K; a++) {
            for (int b = 0; b < K; b++) {
                if (a == b || P[a][b] <= 0) {
                    continue;
                }
                MddLocalMatrix Wa = departure(D1[a], comp[a], lut.get(a), h[a], d[a],
                        servers[a], P[a][b]);
                MddLocalMatrix Wb = arrival(entry[b], comp[b], lut.get(b), h[b], N, d[b]);
                events.add(new MddEvent(a, b, new int[]{a, b},
                        new MddLocalMatrix[]{Wa, Wb}));
            }
        }
        desc.events = events;
        return desc;
    }

    // -----------------------------------------------------------------------
    /** Compositions of 0..N over h phases, in the reference's row order. */
    static int[][] compositions(int h, int N) {
        if (h == 1) {
            int[][] C = new int[N + 1][1];
            for (int n = 0; n <= N; n++) {
                C[n][0] = n;
            }
            return C;
        }
        int[][] sub = compositions(h - 1, N);
        List<int[]> out = new ArrayList<int[]>();
        for (int v1 = 0; v1 <= N; v1++) {
            for (int r = 0; r < sub.length; r++) {
                int s = 0;
                for (int a = 0; a < sub[r].length; a++) {
                    s += sub[r][a];
                }
                if (s > N - v1) {
                    continue;
                }
                int[] row = new int[h];
                row[0] = v1;
                System.arraycopy(sub[r], 0, row, 1, h - 1);
                out.add(row);
            }
        }
        return out.toArray(new int[out.size()][]);
    }

    /** Service-rate scaling: 1/n shared by n jobs at PS, unscaled at IS. */
    private static double share(int n, double srv) {
        if (n == 0) {
            return 0.0;
        }
        if (Double.isInfinite(srv)) {
            return 1.0;
        }
        return 1.0 / n;
    }

    private static MddLocalMatrix internal(double[][] D0i, int[][] C,
                                           Map<CompKey, Integer> lut, int h, int d,
                                           double srv) {
        MddLocalMatrix.Builder bld = new MddLocalMatrix.Builder(d);
        for (int r = 0; r < d; r++) {
            int[] v = C[r];
            int n = 0;
            for (int a = 0; a < h; a++) {
                n += v[a];
            }
            if (n == 0) {
                continue;
            }
            double sc = share(n, srv);
            for (int a = 0; a < h; a++) {
                if (v[a] == 0) {
                    continue;
                }
                for (int b = 0; b < h; b++) {
                    if (a == b || D0i[a][b] == 0) {
                        continue;
                    }
                    int[] w = v.clone();
                    w[a]--;
                    w[b]++;
                    bld.add(r, lut.get(new CompKey(w)).intValue(), v[a] * D0i[a][b] * sc);
                }
            }
        }
        return bld.build();
    }

    private static MddLocalMatrix departure(double[][] D1i, int[][] C,
                                            Map<CompKey, Integer> lut, int h, int d,
                                            double srv, double pr) {
        double[] t = new double[h];
        for (int a = 0; a < h; a++) {
            double s = 0;
            for (int b = 0; b < h; b++) {
                s += D1i[a][b];
            }
            t[a] = s;
        }
        MddLocalMatrix.Builder bld = new MddLocalMatrix.Builder(d);
        for (int r = 0; r < d; r++) {
            int[] v = C[r];
            int n = 0;
            for (int a = 0; a < h; a++) {
                n += v[a];
            }
            if (n == 0) {
                continue;
            }
            double sc = share(n, srv);
            for (int a = 0; a < h; a++) {
                if (v[a] == 0 || t[a] == 0) {
                    continue;
                }
                int[] w = v.clone();
                w[a]--;
                bld.add(r, lut.get(new CompKey(w)).intValue(), v[a] * t[a] * sc * pr);
            }
        }
        return bld.build();
    }

    private static MddLocalMatrix arrival(double[] pieb, int[][] C,
                                          Map<CompKey, Integer> lut, int h, int N, int d) {
        MddLocalMatrix.Builder bld = new MddLocalMatrix.Builder(d);
        for (int r = 0; r < d; r++) {
            int[] v = C[r];
            int n = 0;
            for (int a = 0; a < h; a++) {
                n += v[a];
            }
            if (n >= N) {
                continue;
            }
            for (int b = 0; b < h; b++) {
                if (pieb[b] == 0) {
                    continue;
                }
                int[] w = v.clone();
                w[b]++;
                bld.add(r, lut.get(new CompKey(w)).intValue(), pieb[b]);
            }
        }
        return bld.build();
    }

    private static int[][] successors(int[] s, double[][][] D0, double[][][] D1,
                                      double[][] pie, int[] h, int[][][] comp,
                                      List<Map<CompKey, Integer>> lut, int N, double[][] P) {
        int K = s.length;
        List<int[]> T = new ArrayList<int[]>();
        for (int i = 0; i < K; i++) {
            int[] v = comp[i][s[i]];
            int n = 0;
            for (int a = 0; a < h[i]; a++) {
                n += v[a];
            }
            if (n == 0) {
                continue;
            }
            // internal phase moves
            for (int a = 0; a < h[i]; a++) {
                if (v[a] == 0) {
                    continue;
                }
                for (int b = 0; b < h[i]; b++) {
                    if (a == b || D0[i][a][b] == 0) {
                        continue;
                    }
                    int[] w = v.clone();
                    w[a]--;
                    w[b]++;
                    int[] t = s.clone();
                    t[i] = lut.get(i).get(new CompKey(w)).intValue();
                    T.add(t);
                }
            }
            // completions routed to j
            for (int a = 0; a < h[i]; a++) {
                if (v[a] == 0) {
                    continue;
                }
                double ta = 0;
                for (int b = 0; b < h[i]; b++) {
                    ta += D1[i][a][b];
                }
                if (ta == 0) {
                    continue;
                }
                int[] w = v.clone();
                w[a]--;
                for (int j = 0; j < K; j++) {
                    if (j == i || P[i][j] <= 0) {
                        continue;
                    }
                    int[] vj = comp[j][s[j]];
                    int nj = 0;
                    for (int b = 0; b < h[j]; b++) {
                        nj += vj[b];
                    }
                    if (nj >= N) {
                        continue;
                    }
                    for (int b = 0; b < h[j]; b++) {
                        if (pie[j][b] == 0) {
                            continue;
                        }
                        int[] wj = vj.clone();
                        wj[b]++;
                        int[] t = s.clone();
                        t[i] = lut.get(i).get(new CompKey(w)).intValue();
                        t[j] = lut.get(j).get(new CompKey(wj)).intValue();
                        T.add(t);
                    }
                }
            }
        }
        return T.toArray(new int[T.size()][]);
    }

    /** Hashable wrapper of a per-phase count vector. */
    static final class CompKey {
        private final int[] v;
        private final int hash;

        CompKey(int[] v) {
            this.v = v;
            this.hash = java.util.Arrays.hashCode(v);
        }

        @Override
        public boolean equals(Object o) {
            if (this == o) {
                return true;
            }
            if (!(o instanceof CompKey)) {
                return false;
            }
            return java.util.Arrays.equals(this.v, ((CompKey) o).v);
        }

        @Override
        public int hashCode() {
            return hash;
        }
    }
}
