package jline.api.mdd;

import java.util.ArrayList;
import java.util.List;

/**
 * Kronecker rate descriptor of a single-class closed queueing network.
 *
 * <p>For the Miner-Ciardo-Donatelli approximate-aggregation solver
 * ({@link Mdd_mcd}), after A.S. Miner, G. Ciardo, S. Donatelli, "Using the
 * exact state space of a Markov model to compute approximate stationary
 * measures", SIGMETRICS 2000.</p>
 *
 * <p>The transition rate matrix is expressed compositionally as
 * R = sum_e (kron_k W_k^e) restricted to the reachable set, with
 * W_k^e[i,j] = lambda_k^e[i] * Prob_k^e(i,j) (Eq. 1). Each level k is a station
 * and each event e is a completion at station a routed to b.</p>
 *
 * <p><b>Exponential stations.</b> The local state is the population alone:
 * W_a^e[i,i-1] = mu[a]*min(i,servers[a])*P[a][b] for i &gt;= 1 (departure),
 * W_b^e[i,i+1] = 1 for i &lt;= N-1 (arrival), and the identity elsewhere.</p>
 *
 * <p><b>Phase-type stations.</b> The local state is the PAIR (population, phase
 * of the job in service), encoded in one level rather than two. Splitting them
 * does not work: on a completion routed into station b the phase at b restarts
 * only when b was empty, a joint condition on b's two components, which is not a
 * product of per-level terms. Merging them keeps every event local:</p>
 * <pre>
 *   index 0                   : station empty
 *   index 1 + (n-1)*h + (a-1) : n jobs present, job in service in phase a
 *   domain                    = 1 + N*h    (h = 1 reproduces index = n)
 * </pre>
 * <p>With exit vector t = D1*1 and entry law pie, departure (n,a)-&gt;(n-1,b) at
 * t[a]*P[a][b]*pie[b] for n &gt;= 2 and (1,a)-&gt;0 at t[a]*P[a][b]; arrival
 * 0-&gt;(1,b) at pie[b] and (m,c)-&gt;(m+1,c) at 1 for m &gt;= 1; internal
 * (n,a)-&gt;(n,b) at D0[a][b] for n &gt;= 1, a != b.</p>
 *
 * <p><b>Restrictions.</b> A phase-type station must be single-server: with
 * c &gt; 1 or an infinite server the local state would have to count jobs per
 * phase rather than name one phase, a different and much larger encoding. It
 * must also be NON-preemptive, because the composite level names the phase of
 * the one job in service and restarts it at pie when the next job starts; under
 * preemptive resume an arrival suspends that job and its phase has to be
 * remembered, so the local state would need a stack of phases. That matters for
 * LCFSPR, which is BCMP type 2 and stays product-form under general service:
 * the insensitivity is real but is NOT reachable through this encoding.
 * Exponential service is unaffected, preemption being immaterial by
 * memorylessness. Pass the disciplines to have the case rejected rather than
 * silently modelled as non-preemptive.</p>
 */
public class Mdd_descriptor {

    private Mdd_descriptor() {}

    /** Disciplines whose preemptive-resume semantics the composite encoding cannot carry. */
    private static final String[] PREEMPTIVE_RESUME =
            {"LCFSPR", "FCFSPR", "LCFSPRPRIO", "FCFSPRPRIO"};
    /** Shared-server disciplines, where every job present holds its own phase. */
    private static final String[] SHARED_SERVER = {"PS", "DPS", "GPS"};

    /**
     * Build the descriptor.
     *
     * @param mu station service rates, 1/E[S]; entry i is ignored when station i
     *        is given a phase-type law through proc
     * @param P station-to-station routing matrix, row-stochastic
     * @param servers servers per station, Double.POSITIVE_INFINITY for delay/IS
     * @param N closed population
     * @param proc per-station service law, or null for an exponential station;
     *        the array itself may be null when every station is exponential
     * @param sched per-station discipline names, consulted only to REJECT a
     *        phase-type law at a preemptive-resume or shared-server station; may
     *        be null when every station is non-preemptive
     */
    public static MddDescriptor mdd_descriptor(double[] mu, double[][] P, double[] servers,
                                               int N, MddServiceLaw[] proc, String[] sched) {
        int K = mu.length;
        double[][][] D0 = new double[K][][];
        double[][][] D1 = new double[K][][];
        double[][] entry = new double[K][];
        int[] h = new int[K];

        for (int i = 0; i < K; i++) {
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
            entry[i] = entryLaw(law.pie, D1[i], h[i], i, "mdd_descriptor");
            if (h[i] > 1 && servers[i] != 1) {
                throw new RuntimeException("mdd_descriptor: station " + (i + 1)
                        + " has a phase-type service law and " + servers[i] + " servers; a "
                        + "multi-server or delay station would have to count jobs per phase "
                        + "rather than name the phase of one job in service, which this "
                        + "encoding does not carry");
            }
            if (h[i] > 1 && sched != null && i < sched.length && sched[i] != null) {
                String nm = sched[i].toUpperCase();
                if (contains(PREEMPTIVE_RESUME, nm)) {
                    throw new RuntimeException("mdd_descriptor: station " + (i + 1)
                            + " combines a phase-type service law with a preemptive-resume "
                            + "discipline; the suspended jobs' phases would have to be stacked "
                            + "in the local state, which this encoding does not carry, and the "
                            + "descriptor would silently model the non-preemptive chain instead");
                }
                if (contains(SHARED_SERVER, nm)) {
                    throw new RuntimeException("mdd_descriptor: station " + (i + 1)
                            + " combines a phase-type service law with a shared-server "
                            + "discipline; every job present is in service and holds its own "
                            + "phase, which this encoding does not carry. Use Mdd_ps, whose "
                            + "local state is the per-phase count vector.");
                }
            }
        }

        int[] d = new int[K];
        for (int i = 0; i < K; i++) {
            d[i] = 1 + N * h[i];
        }

        MddDescriptor desc = new MddDescriptor();
        desc.K = K;
        desc.N = N;
        desc.domain = d;
        desc.mu = mu;
        desc.servers = servers;
        desc.P = P;
        desc.nphases = h;

        // ---- index maps
        desc.valuemap = new double[K][];
        for (int i = 0; i < K; i++) {
            double[] vm = new double[d[i]];
            for (int n = 1; n <= N; n++) {
                for (int a = 0; a < h[i]; a++) {
                    vm[1 + (n - 1) * h[i] + a] = n;
                }
            }
            desc.valuemap[i] = vm;
        }

        // ---- initial state: all jobs at station 1, in its entry phase
        int[] init = new int[K];
        int firstPhase = 0;
        for (int a = 0; a < h[0]; a++) {
            if (entry[0][a] > 0) {
                firstPhase = a;
                break;
            }
        }
        init[0] = idx(N, firstPhase + 1, h[0]);
        desc.init = init;

        final double[][][] fD0 = D0;
        final double[][][] fD1 = D1;
        final double[][] fPie = entry;
        final int[] fh = h;
        final int fN = N;
        final double[][] fP = P;
        desc.nextfun = new MddNextState() {
            @Override
            public int[][] next(int[] state) {
                return successors(state, fD0, fD1, fPie, fh, fN, fP);
            }
        };

        // ---- events
        List<MddEvent> events = new ArrayList<MddEvent>();
        // internal phase changes, one event per phase-type station
        for (int i = 0; i < K; i++) {
            if (h[i] == 1) {
                continue;
            }
            MddLocalMatrix Wi = internal(D0[i], h[i], N, d[i]);
            if (Wi.nnz == 0) {
                continue;
            }
            events.add(new MddEvent(i, i, new int[]{i}, new MddLocalMatrix[]{Wi}));
        }
        // completions routed a -> b
        for (int a = 0; a < K; a++) {
            for (int b = 0; b < K; b++) {
                if (a == b || P[a][b] <= 0) {
                    continue;
                }
                MddLocalMatrix Wa = departure(D1[a], entry[a], h[a], N, d[a], mu[a],
                        servers[a], P[a][b]);
                MddLocalMatrix Wb = arrival(entry[b], h[b], N, d[b]);
                events.add(new MddEvent(a, b, new int[]{a, b},
                        new MddLocalMatrix[]{Wa, Wb}));
            }
        }
        desc.events = events;
        return desc;
    }

    // -----------------------------------------------------------------------
    private static boolean contains(String[] names, String nm) {
        for (int i = 0; i < names.length; i++) {
            if (names[i].equals(nm)) {
                return true;
            }
        }
        return false;
    }

    /**
     * Entry law of a phase-type station, taken as given or derived from D1.
     *
     * <p>A {D0,D1} pair carries its own restart law: for a renewal process
     * D1 = t0 * pie, so every row with a positive exit rate is proportional to
     * pie. Deriving it is not optional -- defaulting to e_1 instead silently
     * replaces a hyperexponential (whose D0 is diagonal, so a job entering phase
     * 1 can never leave it) by an exponential at the phase-1 rate.</p>
     */
    static double[] entryLaw(double[] given, double[][] D1, int h, int i, String caller) {
        if (given != null && given.length > 0) {
            return normalize(given.clone());
        }
        double[] t0 = new double[h];
        int live = -1;
        int nlive = 0;
        for (int a = 0; a < h; a++) {
            double s = 0;
            for (int b = 0; b < h; b++) {
                s += D1[a][b];
            }
            t0[a] = s;
            if (s > 0) {
                nlive++;
                if (live < 0) {
                    live = a;
                }
            }
        }
        if (live < 0) {
            double[] pie = new double[h];
            pie[0] = 1.0;
            return pie;
        }
        double[] first = new double[h];
        for (int b = 0; b < h; b++) {
            first[b] = D1[live][b] / t0[live];
        }
        // A non-renewal MAP restarts in a law that depends on the phase it left
        // from, which one entry vector cannot express; the composite level would
        // silently model the renewal process instead.
        if (nlive > 1) {
            for (int a = live + 1; a < h; a++) {
                if (t0[a] <= 0) {
                    continue;
                }
                for (int b = 0; b < h; b++) {
                    if (Math.abs(D1[a][b] / t0[a] - first[b]) > 1e-9) {
                        throw new RuntimeException(caller + ": station " + (i + 1)
                                + " carries a service law whose restart distribution depends on "
                                + "the completing phase (a non-renewal MAP); the local state "
                                + "names one entry law, so this encoding cannot represent it");
                    }
                }
            }
        }
        return first;
    }

    static double[] normalize(double[] v) {
        double s = 0;
        for (int i = 0; i < v.length; i++) {
            s += v[i];
        }
        for (int i = 0; i < v.length; i++) {
            v[i] /= s;
        }
        return v;
    }

    /** Local index of (population n, service phase a); 0 when the station is empty. */
    private static int idx(int n, int a, int h) {
        if (n == 0) {
            return 0;
        }
        return 1 + (n - 1) * h + (a - 1);
    }

    /** Population of a local index; 0 when empty. */
    private static int population(int index, int h) {
        if (index == 0) {
            return 0;
        }
        return (index - 1) / h + 1;
    }

    /** Service phase (1-based) of a local index; 0 when empty. */
    private static int phase(int index, int h) {
        if (index == 0) {
            return 0;
        }
        return (index - 1) % h + 1;
    }

    /** Phase changes that do not complete a service, at any population n &gt;= 1. */
    private static MddLocalMatrix internal(double[][] D0i, int h, int N, int d) {
        MddLocalMatrix.Builder bld = new MddLocalMatrix.Builder(d);
        for (int n = 1; n <= N; n++) {
            for (int a = 1; a <= h; a++) {
                for (int b = 1; b <= h; b++) {
                    if (a == b || D0i[a - 1][b - 1] == 0) {
                        continue;
                    }
                    bld.add(idx(n, a, h), idx(n, b, h), D0i[a - 1][b - 1]);
                }
            }
        }
        return bld.build();
    }

    /** Completion at this station, routed out with probability pr. */
    private static MddLocalMatrix departure(double[][] D1i, double[] piei, int h, int N,
                                            int d, double mui, double srv, double pr) {
        MddLocalMatrix.Builder bld = new MddLocalMatrix.Builder(d);
        if (h == 1) {
            // exponential: the multi-server and delay rate laws live here
            for (int n = 1; n <= N; n++) {
                bld.add(n, n - 1, mui * Math.min(n, srv) * pr);
            }
        } else {
            double[] t = new double[h];
            for (int a = 0; a < h; a++) {
                double s = 0;
                for (int b = 0; b < h; b++) {
                    s += D1i[a][b];
                }
                t[a] = s;
            }
            for (int n = 1; n <= N; n++) {
                for (int a = 1; a <= h; a++) {
                    if (t[a - 1] == 0) {
                        continue;
                    }
                    if (n == 1) {
                        bld.add(idx(1, a, h), 0, t[a - 1] * pr);
                    } else {
                        for (int b = 1; b <= h; b++) {
                            if (piei[b - 1] == 0) {
                                continue;
                            }
                            bld.add(idx(n, a, h), idx(n - 1, b, h),
                                    t[a - 1] * pr * piei[b - 1]);
                        }
                    }
                }
            }
        }
        return bld.build();
    }

    /** An arrival starts service only when the station was empty. */
    private static MddLocalMatrix arrival(double[] pieb, int h, int N, int d) {
        MddLocalMatrix.Builder bld = new MddLocalMatrix.Builder(d);
        if (h == 1) {
            for (int m = 0; m < N; m++) {
                bld.add(m, m + 1, 1.0);
            }
        } else {
            for (int b = 1; b <= h; b++) {
                if (pieb[b - 1] == 0) {
                    continue;
                }
                bld.add(0, idx(1, b, h), pieb[b - 1]);
            }
            for (int m = 1; m < N; m++) {
                for (int c = 1; c <= h; c++) {
                    bld.add(idx(m, c, h), idx(m + 1, c, h), 1.0);
                }
            }
        }
        return bld.build();
    }

    /** Successor local-index vectors of a state, used to generate the reachable set. */
    private static int[][] successors(int[] s, double[][][] D0, double[][][] D1,
                                      double[][] pie, int[] h, int N, double[][] P) {
        int K = s.length;
        List<int[]> T = new ArrayList<int[]>();
        for (int i = 0; i < K; i++) {
            int ni = population(s[i], h[i]);
            int ai = phase(s[i], h[i]);
            if (ni == 0) {
                continue;
            }
            // internal phase change
            if (h[i] > 1) {
                for (int b = 1; b <= h[i]; b++) {
                    if (b == ai || D0[i][ai - 1][b - 1] == 0) {
                        continue;
                    }
                    int[] t = s.clone();
                    t[i] = idx(ni, b, h[i]);
                    T.add(t);
                }
            }
            // completion routed to j
            double exits = 1.0;
            if (h[i] > 1) {
                exits = 0;
                for (int b = 0; b < h[i]; b++) {
                    exits += D1[i][ai - 1][b];
                }
            }
            if (exits == 0) {
                continue;
            }
            for (int j = 0; j < K; j++) {
                if (j == i || P[i][j] <= 0) {
                    continue;
                }
                int nj = population(s[j], h[j]);
                int aj = phase(s[j], h[j]);
                int newi = (h[i] == 1) ? idx(ni - 1, 1, 1) : 0;
                for (int bi = 1; bi <= h[i]; bi++) {
                    if (h[i] > 1) {
                        if (ni == 1) {
                            newi = 0;
                        } else if (pie[i][bi - 1] == 0) {
                            continue;
                        } else {
                            newi = idx(ni - 1, bi, h[i]);
                        }
                    } else if (bi > 1) {
                        continue;
                    }
                    for (int bj = 1; bj <= h[j]; bj++) {
                        int newj;
                        if (nj == 0) {
                            if (pie[j][bj - 1] == 0) {
                                continue;
                            }
                            newj = idx(1, bj, h[j]);
                        } else if (bj > 1) {
                            continue;
                        } else {
                            newj = idx(nj + 1, aj, h[j]);
                        }
                        int[] t = s.clone();
                        t[i] = newi;
                        t[j] = newj;
                        T.add(t);
                    }
                    if (ni == 1 && h[i] > 1) {
                        break;
                    }
                }
            }
        }
        return T.toArray(new int[T.size()][]);
    }
}
