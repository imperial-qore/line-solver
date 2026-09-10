/**
 * CTMC Transient Analysis via Fast Adaptive Uniformization
 *
 * Adaptive uniformization (van Moorsel and Sanders, 1994) draws the uniformization
 * rate of each step from the states the iterate actually occupies rather than from
 * the whole state space, so the subordinating process is a pure birth process
 * instead of a Poisson process. Fast adaptive uniformization (Mateescu, Wolf,
 * Didier and Henzinger, 2010) adds the dropping of states below an occupancy
 * threshold, which is what turns a population cutoff into a numerical one.
 *
 * @since LINE 3.0
 */
package jline.api.mc;

import java.util.Arrays;
import java.util.Iterator;

import org.apache.commons.math3.util.FastMath;

import jline.util.matrix.Matrix;
import jline.util.matrix.MatrixEntry;

public final class Ctmc_fau {
    private Ctmc_fau() {}

    /**
     * Default cap on the number of birth steps, so a pathological horizon fails
     * loudly through {@link CtmcFauResult#truncated} rather than running forever.
     */
    public static final int FAU_MAX_STEPS = 1000000;

    /**
     * Transient distribution of a fast adaptive uniformization sweep, together
     * with the diagnostics that make its error auditable.
     */
    public static final class CtmcFauResult {
        /** Defective distribution at time t, a componentwise lower bound on pi(t). */
        public final Matrix pit;
        /** Number of birth steps K+1 actually taken. */
        public final int steps;
        /** Smallest adaptive rate used. */
        public final double lambdaMin;
        /** Largest adaptive rate used, the Lstar of the weight computation. */
        public final double lambdaMax;
        /** max_i |q_ii|, the rate ordinary uniformization would have used. */
        public final double uniformRate;
        /** Mass reaching the overflow index, that is P{N(t) &gt; K}. */
        public final double weightTail;
        /** Poisson mass outside the Fox-Glynn window of the weights. */
        public final double weightWindow;
        /** Probability removed by the occupancy threshold. */
        public final double droppedMass;
        /** sum(pi0) - sum(pit), which IS the L1 error against pi(t). */
        public final double errorBound;
        /** Largest occupied support over the sweep. */
        public final int supportMax;
        /** Support at the last step. */
        public final int supportFinal;
        /** True if maxsteps stopped the sweep. */
        public final boolean truncated;
        /** True if the support emptied or became absorbing. */
        public final boolean absorbed;

        CtmcFauResult(Matrix pit, int steps, double lambdaMin, double lambdaMax, double uniformRate,
                      double weightTail, double weightWindow, double droppedMass, double errorBound,
                      int supportMax, int supportFinal, boolean truncated, boolean absorbed) {
            this.pit = pit;
            this.steps = steps;
            this.lambdaMin = lambdaMin;
            this.lambdaMax = lambdaMax;
            this.uniformRate = uniformRate;
            this.weightTail = weightTail;
            this.weightWindow = weightWindow;
            this.droppedMass = droppedMass;
            this.errorBound = errorBound;
            this.supportMax = supportMax;
            this.supportFinal = supportFinal;
            this.truncated = truncated;
            this.absorbed = absorbed;
        }
    }

    /**
     * Return the transient distribution of the CTMC at time t by fast adaptive
     * uniformization, at the default tolerances.
     *
     * @param pi0 Initial distribution of the CTMC
     * @param Q   Infinitesimal generator of the CTMC
     * @param t   Transient analysis period boundary [0,t]
     * @return Distribution at time t with its error diagnostics
     */
    public static CtmcFauResult ctmc_fau(Matrix pi0, Matrix Q, double t) {
        return ctmc_fau(pi0, Q, t, 1e-6, 1e-12, -1);
    }

    /**
     * Return the transient distribution of the CTMC at time t by fast adaptive
     * uniformization.
     *
     * <p>Ordinary uniformization fixes one rate q &gt;= max_i |q_ii| over the whole
     * state space and mixes the powers of P = I + Q/q against a Poisson(q*t) law,
     * so its cost is set by the fastest state anywhere, including states that carry
     * no probability at time t. Adaptive uniformization instead picks a rate per
     * step from the states the iterate occupies,</p>
     *
     * <pre>
     *   Lambda_n &gt;= max{|q_ii| : i in supp(u^(n))},  u^(n+1) = u^(n)(I + Q/Lambda_n),
     * </pre>
     *
     * <p>which keeps every entry of u^(n+1) nonnegative. The subordinating process
     * is then the pure birth process N(t) with rates Lambda_0, Lambda_1, ... and
     * pi(t) = sum_n P{N(t) = n} u^(n). The fast variant drops an entry below delta
     * rather than propagating it, so the support tracks the states of
     * non-negligible occupancy instead of the reachable set.</p>
     *
     * <p>Nothing is renormalized anywhere, so the error is not estimated but
     * measured: the birth index truncated at K, the Poisson window of the weight
     * computation and the delta threshold each remove mass and none puts any back,
     * whence 0 &lt;= pi(t) - pit componentwise and
     * |pi(t) - pit|_1 = sum(pi0) - sum(pit) = errorBound.</p>
     *
     * <p>The birth weights are exact rather than quadratured: the rates generate a
     * bidiagonal generator on the birth index plus one absorbing overflow index,
     * and its transient distribution is obtained by uniformizing that scalar chain
     * at Lstar = max_n Lambda_n and applying the Fox-Glynn weights. The sweep runs
     * twice because b_n(t) needs the rates up to n, which are not known before the
     * sweep ends, while u^(n) is needed after them, and storing every iterate would
     * cost K times the support. Stopping is certified by the stochastic domination
     * of the birth epochs by an Erlang, so this method never takes more steps than
     * uniformization at the largest rate it visited.</p>
     *
     * <p>This is a transient method: it produces no stationary distribution.</p>
     *
     * @param pi0      Initial distribution of the CTMC
     * @param Q        Infinitesimal generator of the CTMC
     * @param t        Transient analysis period boundary [0,t]
     * @param epsilon  Birth-process truncation tolerance
     * @param delta    Occupancy threshold below which a state is dropped
     * @param maxsteps Cap on birth steps; nonpositive for the default cap
     * @return Distribution at time t with its error diagnostics
     */
    public static CtmcFauResult ctmc_fau(Matrix pi0, Matrix Q, double t, double epsilon, double delta,
                                         int maxsteps) {
        double eps = (epsilon > 0.0) ? epsilon : 1e-6;
        double drop = (delta > 0.0) ? delta : 0.0;
        int cap = (maxsteps > 0) ? maxsteps : FAU_MAX_STEPS;

        int n = Q.getNumRows();
        if (Q.getNumCols() != n) {
            throw new IllegalArgumentException("Q must be square.");
        }
        double[] p0 = rowVector(pi0, n);
        if (t < 0.0) {
            throw new IllegalArgumentException("t must be nonnegative.");
        }

        double[] d = new double[n];
        double uniformRate = 0.0;
        for (int i = 0; i < n; i++) {
            d[i] = -Q.get(i, i);
            if (d[i] > uniformRate) {
                uniformRate = d[i];
            }
        }
        if (t == 0.0 || n == 0) {
            int support = 0;
            for (int i = 0; i < n; i++) {
                if (p0[i] != 0.0) {
                    support++;
                }
            }
            return new CtmcFauResult(asRow(p0), 1, 0.0, 0.0, uniformRate, 0.0, 0.0, 0.0, 0.0,
                    support, support, false, false);
        }

        Csr csr = new Csr(Q, n);

        // Pass one: the adaptive rate sequence, and where it stops.
        RateSweep sweep = rates(p0, csr, d, t, drop, cap, eps);

        // The birth-process weights of that rate sequence, exactly.
        Weights weights = birthWeights(sweep.lambda, sweep.steps, t, eps);

        // Pass two: the same sweep again, accumulating sum_n b_n u^(n).
        Accumulation acc = accumulate(p0, csr, d, drop, sweep.steps, weights.b);

        double mass0 = 0.0;
        double massT = 0.0;
        for (int i = 0; i < n; i++) {
            mass0 += p0[i];
            massT += acc.pit[i];
        }
        double lambdaMin = 0.0;
        double lambdaMax = 0.0;
        if (sweep.steps > 0) {
            lambdaMin = sweep.lambda[0];
            lambdaMax = sweep.lambda[0];
            for (int m = 1; m < sweep.steps; m++) {
                lambdaMin = FastMath.min(lambdaMin, sweep.lambda[m]);
                lambdaMax = FastMath.max(lambdaMax, sweep.lambda[m]);
            }
        }
        return new CtmcFauResult(asRow(acc.pit), sweep.steps, lambdaMin, lambdaMax, uniformRate,
                weights.tail, weights.window, acc.dropped, mass0 - massT, acc.supportMax,
                acc.supportFinal, sweep.truncated, sweep.absorbed);
    }

    /**
     * Rows of Q in compressed form. Both sweeps read the rows of the occupied
     * states only, which is the whole point of the method, and neither the dense
     * nor the column-compressed layout of {@link Matrix} indexes rows directly.
     */
    private static final class Csr {
        final int[] rowPtr;
        final int[] colIdx;
        final double[] val;

        Csr(Matrix Q, int n) {
            int[] counts = new int[n];
            int nnz = 0;
            Iterator<MatrixEntry> it = Q.nonZeroIterator();
            while (it.hasNext()) {
                MatrixEntry e = it.next();
                counts[e.row]++;
                nnz++;
            }
            rowPtr = new int[n + 1];
            for (int i = 0; i < n; i++) {
                rowPtr[i + 1] = rowPtr[i] + counts[i];
            }
            colIdx = new int[nnz];
            val = new double[nnz];
            int[] fill = new int[n];
            it = Q.nonZeroIterator();
            while (it.hasNext()) {
                MatrixEntry e = it.next();
                int pos = rowPtr[e.row] + fill[e.row];
                colIdx[pos] = e.col;
                val[pos] = e.value;
                fill[e.row]++;
            }
        }
    }

    private static final class RateSweep {
        final double[] lambda;
        final int steps;
        final boolean truncated;
        final boolean absorbed;

        RateSweep(double[] lambda, int steps, boolean truncated, boolean absorbed) {
            this.lambda = lambda;
            this.steps = steps;
            this.truncated = truncated;
            this.absorbed = absorbed;
        }
    }

    private static final class Weights {
        final double[] b;
        final double tail;
        final double window;

        Weights(double[] b, double tail, double window) {
            this.b = b;
            this.tail = tail;
            this.window = window;
        }
    }

    private static final class Accumulation {
        final double[] pit;
        final double dropped;
        final int supportMax;
        final int supportFinal;

        Accumulation(double[] pit, double dropped, int supportMax, int supportFinal) {
            this.pit = pit;
            this.dropped = dropped;
            this.supportMax = supportMax;
            this.supportFinal = supportFinal;
        }
    }

    /**
     * Sweep the iterate to collect the adaptive rates Lambda_0..Lambda_K, stopping
     * when the Poisson-dominance bound on P{N(t) &gt; K} falls to epsilon.
     */
    private static RateSweep rates(double[] pi0, Csr csr, double[] d, double t, double delta,
                                   int maxsteps, double epsilon) {
        int n = pi0.length;
        double[] u = pi0.clone();
        int[] act = support(u);
        double[] lambda = new double[16];
        int steps = 0;
        boolean truncated = false;
        boolean absorbed = false;
        double lstar = 0.0;
        double[] scratch = new double[n];
        boolean[] touchedFlag = new boolean[n];
        int[] touched = new int[n];
        while (true) {
            if (act.length == 0) {
                absorbed = true;
                break;
            }
            double L = maxExit(d, act);
            if (steps == lambda.length) {
                lambda = Arrays.copyOf(lambda, 2 * lambda.length);
            }
            lambda[steps++] = L;
            if (L <= 0.0) {
                // Every occupied state is absorbing: the birth process stops here
                // and the remaining weight falls entirely on this iterate.
                absorbed = true;
                break;
            }
            lstar = FastMath.max(lstar, L);
            if (tailBound(lstar, t, steps) <= epsilon) {
                break;
            }
            if (steps >= maxsteps) {
                truncated = true;
                break;
            }
            act = step(u, act, csr, L, delta, scratch, touchedFlag, touched, null);
        }
        return new RateSweep(lambda, steps, truncated, absorbed);
    }

    /**
     * Replay the sweep of {@link #rates}, accumulating sum_n b_n u^(n). The
     * arithmetic is identical, so the rates and the drops reproduce those of the
     * first pass.
     */
    private static Accumulation accumulate(double[] pi0, Csr csr, double[] d, double delta, int nsteps,
                                           double[] b) {
        int n = pi0.length;
        double[] u = pi0.clone();
        double[] pit = new double[n];
        int[] act = support(u);
        int supportMax = act.length;
        int supportFinal = act.length;
        double[] dropped = new double[1];
        double[] scratch = new double[n];
        boolean[] touchedFlag = new boolean[n];
        int[] touched = new int[n];
        for (int m = 0; m < nsteps; m++) {
            if (act.length == 0) {
                break;
            }
            supportMax = FastMath.max(supportMax, act.length);
            supportFinal = act.length;
            for (int k = 0; k < act.length; k++) {
                pit[act[k]] += b[m] * u[act[k]];
            }
            if (m < nsteps - 1) {
                double L = maxExit(d, act);
                if (L <= 0.0) {
                    break;
                }
                act = step(u, act, csr, L, delta, scratch, touchedFlag, touched, dropped);
            }
        }
        return new Accumulation(pit, dropped[0], supportMax, supportFinal);
    }

    /**
     * One adaptive uniformization step u &lt;- u(I + Q/L), touching only the rows of
     * Q in the current support, followed by the drop rule. A state with a zero exit
     * rate is absorbing: its row of Q is empty, so it holds its mass and stays in
     * the support.
     */
    private static int[] step(double[] u, int[] act, Csr csr, double L, double delta, double[] scratch,
                              boolean[] touchedFlag, int[] touched, double[] dropped) {
        int ntouched = 0;
        for (int k = 0; k < act.length; k++) {
            int i = act[k];
            double ui = u[i];
            for (int pos = csr.rowPtr[i]; pos < csr.rowPtr[i + 1]; pos++) {
                int j = csr.colIdx[pos];
                if (!touchedFlag[j]) {
                    touchedFlag[j] = true;
                    touched[ntouched++] = j;
                }
                scratch[j] += ui * csr.val[pos];
            }
        }
        if (ntouched == 0) {
            return act;
        }
        Arrays.sort(touched, 0, ntouched);
        int kept = 0;
        for (int k = 0; k < ntouched; k++) {
            int j = touched[k];
            double contrib = scratch[j];
            scratch[j] = 0.0;
            touchedFlag[j] = false;
            if (contrib == 0.0) {
                // A row that cancels exactly leaves its state untouched, and the
                // union below carries the surviving part of the old support anyway.
                touched[k] = -1;
                continue;
            }
            double v = u[j] + contrib / L;
            if (v < delta) {
                if (dropped != null && v > 0.0) {
                    dropped[0] += v;
                }
                v = 0.0;
            }
            u[j] = v;
            if (v > 0.0) {
                touched[kept++] = j;
            } else {
                touched[k] = -1;
            }
        }
        // Sorted union of the surviving old support with the states just written.
        int[] survivors = new int[act.length];
        int nsurv = 0;
        for (int k = 0; k < act.length; k++) {
            if (u[act[k]] > 0.0) {
                survivors[nsurv++] = act[k];
            }
        }
        return union(survivors, nsurv, touched, kept);
    }

    private static int[] union(int[] a, int na, int[] b, int nb) {
        int[] out = new int[na + nb];
        int i = 0;
        int j = 0;
        int m = 0;
        while (i < na && j < nb) {
            if (a[i] == b[j]) {
                out[m++] = a[i++];
                j++;
            } else if (a[i] < b[j]) {
                out[m++] = a[i++];
            } else {
                out[m++] = b[j++];
            }
        }
        while (i < na) {
            out[m++] = a[i++];
        }
        while (j < nb) {
            out[m++] = b[j++];
        }
        return Arrays.copyOf(out, m);
    }

    /**
     * Transient distribution of the pure birth process with the given rates at time
     * t, that is b_n = P{N(t) = n} for n = 0..K, plus the mass that reached the
     * absorbing overflow index K+1 and therefore measures P{N(t) &gt; K}.
     *
     * The chain is uniformized at Lstar = max(lambda) and mixed against Fox-Glynn
     * Poisson weights, so the kernel entries 1 - Lambda_n/Lstar and Lambda_n/Lstar
     * are probabilities and nothing cancels. The weights are taken UNNORMALIZED, so
     * the Poisson mass outside the window is missing from b rather than
     * redistributed over it: b is then a sub-distribution, every term of the
     * mixture is an underestimate, and the error stays measurable as missing mass.
     */
    private static Weights birthWeights(double[] lambda, int k1, double t, double tol) {
        if (k1 == 0) {
            return new Weights(new double[0], 0.0, 0.0);
        }
        double lstar = 0.0;
        for (int m = 0; m < k1; m++) {
            lstar = FastMath.max(lstar, lambda[m]);
        }
        double[] b = new double[k1];
        if (lstar <= 0.0 || t <= 0.0) {
            b[0] = 1.0;
            return new Weights(b, 0.0, 0.0);
        }
        Ctmc_foxglynn.FoxGlynnWeights fg = Ctmc_foxglynn.ctmc_foxglynn_weights(lstar * t, tol, -1, false);
        double wsum = 0.0;
        for (int i = 0; i < fg.w.length; i++) {
            wsum += fg.w[i];
        }
        double window = FastMath.max(1.0 - wsum, 0.0);

        double[] v = new double[k1 + 1];
        double[] acc = new double[k1 + 1];
        v[0] = 1.0;
        for (int k = 0; k <= fg.right; k++) {
            if (k >= fg.left) {
                double wk = fg.w[k - fg.left];
                for (int m = 0; m <= k1; m++) {
                    acc[m] += wk * v[m];
                }
            }
            if (k < fg.right) {
                for (int m = k1; m >= 1; m--) {
                    double forward = (m - 1 < k1) ? v[m - 1] * lambda[m - 1] / lstar : 0.0;
                    double stay = (m < k1) ? v[m] * (1.0 - lambda[m] / lstar) : v[m];
                    v[m] = stay + forward;
                }
                v[0] = v[0] * (1.0 - lambda[0] / lstar);
            }
        }
        System.arraycopy(acc, 0, b, 0, k1);
        return new Weights(b, acc[k1], window);
    }

    /**
     * Upper bound on P{N(t) &gt;= k} for the birth process, through the stochastic
     * domination of its k-th jump epoch by an Erlang(k, lstar): the bound is the
     * Poisson(lstar*t) upper tail P{X &gt;= k}, taken at its Chernoff exponent
     * lambda*h(k/lambda) with h(u) = u*log(u) - u + 1. That exponent bounds the
     * upper tail only above the mean, so below it the bound is left vacuous.
     */
    private static double tailBound(double lstar, double t, int k) {
        double lambda = lstar * t;
        if (lambda <= 0.0 || k <= lambda) {
            return 1.0;
        }
        return FastMath.exp(-(lambda - k + k * FastMath.log(k / lambda)));
    }

    private static double maxExit(double[] d, int[] act) {
        double L = 0.0;
        for (int k = 0; k < act.length; k++) {
            if (d[act[k]] > L) {
                L = d[act[k]];
            }
        }
        return L;
    }

    private static int[] support(double[] u) {
        int count = 0;
        for (int i = 0; i < u.length; i++) {
            if (u[i] > 0.0) {
                count++;
            }
        }
        int[] act = new int[count];
        int m = 0;
        for (int i = 0; i < u.length; i++) {
            if (u[i] > 0.0) {
                act[m++] = i;
            }
        }
        return act;
    }

    private static double[] rowVector(Matrix pi0, int n) {
        if (pi0.getNumElements() != n) {
            throw new IllegalArgumentException("pi0 and Q have inconsistent sizes.");
        }
        double[] out = new double[n];
        for (int i = 0; i < n; i++) {
            out[i] = pi0.get(i);
        }
        return out;
    }

    private static Matrix asRow(double[] v) {
        Matrix out = new Matrix(1, v.length);
        for (int i = 0; i < v.length; i++) {
            out.set(0, i, v[i]);
        }
        return out;
    }
}
