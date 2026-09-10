/**
 * @file Chen-O'Cinneide regularization: Markov chain Monte Carlo estimator of the
 *       class throughputs and queue lengths of a closed multiclass product-form network.
 *
 * @since LINE 3.0
 */
package jline.api.pfqn.nc;

import java.util.Random;

import jline.GlobalConstants;
import jline.io.Ret;
import jline.solvers.SolverOptions;
import jline.util.matrix.Matrix;

/**
 * Markov chain Monte Carlo estimator of the class throughputs
 * {@code X(r) = G(N-e_r)/G(N)} and of the mean queue lengths {@code Q(i,r)} of a
 * CLOSED multiclass product-form (BCMP, no type changes) network, by the
 * REGULARIZATION algorithm of
 *
 * <p>W. Chen, C. A. O'Cinneide, "Towards a Polynomial-Time Randomized Algorithm for
 * Closed Product-Form Networks", ACM TOMACS 8(3):227-253, 1998.</p>
 *
 * <p>The three steps of the paper are:</p>
 *
 * <p>I. CONSTRUCT THE REGULARIZED NETWORK. Write {@code rho(i,r)} for the surrogate
 * traffic intensity of class r at station i -- here the service demand, since
 * {@code rho = lambda/mu} is a visit ratio over a service rate -- and
 * {@code rho(r) = sum_i rho(i,r)}. The regularized network has the same stations,
 * classes and populations, UNIT service rates at every station, the processor-sharing
 * discipline, and a routing matrix that depends on the destination only,
 * {@code P*(i->m | class r) = rho(m,r)/rho(r)}. By Theorem 2.1 it is a REVERSIBLE chain
 * with the SAME steady-state distribution as the original network, and its throughputs
 * satisfy {@code Theta*(r) = rho(r)*Theta(r)}.</p>
 *
 * <p>II. SIMULATE IT at service-completion epochs. With {@code Y(i,r)} the number of
 * class-r jobs at station i, {@code Y(i)} their total and {@code Psi_i(k)=min(s_i,k)}
 * the number of busy servers,</p>
 *
 * <pre>
 *     r(i,r) = Y(i,r)/Y(i) * Psi_i(Y(i)),   r(r) = sum_i r(i,r),
 *     r      = sum_i Psi_i(Y(i)),
 * </pre>
 *
 * <p>the next completion is of class r at station i with probability {@code r(i,r)/r},
 * and the conditional expected time to it is {@code 1/r}. Equation (10) of the paper is
 * the holding-time weighted ratio estimator
 * {@code Theta*(r) = sum_t r(r,t)/r(t) / sum_t 1/r(t)}, and the same weights give the
 * time-average queue lengths, which need no transformation at all because the two
 * networks share their steady state.</p>
 *
 * <p>III. TRANSFORM BACK: {@code X(r) = Theta*(r)/rho(r)}.</p>
 *
 * <p>Because P* forgets the station of origin and every station serves at unit rate, the
 * regularized chain has neither the slowly mixing routing chain nor the customer-trapping
 * slow station that make the original chain converge slowly. The paper proves
 * {@code O(N^2*M^3)} mixing in two special cases (Section 4) and reports the general
 * behaviour experimentally (Section 5).</p>
 *
 * <p>Delay (infinite-server) demand enters as ONE extra station with {@code s = Inf} and
 * demand Z. Aggregating infinite-server stations that way is exact in the product form,
 * since their joint term is multinomial in the per-class totals.</p>
 *
 * <p>Confidence: the run is split into non-overlapping batches (Schmeiser 1982, 30 by
 * default, the count used in the tables of the paper), the batch means of the ratio
 * estimator give a standard error, and the result reports the paper's two-sigma interval.
 * The estimator is a ratio of correlated averages, so it carries an {@code O(1/samples)}
 * bias on top of the initialization bias; the paper ignores both, this implementation
 * additionally discards a warm-up fraction (10% by default).</p>
 *
 * @see Pfqn_nc
 */
public final class Pfqn_mcmc {
    private Pfqn_mcmc() {}

    /** Default run length when the options carry none, matching pfqn_mcmc.m. */
    private static final int DEFAULT_SAMPLES = 100000;
    /** Schmeiser (1982), the batch count used in the tables of the paper. */
    private static final int DEFAULT_BATCHES = 30;
    /** Warm-up fraction discarded before accumulation starts. */
    private static final double DEFAULT_BURNIN = 0.1;

    /**
     * Single-server form: every queueing station serves one job at a time.
     *
     * @param L       (M x R) per-class service demands at the M queueing stations
     * @param N       (1 x R) closed population vector; finite and integer
     * @param Z       (1 x R) aggregated think times, or a matrix summed over its rows
     * @param options solver options; may be null
     * @return the throughput and queue-length estimates with their two-sigma intervals
     */
    public static Ret.pfqnMcmc pfqn_mcmc(Matrix L, Matrix N, Matrix Z, SolverOptions options) {
        return pfqn_mcmc(L, N, Z, null, options);
    }

    /**
     * Multiserver form, the paper's own selling point (its Tables IV and V).
     *
     * @param L       (M x R) per-class service demands at the M queueing stations
     * @param N       (1 x R) closed population vector; finite and integer
     * @param Z       (1 x R) aggregated think times, or a matrix summed over its rows
     * @param s       (M x 1) servers per station, infinite for an infinite server; null
     *                means all stations single-server
     * @param options solver options; may be null
     * @return the throughput and queue-length estimates with their two-sigma intervals
     */
    public static Ret.pfqnMcmc pfqn_mcmc(Matrix L, Matrix N, Matrix Z, Matrix s,
                                         SolverOptions options) {
        final int M = L.getNumRows();
        final int R = L.getNumCols();

        // Populations. The chain lives on the integer lattice sum_i Y(i,r) = N(r), so a
        // fractional population has no state space at all; this is not a matter of
        // accuracy and must not be rounded away silently.
        int[] pop = new int[R];
        for (int r = 0; r < R; r++) {
            double nr = N.get(r);
            if (Double.isInfinite(nr) || Double.isNaN(nr)) {
                throw new RuntimeException("pfqn_mcmc requires a closed model, but the "
                        + "population vector has an infinite entry.");
            }
            if (Math.abs(nr - Math.round(nr)) > GlobalConstants.FineTol) {
                throw new RuntimeException("pfqn_mcmc simulates a state space of integer "
                        + "populations, but N(" + (r + 1) + ") = " + nr + " is fractional. "
                        + "Use an asymptotic method ('le', 'ble', 'kt') on fractional "
                        + "populations.");
            }
            pop[r] = (int) Math.round(nr);
        }

        Matrix X = new Matrix(1, R);
        Matrix Q = new Matrix(M, R);
        Matrix Xse = new Matrix(1, R);
        Matrix Qse = new Matrix(M, R);
        int totalPop = 0;
        for (int r = 0; r < R; r++) totalPop += pop[r];
        if (totalPop == 0) {
            return new Ret.pfqnMcmc(X, Q, Xse, X.copy(), X.copy(), Qse, Q.copy(), Q.copy(),
                    0, 0L, 0L);
        }

        // Think times, summed over their rows as the reference does.
        double[] Zr = new double[R];
        if (Z != null && !Z.isEmpty()) {
            for (int r = 0; r < R && r < Z.getNumCols(); r++) {
                double acc = 0.0;
                for (int k = 0; k < Z.getNumRows(); k++) {
                    double v = Z.get(k, r);
                    if (!Double.isInfinite(v) && !Double.isNaN(v)) acc += v;
                }
                Zr[r] = acc;
            }
        }

        // Server counts.
        double[] svec0 = new double[M];
        if (s == null || s.isEmpty()) {
            for (int i = 0; i < M; i++) svec0[i] = 1.0;
        } else {
            if (s.length() != M) {
                throw new RuntimeException("pfqn_mcmc: the server count vector has "
                        + s.length() + " entries but L has " + M + " stations.");
            }
            for (int i = 0; i < M; i++) svec0[i] = s.get(i);
        }

        // ---- Step I: the regularized network ---------------------------------------
        // Only the surrogate traffic intensities rho(i,r) enter the product form, and
        // scaling a whole class column by a constant leaves the steady-state distribution
        // unchanged, so the demands are used as they are.
        boolean hasDelay = false;
        for (int r = 0; r < R; r++) {
            if (Zr[r] > 0) { hasDelay = true; break; }
        }
        final int Mx = hasDelay ? M + 1 : M;
        double[][] rho = new double[Mx][R];
        double[] svec = new double[Mx];
        for (int i = 0; i < M; i++) {
            svec[i] = svec0[i];
            for (int r = 0; r < R; r++) {
                double v = L.get(i, r);
                if (Double.isInfinite(v) || Double.isNaN(v) || v < 0) v = 0.0;
                rho[i][r] = v;
            }
        }
        if (hasDelay) {
            svec[M] = Double.POSITIVE_INFINITY;
            for (int r = 0; r < R; r++) rho[M][r] = Zr[r];
        }

        double[] rhoTot = new double[R];
        for (int r = 0; r < R; r++) {
            for (int i = 0; i < Mx; i++) rhoTot[r] += rho[i][r];
            if (pop[r] > 0 && rhoTot[r] <= 0) {
                throw new RuntimeException("pfqn_mcmc: class " + (r + 1) + " has a positive "
                        + "population but no demand anywhere in the network.");
            }
        }

        // Routing of the regularized network, P*(m|r) = rho(m,r)/rho(r), held as one
        // column of cumulative probabilities per class.
        double[][] Pstar = new double[Mx][R];
        double[][] cumP = new double[Mx][R];
        for (int r = 0; r < R; r++) {
            double den = Math.max(rhoTot[r], Double.MIN_VALUE);
            double acc = 0.0;
            for (int i = 0; i < Mx; i++) {
                Pstar[i][r] = rho[i][r] / den;
                acc += Pstar[i][r];
                cumP[i][r] = acc;
            }
            cumP[Mx - 1][r] = 1.0; // guard the last bin against a floating-point shortfall
        }

        // Run length, batching and warm-up.
        int samples = DEFAULT_SAMPLES;
        int nbatches = DEFAULT_BATCHES;
        double burninFrac = DEFAULT_BURNIN;
        long seed = 0L;
        boolean seeded = false;
        if (options != null) {
            if (options.samples > 0) samples = Math.max(1, options.samples);
            seed = options.seed;
            seeded = true;
            if (options.config != null) {
                if (options.config.mcmc_batches != null) {
                    nbatches = Math.max(1, options.config.mcmc_batches.intValue());
                }
                if (options.config.mcmc_burnin != null) {
                    burninFrac = Math.min(0.9, Math.max(0.0, options.config.mcmc_burnin));
                }
            }
        }
        final int batchLen = Math.max(1, samples / nbatches);
        samples = batchLen * nbatches;
        final long nburn = Math.round(burninFrac * samples);
        Random rand = seeded ? new Random(seed) : new Random();

        // Initial state: spread each class over the stations it can occupy in the
        // proportions P*(.|r), by largest remainder. That is the marginal the regularized
        // network would have with no queueing, so it costs nothing and starts the chain
        // far closer to stationarity than a single-station state.
        int[][] Y = new int[Mx][R];
        double[] Ytot = new double[Mx];
        for (int r = 0; r < R; r++) {
            if (pop[r] == 0) continue;
            double[] target = new double[Mx];
            int placed = 0;
            for (int i = 0; i < Mx; i++) {
                target[i] = pop[r] * Pstar[i][r];
                Y[i][r] = (int) Math.floor(target[i]);
                placed += Y[i][r];
            }
            int shortfall = pop[r] - placed;
            // Largest remainder: hand each remaining unit to the biggest fraction still
            // outstanding. Once a station is topped up its residue target-Y drops by one,
            // so it cannot win twice, which is the 'descend' sort of the reference.
            for (int k = 0; k < shortfall; k++) {
                int best = 0;
                double bestRem = Double.NEGATIVE_INFINITY;
                for (int i = 0; i < Mx; i++) {
                    double rem = target[i] - Y[i][r];
                    if (rem > bestRem) { bestRem = rem; best = i; }
                }
                Y[best][r]++;
            }
        }
        for (int i = 0; i < Mx; i++) {
            double tot = 0.0;
            for (int r = 0; r < R; r++) tot += Y[i][r];
            Ytot[i] = tot;
        }

        // ---- Step II: simulate at service-completion epochs -------------------------
        double[][] xnum = new double[nbatches][R];      // sum_t r(r,t)/r(t) within the batch
        double[][][] qnum = new double[nbatches][Mx][R];// sum_t Y(t)/r(t)   within the batch
        double[] den = new double[nbatches];            // sum_t 1/r(t)      within the batch
        double[] Psi = new double[Mx];
        double[] cPsi = new double[Mx];
        double[] rvec = new double[R];
        double[] cYv = new double[R];
        double w = 0.0;
        boolean stale = true;
        final long horizon = nburn + samples;
        for (long t = 1; t <= horizon; t++) {
            if (stale) {
                // (8)-(9): busy servers, per-class completion rates, total rate
                double totPsi = 0.0;
                for (int i = 0; i < Mx; i++) {
                    Psi[i] = Math.min(svec[i], Ytot[i]);
                    totPsi += Psi[i];
                    cPsi[i] = totPsi;
                }
                for (int r = 0; r < R; r++) rvec[r] = 0.0;
                for (int i = 0; i < Mx; i++) {
                    if (Ytot[i] <= 0) continue;
                    double rw = Psi[i] / Ytot[i];
                    for (int r = 0; r < R; r++) {
                        if (Y[i][r] != 0) rvec[r] += rw * Y[i][r];
                    }
                }
                w = 1.0 / totPsi;
                stale = false;
            }
            if (t > nburn) {
                int b = (int) ((t - nburn - 1) / batchLen);
                den[b] += w;
                for (int r = 0; r < R; r++) xnum[b][r] += w * rvec[r];
                for (int i = 0; i < Mx; i++) {
                    for (int r = 0; r < R; r++) {
                        if (Y[i][r] != 0) qnum[b][i][r] += w * Y[i][r];
                    }
                }
            }
            // Pick the completing station with probability Psi(i)/r, then the completing
            // class within it with probability Y(i,r)/Y(i); the product is the r(i,r)/r of
            // the paper, since sum_r Y(i,r)/Y(i)*Psi(i) = Psi(i).
            int i = pickCumulative(cPsi, Mx, rand.nextDouble() * cPsi[Mx - 1]);
            if (i < 0) {
                for (int k = Mx - 1; k >= 0; k--) {
                    if (Psi[k] > 0) { i = k; break; }
                }
            }
            double cY = 0.0;
            for (int r = 0; r < R; r++) { cY += Y[i][r]; cYv[r] = cY; }
            int cls = pickCumulative(cYv, R, rand.nextDouble() * cYv[R - 1]);
            if (cls < 0) {
                for (int k = R - 1; k >= 0; k--) {
                    if (Y[i][k] > 0) { cls = k; break; }
                }
            }
            // Route it. A self-transition leaves the state, hence the rates and the
            // weight, unchanged: skipping the recomputation is the saving described at
            // the end of Section 2 of the paper.
            double u = rand.nextDouble();
            int m = -1;
            for (int k = 0; k < Mx; k++) {
                if (cumP[k][cls] >= u) { m = k; break; }
            }
            if (m >= 0 && m != i) {
                Y[i][cls]--;
                Y[m][cls]++;
                Ytot[i] -= 1.0;
                Ytot[m] += 1.0;
                stale = true;
            }
        }

        // ---- Step III: back to the original network ---------------------------------
        // Theta(r) = Theta*(r)/rho(r) by (7) and (11); the queue lengths transfer
        // unchanged, the two networks sharing their steady-state distribution.
        double denTot = 0.0;
        for (int b = 0; b < nbatches; b++) denTot += den[b];
        for (int r = 0; r < R; r++) {
            double num = 0.0;
            for (int b = 0; b < nbatches; b++) num += xnum[b][r];
            X.set(0, r, (num / denTot) / rhoTot[r]);
        }
        for (int i = 0; i < M; i++) {
            for (int r = 0; r < R; r++) {
                double num = 0.0;
                for (int b = 0; b < nbatches; b++) num += qnum[b][i][r];
                Q.set(i, r, num / denTot);
            }
        }

        // Batch-means standard error and the two-sigma interval of the paper.
        if (nbatches > 1) {
            for (int r = 0; r < R; r++) {
                double[] v = new double[nbatches];
                for (int b = 0; b < nbatches; b++) v[b] = (xnum[b][r] / den[b]) / rhoTot[r];
                Xse.set(0, r, sampleStd(v) / Math.sqrt(nbatches));
            }
            for (int i = 0; i < M; i++) {
                for (int r = 0; r < R; r++) {
                    double[] v = new double[nbatches];
                    for (int b = 0; b < nbatches; b++) v[b] = qnum[b][i][r] / den[b];
                    Qse.set(i, r, sampleStd(v) / Math.sqrt(nbatches));
                }
            }
        }
        Matrix Xlo = new Matrix(1, R);
        Matrix Xhi = new Matrix(1, R);
        for (int r = 0; r < R; r++) {
            Xlo.set(0, r, X.get(0, r) - 2 * Xse.get(0, r));
            Xhi.set(0, r, X.get(0, r) + 2 * Xse.get(0, r));
        }
        Matrix Qlo = new Matrix(M, R);
        Matrix Qhi = new Matrix(M, R);
        for (int i = 0; i < M; i++) {
            for (int r = 0; r < R; r++) {
                Qlo.set(i, r, Q.get(i, r) - 2 * Qse.get(i, r));
                Qhi.set(i, r, Q.get(i, r) + 2 * Qse.get(i, r));
            }
        }
        return new Ret.pfqnMcmc(X, Q, Xse, Xlo, Xhi, Qse, Qlo, Qhi, nbatches, samples, nburn);
    }

    /** First index whose cumulative value reaches {@code u}, or -1 when none does. */
    private static int pickCumulative(double[] cum, int n, double u) {
        for (int i = 0; i < n; i++) {
            if (cum[i] >= u) return i;
        }
        return -1;
    }

    /** Sample standard deviation, normalized by n-1 as MATLAB's std is. */
    private static double sampleStd(double[] v) {
        int n = v.length;
        if (n < 2) return 0.0;
        double mean = 0.0;
        for (int i = 0; i < n; i++) mean += v[i];
        mean /= n;
        double acc = 0.0;
        for (int i = 0; i < n; i++) {
            double d = v[i] - mean;
            acc += d * d;
        }
        return Math.sqrt(acc / (n - 1));
    }
}
