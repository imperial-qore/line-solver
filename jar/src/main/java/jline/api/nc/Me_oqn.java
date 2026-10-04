/**
 * Maximum Entropy algorithm for Open Queueing Networks.
 *
 * Implements the ME algorithm from Kouvatsos (1994) "Entropy Maximisation
 * and Queueing Network Models", Section 3.2, with the GE/GE/c building
 * block of Section 3.4 (eq. 3.9) and the GE/GE/inf building block.
 *
 * @since LINE 3.0
 */
package jline.api.nc;

import jline.io.InputOutput;
import jline.util.matrix.Matrix;

public final class Me_oqn {
    private Me_oqn() {}

    public static MeOqnResult me_oqn(int M, int R, Matrix lambda0, Matrix Ca0,
                                      Matrix mu, Matrix Cs, Matrix[][] P) {
        return me_oqn(M, R, lambda0, Ca0, mu, Cs, P, Matrix.ones(M, 1), new MeOqnOptions());
    }

    public static MeOqnResult me_oqn(int M, int R, Matrix lambda0, Matrix Ca0,
                                      Matrix mu, Matrix Cs, Matrix[][] P,
                                      MeOqnOptions options) {
        return me_oqn(M, R, lambda0, Ca0, mu, Cs, P, Matrix.ones(M, 1), options);
    }

    /**
     * Maximum Entropy algorithm for Open Queueing Networks.
     *
     * @param M       number of queues (stations)
     * @param R       number of job classes
     * @param lambda0 external arrival rates [M x R]
     * @param Ca0     external arrival scvs [M x R]
     * @param mu      service rates [M x R]
     * @param Cs      service scvs [M x R]
     * @param P       routing probabilities, P[j][i].get(r,0) = p_ji,r
     * @param c       servers per queue [M x 1]; Double.POSITIVE_INFINITY marks
     *                an infinite-server (IS) queue
     * @param options algorithm options (tolerance, max iterations, verbosity)
     * @return mean queue lengths, response times, arrival/departure scvs,
     *         arrival rates (inclusive of self-loop revisits), utilizations
     *         and the iteration count
     */
    public static MeOqnResult me_oqn(int M, int R, Matrix lambda0, Matrix Ca0,
                                      Matrix mu, Matrix Cs, Matrix[][] P,
                                      Matrix c, MeOqnOptions options) {
        return me_oqn(M, R, lambda0, Ca0, mu, Cs, P, c, null, options);
    }

    /**
     * Maximum Entropy algorithm for Open Queueing Networks with
     * discipline-aware building blocks: insens[i] marks a station with an
     * insensitive scheduling discipline (PS, LCFS-PR), solved with the
     * product-form mean queue length L_r = rho_r/(1-rho) instead of the
     * FCFS GE formula.
     */
    public static MeOqnResult me_oqn(int M, int R, Matrix lambda0, Matrix Ca0,
                                      Matrix mu, Matrix Cs, Matrix[][] P,
                                      Matrix c, boolean[] insens, MeOqnOptions options) {
        if (insens == null) {
            insens = new boolean[M];
        }
        // see _kb/03-api-layer.md for rationale
        Matrix[][] P_eff = new Matrix[M][M];
        for (int j = 0; j < M; j++) {
            for (int i = 0; i < M; i++) {
                P_eff[j][i] = P[j][i].copy();
            }
        }
        Matrix mu_eff = mu.copy();
        Matrix Cs_eff = Cs.copy();

        for (int i = 0; i < M; i++) {
            for (int r = 0; r < R; r++) {
                double pii = P[i][i].get(r, 0);
                if (pii > 0) {
                    mu_eff.set(i, r, mu.get(i, r) * (1 - pii));
                    Cs_eff.set(i, r, pii + (1 - pii) * Cs.get(i, r));
                    for (int dest = 0; dest < M; dest++) {
                        P_eff[i][dest].set(r, 0, P[i][dest].get(r, 0) / (1 - pii));
                    }
                    P_eff[i][i].set(r, 0, 0.0);
                }
            }
        }

        // see _kb/03-api-layer.md for rationale
        Matrix lambda = Matrix.zeros(M, R);
        Matrix lambda_eff = Matrix.zeros(M, R);
        for (int r = 0; r < R; r++) {
            Matrix Pr = Matrix.zeros(M, M);
            for (int j = 0; j < M; j++) {
                for (int i = 0; i < M; i++) {
                    Pr.set(j, i, P[j][i].get(r, 0));
                }
            }
            Matrix PrT = Pr.transpose();
            Matrix A = Matrix.eye(M).sub(PrT);
            Matrix lambda0_r = Matrix.extractColumn(lambda0, r, null);
            Matrix lambda_r = Matrix.zeros(M, 1);
            Matrix.solve(A, lambda0_r, lambda_r);
            for (int i = 0; i < M; i++) {
                lambda.set(i, r, lambda_r.get(i, 0));
                lambda_eff.set(i, r, lambda_r.get(i, 0) * (1 - P[i][i].get(r, 0)));
            }
        }

        // see _kb/03-api-layer.md for rationale
        Matrix rho = Matrix.zeros(M, R);
        for (int i = 0; i < M; i++) {
            for (int r = 0; r < R; r++) {
                if (mu.get(i, r) > 0) {
                    if (Double.isInfinite(c.get(i, 0))) {
                        rho.set(i, r, lambda_eff.get(i, r) / mu_eff.get(i, r));
                    } else {
                        rho.set(i, r, lambda.get(i, r) / (c.get(i, 0) * mu.get(i, r)));
                    }
                }
            }
        }

        // Stability check (finite-server queues only)
        boolean[] unstable = new boolean[M];
        for (int i = 0; i < M; i++) {
            double sum = 0.0;
            for (int r = 0; r < R; r++) {
                sum += rho.get(i, r);
            }
            if (!Double.isInfinite(c.get(i, 0)) && sum >= 1.0) {
                unstable[i] = true;
                InputOutput.line_warning("me_oqn",
                        "Network is unstable (utilization >= 1 at queue %d)", i);
            }
        }

        // Step 2: Initialize arrival scvs
        Matrix Ca = Matrix.ones(M, R);
        Matrix Cd = Matrix.ones(M, R);
        Matrix L = Matrix.zeros(M, R);

        // Steps 4-5: fixed-point iteration on the arrival scvs
        double delta = Double.POSITIVE_INFINITY;
        int iters = 0;
        for (int iter = 1; iter <= options.getMaxIter(); iter++) {
            iters = iter;
            Matrix Ca_old = Ca.copy();

            // Step 4: GE-type mean queue length formulae
            for (int i = 0; i < M; i++) {
                double rho_i = 0.0;
                for (int r = 0; r < R; r++) {
                    rho_i += rho.get(i, r);
                }
                double ci = c.get(i, 0);
                if (Double.isInfinite(ci)) {
                    // GE/GE/inf queue: L = lambda/mu
                    for (int r = 0; r < R; r++) {
                        if (lambda_eff.get(i, r) > 0 && mu_eff.get(i, r) > 0) {
                            L.set(i, r, lambda_eff.get(i, r) / mu_eff.get(i, r));
                        }
                    }
                } else if (unstable[i]) {
                    for (int r = 0; r < R; r++) {
                        if (lambda_eff.get(i, r) > 0) {
                            L.set(i, r, Double.POSITIVE_INFINITY);
                        }
                    }
                } else if (ci == 1.0) {
                    if (insens[i]) {
                        // Insensitive disciplines (PS, LCFS-PR): product-form
                        // mql, exact irrespective of the service distribution
                        for (int r = 0; r < R; r++) {
                            if (lambda_eff.get(i, r) > 0 && mu_eff.get(i, r) > 0) {
                                L.set(i, r, rho.get(i, r) / (1 - rho_i));
                            }
                        }
                    } else {
                        // see _kb/03-api-layer.md for rationale
                        double resid = 0.0;
                        for (int u = 0; u < R; u++) {
                            if (lambda_eff.get(i, u) > 0 && mu_eff.get(i, u) > 0) {
                                resid += lambda_eff.get(i, u) * (Cs_eff.get(i, u) + Ca.get(i, u))
                                        / (mu_eff.get(i, u) * mu_eff.get(i, u));
                            }
                        }
                        for (int r = 0; r < R; r++) {
                            if (lambda_eff.get(i, r) > 0 && mu_eff.get(i, r) > 0) {
                                L.set(i, r, rho.get(i, r) * (Ca.get(i, r) + 1) / 2
                                        + lambda_eff.get(i, r) * resid / (2 * (1 - rho_i)));
                            }
                        }
                    }
                } else {
                    // see _kb/03-api-layer.md for rationale
                    double lam_a = 0.0;
                    for (int u = 0; u < R; u++) {
                        if (lambda_eff.get(i, u) > 0 && mu_eff.get(i, u) > 0) {
                            lam_a += lambda_eff.get(i, u);
                        }
                    }
                    if (lam_a > 0) {
                        double inv_a = 0.0;
                        double ES = 0.0;
                        double ES2 = 0.0;
                        for (int u = 0; u < R; u++) {
                            if (lambda_eff.get(i, u) > 0 && mu_eff.get(i, u) > 0) {
                                double wu = lambda_eff.get(i, u) / lam_a;
                                inv_a += wu / (Ca.get(i, u) + 1);
                                ES += wu / mu_eff.get(i, u);
                                ES2 += wu * (Cs_eff.get(i, u) + 1)
                                        / (mu_eff.get(i, u) * mu_eff.get(i, u));
                            }
                        }
                        double Ca_a = -1 + 1 / inv_a;
                        double Cs_a = ES2 / (ES * ES) - 1;
                        double L_a = geGecMql(lam_a, Ca_a, 1 / ES, Cs_a, (int) ci);
                        double Lq_a = L_a - lam_a * ES; // mean waiting-line length
                        for (int r = 0; r < R; r++) {
                            if (lambda_eff.get(i, r) > 0 && mu_eff.get(i, r) > 0) {
                                L.set(i, r, ci * rho.get(i, r)
                                        + (lambda_eff.get(i, r) / lam_a) * Lq_a);
                            }
                        }
                    }
                }
            }

            // Step 5a: departure scvs
            for (int j = 0; j < M; j++) {
                double rho_j = 0.0;
                for (int r = 0; r < R; r++) {
                    rho_j += rho.get(j, r);
                }
                double cj = c.get(j, 0);
                for (int r = 0; r < R; r++) {
                    if (lambda_eff.get(j, r) > 0) {
                        if (Double.isInfinite(cj)) {
                            // GE/GE/inf: interdeparture scv = interarrival scv
                            Cd.set(j, r, Ca.get(j, r));
                        } else if (unstable[j]) {
                            // Saturated server: departures follow the service process
                            Cd.set(j, r, Cs_eff.get(j, r));
                        } else if (cj == 1.0) {
                            // Eq. (3.6) on the class-r virtual queue, with the
                            // marginal utilization rhohat_r of eq. (3.3)
                            double rhohat = rho.get(j, r) * L.get(j, r)
                                    / (L.get(j, r) + rho_j - rho.get(j, r));
                            Cd.set(j, r, 2 * L.get(j, r) * (1 - rhohat)
                                    + Ca.get(j, r) * (1 - 2 * rhohat));
                        } else {
                            // GE/GE/c interdeparture scv (Section 4.2)
                            Cd.set(j, r, rho_j * (1 - rho_j) + (1 - rho_j) * Ca.get(j, r)
                                    + rho_j * rho_j * Cs_eff.get(j, r));
                        }
                    }
                }
            }

            // see _kb/03-api-layer.md for rationale
            for (int i = 0; i < M; i++) {
                for (int r = 0; r < R; r++) {
                    if (lambda_eff.get(i, r) > 0) {
                        double sumInv = 0.0;
                        for (int j = 0; j < M; j++) {
                            double pji = P_eff[j][i].get(r, 0);
                            if (pji > 0 && lambda_eff.get(j, r) > 0) {
                                double Cdji = 1.0 + pji * (Cd.get(j, r) - 1.0);
                                double weight = (lambda_eff.get(j, r) * pji) / lambda_eff.get(i, r);
                                sumInv += weight / (Cdji + 1.0);
                            }
                        }
                        if (lambda0.get(i, r) > 0) {
                            double weight0 = lambda0.get(i, r) / lambda_eff.get(i, r);
                            sumInv += weight0 / (Ca0.get(i, r) + 1.0);
                        }
                        if (sumInv > 0) {
                            Ca.set(i, r, -1.0 + 1.0 / sumInv);
                        }
                    }
                }
            }

            // Check convergence
            delta = 0.0;
            for (int i = 0; i < M; i++) {
                for (int r = 0; r < R; r++) {
                    delta = Math.max(delta, Math.abs(Ca.get(i, r) - Ca_old.get(i, r)));
                }
            }

            if (options.getVerbose()) {
                System.out.println("Iteration " + iter + ": max delta = " + delta);
            }

            if (delta < options.getTol()) {
                if (options.getVerbose()) {
                    System.out.println("Converged after " + iter + " iterations");
                }
                break;
            }
        }

        if (iters == options.getMaxIter() && delta >= options.getTol()) {
            InputOutput.line_warning("me_oqn",
                    "Did not converge within %d iterations (delta=%f)",
                    options.getMaxIter(), delta);
        }

        // Step 6: response times by Little's law on the reported arrival rates
        Matrix W = Matrix.zeros(M, R);
        for (int i = 0; i < M; i++) {
            for (int r = 0; r < R; r++) {
                if (lambda.get(i, r) > 0) {
                    W.set(i, r, L.get(i, r) / lambda.get(i, r));
                }
            }
        }

        return new MeOqnResult(L, W, Ca, Cd, lambda, rho, iters);
    }

    /**
     * Mean queue length of a stable GE/GE/c/FCFS queue via the exact ME
     * solution of Kouvatsos (1994), eq. (3.9).
     */
    /** Mean queue length of a stable infinite-capacity GE/GE/c/FCFS queue, eq. (3.9). Shared with {@link Me_oqn_blk}. */
    static double geGecMql(double lambda, double Ca, double mu, double Cs, int c) {
        double alpha2 = 2.0 / (Cs + 1.0);
        double alpha1 = 1.0 - alpha2;
        double beta2 = 2.0 / (Ca + 1.0);
        double beta1 = 1.0 - beta2;
        double lambda2 = beta2 * lambda;
        double mu2 = alpha2 * mu;
        double[] g = new double[c];
        for (int j = 1; j < c; j++) {
            g[j - 1] = (lambda2 + (j - 1) * mu2 * beta1) * alpha2
                    / (j * mu2 * (1.0 - alpha1 * beta1));
        }
        g[c - 1] = (lambda2 + (c - 1) * mu2 * beta1) * alpha2
                / (lambda2 * alpha1 + c * mu2);
        double x = (lambda2 + c * mu2 * beta1) / (lambda2 * alpha1 + c * mu2);
        double[] Gn = new double[c];
        double prod = 1.0;
        for (int n = 0; n < c; n++) {
            prod *= g[n];
            Gn[n] = prod;
        }
        double Z = 1.0 + Gn[c - 1] / (1.0 - x);
        double S1 = 0.0;
        for (int n = 1; n <= c - 1; n++) {
            Z += Gn[n - 1];
            S1 += n * Gn[n - 1];
        }
        double S2 = Gn[c - 1] * (c / (1.0 - x) + x / ((1.0 - x) * (1.0 - x)));
        return (S1 + S2) / Z;
    }
}
