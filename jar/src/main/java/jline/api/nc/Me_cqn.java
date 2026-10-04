/**
 * Maximum Entropy algorithm for Closed Queueing Networks.
 *
 * Implements the two-stage ME algorithm from Kouvatsos (1994) "Entropy
 * Maximisation and Queueing Network Models", Section 3.3, for closed
 * multiclass networks of G/G/1 and G/G/inf queues.
 *
 * @since LINE 3.0
 */
package jline.api.nc;

import jline.io.InputOutput;
import jline.util.matrix.Matrix;

public final class Me_cqn {
    private Me_cqn() {}

    public static MeCqnResult me_cqn(int M, int R, Matrix N, Matrix mu, Matrix Cs,
                                      Matrix[][] P, Matrix c, int[] refstat) {
        return me_cqn(M, R, N, mu, Cs, P, c, refstat, new MeOqnOptions());
    }

    /**
     * Maximum Entropy algorithm for Closed Queueing Networks.
     *
     * Stage 1 solves a pseudo-open network (no external arrivals) subject
     * to job flow conservation and the population constraints
     * sum_i L(i,r)=N(r), using the GE-type fixed point of the open
     * algorithm (Section 3.2) on class-composed flows. Stage 2 builds the
     * ME product-form solution (3.8) from the Lagrangian coefficients of
     * Stage 1, computes the normalising constant Z(N) by a multiclass
     * convolution, and iterates the flow (work rate) equations until the
     * class throughputs implied by the closed ME solution agree with
     * those used to parametrise the building blocks.
     *
     * @param M       number of queues (stations)
     * @param R       number of job classes
     * @param N       class populations [1 x R]
     * @param mu      service rates [M x R]
     * @param Cs      service scvs [M x R]
     * @param P       routing probabilities, P[j][i].get(r,0) = p_ji,r
     * @param c       servers per queue [M x 1]; Double.POSITIVE_INFINITY
     *                marks an IS queue; finite values must be 1
     * @param refstat reference station per class [R], 0-based; negative
     *                entries select the first station visited by the class
     * @param options algorithm options (tolerance, max iterations, verbosity)
     * @return closed mean queue lengths, response times, pseudo-open flow
     *         scvs, per-station throughputs, closed utilizations, class
     *         throughputs at the reference stations and iteration count
     */
    public static MeCqnResult me_cqn(int M, int R, Matrix N, Matrix mu, Matrix Cs,
                                      Matrix[][] P, Matrix c, int[] refstat,
                                      MeOqnOptions options) {
        return me_cqn(M, R, N, mu, Cs, P, c, refstat, null, options);
    }

    /**
     * Maximum Entropy algorithm for Closed Queueing Networks with
     * discipline-aware building blocks: insens[i] marks a station with an
     * insensitive scheduling discipline (PS, LCFS-PR), solved with the
     * product-form mean queue length instead of the FCFS GE formula.
     */
    public static MeCqnResult me_cqn(int M, int R, Matrix N, Matrix mu, Matrix Cs,
                                      Matrix[][] P, Matrix c, int[] refstat,
                                      boolean[] insens, MeOqnOptions options) {
        if (insens == null) {
            insens = new boolean[M];
        }
        double tol = options.getTol();
        int maxiter = options.getMaxIter();
        boolean verbose = options.getVerbose();

        // Feedback correction (as in the open algorithm)
        Matrix[][] P_eff = new Matrix[M][M];
        for (int j = 0; j < M; j++) {
            for (int i = 0; i < M; i++) {
                P_eff[j][i] = P[j][i].copy();
            }
        }
        Matrix mu_eff = mu.copy();
        Matrix Cs_eff = Cs.copy();
        Matrix selfp = new Matrix(M, R);
        for (int i = 0; i < M; i++) {
            for (int r = 0; r < R; r++) {
                double pii = P[i][i].get(r, 0);
                if (pii > 0) {
                    selfp.set(i, r, pii);
                    mu_eff.set(i, r, mu.get(i, r) * (1 - pii));
                    Cs_eff.set(i, r, pii + (1 - pii) * Cs.get(i, r));
                    for (int dest = 0; dest < M; dest++) {
                        P_eff[i][dest].set(r, 0, P[i][dest].get(r, 0) / (1 - pii));
                    }
                    P_eff[i][i].set(r, 0, 0.0);
                }
            }
        }

        // Visit ratios from the original routing (visit-inclusive),
        // normalised at the reference station of each class
        Matrix V = new Matrix(M, R);
        for (int r = 0; r < R; r++) {
            Matrix A = Matrix.eye(M);
            for (int j = 0; j < M; j++) {
                for (int i = 0; i < M; i++) {
                    A.set(j, i, A.get(j, i) - P[i][j].get(r, 0)); // I - P'
                }
            }
            int ref = refstat[r];
            if (ref < 0) {
                for (int i = 0; i < M; i++) {
                    if (mu.get(i, r) > 0) {
                        ref = i;
                        break;
                    }
                }
            }
            for (int i = 0; i < M; i++) {
                A.set(ref, i, 0.0);
            }
            A.set(ref, ref, 1.0);
            Matrix b = new Matrix(M, 1);
            b.set(ref, 0, 1.0);
            Matrix v = new Matrix(M, 1);
            Matrix.solve(A, b, v);
            for (int i = 0; i < M; i++) {
                double vi = v.get(i, 0);
                V.set(i, r, Math.abs(vi) < 1e-14 ? 0.0 : vi);
            }
            refstat[r] = ref;
        }

        // Stage 1: pseudo-open network, find X such that sum_i L(i,r)=N(r)
        double[] X = new double[R];
        for (int r = 0; r < R; r++) {
            double capr = Double.POSITIVE_INFINITY;
            for (int i = 0; i < M; i++) {
                if (!Double.isInfinite(c.get(i, 0)) && V.get(i, r) > 0 && mu.get(i, r) > 0) {
                    capr = Math.min(capr, mu.get(i, r) / V.get(i, r));
                }
            }
            if (Double.isInfinite(capr)) { // IS-only class
                capr = 1.0;
            }
            X[r] = 0.5 * capr / R;
        }

        Matrix Ca = Matrix.ones(M, R);
        Matrix Cd = Matrix.ones(M, R);
        Matrix Lpo = new Matrix(M, R);
        Matrix rho_po = new Matrix(M, R);
        Matrix lambda = new Matrix(M, R);
        int iters = 0;
        // see _kb/03-api-layer.md for rationale
        int maxit1 = Math.min(maxiter, 100);
        for (int it1 = 1; it1 <= maxit1; it1++) {
            iters++;
            capacityCap(X, V, mu, c, M, R);
            for (int i = 0; i < M; i++) {
                for (int r = 0; r < R; r++) {
                    lambda.set(i, r, V.get(i, r) * X[r]);
                }
            }
            pseudoOpen(M, R, lambda, mu, mu_eff, Cs_eff, P_eff, selfp, c, insens, Ca, Cd, Lpo, rho_po, tol, maxiter);
            double err1 = 0.0;
            double[] Ltot = new double[R];
            for (int r = 0; r < R; r++) {
                for (int i = 0; i < M; i++) {
                    Ltot[r] += Lpo.get(i, r);
                }
                if (N.get(0, r) > 0 && Ltot[r] > 0) {
                    err1 = Math.max(err1, Math.abs(Ltot[r] - N.get(0, r)) / N.get(0, r));
                }
            }
            if (verbose) {
                System.out.println("Stage 1 iteration " + it1 + ": max population error = " + err1);
            }
            if (err1 < tol) {
                break;
            }
            double[] Xold = X.clone();
            for (int r = 0; r < R; r++) {
                if (Ltot[r] > 0) {
                    double fac = Math.min(Math.max(Math.sqrt(N.get(0, r) / Ltot[r]), 0.25), 4.0);
                    X[r] = 0.5 * X[r] + 0.5 * X[r] * fac; // damped, clamped update
                }
            }
            // Stall guard: the stability cap can bind before the population
            // target is met; stop when X no longer moves
            double[] Xcap = X.clone();
            capacityCap(Xcap, V, mu, c, M, R);
            double move = 0.0;
            for (int r = 0; r < R; r++) {
                move = Math.max(move, Math.abs(Xcap[r] - Xold[r]) / Math.max(Xold[r], 1e-12));
            }
            if (move < tol) {
                break;
            }
        }

        // Stage 2: closed ME solution by convolution, iterated on the flow
        // (work rate) equations
        int[] sz = new int[R];
        int PIdx = 1;
        for (int r = 0; r < R; r++) {
            sz[r] = (int) Math.round(N.get(0, r)) + 1;
            PIdx *= sz[r];
        }
        int[][] Dec = new int[PIdx][R];
        for (int p = 0; p < PIdx; p++) {
            int q = p;
            for (int r = 0; r < R; r++) {
                Dec[p][r] = q % sz[r];
                q /= sz[r];
            }
        }
        int[] rad = new int[R];
        rad[0] = 1;
        for (int r = 1; r < R; r++) {
            rad[r] = rad[r - 1] * sz[r - 1];
        }
        int[] Nint = new int[R];
        for (int r = 0; r < R; r++) {
            Nint[r] = sz[r] - 1;
        }
        Matrix L = Lpo.copy();
        Matrix rho = rho_po.copy();
        double err2 = Double.POSITIVE_INFINITY;
        int it2 = 0;
        for (it2 = 1; it2 <= maxiter; it2++) {
            iters++;
            double[][] F = coefficients(M, R, PIdx, Dec, sz, Lpo, rho_po, lambda, mu_eff, Cs_eff, Ca, c, selfp);
            double[] U = new double[M];
            convolve(M, R, Nint, PIdx, Dec, rad, F, L, U);
            // Utilization split by pseudo-open per-class load; implied
            // throughputs from the work rate theorem with visit weights
            rho = new Matrix(M, R);
            double[] Xhat = new double[R];
            for (int r = 0; r < R; r++) {
                double num = 0.0;
                double den = 0.0;
                for (int i = 0; i < M; i++) {
                    if (lambda.get(i, r) > 0) {
                        if (Double.isInfinite(c.get(i, 0))) {
                            rho.set(i, r, L.get(i, r));
                            // IS work rate: lambda_eff = L*mu_eff, revisits add 1/(1-p)
                            num += L.get(i, r) * mu_eff.get(i, r) / (1 - selfp.get(i, r));
                        } else {
                            double rho_i = 0.0;
                            for (int u = 0; u < R; u++) {
                                rho_i += rho_po.get(i, u);
                            }
                            if (rho_i > 0) {
                                rho.set(i, r, U[i] * rho_po.get(i, r) / rho_i);
                            }
                            num += rho.get(i, r) * mu.get(i, r);
                        }
                        den += V.get(i, r);
                    }
                }
                if (den > 0) {
                    Xhat[r] = num / den;
                }
            }
            err2 = 0.0;
            for (int r = 0; r < R; r++) {
                if (X[r] > 0) {
                    err2 = Math.max(err2, Math.abs(Xhat[r] - X[r]) / X[r]);
                }
            }
            if (verbose) {
                System.out.println("Stage 2 iteration " + it2 + ": max flow error = " + err2);
            }
            if (err2 < tol) {
                break;
            }
            double[] Xold = X.clone();
            for (int r = 0; r < R; r++) {
                X[r] = 0.5 * X[r] + 0.5 * Xhat[r];
            }
            capacityCap(X, V, mu, c, M, R);
            // Stall guard: when the stability cap binds, X stops moving even
            // though the residual flow error stays above tolerance
            double move = 0.0;
            for (int r = 0; r < R; r++) {
                move = Math.max(move, Math.abs(X[r] - Xold[r]) / Math.max(Xold[r], 1e-12));
            }
            if (move < tol) {
                break;
            }
            for (int i = 0; i < M; i++) {
                for (int r = 0; r < R; r++) {
                    lambda.set(i, r, V.get(i, r) * X[r]);
                }
            }
            pseudoOpen(M, R, lambda, mu, mu_eff, Cs_eff, P_eff, selfp, c, insens, Ca, Cd, Lpo, rho_po, tol, maxiter);
        }

        if (it2 >= maxiter && err2 >= tol) {
            InputOutput.line_warning("me_cqn",
                    "Did not converge within %d iterations (flow error=%f)", maxiter, err2);
        }

        // Response times by Little's law on the visit-inclusive throughputs
        Matrix W = new Matrix(M, R);
        Matrix Xout = new Matrix(1, R);
        for (int r = 0; r < R; r++) {
            Xout.set(0, r, X[r]);
            for (int i = 0; i < M; i++) {
                double lir = V.get(i, r) * X[r];
                lambda.set(i, r, lir);
                if (lir > 0) {
                    W.set(i, r, L.get(i, r) / lir);
                }
            }
        }

        return new MeCqnResult(L, W, Ca, Cd, lambda, rho, Xout, iters);
    }

    /**
     * Scales the class throughputs uniformly so that every single-server
     * queue in the pseudo-open network remains stable.
     */
    private static void capacityCap(double[] X, Matrix V, Matrix mu, Matrix c, int M, int R) {
        double maxrho = 0.0;
        for (int i = 0; i < M; i++) {
            if (!Double.isInfinite(c.get(i, 0))) {
                double rho_i = 0.0;
                for (int r = 0; r < R; r++) {
                    if (V.get(i, r) > 0 && mu.get(i, r) > 0) {
                        rho_i += X[r] * V.get(i, r) / mu.get(i, r);
                    }
                }
                maxrho = Math.max(maxrho, rho_i);
            }
        }
        if (maxrho >= 0.999) {
            for (int r = 0; r < R; r++) {
                X[r] = X[r] * (0.999 / maxrho);
            }
        }
    }

    /**
     * GE-type fixed point of the open algorithm (Section 3.2) on the
     * pseudo-open network: no external arrivals, flows given by lambda.
     * The flow scvs are computed on the class-composed (aggregate) streams
     * and disaggregated per class by thinning, following the class
     * composition and disaggregation principle of the closed ME algorithm.
     */
    private static void pseudoOpen(int M, int R, Matrix lambda, Matrix mu, Matrix mu_eff,
                                   Matrix Cs_eff, Matrix[][] P_eff, Matrix selfp, Matrix c,
                                   boolean[] insens, Matrix Ca, Matrix Cd, Matrix L, Matrix rho,
                                   double tol, int maxiter) {
        Matrix lambda_eff = new Matrix(M, R);
        for (int i = 0; i < M; i++) {
            for (int r = 0; r < R; r++) {
                lambda_eff.set(i, r, lambda.get(i, r) * (1 - selfp.get(i, r)));
                rho.set(i, r, 0.0);
                if (mu.get(i, r) > 0) {
                    if (Double.isInfinite(c.get(i, 0))) {
                        rho.set(i, r, lambda_eff.get(i, r) / mu_eff.get(i, r));
                    } else {
                        rho.set(i, r, lambda.get(i, r) / mu.get(i, r));
                    }
                }
            }
        }
        // Class composition per station: aggregate flow, service process
        // moments and flow-weighted aggregate routing
        double[] lam_a = new double[M];
        double[] mu_a = new double[M];
        double[] Cs_a = new double[M];
        for (int i = 0; i < M; i++) {
            Cs_a[i] = 1.0;
            for (int r = 0; r < R; r++) {
                lam_a[i] += lambda_eff.get(i, r);
            }
            if (lam_a[i] > 0) {
                double ES = 0.0;
                double ES2 = 0.0;
                for (int u = 0; u < R; u++) {
                    if (lambda_eff.get(i, u) > 0 && mu_eff.get(i, u) > 0) {
                        double wu = lambda_eff.get(i, u) / lam_a[i];
                        ES += wu / mu_eff.get(i, u);
                        ES2 += wu * (Cs_eff.get(i, u) + 1) / (mu_eff.get(i, u) * mu_eff.get(i, u));
                    }
                }
                if (ES > 0) {
                    mu_a[i] = 1.0 / ES;
                    Cs_a[i] = ES2 / (ES * ES) - 1.0;
                }
            }
        }
        double[][] Pa = new double[M][M];
        for (int j = 0; j < M; j++) {
            if (lam_a[j] > 0) {
                for (int i = 0; i < M; i++) {
                    double num = 0.0;
                    for (int r = 0; r < R; r++) {
                        if (lambda_eff.get(j, r) > 0) {
                            num += lambda_eff.get(j, r) * P_eff[j][i].get(r, 0);
                        }
                    }
                    Pa[j][i] = num / lam_a[j];
                }
            }
        }
        // Fixed point on the aggregate arrival scvs
        double[] Ca_a = new double[M];
        for (int i = 0; i < M; i++) {
            Ca_a[i] = 1.0;
            if (lam_a[i] > 0) {
                for (int r = 0; r < R; r++) {
                    if (lambda_eff.get(i, r) > 0) {
                        Ca_a[i] = 1.0 + (Ca.get(i, r) - 1.0) * lam_a[i]
                                / Math.max(lambda_eff.get(i, r), Double.MIN_NORMAL); // warm start
                        break;
                    }
                }
            }
        }
        double[] Cd_a = new double[M];
        double[] L_a = new double[M];
        for (int i = 0; i < M; i++) {
            Cd_a[i] = 1.0;
        }
        for (int it = 1; it <= maxiter; it++) {
            double delta = 0.0;
            double[] Ca_old = Ca_a.clone();
            for (int i = 0; i < M; i++) {
                if (lam_a[i] <= 0) {
                    continue;
                }
                double rho_i = 0.0;
                for (int r = 0; r < R; r++) {
                    rho_i += rho.get(i, r);
                }
                if (Double.isInfinite(c.get(i, 0))) {
                    // GE/GE/inf: L = lambda/mu, departures inherit the arrival scv
                    L_a[i] = lam_a[i] / mu_a[i];
                    Cd_a[i] = Ca_a[i];
                } else if (rho_i < 1) {
                    if (insens[i]) {
                        // Insensitive disciplines (PS, LCFS-PR): product-form mql
                        L_a[i] = rho_i / (1 - rho_i);
                    } else {
                        // Single-class GE/GE/1 mql, eq. (3.6)
                        L_a[i] = rho_i * (Ca_a[i] + 1) / 2 + rho_i * rho_i * (Ca_a[i] + Cs_a[i]) / (2 * (1 - rho_i));
                    }
                    Cd_a[i] = 2 * L_a[i] * (1 - rho_i) + Ca_a[i] * (1 - 2 * rho_i);
                }
            }
            // GE-type merging, eq. (3.7) with lambda_o = 0, on aggregate flows
            for (int i = 0; i < M; i++) {
                if (lam_a[i] > 0) {
                    double sumInv = 0.0;
                    for (int j = 0; j < M; j++) {
                        if (Pa[j][i] > 0 && lam_a[j] > 0) {
                            double Cdji = 1.0 + Pa[j][i] * (Cd_a[j] - 1.0);
                            sumInv += (lam_a[j] * Pa[j][i] / lam_a[i]) / (Cdji + 1.0);
                        }
                    }
                    if (sumInv > 0) {
                        Ca_a[i] = -1.0 + 1.0 / sumInv;
                    }
                }
            }
            for (int i = 0; i < M; i++) {
                delta = Math.max(delta, Math.abs(Ca_a[i] - Ca_old[i]));
            }
            if (delta < tol) {
                break;
            }
        }
        // Disaggregation: per-class arrival scvs by thinning of the composed
        // stream, then per-class mean queue lengths (Section 3.1.1)
        for (int i = 0; i < M; i++) {
            double rho_i = 0.0;
            for (int r = 0; r < R; r++) {
                rho_i += rho.get(i, r);
                L.set(i, r, 0.0);
                Cd.set(i, r, 1.0);
            }
            for (int r = 0; r < R; r++) {
                if (lambda_eff.get(i, r) > 0) {
                    double pr = lambda_eff.get(i, r) / lam_a[i];
                    Ca.set(i, r, 1.0 + pr * (Ca_a[i] - 1.0));
                    Cd.set(i, r, 1.0 + pr * (Cd_a[i] - 1.0));
                }
            }
            if (Double.isInfinite(c.get(i, 0))) {
                for (int r = 0; r < R; r++) {
                    if (lambda_eff.get(i, r) > 0 && mu_eff.get(i, r) > 0) {
                        L.set(i, r, lambda_eff.get(i, r) / mu_eff.get(i, r));
                    }
                }
            } else if (rho_i < 1) {
                if (insens[i]) {
                    // Insensitive disciplines (PS, LCFS-PR): product-form mql
                    for (int r = 0; r < R; r++) {
                        if (lambda_eff.get(i, r) > 0 && mu_eff.get(i, r) > 0) {
                            L.set(i, r, rho.get(i, r) / (1 - rho_i));
                        }
                    }
                } else {
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
            }
        }
    }

    /**
     * Auxiliary functions f_i(n) of the ME solution (3.8): the right-hand
     * sides of (3.2) and (3.4) with the (1-rho) factor removed, evaluated
     * from the Stage 1 Lagrangian coefficients. Each f_i is rescaled by its
     * maximum for numerical stability (per-station constants cancel in the
     * marginal probabilities).
     */
    private static double[][] coefficients(int M, int R, int PIdx, int[][] Dec, int[] sz,
                                           Matrix Lpo, Matrix rho_po, Matrix lambda,
                                           Matrix mu_eff, Matrix Cs_eff, Matrix Ca,
                                           Matrix c, Matrix selfp) {
        double[][] F = new double[PIdx][M];
        Matrix lambda_eff = new Matrix(M, R);
        for (int i = 0; i < M; i++) {
            for (int r = 0; r < R; r++) {
                lambda_eff.set(i, r, lambda.get(i, r) * (1 - selfp.get(i, r)));
            }
        }
        for (int i = 0; i < M; i++) {
            if (Double.isInfinite(c.get(i, 0))) {
                // GE/GE/inf: f(n) = prod_r prod_{k=1}^{n_r} g_r(k)
                double[][] logg = new double[R][];
                for (int r = 0; r < R; r++) {
                    logg[r] = new double[Math.max(sz[r] - 1, 0)];
                    for (int j = 1; j < sz[r]; j++) {
                        if (lambda_eff.get(i, r) > 0 && mu_eff.get(i, r) > 0) {
                            double gj = (lambda_eff.get(i, r) * (1 + Cs_eff.get(i, r))
                                    + (j - 1) * mu_eff.get(i, r) * (Ca.get(i, r) - 1))
                                    / (j * mu_eff.get(i, r) * (Ca.get(i, r) + Cs_eff.get(i, r)));
                            logg[r][j - 1] = gj > 0 ? Math.log(gj) : Double.NEGATIVE_INFINITY;
                        } else {
                            logg[r][j - 1] = Double.NEGATIVE_INFINITY;
                        }
                    }
                }
                for (int p = 0; p < PIdx; p++) {
                    double val = 0.0;
                    for (int r = 0; r < R; r++) {
                        for (int j = 1; j <= Dec[p][r]; j++) {
                            val += logg[r][j - 1];
                        }
                    }
                    F[p][i] = Double.isInfinite(val) ? 0.0 : Math.exp(val);
                }
                F[0][i] = 1.0;
            } else {
                // see _kb/03-api-layer.md for rationale
                double rho_i = 0.0;
                double Li = 0.0;
                for (int r = 0; r < R; r++) {
                    rho_i += rho_po.get(i, r);
                    Li += Lpo.get(i, r);
                }
                double[] x = new double[R];
                double[] gx = new double[R];
                if (Li > 0 && rho_i < 1) {
                    for (int r = 0; r < R; r++) {
                        if (lambda.get(i, r) > 0) {
                            x[r] = Math.max(Lpo.get(i, r) - rho_po.get(i, r), 0.0) / Li;
                            gx[r] = rho_po.get(i, r) * rho_i / ((1 - rho_i) * Li);
                        }
                    }
                }
                for (int p = 0; p < PIdx; p++) {
                    int ntot = 0;
                    boolean absent = false;
                    for (int r = 0; r < R; r++) {
                        ntot += Dec[p][r];
                        if (Dec[p][r] > 0 && lambda.get(i, r) <= 0) {
                            absent = true; // class not visiting this station
                        }
                    }
                    if (ntot == 0) {
                        F[p][i] = 1.0;
                        continue;
                    }
                    if (absent) {
                        F[p][i] = 0.0;
                        continue;
                    }
                    double logmult = lgamma(ntot);
                    for (int r = 0; r < R; r++) {
                        logmult -= lgamma(Dec[p][r] + 1);
                    }
                    double tot = 0.0;
                    for (int r = 0; r < R; r++) {
                        if (Dec[p][r] > 0 && gx[r] > 0) {
                            double lterm = Math.log(Dec[p][r]) + Math.log(gx[r]);
                            boolean ok = true;
                            for (int s = 0; s < R; s++) {
                                int es = (s == r) ? Dec[p][s] - 1 : Dec[p][s];
                                if (es > 0) {
                                    if (x[s] > 0) {
                                        lterm += es * Math.log(x[s]);
                                    } else {
                                        ok = false;
                                        break;
                                    }
                                }
                            }
                            if (ok) {
                                tot += Math.exp(logmult + lterm);
                            }
                        }
                    }
                    F[p][i] = tot;
                }
            }
            double fmax = 0.0;
            for (int p = 0; p < PIdx; p++) {
                fmax = Math.max(fmax, F[p][i]);
            }
            if (fmax > 0) {
                for (int p = 0; p < PIdx; p++) {
                    F[p][i] /= fmax;
                }
            }
            F[0][i] = Math.max(F[0][i], Double.MIN_NORMAL);
        }
        return F;
    }

    /** log((k-1)!) via the log-gamma function. */
    private static double lgamma(int k) {
        double s = 0.0;
        for (int j = 2; j < k; j++) {
            s += Math.log(j);
        }
        return s;
    }

    /**
     * Convolution of a partial normalising constant with one station term
     * over the population lattice.
     */
    private static double[] convPair(double[] G, double[] f, int PIdx, int[][] Dec,
                                     int[] rad, int[] N, int R) {
        double[] G2 = new double[PIdx];
        for (int p = 0; p < PIdx; p++) {
            if (f[p] == 0) {
                continue;
            }
            for (int q = 0; q < PIdx; q++) {
                if (G[q] == 0) {
                    continue;
                }
                int idx = 0;
                boolean ok = true;
                for (int r = 0; r < R; r++) {
                    int t = Dec[p][r] + Dec[q][r];
                    if (t > N[r]) {
                        ok = false;
                        break;
                    }
                    idx += t * rad[r];
                }
                if (ok) {
                    G2[idx] += f[p] * G[q];
                }
            }
        }
        return G2;
    }

    /**
     * Computes the normalising constant by convolving the f_i over the
     * population lattice, and the per-station marginals by prefix/suffix
     * convolutions; fills the closed mean queue lengths and the per-station
     * busy probabilities.
     */
    private static void convolve(int M, int R, int[] N, int PIdx, int[][] Dec, int[] rad,
                                 double[][] F, Matrix L, double[] U) {
        double[] G0 = new double[PIdx];
        G0[0] = 1.0;
        double[][] Gpre = new double[M + 1][];
        Gpre[0] = G0;
        double[][] Fi = new double[M][PIdx];
        for (int k = 0; k < M; k++) {
            for (int p = 0; p < PIdx; p++) {
                Fi[k][p] = F[p][k];
            }
            Gpre[k + 1] = convPair(Gpre[k], Fi[k], PIdx, Dec, rad, N, R);
        }
        double[][] Gsuf = new double[M + 1][];
        Gsuf[M] = G0;
        for (int k = M - 1; k >= 0; k--) {
            Gsuf[k] = convPair(Gsuf[k + 1], Fi[k], PIdx, Dec, rad, N, R);
        }
        double Z = Gpre[M][PIdx - 1];
        for (int i = 0; i < M; i++) {
            U[i] = 0.0;
            for (int r = 0; r < R; r++) {
                L.set(i, r, 0.0);
            }
            double[] Grest = convPair(Gpre[i], Gsuf[i + 1], PIdx, Dec, rad, N, R);
            for (int p = 0; p < PIdx; p++) {
                if (F[p][i] > 0) {
                    int q = 0;
                    for (int r = 0; r < R; r++) {
                        q += (N[r] - Dec[p][r]) * rad[r];
                    }
                    double pin = F[p][i] * Grest[q] / Z;
                    if (p > 0) {
                        U[i] += pin;
                    }
                    for (int r = 0; r < R; r++) {
                        if (Dec[p][r] > 0) {
                            L.set(i, r, L.get(i, r) + Dec[p][r] * pin);
                        }
                    }
                }
            }
        }
    }
}
