/*
 * Copyright (c) 2012-2026, QORE Lab, Imperial College London
 * All rights reserved.
 */

package jline.lib.rmf;

import jline.io.Ret;
import jline.util.matrix.Matrix;
import odesolver.LSODA;
import org.apache.commons.math3.ode.FirstOrderDifferentialEquations;

/**
 * Multi-list cache with RANDOM(m) replacement as a Density-Dependent Population Process (DDPP).
 *
 * <p>Models a cache with h lists of capacities m = [m_1, ..., m_h] and n items
 * with popularity distribution p = [p_1, ..., p_n]. The state space tracks
 * in which list each item resides (or if it is outside the cache).</p>
 *
 * <p>State vector layout:
 *   x[i + k * n] = density of item i in list k
 *   where k=0 means "outside cache", k=1..h are the cache lists</p>
 *
 * <p>Reference:
 *   N. Gast, "Expected Values Estimated via Mean-Field Approximation are
 *   1/N-Accurate", Proc. ACM Meas. Anal. Comput. Syst., 2017.</p>
 */
public class CacheRMF {

    /** Item popularity distribution (normalized). */
    private final double[] p;

    /** List capacities. m[0] is outermost, m[h-1] is innermost. */
    private final int[] m;

    /** Number of items in the catalog (n). */
    private final int numberOfItems;

    /** Number of cache lists (h). */
    private final int numberOfLists;

    /** Total state dimension: n * (h + 1). */
    private final int modelDimension;

    /** Initial state vector. */
    private final double[] x0;

    /** Optional per-item access graph [item][h+1][h+1]; null = standard linear chain. */
    private double[][][] itemGraph = null;

    /** Set the per-item access graph enabling the general (accost-aware) drift. */
    public void setItemGraph(double[][][] g) {
        this.itemGraph = g;
    }

    /**
     * Result of dimension reduction.
     */
    public static class DimensionReduction {
        public final double[][] P;
        public final double[][] Pinv;
        public final int rank;

        public DimensionReduction(double[][] P, double[][] Pinv, int rank) {
            this.P = P;
            this.Pinv = Pinv;
            this.rank = rank;
        }
    }

    /**
     * Result of reduced system computation.
     */
    public static class ReducedSystem {
        public final double[][] Fp_r;
        public final double[][][] Fpp_r;
        public final double[][] Q_r;
        public final double[][] P;
        public final double[][] Pinv;
        public final int rank;

        public ReducedSystem(double[][] Fp_r, double[][][] Fpp_r, double[][] Q_r,
                             double[][] P, double[][] Pinv, int rank) {
            this.Fp_r = Fp_r;
            this.Fpp_r = Fpp_r;
            this.Q_r = Q_r;
            this.P = P;
            this.Pinv = Pinv;
            this.rank = rank;
        }
    }

    /**
     * Constructs a CacheRMF model.
     *
     * @param p Item request probabilities (popularity distribution), normalized.
     * @param m Capacity of each cache list.
     */
    public CacheRMF(double[] p, int[] m) {
        this.numberOfItems = p.length;
        this.numberOfLists = m.length;
        this.modelDimension = numberOfItems * (numberOfLists + 1);
        this.p = new double[p.length];
        System.arraycopy(p, 0, this.p, 0, p.length);
        this.m = new int[m.length];
        System.arraycopy(m, 0, this.m, 0, m.length);

        // Build initial state: first m[0] items in list 1, next m[1] in list 2, etc.
        this.x0 = new double[modelDimension];
        int objIdx = 0;
        for (int k = 0; k < numberOfLists; k++) {
            for (int cnt = 0; cnt < m[k]; cnt++) {
                if (objIdx < numberOfItems) {
                    x0[index(objIdx, k + 1)] = 1.0;
                    objIdx++;
                }
            }
        }
        for (int i = objIdx; i < numberOfItems; i++) {
            x0[index(i, 0)] = 1.0;
        }
    }

    /**
     * Map (item i, list k) to flat state index.
     *
     * @param i Item index (0-based).
     * @param k List index (0 = outside cache, 1..h = cache lists).
     * @return Flat index into state vector.
     */
    public int index(int i, int k) {
        return i + k * numberOfItems;
    }

    /**
     * Compute hit rate contribution from a specific list.
     *
     * @param x    State vector of dimension modelDimension.
     * @param listNumber List index (0 = outside cache, 1..h = cache lists).
     * @return Sum of p[i] * x[index(i, listNumber)] over all items i.
     */
    public double hitRate(double[] x, int listNumber) {
        double sum = 0.0;
        for (int i = 0; i < numberOfItems; i++) {
            sum += p[i] * x[index(i, listNumber)];
        }
        return sum;
    }

    /**
     * Compute the mean field drift F(x).
     *
     * <p>The drift for RANDOM(m) replacement is:
     *   dx[i,k]/dt = -p_i * x[i,k] + hitRate[k] * x[i,k+1] / m[k]
     *   dx[i,k+1]/dt = p_i * x[i,k] - hitRate[k] * x[i,k+1] / m[k]</p>
     *
     * @param x State vector of dimension modelDimension.
     * @return Drift vector dx/dt of same dimension.
     */
    public double[] drift(double[] x) {
        double[] hitRates = new double[numberOfLists + 1];
        for (int k = 0; k <= numberOfLists; k++) {
            hitRates[k] = hitRate(x, k);
        }
        double[] dX = new double[modelDimension];
        for (int i = 0; i < numberOfItems; i++) {
            for (int k = 0; k < numberOfLists; k++) {
                double flow = p[i] * x[index(i, k)] - hitRates[k] * x[index(i, k + 1)] / m[k];
                dX[index(i, k)] -= flow;
                dX[index(i, k + 1)] += flow;
            }
        }
        return dX;
    }

    /**
     * General RANDOM(m) mean-field drift honouring the per-item access graph
     * {@link #itemGraph}. Miss admission is weighted by row 0 of each item's
     * graph, hit promotion by row 1+i, and a uniformly random occupant of the
     * target list is displaced (evicted out on a miss, swapped to the source
     * list on a hit) -- the exact RR sample-path semantics
     * (State.afterEventCache). Reduces to {@link #drift(double[])} for the
     * linear chain.
     *
     * @param x state vector of dimension modelDimension.
     * @return drift vector dx/dt of the same dimension.
     */
    public double[] driftGraph(double[] x) {
        int n = numberOfItems, h = numberOfLists;
        double[] xc = new double[modelDimension];
        for (int a = 0; a < modelDimension; a++) {
            double v = x[a];
            xc[a] = v < 0.0 ? 0.0 : (v > 1.0 ? 1.0 : v);
        }
        // A[s][i]: aggregate insertion/promotion into list i from source s (s=0 miss)
        double[][] A = new double[h + 1][h + 1];
        for (int s = 0; s <= h; s++) {
            for (int j = 0; j < n; j++) {
                double xjs = xc[index(j, s)];
                if (xjs == 0.0) {
                    continue;
                }
                double[][] gj = itemGraph[j];
                for (int i = 1; i <= h; i++) {
                    A[s][i] += p[j] * xjs * gj[s][i];
                }
            }
        }
        double[] dX = new double[modelDimension];
        for (int k = 0; k < n; k++) {
            double outk = xc[index(k, 0)];
            double[][] gk = itemGraph[k];
            for (int i = 1; i <= h; i++) {
                double xki = xc[index(k, i)];
                double infl = p[k] * outk * gk[0][i];
                for (int s = 1; s < i; s++) {
                    infl += p[k] * xc[index(k, s)] * gk[s][i];
                }
                for (int b = i + 1; b <= h; b++) {
                    infl += A[i][b] * xc[index(k, b)] / m[b - 1];
                }
                double outfl = p[k] * xki * (1.0 - gk[i][i]);
                double disp = 0.0;
                for (int s = 0; s < i; s++) {
                    disp += A[s][i];
                }
                outfl += disp * xki / m[i - 1];
                dX[index(k, i)] += infl - outfl;
            }
            double acc = 0.0;
            for (int i = 1; i <= h; i++) {
                acc += dX[index(k, i)];
            }
            dX[index(k, 0)] = -acc;
        }
        return dX;
    }

    /**
     * Plain mean-field fixed point of the general (accost-aware) drift by LSODA
     * integration. Requires {@link #setItemGraph(double[][][])} to have been set.
     *
     * @param tmax integration horizon.
     * @return fixed-point state vector of dimension modelDimension.
     */
    public double[] fixedPointGraph(double tmax) {
        FirstOrderDifferentialEquations ode = new FirstOrderDifferentialEquations() {
            @Override
            public int getDimension() {
                return modelDimension;
            }

            @Override
            public void computeDerivatives(double t, double[] y, double[] yDot) {
                double[] d = driftGraph(y);
                System.arraycopy(d, 0, yDot, 0, d.length);
            }
        };
        double[] result = new double[modelDimension];
        LSODA lsoda = new LSODA(1e-12, 1.0, 1e-8, 1e-10, 12, 5);
        lsoda.integrate(ode, 0.0, x0, tmax, result);
        return result;
    }

    public double[] fixedPointGraph() {
        return fixedPointGraph(20000.0);
    }

    /**
     * Compute Jacobian dF/dx at state x.
     *
     * @param x State vector of dimension modelDimension.
     * @return Jacobian matrix of shape (modelDimension, modelDimension).
     */
    public double[][] jacobian(double[] x) {
        double[] hitRates = new double[numberOfLists + 1];
        for (int k = 0; k <= numberOfLists; k++) {
            hitRates[k] = hitRate(x, k);
        }
        int dim = modelDimension;
        double[][] Fp = new double[dim][dim];

        for (int i = 0; i < numberOfItems; i++) {
            for (int k = 0; k < numberOfLists; k++) {
                int ik = index(i, k);
                int ik1 = index(i, k + 1);

                // Direct rate terms
                Fp[ik][ik] -= p[i];
                Fp[ik1][ik] += p[i];
                Fp[ik][ik1] += hitRates[k] / m[k];
                Fp[ik1][ik1] -= hitRates[k] / m[k];

                // Indirect terms via hit rate dependence on x[j,k]
                for (int j = 0; j < numberOfItems; j++) {
                    int jk = index(j, k);
                    int jk1 = index(j, k + 1);
                    // d(hitRate[k])/d(x[j,k]) = p[j], affects x[i,k+1] terms
                    Fp[ik][jk1] -= p[i] * x[ik] / m[k];
                    Fp[ik1][jk1] += p[i] * x[ik] / m[k];
                    Fp[ik][jk] += p[j] * x[ik1] / m[k];
                    Fp[ik1][jk] -= p[j] * x[ik1] / m[k];
                }
            }
        }
        return Fp;
    }

    /**
     * Compute Hessian d^2F/dx^2.
     *
     * <p>The Hessian is constant (drift is quadratic in x), so the x argument
     * is unused but kept for interface consistency.</p>
     *
     * @param x State vector (unused - Hessian is constant for this model).
     * @return Hessian tensor of shape (modelDimension, modelDimension, modelDimension).
     */
    public double[][][] hessian(double[] x) {
        int dim = modelDimension;
        double[][][] Fpp = new double[dim][dim][dim];

        for (int i = 0; i < numberOfItems; i++) {
            for (int k = 0; k < numberOfLists; k++) {
                int ik = index(i, k);
                int ik1 = index(i, k + 1);
                for (int j = 0; j < numberOfItems; j++) {
                    if (j != i) {
                        int jk = index(j, k);
                        int jk1 = index(j, k + 1);
                        double pjOverMk = p[j] / m[k];
                        double piOverMk = p[i] / m[k];

                        // d^2 F[ik] / (d x[jk] d x[ik1]) = p[j]/m[k]
                        Fpp[ik][jk][ik1] += pjOverMk;
                        Fpp[ik][ik1][jk] += pjOverMk;
                        // d^2 F[ik] / (d x[jk1] d x[ik]) = -p[i]/m[k]
                        Fpp[ik][jk1][ik] += -piOverMk;
                        Fpp[ik][ik][jk1] += -piOverMk;
                        // Symmetric for ik1
                        Fpp[ik1][jk][ik1] -= pjOverMk;
                        Fpp[ik1][ik1][jk] -= pjOverMk;
                        Fpp[ik1][jk1][ik] -= -piOverMk;
                        Fpp[ik1][ik][jk1] -= -piOverMk;
                    }
                }
            }
        }
        return Fpp;
    }

    /**
     * Compute noise intensity matrix Q(x) for the DDPP.
     *
     * <p>Q[a,b] = sum_ell ell[a] * ell[b] * beta_ell(x) where each transition ell
     * is a swap between items i and j across lists k and k+1.</p>
     *
     * @param x State vector of dimension modelDimension.
     * @return Noise matrix of shape (modelDimension, modelDimension).
     */
    public double[][] noiseMatrix(double[] x) {
        int dim = modelDimension;
        double[][] Q = new double[dim][dim];

        int[] indices = new int[4];
        int[] signs = {-1, 1, 1, -1};

        for (int i = 0; i < numberOfItems; i++) {
            for (int k = 0; k < numberOfLists; k++) {
                for (int j = 0; j < numberOfItems; j++) {
                    double rate = p[i] * x[index(i, k)] * x[index(j, k + 1)] / m[k];
                    indices[0] = index(i, k);
                    indices[1] = index(j, k);
                    indices[2] = index(i, k + 1);
                    indices[3] = index(j, k + 1);
                    for (int ia = 0; ia < 4; ia++) {
                        for (int ib = 0; ib < 4; ib++) {
                            Q[indices[ia]][indices[ib]] += rate * signs[ia] * signs[ib];
                        }
                    }
                }
            }
        }
        return Q;
    }

    /**
     * Compute mean field fixed point by ODE integration.
     *
     * <p>Integrates dx/dt = F(x) until steady state using LSODA.</p>
     *
     * @param tmax Maximum integration time (should be large enough for convergence).
     * @return Fixed point state vector of dimension modelDimension.
     */
    public double[] fixedPoint(double tmax) {
        FirstOrderDifferentialEquations ode = new FirstOrderDifferentialEquations() {
            @Override
            public int getDimension() {
                return modelDimension;
            }

            @Override
            public void computeDerivatives(double t, double[] y, double[] yDot) {
                double[] d = drift(y);
                System.arraycopy(d, 0, yDot, 0, d.length);
            }
        };

        double[] result = new double[modelDimension];
        LSODA lsoda = new LSODA(1e-12, 1.0, 1e-8, 1e-10, 12, 5);
        lsoda.integrate(ode, 0.0, x0, tmax, result);
        return result;
    }

    /**
     * Compute mean field fixed point with default tmax=10000.
     *
     * @return Fixed point state vector.
     */
    public double[] fixedPoint() {
        return fixedPoint(10000.0);
    }

    /**
     * Stationary covariance of the occupancy process under the linear noise
     * approximation: the solution of F'(x) W + W F'(x)' + Q(x) = 0 with the same
     * {@link #jacobian} and {@link #noiseMatrix} the 1/N correction is built
     * from, so mean and covariance linearise about the identical drift.
     *
     * <p>THE SUBSPACE IS THE POINT. The Jacobian is singular twice over, because
     * the cache conserves two things: each item is in exactly one list
     * (sum_k x[i,k] = 1) and each list holds exactly its capacity
     * (sum_i x[i,k] = m[k]). Every jump is a SWAP,
     * l(i,j,k) = (e_i - e_j) tensor (e_{k+1} - e_k), so the fluctuation lives on
     * the tensor product of the zero-sum item space with the zero-sum list
     * space: the double-centred subspace, of dimension (n-1)*h. Restricting to
     * an orthonormal basis of it is exact, and it is what makes the covariance
     * of a deterministic total come out as zero. The reduction used by
     * {@link #reduceFpFppQ} does not: it drops the last item's rows and pads
     * with null vectors, which leaves the miss-indicator covariance summing to a
     * nonzero number.</p>
     *
     * @param x occupancy at which to linearise
     * @return the modelDimension-square covariance, or null when the fixed point
     *         is not exponentially stable on the reachable subspace
     */
    public double[][] lnaCovariance(double[] x) {
        double[][] Fp = jacobian(x);
        double[][] Q = noiseMatrix(x);

        double[][] Ui = centeredBasis(numberOfItems);
        double[][] Ul = centeredBasis(numberOfLists + 1);
        int nv = Ui[0].length * Ul[0].length;
        if (nv == 0) {
            return new double[modelDimension][modelDimension];
        }
        // V = kron(Ul, Ui) on the item-major flat index i + k*n
        double[][] V = new double[modelDimension][nv];
        for (int k = 0; k <= numberOfLists; k++) {
            for (int i = 0; i < numberOfItems; i++) {
                int row = index(i, k);
                for (int b = 0; b < Ul[0].length; b++) {
                    for (int a = 0; a < Ui[0].length; a++) {
                        V[row][b * Ui[0].length + a] = Ul[k][b] * Ui[i][a];
                    }
                }
            }
        }

        Matrix Vm = arrayToMatrix(V);
        Matrix Vt = Vm.transpose();
        Matrix Ar = Vt.mult(arrayToMatrix(Fp)).mult(Vm);
        Matrix Qr = Vt.mult(arrayToMatrix(symmetrize(Q))).mult(Vm);
        Qr = arrayToMatrix(symmetrize(matrixToArray(Qr)));

        // the LNA has a stationary covariance only at an exponentially stable
        // fixed point
        java.util.List<org.apache.commons.math3.complex.Complex> ev = Ar.eig();
        double maxRe = Double.NEGATIVE_INFINITY;
        for (int i = 0; i < ev.size(); i++) {
            maxRe = Math.max(maxRe, ev.get(i).getReal());
        }
        if (!(maxRe < -Math.sqrt(2.220446049250313e-16))) {
            return null;
        }

        Matrix Wr = Matrix.sylv(Ar, Ar.transpose(), Qr);
        double[][] WrArr = symmetrize(matrixToArray(Wr));
        Matrix W = Vm.mult(arrayToMatrix(WrArr)).mult(Vt);
        return symmetrize(matrixToArray(W));
    }

    /** Orthonormal basis of {u in R^n : sum(u) = 0}, n-by-(n-1), by Gram-Schmidt. */
    private static double[][] centeredBasis(int n) {
        if (n <= 1) {
            return new double[Math.max(n, 1)][0];
        }
        double[][] U = new double[n][n - 1];
        // Helmert basis: column j has j entries 1/sqrt(j(j+1)) and one -j/sqrt(j(j+1))
        for (int j = 0; j < n - 1; j++) {
            double d = Math.sqrt((j + 1.0) * (j + 2.0));
            for (int i = 0; i <= j; i++) {
                U[i][j] = 1.0 / d;
            }
            U[j + 1][j] = -(j + 1.0) / d;
        }
        return U;
    }

    /** (M + M')/2. */
    private static double[][] symmetrize(double[][] M) {
        int n = M.length;
        double[][] S = new double[n][n];
        for (int i = 0; i < n; i++) {
            for (int j = 0; j < n; j++) {
                S[i][j] = 0.5 * (M[i][j] + M[j][i]);
            }
        }
        return S;
    }

    /**
     * Integrate the plain mean-field drift over a finite window on a uniform
     * time grid, from a supplied initial occupancy.
     *
     * <p>This is the transient counterpart of {@link #fixedPoint()}: the fixed
     * point drives the same drift to steady state, whereas this returns the full
     * trajectory. It is the order-0 (no 1/N correction) transient used to carry
     * the cache mean occupancy across environment switches. Mirrors the MATLAB
     * {@code cache_miss_rmf} tspan path.</p>
     *
     * @param time    end time of the window (start is 0).
     * @param nPoints number of uniform grid points (>= 2).
     * @param xinit   initial occupancy of dimension modelDimension; when null,
     *                the default first-m-in-list initial state is used.
     * @return {@code Object[]{ T (double[nPoints]), X (double[nPoints][dim]) }}.
     */
    public Object[] driftTrajectory(double time, int nPoints, double[] xinit) {
        final int dim = modelDimension;
        double[] T = new double[nPoints];
        for (int i = 0; i < nPoints; i++) {
            T[i] = time * i / (nPoints - 1);
        }

        FirstOrderDifferentialEquations ode = new FirstOrderDifferentialEquations() {
            @Override
            public int getDimension() {
                return dim;
            }

            @Override
            public void computeDerivatives(double t, double[] y, double[] yDot) {
                double[] d = drift(y);
                System.arraycopy(d, 0, yDot, 0, d.length);
            }
        };

        double[][] X = new double[nPoints][dim];
        double[] state = new double[dim];
        if (xinit != null) {
            System.arraycopy(xinit, 0, state, 0, dim);
        } else {
            System.arraycopy(x0, 0, state, 0, dim);
        }
        System.arraycopy(state, 0, X[0], 0, dim);
        for (int i = 1; i < nPoints; i++) {
            double[] result = new double[dim];
            LSODA lsoda = new LSODA(1e-12, 1.0, 1e-8, 1e-10, 12, 5);
            lsoda.integrate(ode, T[i - 1], state, T[i], result);
            System.arraycopy(result, 0, X[i], 0, dim);
            System.arraycopy(result, 0, state, 0, dim);
        }
        return new Object[]{T, X};
    }

    /**
     * Compute dimension reduction matrices for the singular Jacobian.
     *
     * <p>The Jacobian is singular because item populations are conserved
     * (sum over lists for each item = 1). This method finds a change of
     * basis that separates the rank-deficient directions.</p>
     *
     * @param Fp Jacobian matrix at fixed point.
     * @return DimensionReduction containing P, Pinv, rank.
     */
    public DimensionReduction dimensionReduction(double[][] Fp) {
        int dim = modelDimension;

        // Compute rank of Fp using Matrix
        Matrix FpMat = arrayToMatrix(Fp);
        int rank = FpMat.rank();

        // Build change-of-basis: first rank rows are independent coordinates,
        // remaining rows span the null space of Fp
        double[][] C = new double[dim][dim];
        int d = 0;
        for (int lIdx = 0; lIdx <= numberOfLists; lIdx++) {
            for (int i = 0; i < numberOfItems - 1; i++) {
                C[d][index(i, lIdx)] = 1.0;
                d++;
            }
        }

        // SVD of Fp
        Ret.SVD svdResult = FpMat.svd();
        Matrix U = svdResult.u;

        // C[rank:, :] = U.T[rank:, :]
        for (int row = rank; row < dim; row++) {
            for (int col = 0; col < dim; col++) {
                C[row][col] = U.get(col, row); // U.T[row][col] = U[col][row]
            }
        }

        // Cinv = inv(C)
        Matrix CMat = arrayToMatrix(C);
        Matrix CinvMat = CMat.inv();
        double[][] Cinv = matrixToArray(CinvMat);

        return new DimensionReduction(C, Cinv, rank);
    }

    /**
     * Apply dimension reduction to Fp, Fpp, Q.
     *
     * <p>Projects the Jacobian, Hessian, and noise matrix onto the
     * non-singular subspace identified by dimension reduction.</p>
     *
     * @param Fp  Jacobian matrix.
     * @param Fpp Hessian tensor.
     * @param Q   Noise matrix.
     * @return ReducedSystem with reduced Fp, Fpp, Q, and basis matrices.
     */
    public ReducedSystem reduceFpFppQ(double[][] Fp, double[][][] Fpp, double[][] Q) {
        DimensionReduction dr = dimensionReduction(Fp);
        double[][] P = dr.P;
        double[][] Pinv = dr.Pinv;
        int rank = dr.rank;
        int dim = modelDimension;

        // Fp_r = (P @ Fp @ Pinv)[:rank, :rank]
        double[][] PFp = matMul(P, Fp);
        double[][] PFpPinv = matMul(PFp, Pinv);
        double[][] Fp_r = subMatrix(PFpPinv, rank, rank);

        // Fpp_r[a,b,c] = sum_{i,j,k} P[a,i] * Fpp[i,j,k] * Pinv[j,b] * Pinv[k,c]
        // Take only [:rank, :rank, :rank]
        double[][][] Fpp_r = new double[rank][rank][rank];
        // First contract: tmp1[a,j,k] = sum_i P[a,i] * Fpp[i,j,k]
        // Then: tmp2[a,b,k] = sum_j tmp1[a,j,k] * Pinv[j,b]
        // Then: Fpp_r[a,b,c] = sum_k tmp2[a,b,k] * Pinv[k,c]
        for (int a = 0; a < rank; a++) {
            for (int b = 0; b < rank; b++) {
                for (int c = 0; c < rank; c++) {
                    double val = 0.0;
                    for (int ii = 0; ii < dim; ii++) {
                        double Pa_ii = P[a][ii];
                        if (Pa_ii == 0.0) continue;
                        for (int jj = 0; jj < dim; jj++) {
                            double Pinv_jj_b = Pinv[jj][b];
                            if (Pinv_jj_b == 0.0) continue;
                            double Pa_Pinv = Pa_ii * Pinv_jj_b;
                            for (int kk = 0; kk < dim; kk++) {
                                double Fpp_ijk = Fpp[ii][jj][kk];
                                if (Fpp_ijk == 0.0) continue;
                                val += Pa_Pinv * Fpp_ijk * Pinv[kk][c];
                            }
                        }
                    }
                    Fpp_r[a][b][c] = val;
                }
            }
        }

        // Q_r = (P @ Q @ P.T)[:rank, :rank]
        double[][] PQ = matMul(P, Q);
        double[][] PT = transpose(P);
        double[][] PQPT = matMul(PQ, PT);
        double[][] Q_r = subMatrix(PQPT, rank, rank);

        return new ReducedSystem(Fp_r, Fpp_r, Q_r, P, Pinv, rank);
    }

    /**
     * Expand reduced V back to full dimension.
     *
     * @param V_r  Reduced correction vector (rank).
     * @param Pinv Inverse basis matrix (dim x dim).
     * @param rank Rank of the non-singular subspace.
     * @return Full-dimension correction vector (modelDimension).
     */
    public double[] expandV(double[] V_r, double[][] Pinv, int rank) {
        int dim = modelDimension;
        double[] V = new double[dim];
        // V = Pinv[:, :rank] @ V_r
        for (int i = 0; i < dim; i++) {
            double val = 0.0;
            for (int j = 0; j < rank; j++) {
                val += Pinv[i][j] * V_r[j];
            }
            V[i] = val;
        }
        return V;
    }

    /**
     * Expand reduced W back to full dimension.
     *
     * @param W_r  Reduced covariance matrix (rank x rank).
     * @param Pinv Inverse basis matrix (dim x dim).
     * @param rank Rank of the non-singular subspace.
     * @return Full-dimension covariance matrix (modelDimension x modelDimension).
     */
    public double[][] expandW(double[][] W_r, double[][] Pinv, int rank) {
        int dim = modelDimension;
        // W = Pinv[:, :rank] @ W_r @ Pinv.T[:rank, :]
        // = Pinv[:, :rank] @ W_r @ (Pinv[:, :rank]).T
        // Step 1: tmp[i][j] = sum_k Pinv[i][k] * W_r[k][j] for k in 0..rank-1
        double[][] tmp = new double[dim][rank];
        for (int i = 0; i < dim; i++) {
            for (int j = 0; j < rank; j++) {
                double val = 0.0;
                for (int k = 0; k < rank; k++) {
                    val += Pinv[i][k] * W_r[k][j];
                }
                tmp[i][j] = val;
            }
        }
        // Step 2: W[i][j] = sum_k tmp[i][k] * Pinv[j][k] for k in 0..rank-1
        // (Pinv.T[:rank, :] has element [k][j] = Pinv[j][k])
        double[][] W = new double[dim][dim];
        for (int i = 0; i < dim; i++) {
            for (int j = 0; j < dim; j++) {
                double val = 0.0;
                for (int k = 0; k < rank; k++) {
                    val += tmp[i][k] * Pinv[j][k];
                }
                W[i][j] = val;
            }
        }
        return W;
    }

    /**
     * Compute refined mean field steady-state expansion.
     *
     * <p>Computes the mean field fixed point pi and the 1/N correction V
     * using the Lyapunov equation approach with dimension reduction to
     * handle the singular Jacobian (due to per-item conservation constraints).</p>
     *
     * <p>The refined approximation for a system of N items is:
     *   E[X] ~ pi + V/N + O(1/N^2)</p>
     *
     * @param order Expansion order (0 = plain mean field, 1 = with 1/N correction).
     * @return Object array {pi (double[]), V (double[]), W (double[][])}.
     */
    public Object[] meanFieldExpansionSteadyState(int order) {
        double[] pi = fixedPoint();

        if (order == 0) {
            double[] V = new double[modelDimension];
            double[][] W = new double[modelDimension][modelDimension];
            return new Object[]{pi, V, W};
        }

        double[][] Fp = jacobian(pi);
        double[][][] Fpp = hessian(pi);
        double[][] Q = noiseMatrix(pi);

        // Dimension reduction: project onto non-singular subspace
        ReducedSystem rs = reduceFpFppQ(Fp, Fpp, Q);
        double[][] Fp_r = rs.Fp_r;
        double[][][] Fpp_r = rs.Fpp_r;
        double[][] Q_r = rs.Q_r;
        double[][] P_basis = rs.P;
        double[][] Pinv = rs.Pinv;
        int rank = rs.rank;

        // Solve Lyapunov equation in reduced space: Fp_r @ W_r + W_r @ Fp_r.T + Q_r = 0,
        // i.e. Fp_r * W_r + W_r * Fp_r^T = -Q_r. Matrix.sylv(A, B, C) solves
        // A*X + X*B = -C, so pass C = Q_r directly (it is negated internally);
        // pre-negating here would double-negate and flip the sign of the 1/N
        // correction. Mirrors Python scipy.linalg.solve_lyapunov(Fp_r, -Q_r).
        Matrix Fp_r_mat = arrayToMatrix(Fp_r);
        Matrix Q_r_mat = arrayToMatrix(Q_r);
        Matrix W_r_mat = Matrix.sylv(Fp_r_mat, Fp_r_mat.transpose(), Q_r_mat);
        double[][] W_r = matrixToArray(W_r_mat);

        // Hessian contraction: C_r[a] = sum_{b,c} Fpp_r[a,b,c] * W_r[b,c]
        double[] C_r = new double[rank];
        for (int a = 0; a < rank; a++) {
            double val = 0.0;
            for (int b = 0; b < rank; b++) {
                for (int c = 0; c < rank; c++) {
                    val += Fpp_r[a][b][c] * W_r[b][c];
                }
            }
            C_r[a] = val;
        }

        // First-order correction: V_r = -solve(Fp_r, C_r / 2)
        Matrix C_r_half = new Matrix(rank, 1);
        for (int i = 0; i < rank; i++) {
            C_r_half.set(i, 0, C_r[i] / 2.0);
        }
        Matrix V_r_mat = new Matrix(rank, 1);
        Matrix.solveSafe(Fp_r_mat, C_r_half, V_r_mat);
        // Negate
        double[] V_r = new double[rank];
        for (int i = 0; i < rank; i++) {
            V_r[i] = -V_r_mat.get(i, 0);
        }

        // Expand back to full dimension
        double[] V = expandV(V_r, Pinv, rank);
        double[][] W = expandW(W_r, Pinv, rank);

        return new Object[]{pi, V, W};
    }

    /**
     * Compute refined mean field steady-state expansion with default order=1.
     *
     * @return Object array {pi (double[]), V (double[]), W (double[][])}.
     */
    public Object[] meanFieldExpansionSteadyState() {
        return meanFieldExpansionSteadyState(1);
    }

    /**
     * Compute refined mean field transient expansion.
     *
     * <p>Integrates the coupled ODE system for (X, V, W) where:
     * <ul>
     *   <li>X(t): mean field trajectory</li>
     *   <li>V(t): 1/N correction trajectory</li>
     *   <li>W(t): covariance trajectory</li>
     * </ul>
     *
     * @param time     Maximum integration time.
     * @param nPoints  Number of output time points.
     * @param order    Expansion order (0 or 1).
     * @return Object array {T (double[]), X (double[][]), V (double[][]), W (double[][][])}.
     */
    public Object[] meanFieldExpansionTransient(double time, int nPoints, int order) {
        int dim = modelDimension;
        double[] T = new double[nPoints];
        for (int i = 0; i < nPoints; i++) {
            T[i] = time * i / (nPoints - 1);
        }

        if (order == 0) {
            // Simple mean field ODE
            FirstOrderDifferentialEquations ode = new FirstOrderDifferentialEquations() {
                @Override
                public int getDimension() {
                    return dim;
                }

                @Override
                public void computeDerivatives(double t, double[] y, double[] yDot) {
                    double[] d = drift(y);
                    System.arraycopy(d, 0, yDot, 0, d.length);
                }
            };

            double[][] X = new double[nPoints][dim];
            double[][] V = new double[nPoints][dim];
            double[][][] W = new double[nPoints][dim][dim];

            // Integrate to each time point
            double[] state = new double[dim];
            System.arraycopy(x0, 0, state, 0, dim);
            X[0] = new double[dim];
            System.arraycopy(x0, 0, X[0], 0, dim);

            for (int i = 1; i < nPoints; i++) {
                double[] result = new double[dim];
                LSODA lsoda = new LSODA(1e-12, 1.0, 1e-6, 1e-10, 12, 5);
                lsoda.integrate(ode, T[i - 1], state, T[i], result);
                X[i] = new double[dim];
                System.arraycopy(result, 0, X[i], 0, dim);
                System.arraycopy(result, 0, state, 0, dim);
            }
            return new Object[]{T, X, V, W};
        }

        // Coupled ODE: state = [X (dim), V (dim), W_flat (dim^2)]
        int totalDim = dim + dim + dim * dim;
        FirstOrderDifferentialEquations coupledOde = new FirstOrderDifferentialEquations() {
            @Override
            public int getDimension() {
                return totalDim;
            }

            @Override
            public void computeDerivatives(double t, double[] y, double[] yDot) {
                // Extract X, V, W from y
                double[] xState = new double[dim];
                System.arraycopy(y, 0, xState, 0, dim);
                double[] vState = new double[dim];
                System.arraycopy(y, dim, vState, 0, dim);
                double[][] wState = new double[dim][dim];
                for (int ii = 0; ii < dim; ii++) {
                    System.arraycopy(y, 2 * dim + ii * dim, wState[ii], 0, dim);
                }

                double[] F = drift(xState);
                double[][] Fp = jacobian(xState);
                double[][][] Fpp = hessian(xState);
                double[][] Qmat = noiseMatrix(xState);

                // dx = F
                System.arraycopy(F, 0, yDot, 0, dim);

                // dv = Fp @ v + 0.5 * tensordot(Fpp, w)
                for (int a = 0; a < dim; a++) {
                    double dvA = 0.0;
                    for (int b = 0; b < dim; b++) {
                        dvA += Fp[a][b] * vState[b];
                    }
                    // Hessian contraction: 0.5 * sum_{b,c} Fpp[a,b,c] * W[b,c]
                    double contraction = 0.0;
                    for (int b = 0; b < dim; b++) {
                        for (int c = 0; c < dim; c++) {
                            contraction += Fpp[a][b][c] * wState[b][c];
                        }
                    }
                    yDot[dim + a] = dvA + 0.5 * contraction;
                }

                // dw = Fp @ W + W @ Fp' + Q
                for (int ii = 0; ii < dim; ii++) {
                    for (int jj = 0; jj < dim; jj++) {
                        double val = Qmat[ii][jj];
                        for (int kk = 0; kk < dim; kk++) {
                            val += Fp[ii][kk] * wState[kk][jj] + wState[ii][kk] * Fp[jj][kk];
                        }
                        yDot[2 * dim + ii * dim + jj] = val;
                    }
                }
            }
        };

        // Initial condition: X = x0, V = 0, W = 0
        double[] y0 = new double[totalDim];
        System.arraycopy(x0, 0, y0, 0, dim);

        // Integrate step by step
        double[][] X = new double[nPoints][dim];
        double[][] Vout = new double[nPoints][dim];
        double[][][] Wout = new double[nPoints][dim][dim];

        System.arraycopy(x0, 0, X[0], 0, dim);

        double[] state = new double[totalDim];
        System.arraycopy(y0, 0, state, 0, totalDim);

        for (int i = 1; i < nPoints; i++) {
            double[] result = new double[totalDim];
            LSODA lsoda = new LSODA(1e-12, 1.0, 1e-6, 1e-10, 12, 5);
            lsoda.integrate(coupledOde, T[i - 1], state, T[i], result);

            X[i] = new double[dim];
            System.arraycopy(result, 0, X[i], 0, dim);
            Vout[i] = new double[dim];
            System.arraycopy(result, dim, Vout[i], 0, dim);
            for (int ii = 0; ii < dim; ii++) {
                System.arraycopy(result, 2 * dim + ii * dim, Wout[i][ii], 0, dim);
            }
            System.arraycopy(result, 0, state, 0, totalDim);
        }

        return new Object[]{T, X, Vout, Wout};
    }

    /**
     * Compute refined mean field transient expansion with default parameters.
     *
     * @return Object array {T, X, V, W}.
     */
    public Object[] meanFieldExpansionTransient() {
        return meanFieldExpansionTransient(50.0, 200, 1);
    }

    /**
     * Compute hit rates for all lists.
     *
     * @param x State vector of dimension modelDimension.
     * @return Array of shape (numberOfLists + 1) with hit rate per list.
     */
    public double[] hitRatesAll(double[] x) {
        double[] rates = new double[numberOfLists + 1];
        for (int k = 0; k <= numberOfLists; k++) {
            rates[k] = hitRate(x, k);
        }
        return rates;
    }

    // ---- Getters ----

    public int getNumberOfItems() {
        return numberOfItems;
    }

    public int getNumberOfLists() {
        return numberOfLists;
    }

    public int getModelDimension() {
        return modelDimension;
    }

    public double[] getP() {
        double[] copy = new double[p.length];
        System.arraycopy(p, 0, copy, 0, p.length);
        return copy;
    }

    public int[] getM() {
        int[] copy = new int[m.length];
        System.arraycopy(m, 0, copy, 0, m.length);
        return copy;
    }

    public double[] getX0() {
        double[] copy = new double[x0.length];
        System.arraycopy(x0, 0, copy, 0, x0.length);
        return copy;
    }

    // ---- Private helpers ----

    /**
     * Convert a 2D double array to a Matrix.
     */
    private static Matrix arrayToMatrix(double[][] arr) {
        int rows = arr.length;
        int cols = arr[0].length;
        Matrix mat = new Matrix(rows, cols);
        for (int i = 0; i < rows; i++) {
            for (int j = 0; j < cols; j++) {
                mat.set(i, j, arr[i][j]);
            }
        }
        return mat;
    }

    /**
     * Convert a Matrix to a 2D double array.
     */
    private static double[][] matrixToArray(Matrix mat) {
        int rows = mat.getNumRows();
        int cols = mat.getNumCols();
        double[][] arr = new double[rows][cols];
        for (int i = 0; i < rows; i++) {
            for (int j = 0; j < cols; j++) {
                arr[i][j] = mat.get(i, j);
            }
        }
        return arr;
    }

    /**
     * Matrix multiplication of two 2D arrays.
     */
    private static double[][] matMul(double[][] A, double[][] B) {
        int m = A.length;
        int n = B[0].length;
        int p = B.length;
        double[][] C = new double[m][n];
        for (int i = 0; i < m; i++) {
            for (int k = 0; k < p; k++) {
                double aik = A[i][k];
                if (aik == 0.0) continue;
                for (int j = 0; j < n; j++) {
                    C[i][j] += aik * B[k][j];
                }
            }
        }
        return C;
    }

    /**
     * Transpose a 2D array.
     */
    private static double[][] transpose(double[][] A) {
        int rows = A.length;
        int cols = A[0].length;
        double[][] T = new double[cols][rows];
        for (int i = 0; i < rows; i++) {
            for (int j = 0; j < cols; j++) {
                T[j][i] = A[i][j];
            }
        }
        return T;
    }

    /**
     * Extract top-left submatrix of size rows x cols.
     */
    private static double[][] subMatrix(double[][] A, int rows, int cols) {
        double[][] sub = new double[rows][cols];
        for (int i = 0; i < rows; i++) {
            System.arraycopy(A[i], 0, sub[i], 0, cols);
        }
        return sub;
    }
}
