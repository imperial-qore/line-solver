package jline.solvers.mam.handlers;

import java.util.ArrayList;
import java.util.Collections;
import java.util.List;

import jline.api.da.Da_fpi;
import jline.api.mc.Ctmc_makeinfgen;
import jline.api.mc.Ctmc_solve;
import jline.util.Pair;
import jline.util.matrix.Matrix;

public final class Solver_mam_inap {
    private Solver_mam_inap() {}

    public static INAPResult solver_mam_inap(RCATModel rcat, double tol, int maxiter, int seed, String method) {
        Matrix[][] R = rcat.R;
        Matrix AP = rcat.AP;
        int[] N = rcat.N;

        int numActions = rcat.actionMap.size();
        int numProcesses = N.length;

        if (numProcesses == 0) {
            return new INAPResult(new Matrix(0, 1), Collections.<Matrix>emptyList(), Collections.<Matrix>emptyList(), 0);
        }

        Matrix[] Aa = new Matrix[numActions];
        Matrix[] Pb = new Matrix[numActions];
        for (int a = 0; a < numActions; a++) {
            Aa[a] = R[a][0];
            Pb[a] = R[a][1];
        }
        Matrix[] L = new Matrix[numProcesses];
        for (int k = 0; k < numProcesses; k++) {
            L[k] = R[numActions][k];
        }

        int[] ACT = new int[numActions];
        int[] PSV = new int[numActions];
        for (int a = 0; a < numActions; a++) {
            ACT[a] = (int) AP.get(a, 0);
            PSV[a] = (int) AP.get(a, 1);
        }

        // Processes whose local dynamics are not birth-death. A catastrophe or
        // a batch-removal signal moves the process down by more than one
        // level, so the reversed rate of an action of that process is
        // state-dependent and the mean-of-ratios estimator of INAP has no
        // fixed point of RCAT type (it overestimates the departure rate, e.g.
        // 2x on a tandem G-network with catastrophes). For those processes the
        // action rate is instead set by rate conservation,
        // x(a) = sum_ij Aa(i,j)*pi(i), which is the actual departure rate of
        // the active process (the INAP+ estimator). Birth-death processes,
        // including the classic single-removal negative customer, keep the
        // standard INAP estimator, which is exact for them.
        final boolean[] notBirthDeath = new boolean[numProcesses];
        for (int k = 0; k < numProcesses; k++) {
            if (L[k] == null) continue;
            for (int n = 0; n < N[k]; n++) {
                for (int m = 0; m < N[k]; m++) {
                    if (L[k].get(n, m) > 0 && (m < n - 1 || m > n + 1)) {
                        notBirthDeath[k] = true;
                    }
                }
            }
        }

        // see _kb/06-solver-catalog.md for rationale
        Matrix x = new Matrix(numActions, 1);
        for (int a = 0; a < numActions; a++) {
            x.set(a, 0, (a + 1.0) / (numActions + 1.0));
        }

        Pair<List<Matrix>, List<Matrix>> initial = computeEquilibrium(x, Aa, Pb, L, ACT, PSV, numProcesses, numActions, N);
        List<Matrix> pi = initial.getLeft();
        List<Matrix> Q = initial.getRight();

        // reversed-rate fixed point on the isolated-component equilibria,
        // driven by the generic DA driver
        final List<List<Matrix>> QBox = new ArrayList<List<Matrix>>();
        QBox.add(Q);
        Da_fpi.Options<List<Matrix>> fpopts = new Da_fpi.Options<List<Matrix>>(
                maxiter, tol, new Da_fpi.Norm<List<Matrix>>() {
            @Override
            public double eval(List<Matrix> xnew, List<Matrix> xref) {
                double maxErr = 0.0;
                for (int k = 0; k < numProcesses; k++) {
                    double errK = 0.0;
                    for (int j = 0; j < N[k]; j++) {
                        errK += Math.abs(xnew.get(k).get(0, j) - xref.get(k).get(0, j));
                    }
                    if (errK > maxErr) maxErr = errK;
                }
                return maxErr;
            }
        });
        Da_fpi.Result<List<Matrix>> fpres = Da_fpi.run(new Da_fpi.Sweep<List<Matrix>>() {
            @Override
            public Da_fpi.SweepResult<List<Matrix>> sweep(List<Matrix> piprev, int iter) {
            for (int a = 0; a < numActions; a++) {
                int k = ACT[a];
                Matrix activeMatrix = Aa[a];
                if (activeMatrix == null) continue;
                Matrix piK = piprev.get(k);

                if ("inapplus".equals(method) || notBirthDeath[k]) {
                    // inapplus: LAMBDA(i,j) = Aa(i,j) * pi(i)
                    double lambdaSum = 0.0;
                    for (int i = 0; i < N[k]; i++) {
                        for (int j = 0; j < N[k]; j++) {
                            double aij = activeMatrix.get(i, j);
                            if (aij > 0) {
                                double piI = piK.get(0, i);
                                lambdaSum += aij * piI;
                            }
                        }
                    }
                    if (lambdaSum > 0) {
                        x.set(a, 0, lambdaSum);
                    }
                } else {
                    // inap: LAMBDA(i,j) = Aa(i,j) * pi(i) / pi(j)
                    List<Double> lambdaVec = new ArrayList<Double>();
                    for (int i = 0; i < N[k]; i++) {
                        for (int j = 0; j < N[k]; j++) {
                            double aij = activeMatrix.get(i, j);
                            double piJ = piK.get(0, j);
                            if (aij > 0 && piJ > 0) {
                                double piI = piK.get(0, i);
                                lambdaVec.add(aij * piI / piJ);
                            }
                        }
                    }
                    if (!lambdaVec.isEmpty()) {
                        double sum = 0.0;
                        for (double v : lambdaVec) sum += v;
                        x.set(a, 0, sum / lambdaVec.size());
                    }
                }
            }
                Pair<List<Matrix>, List<Matrix>> result = computeEquilibrium(x, Aa, Pb, L, ACT, PSV, numProcesses, numActions, N);
                QBox.set(0, result.getRight());
                return new Da_fpi.SweepResult<List<Matrix>>(result.getLeft(), piprev);
            }
        }, pi, fpopts);
        pi = fpres.x;
        Q = QBox.get(0);
        int iter = fpres.it;
        if (!fpres.converged) iter++;
        return new INAPResult(x, pi, Q, iter);
    }

    public static INAPResult solver_mam_inap(RCATModel rcat, double tol, int maxiter, int seed) {
        return solver_mam_inap(rcat, tol, maxiter, seed, "inap");
    }

    public static INAPResult solver_mam_inap(RCATModel rcat, double tol, int maxiter) {
        return solver_mam_inap(rcat, tol, maxiter, 0, "inap");
    }

    public static INAPResult solver_mam_inap(RCATModel rcat) {
        return solver_mam_inap(rcat, 1e-6, 1000, 0, "inap");
    }

    /**
     * Matrix-geometric INAP ('inapinf'): the isolated OPEN components are solved
     * directly on their infinite state space by a scalar matrix-geometric
     * (QBD / catastrophe) decomposition -- the marginal is geometric
     * pi_n = (1-rho) rho^n with rho the sub-unit root of the QBD characteristic
     * equation, and any catastrophe drain to the empty state folds into the
     * local outflow. Closed components remain finite. The reversed rate is
     * updated by the weighted-mean formula Eq. (4) evaluated in closed form on
     * the geometric tail, and the RCAT product-form residual (Remark 2) is
     * returned. Mirrors solver_mam_ag.m (inap_inf).
     *
     * Reference: A. Marin, S. Rota Bulo, S. Balsamo, "A Numerical Algorithm for
     * the Decomposition of Cooperating Structured Markov Processes", MASCOTS 2012.
     */
    public static INAPResult solver_mam_inapinf(RCATModel rcat, boolean[] isOpenProc,
                                                double tol, int maxiter) {
        Matrix[][] R = rcat.R;
        Matrix AP = rcat.AP;
        int[] N = rcat.N;

        int numActions = rcat.actionMap.size();
        int numProcesses = N.length;

        if (numProcesses == 0) {
            return new INAPResult(new Matrix(0, 1), Collections.<Matrix>emptyList(),
                    Collections.<Matrix>emptyList(), 0,
                    new double[0], new boolean[0], 0.0);
        }

        Matrix[] Aa = new Matrix[numActions];
        Matrix[] Pb = new Matrix[numActions];
        for (int a = 0; a < numActions; a++) {
            Aa[a] = R[a][0];
            Pb[a] = R[a][1];
        }
        Matrix[] L = new Matrix[numProcesses];
        for (int k = 0; k < numProcesses; k++) {
            L[k] = R[numActions][k];
        }

        int[] ACT = new int[numActions];
        int[] PSV = new int[numActions];
        for (int a = 0; a < numActions; a++) {
            ACT[a] = (int) AP.get(a, 0);
            PSV[a] = (int) AP.get(a, 1);
        }

        // Active-transition row sums (rate of the active label out of each state)
        double[][] aRowSum = new double[numActions][];
        for (int a = 0; a < numActions; a++) {
            if (Aa[a] == null) { aRowSum[a] = new double[0]; continue; }
            int rows = Aa[a].getNumRows();
            int cols = Aa[a].getNumCols();
            double[] rs = new double[rows];
            for (int i = 0; i < rows; i++) {
                double s = 0.0;
                for (int j = 0; j < cols; j++) s += Aa[a].get(i, j);
                rs[i] = s;
            }
            aRowSum[a] = rs;
        }

        // Deterministic initial guess (reproducible across back-ends).
        Matrix x = new Matrix(numActions, 1);
        for (int a = 0; a < numActions; a++) {
            x.set(a, 0, (a + 1.0) / (numActions + 1.0));
        }

        QbdEquilibrium eq = computeEquilibriumQbd(x, Aa, Pb, L, ACT, PSV, numProcesses, numActions, N, isOpenProc);
        List<Matrix> pi = eq.pi;
        List<Matrix> Q = eq.Q;
        double[] rhoProc = eq.rhoProc;
        boolean[] isGeomProc = eq.isGeomProc;

        // reversed-rate fixed point on the isolated-component equilibria
        // (matrix-geometric variant), driven by the generic DA driver
        final List<List<Matrix>> QBox = new ArrayList<List<Matrix>>();
        QBox.add(Q);
        final double[][] rhoBox = new double[][]{rhoProc};
        final boolean[][] geomBox = new boolean[][]{isGeomProc};
        Da_fpi.Options<List<Matrix>> fpopts = new Da_fpi.Options<List<Matrix>>(
                maxiter, tol, new Da_fpi.Norm<List<Matrix>>() {
            @Override
            public double eval(List<Matrix> xnew, List<Matrix> xref) {
                double maxErr = 0.0;
                for (int k = 0; k < numProcesses; k++) {
                    int m = Math.min(xnew.get(k).getNumCols(), xref.get(k).getNumCols());
                    double errK = 0.0;
                    for (int j = 0; j < m; j++) errK += Math.abs(xnew.get(k).get(0, j) - xref.get(k).get(0, j));
                    if (errK > maxErr) maxErr = errK;
                }
                return maxErr;
            }
        });
        Da_fpi.Result<List<Matrix>> fpres = Da_fpi.run(new Da_fpi.Sweep<List<Matrix>>() {
            @Override
            public Da_fpi.SweepResult<List<Matrix>> sweep(List<Matrix> piprev, int iter) {
            // Reversed-rate update, Eq. (4): x_l = pi^(alpha_l) T^(l) e.
            for (int a = 0; a < numActions; a++) {
                int k = ACT[a];
                if (Aa[a] == null) continue;
                if (geomBox[0][k]) {
                    int occIdx = Math.min(1, N[k] - 1);
                    x.set(a, 0, aRowSum[a][occIdx] * rhoBox[0][k]);
                } else {
                    Matrix piK = piprev.get(k);
                    double v = 0.0;
                    for (int i = 0; i < N[k]; i++) v += piK.get(0, i) * aRowSum[a][i];
                    x.set(a, 0, v);
                }
            }
                QbdEquilibrium eq = computeEquilibriumQbd(x, Aa, Pb, L, ACT, PSV, numProcesses, numActions, N, isOpenProc);
                QBox.set(0, eq.Q);
                rhoBox[0] = eq.rhoProc;
                geomBox[0] = eq.isGeomProc;
                return new Da_fpi.SweepResult<List<Matrix>>(eq.pi, piprev);
            }
        }, pi, fpopts);
        pi = fpres.x;
        Q = QBox.get(0);
        rhoProc = rhoBox[0];
        isGeomProc = geomBox[0];
        int iter = fpres.it;
        if (!fpres.converged) iter++;

        // RCAT product-form residual (Remark 2):
        // max_l || pi^(alpha_l) (x_l I - T^(l)) ||, zero iff exact product form.
        double rcatRes = 0.0;
        for (int a = 0; a < numActions; a++) {
            if (Aa[a] == null) continue;
            int k = ACT[a];
            Matrix piK = pi.get(k);
            double xa = x.get(a, 0);
            int cols = Aa[a].getNumCols();
            double nrm2 = 0.0;
            for (int j = 0; j < cols; j++) {
                double comp = xa * piK.get(0, j);
                for (int i = 0; i < N[k]; i++) comp -= piK.get(0, i) * Aa[a].get(i, j);
                nrm2 += comp * comp;
            }
            double nrm = Math.sqrt(nrm2);
            if (nrm > rcatRes) rcatRes = nrm;
        }

        return new INAPResult(x, pi, Q, iter, rhoProc, isGeomProc, rcatRes);
    }

    /** Holder for the matrix-geometric equilibrium of all components. */
    private static final class QbdEquilibrium {
        final List<Matrix> pi;
        final List<Matrix> Q;
        final double[] rhoProc;
        final boolean[] isGeomProc;
        QbdEquilibrium(List<Matrix> pi, List<Matrix> Q, double[] rhoProc, boolean[] isGeomProc) {
            this.pi = pi; this.Q = Q; this.rhoProc = rhoProc; this.isGeomProc = isGeomProc;
        }
    }

    private static QbdEquilibrium computeEquilibriumQbd(Matrix x, Matrix[] Aa, Matrix[] Pb, Matrix[] L,
                                                        int[] ACT, int[] PSV, int numProcesses,
                                                        int numActions, int[] N, boolean[] isOpenProc) {
        List<Matrix> Q = new ArrayList<Matrix>();
        List<Matrix> pi = new ArrayList<Matrix>();
        double[] rhoProc = new double[numProcesses];
        boolean[] isGeomProc = new boolean[numProcesses];

        for (int k = 0; k < numProcesses; k++) {
            int size = N[k];

            // Assemble strictly off-diagonal rate matrix Off for component k
            Matrix Off = (L[k] != null) ? L[k].copy() : new Matrix(size, size);
            for (int i = 0; i < size; i++) Off.set(i, i, 0.0);
            for (int c = 0; c < numActions; c++) {
                if (PSV[c] == k && Pb[c] != null) {
                    double xc = x.get(c, 0);
                    for (int i = 0; i < size; i++)
                        for (int j = 0; j < size; j++)
                            if (i != j) Off.set(i, j, Off.get(i, j) + xc * Pb[c].get(i, j));
                } else if (ACT[c] == k && Aa[c] != null) {
                    for (int i = 0; i < size; i++)
                        for (int j = 0; j < size; j++)
                            if (i != j) Off.set(i, j, Off.get(i, j) + Aa[c].get(i, j));
                }
            }

            Matrix Qk = Off.copy();
            for (int i = 0; i < size; i++) {
                double rowSum = 0.0;
                for (int j = 0; j < size; j++) if (i != j) rowSum += Off.get(i, j);
                Qk.set(i, i, -rowSum);
            }
            Qk = Ctmc_makeinfgen.ctmc_makeinfgen(Qk);
            Q.add(Qk);

            boolean solvedGeom = false;
            if (isOpenProc != null && k < isOpenProc.length && isOpenProc[k] && size >= 5) {
                int s0 = size - 2;                 // interior state (level s0), 0-indexed
                double f = Off.get(s0, s0 + 1);    // up-1 rate (arrival)
                double b = Off.get(s0, s0 - 1);    // down-1 rate (service + single removal)
                double g0 = Off.get(s0, 0);        // drain to empty state (catastrophe)
                double interDown = 0.0;            // batch to strictly-interior lower levels
                for (int j = 1; j < s0 - 1; j++) interDown += Off.get(s0, j);
                if (interDown <= 1e-11 && f > 0) {
                    double rho = qbdScalarRho(f, b, g0);
                    if (!Double.isInfinite(rho) && !Double.isNaN(rho) && rho > 0 && rho < 1 - 1e-12) {
                        rhoProc[k] = rho;
                        isGeomProc[k] = true;
                        Matrix piK = new Matrix(1, size);
                        double coeff = 1.0 - rho;
                        double p = coeff;
                        for (int n = 0; n < size; n++) {
                            piK.set(0, n, p);
                            p *= rho;
                        }
                        pi.add(piK);
                        solvedGeom = true;
                    }
                }
            }

            if (!solvedGeom) {
                pi.add(isTridiagonal(Qk) ? birthDeathSolve(Qk) : Ctmc_solve.ctmc_solve(Qk));
            }
        }

        return new QbdEquilibrium(pi, Q, rhoProc, isGeomProc);
    }

    /**
     * Sub-unit root rho of the scalar QBD characteristic equation
     *   b*rho^2 - (f+b+g)*rho + f = 0,
     * where f is the up-1 rate, b the down-1 rate and g the extra local outflow
     * (catastrophe drain to the empty state). Block-size-1 instance of Neuts' R.
     */
    private static double qbdScalarRho(double f, double b, double g) {
        if (b <= 1e-14) {
            return (f + g <= 0) ? Double.POSITIVE_INFINITY : f / (f + g);
        }
        double c1 = -(f + b + g);
        double disc = c1 * c1 - 4.0 * b * f;
        if (disc < 0) {
            return Double.POSITIVE_INFINITY;
        }
        double sq = Math.sqrt(disc);
        double r1 = (-c1 - sq) / (2.0 * b);
        double r2 = (-c1 + sq) / (2.0 * b);
        double lo = Math.min(r1, r2);
        double hi = Math.max(r1, r2);
        return (lo > 0) ? lo : hi;
    }

    private static Pair<List<Matrix>, List<Matrix>> computeEquilibrium(Matrix x, Matrix[] Aa, Matrix[] Pb, Matrix[] L,
                                                                       int[] ACT, int[] PSV, int numProcesses,
                                                                       int numActions, int[] N) {
        List<Matrix> Q = new ArrayList<Matrix>();
        List<Matrix> pi = new ArrayList<Matrix>();

        for (int k = 0; k < numProcesses; k++) {
            int size = N[k];

            Matrix Lk = L[k];
            Matrix Qk = (Lk != null) ? Lk.copy() : new Matrix(size, size);

            Matrix diagVec = new Matrix(1, size);
            for (int i = 0; i < size; i++) {
                double rowSum = 0.0;
                for (int j = 0; j < size; j++) {
                    if (i != j) {
                        rowSum += Qk.get(i, j);
                    }
                }
                diagVec.set(0, i, -rowSum);
            }
            for (int i = 0; i < size; i++) {
                Qk.set(i, i, diagVec.get(0, i));
            }

            for (int c = 0; c < numActions; c++) {
                if (PSV[c] == k) {
                    Matrix pbMatrix = Pb[c];
                    if (pbMatrix == null) continue;
                    double xc = x.get(c, 0);
                    for (int i = 0; i < size; i++) {
                        for (int j = 0; j < size; j++) {
                            if (i != j) {
                                Qk.set(i, j, Qk.get(i, j) + xc * pbMatrix.get(i, j));
                            }
                        }
                    }
                } else if (ACT[c] == k) {
                    Matrix aaMatrix = Aa[c];
                    if (aaMatrix == null) continue;
                    for (int i = 0; i < size; i++) {
                        for (int j = 0; j < size; j++) {
                            if (i != j) {
                                Qk.set(i, j, Qk.get(i, j) + aaMatrix.get(i, j));
                            }
                        }
                    }
                }
            }

            Qk = Ctmc_makeinfgen.ctmc_makeinfgen(Qk);
            Q.add(Qk);

            Matrix piK = isTridiagonal(Qk) ? birthDeathSolve(Qk) : Ctmc_solve.ctmc_solve(Qk);
            pi.add(piK);
        }

        return new Pair<List<Matrix>, List<Matrix>>(pi, Q);
    }

    /**
     * Check if a matrix is tridiagonal.
     */
    private static boolean isTridiagonal(Matrix Q) {
        int n = Q.getNumRows();
        for (int i = 0; i < n; i++) {
            for (int j = 0; j < n; j++) {
                if (Math.abs(i - j) > 1 && Math.abs(Q.get(i, j)) > 1e-14) {
                    return false;
                }
            }
        }
        return true;
    }

    /**
     * Solve equilibrium distribution of a birth-death (tridiagonal) CTMC.
     */
    private static Matrix birthDeathSolve(Matrix Q) {
        int n = Q.getNumRows();
        if (n <= 1) {
            Matrix pi = new Matrix(1, 1);
            pi.set(0, 0, 1.0);
            return pi;
        }
        double[] piArr = new double[n];
        piArr[0] = 1.0;
        for (int i = 1; i < n; i++) {
            double birthRate = Q.get(i - 1, i);
            double deathRate = Q.get(i, i - 1);
            piArr[i] = (deathRate > 0) ? piArr[i - 1] * birthRate / deathRate : 0.0;
        }
        double total = 0.0;
        for (double v : piArr) total += v;
        Matrix pi = new Matrix(1, n);
        if (total > 0) {
            for (int i = 0; i < n; i++) pi.set(0, i, piArr[i] / total);
        } else {
            for (int i = 0; i < n; i++) pi.set(0, i, 1.0 / n);
        }
        return pi;
    }
}
