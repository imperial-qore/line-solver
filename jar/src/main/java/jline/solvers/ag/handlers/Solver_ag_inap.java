package jline.solvers.ag.handlers;

import jline.solvers.ag.AgExec;

import java.util.ArrayList;
import java.util.Collections;
import java.util.List;

import jline.api.da.Da_fpi;
import jline.api.mam.Qbd_R_logred;
import jline.api.mc.Ctmc_makeinfgen;
import jline.api.mc.Ctmc_solve;
import jline.util.Pair;
import jline.util.matrix.Matrix;

public final class Solver_ag_inap {
    private Solver_ag_inap() {}

    public static INAPResult solver_ag_inap(RCATModel rcat, double tol, int maxiter, int seed, String method) {
        return solver_ag_inap(rcat, tol, maxiter, seed, method, null);
    }

    /**
     * The reversed-rate fixed point, with the agents of each sweep evaluated by
     * {@code exec}. A null backend is the serial loop and is the reference the
     * others are asserted against.
     */
    public static INAPResult solver_ag_inap(RCATModel rcat, double tol, int maxiter, int seed,
                                            String method, AgExec exec) {
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
        //
        // A PHASE-EXPANDED component takes the same estimator, for the same
        // reason. On a birth-death chain every state-wise reversed rate equals
        // lambda, so their mean is exact; with a phase block per level they do
        // not, the deep truncation levels dominate the unweighted mean, and the
        // mean-of-ratios overestimates the departure rate exactly as it does on
        // a catastrophe (measured on a tandem with Erlang(2) service at Q1: the
        // reversed rate came out 1.27 against the exact 0.5, so flow was not
        // conserved). Rate conservation has no such failure mode.
        //
        // The level distance, not the state distance, is what makes a component
        // non-birth-death: the within-level phase transitions of a PH sit far
        // off the diagonal and are not a departure from that structure.
        final boolean[] notBirthDeath = new boolean[numProcesses];
        for (int k = 0; k < numProcesses; k++) {
            if (rcat.mph[k] > 1) notBirthDeath[k] = true;
            if (L[k] == null) continue;
            int[] lvl = rcat.level[k];
            for (int n = 0; n < N[k]; n++) {
                for (int m = 0; m < N[k]; m++) {
                    if (L[k].get(n, m) > 0 && Math.abs(lvl[n] - lvl[m]) > 1) {
                        notBirthDeath[k] = true;
                    }
                }
            }
        }

        // Columns of each active matrix that carry any rate. Aa does not depend
        // on x, so this is fixed for the whole fixed point.
        final int[][] activeCols = new int[numActions][];
        for (int a = 0; a < numActions; a++) {
            activeCols[a] = columnsWithRate(Aa[a], N[ACT[a]]);
        }

        // see _kb/06-solver-catalog.md for rationale
        Matrix x = new Matrix(numActions, 1);
        for (int a = 0; a < numActions; a++) {
            x.set(a, 0, (a + 1.0) / (numActions + 1.0));
        }

        Pair<List<Matrix>, List<Matrix>> initial = computeEquilibrium(x, Aa, Pb, L, ACT, PSV, numProcesses, numActions, N, rcat, exec);
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
                double[] v = leftMultiply(piK, activeMatrix, N[k]);

                if ("inapplus".equals(method) || notBirthDeath[k]) {
                    // inapplus: x(a) = sum_ij Aa(i,j) pi(i), the departure rate
                    // of the active component.
                    double lambdaSum = 0.0;
                    for (int j = 0; j < N[k]; j++) lambdaSum += v[j];
                    if (lambdaSum > 0) {
                        x.set(a, 0, lambdaSum);
                    }
                } else {
                    // inap: x(a) = mean over the support of the STATE-WISE
                    // reversed rate (pi Aa)_j / pi_j, which RCAT requires to be
                    // independent of j. On a birth-death component every column
                    // of Aa holds one entry, so this is the reference's
                    // entrywise mean of Aa(i,j) pi(i) / pi(j) term for term.
                    double sum = 0.0;
                    int cnt = 0;
                    for (int t = 0; t < activeCols[a].length; t++) {
                        int j = activeCols[a][t];
                        double piJ = piK.get(0, j);
                        if (piJ > 0 && v[j] > 0) {
                            double ratio = v[j] / piJ;
                            if (!Double.isNaN(ratio) && !Double.isInfinite(ratio)) {
                                sum += ratio;
                                cnt++;
                            }
                        }
                    }
                    if (cnt > 0) {
                        x.set(a, 0, sum / cnt);
                    }
                }
            }
                Pair<List<Matrix>, List<Matrix>> result = computeEquilibrium(x, Aa, Pb, L, ACT, PSV, numProcesses, numActions, N, rcat, exec);
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

    public static INAPResult solver_ag_inap(RCATModel rcat, double tol, int maxiter, int seed) {
        return solver_ag_inap(rcat, tol, maxiter, seed, "inap");
    }

    public static INAPResult solver_ag_inap(RCATModel rcat, double tol, int maxiter) {
        return solver_ag_inap(rcat, tol, maxiter, 0, "inap");
    }

    public static INAPResult solver_ag_inap(RCATModel rcat) {
        return solver_ag_inap(rcat, 1e-6, 1000, 0, "inap");
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
     * returned. Mirrors solver_ag.m (inap_inf).
     *
     * Reference: A. Marin, S. Rota Bulo, S. Balsamo, "A Numerical Algorithm for
     * the Decomposition of Cooperating Structured Markov Processes", MASCOTS 2012.
     */
    public static INAPResult solver_ag_inapinf(RCATModel rcat, boolean[] isOpenProc,
                                                double tol, int maxiter) {
        Matrix[][] R = rcat.R;
        Matrix AP = rcat.AP;
        int[] N = rcat.N;

        int numActions = rcat.actionMap.size();
        int numProcesses = N.length;

        if (numProcesses == 0) {
            return new INAPResult(new Matrix(0, 1), Collections.<Matrix>emptyList(),
                    Collections.<Matrix>emptyList(), 0,
                    new double[0], new boolean[0], new QbdTail[0], 0.0);
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

        QbdEquilibrium eq = computeEquilibriumQbd(x, Aa, Pb, L, ACT, PSV, numProcesses, numActions, N, isOpenProc, rcat);
        List<Matrix> pi = eq.pi;
        List<Matrix> Q = eq.Q;
        double[] rhoProc = eq.rhoProc;
        boolean[] isGeomProc = eq.isGeomProc;
        QbdTail[] geomData = eq.geomData;

        // reversed-rate fixed point on the isolated-component equilibria
        // (matrix-geometric variant), driven by the generic DA driver
        final List<List<Matrix>> QBox = new ArrayList<List<Matrix>>();
        QBox.add(Q);
        final double[][] rhoBox = new double[][]{rhoProc};
        final boolean[][] geomBox = new boolean[][]{isGeomProc};
        final QbdTail[][] tailBox = new QbdTail[][]{geomData};
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
                    if (rcat.mph[k] == 1) {
                        // Geometric tail: the active label fires only in
                        // occupied states, so x = (per-state rate) * P(occupied).
                        int occIdx = Math.min(1, N[k] - 1);
                        x.set(a, 0, aRowSum[a][occIdx] * rhoBox[0][k]);
                    } else {
                        // Matrix-geometric tail: sum_{n>=1} pi_n = pi_1 (I-R)^-1,
                        // and the active label has the same row sums at every
                        // busy level.
                        int mph = rcat.mph[k];
                        double acc = 0.0;
                        for (int i = 0; i < mph; i++) {
                            acc += tailBox[0][k].busy.get(0, i) * aRowSum[a][mph + i];
                        }
                        x.set(a, 0, acc);
                    }
                } else {
                    Matrix piK = piprev.get(k);
                    double v = 0.0;
                    for (int i = 0; i < N[k]; i++) v += piK.get(0, i) * aRowSum[a][i];
                    x.set(a, 0, v);
                }
            }
                QbdEquilibrium eq = computeEquilibriumQbd(x, Aa, Pb, L, ACT, PSV, numProcesses, numActions, N, isOpenProc, rcat);
                QBox.set(0, eq.Q);
                rhoBox[0] = eq.rhoProc;
                geomBox[0] = eq.isGeomProc;
                tailBox[0] = eq.geomData;
                return new Da_fpi.SweepResult<List<Matrix>>(eq.pi, piprev);
            }
        }, pi, fpopts);
        pi = fpres.x;
        Q = QBox.get(0);
        rhoProc = rhoBox[0];
        isGeomProc = geomBox[0];
        geomData = tailBox[0];
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

        return new INAPResult(x, pi, Q, iter, rhoProc, isGeomProc, geomData, rcatRes);
    }

    /** Holder for the matrix-geometric equilibrium of all components. */
    private static final class QbdEquilibrium {
        final List<Matrix> pi;
        final List<Matrix> Q;
        final double[] rhoProc;
        final boolean[] isGeomProc;
        final QbdTail[] geomData;
        QbdEquilibrium(List<Matrix> pi, List<Matrix> Q, double[] rhoProc, boolean[] isGeomProc,
                       QbdTail[] geomData) {
            this.pi = pi; this.Q = Q; this.rhoProc = rhoProc; this.isGeomProc = isGeomProc;
            this.geomData = geomData;
        }
    }

    private static QbdEquilibrium computeEquilibriumQbd(Matrix x, Matrix[] Aa, Matrix[] Pb, Matrix[] L,
                                                        int[] ACT, int[] PSV, int numProcesses,
                                                        int numActions, int[] N, boolean[] isOpenProc,
                                                        RCATModel rcat) {
        List<Matrix> Q = new ArrayList<Matrix>();
        List<Matrix> pi = new ArrayList<Matrix>();
        double[] rhoProc = new double[numProcesses];
        boolean[] isGeomProc = new boolean[numProcesses];
        QbdTail[] geomData = new QbdTail[numProcesses];

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

            int mph = rcat.mph[k];
            int nlev = rcat.nlev[k];
            boolean solvedGeom = false;
            if (isOpenProc != null && k < isOpenProc.length && isOpenProc[k] && nlev >= 5) {
                if (mph == 1) {
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
                            double p = 1.0 - rho;
                            for (int n = 0; n < size; n++) {
                                piK.set(0, n, p);
                                p *= rho;
                            }
                            pi.add(piK);
                            solvedGeom = true;
                        }
                    }
                } else if (isBlockTridiagonal(Qk, rcat.level[k])) {
                    // Read the homogeneous interior blocks one level below the
                    // truncation boundary, for the same reason the scalar branch
                    // reads the interior row there.
                    int s0 = nlev - 2;
                    Matrix A0 = levelBlock(Qk, s0, s0 + 1, mph);   // up: arrival
                    Matrix A1 = levelBlock(Qk, s0, s0, mph);       // local, with diagonal
                    Matrix A2 = levelBlock(Qk, s0, s0 - 1, mph);   // down: departure
                    boolean anyUp = false;
                    for (int i = 0; i < mph && !anyUp; i++) {
                        double rs = 0.0;
                        for (int j = 0; j < mph; j++) rs += A0.get(i, j);
                        if (rs > 0) anyUp = true;
                    }
                    if (anyUp) {
                        QbdTail g = qbdMatrixTail(Qk, A0, A1, A2, mph);
                        if (g != null) {
                            isGeomProc[k] = true;
                            geomData[k] = g;
                            pi.add(qbdTailExpand(g, nlev, mph));
                            solvedGeom = true;
                        }
                    }
                }
            }

            if (!solvedGeom) {
                pi.add(solveComponent(Qk, mph, nlev, rcat.level[k]));
            }
        }

        return new QbdEquilibrium(pi, Q, rhoProc, isGeomProc, geomData);
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

    /**
     * One sweep of the fixed point: every agent's generator and stationary
     * vector at the current reversed rates.
     *
     * <p>Agent k is solved in ISOLATION -- its generator reads the rest of the
     * model only through the scalar reversed rates x, and it writes only its own
     * slot -- so this is a fan-out and not a recurrence. That is what lets
     * {@link jline.solvers.ag.AgExec} evaluate the agents on a thread pool or on
     * remote workers and still walk the same iterates as the serial loop. Any
     * cross-agent read added here would silently make those backends race.</p>
     */
    private static Pair<List<Matrix>, List<Matrix>> computeEquilibrium(Matrix x, Matrix[] Aa, Matrix[] Pb, Matrix[] L,
                                                                       int[] ACT, int[] PSV, int numProcesses,
                                                                       int numActions, int[] N, RCATModel rcat) {
        return computeEquilibrium(x, Aa, Pb, L, ACT, PSV, numProcesses, numActions, N, rcat, null);
    }

    static Pair<List<Matrix>, List<Matrix>> computeEquilibrium(Matrix x, Matrix[] Aa, Matrix[] Pb, Matrix[] L,
                                                               int[] ACT, int[] PSV, int numProcesses,
                                                               int numActions, int[] N, RCATModel rcat,
                                                               AgExec exec) {
        Matrix[] Qs = new Matrix[numProcesses];
        Matrix[] pis = new Matrix[numProcesses];

        if (exec == null) {
            for (int k = 0; k < numProcesses; k++) {
                Qs[k] = agentGenerator(k, x, Aa, Pb, L, ACT, PSV, numActions, N);
                pis[k] = solveComponent(Qs[k], rcat.mph[k], rcat.nlev[k], rcat.level[k]);
            }
        } else {
            exec.sweep(x, Aa, Pb, L, ACT, PSV, numProcesses, numActions, N, rcat, Qs, pis);
        }

        List<Matrix> Q = new ArrayList<Matrix>(numProcesses);
        List<Matrix> pi = new ArrayList<Matrix>(numProcesses);
        for (int k = 0; k < numProcesses; k++) {
            Q.add(Qs[k]);
            pi.add(pis[k]);
        }
        return new Pair<List<Matrix>, List<Matrix>>(pi, Q);
    }

    /**
     * Agent k's generator at the current reversed rates.
     *
     * <p>Split from the stationary solve because the halves cost different
     * orders: assembling the generator is O(N^2) and solving it is O(N^3). The
     * cluster backend therefore ships only the stationary vector back and
     * rebuilds the generator on the coordinator -- sending an N-by-N matrix per
     * agent per sweep to save the cheaper half would cost more on the wire than
     * it saves on the worker.</p>
     */
    public static Matrix agentGenerator(int k, Matrix x, Matrix[] Aa, Matrix[] Pb, Matrix[] L,
                                 int[] ACT, int[] PSV, int numActions, int[] N) {
        {
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

            return Ctmc_makeinfgen.ctmc_makeinfgen(Qk);
        }
    }

    /** Agent k's stationary vector, given its generator. Shared by every backend. */
    public static Matrix agentStationary(Matrix Qk, RCATModel rcat, int k) {
        return solveComponent(Qk, rcat.mph[k], rcat.nlev[k], rcat.level[k]);
    }

    /**
     * The same solve from the agent's layout given explicitly, for a worker that
     * holds one agent rather than the whole RCATModel.
     */
    public static Matrix agentStationaryOf(Matrix Qk, int mph, int nlev, int[] level) {
        return solveComponent(Qk, mph, nlev, level);
    }

    /** The (LI, LJ) level block of Q, MPH phases per level. */
    private static Matrix levelBlock(Matrix Q, int li, int lj, int mph) {
        Matrix out = new Matrix(mph, mph);
        int r0 = li * mph;
        int c0 = lj * mph;
        for (int i = 0; i < mph; i++) {
            for (int j = 0; j < mph; j++) out.set(i, j, Q.get(r0 + i, c0 + j));
        }
        return out;
    }

    /** The row vector pi*A, for a 1 x N pi and an N x N A. */
    private static double[] leftMultiply(Matrix pi, Matrix A, int n) {
        double[] out = new double[n];
        for (int i = 0; i < n; i++) {
            double p = pi.get(0, i);
            if (p == 0.0) continue;
            for (int j = 0; j < n; j++) {
                double aij = A.get(i, j);
                if (aij != 0.0) out[j] += p * aij;
            }
        }
        return out;
    }

    /** Indices of the columns of A that carry any rate. */
    private static int[] columnsWithRate(Matrix A, int n) {
        if (A == null) return new int[0];
        boolean[] hit = new boolean[n];
        int cnt = 0;
        for (int i = 0; i < n; i++) {
            for (int j = 0; j < n; j++) {
                if (!hit[j] && A.get(i, j) > 0) {
                    hit[j] = true;
                    cnt++;
                }
            }
        }
        int[] out = new int[cnt];
        int t = 0;
        for (int j = 0; j < n; j++) if (hit[j]) out[t++] = j;
        return out;
    }

    /**
     * Stationary vector of one isolated component.
     *
     * A component with a single phase per level is the birth-death chain the
     * analyzer has always built, and the ratio recursion is both exact and
     * stable there; a phase-expanded component is block tridiagonal instead, and
     * the matrix analogue of that recursion (linear level reduction) keeps the
     * same stability at the 100-level truncation, where a null-space solve is
     * already ill-conditioned. Anything that reaches beyond the neighbouring
     * level -- a catastrophe, a batch removal -- is neither, and falls back to
     * ctmc_solve.
     */
    private static Matrix solveComponent(Matrix Qk, int mph, int nlev, int[] lvl) {
        if (mph == 1) {
            if (isTridiagonal(Qk)) return birthDeathSolve(Qk);
        } else if (isBlockTridiagonal(Qk, lvl)) {
            return qbdFiniteSolve(Qk, mph, nlev);
        }
        return Ctmc_solve.ctmc_solve(Qk);
    }

    /**
     * True when every transition of Q stays within the neighbouring level, LVL
     * being the level index of each state.
     */
    private static boolean isBlockTridiagonal(Matrix Q, int[] lvl) {
        int n = Q.getNumRows();
        for (int i = 0; i < n; i++) {
            for (int j = 0; j < n; j++) {
                if (Math.abs(lvl[i] - lvl[j]) > 1 && Math.abs(Q.get(i, j)) > 1e-14) return false;
            }
        }
        return true;
    }

    /**
     * Stationary vector of a finite block-tridiagonal generator.
     *
     * Linear level reduction: censor the chain level by level from the top,
     *   C(nlev-1) = B(nlev-1),  C(n) = B(n) + F(n) (-C(n+1))^-1 D(n+1),
     * with B, F and D the diagonal, up and down blocks. C(0) is the generator of
     * the chain censored on level 0, so pi_0 is its stationary vector and the
     * rest follows from pi_(n+1) = pi_n F(n) (-C(n+1))^-1. This is the block form
     * of birthDeathSolve and reduces to it entry for entry when mph == 1.
     */
    private static Matrix qbdFiniteSolve(Matrix Q, int m, int nlev) {
        if (nlev <= 1) return Ctmc_solve.ctmc_solve(Q);

        Matrix[] C = new Matrix[nlev];
        C[nlev - 1] = levelBlock(Q, nlev - 1, nlev - 1, m);
        for (int n = nlev - 2; n >= 0; n--) {
            Matrix F = levelBlock(Q, n, n + 1, m);
            Matrix D = levelBlock(Q, n + 1, n, m);
            Matrix negC = C[n + 1].scale(-1.0);
            C[n] = levelBlock(Q, n, n, m).add(1.0, F.mult(negC.inv()).mult(D));
        }

        Matrix pi = new Matrix(1, nlev * m);
        double[] p0 = statVector(C[0]);
        for (int i = 0; i < m; i++) pi.set(0, i, p0[i]);
        for (int n = 0; n < nlev - 1; n++) {
            Matrix F = levelBlock(Q, n, n + 1, m);
            Matrix step = F.mult(C[n + 1].scale(-1.0).inv());
            for (int j = 0; j < m; j++) {
                double acc = 0.0;
                for (int i = 0; i < m; i++) acc += pi.get(0, n * m + i) * step.get(i, j);
                pi.set(0, (n + 1) * m + j, acc);
            }
        }

        double total = 0.0;
        for (int i = 0; i < nlev * m; i++) total += pi.get(0, i);
        if (total > 0) {
            for (int i = 0; i < nlev * m; i++) pi.set(0, i, pi.get(0, i) / total);
        } else {
            for (int i = 0; i < nlev * m; i++) pi.set(0, i, 1.0 / (nlev * m));
        }
        return pi;
    }

    /**
     * Stationary vector of a generator C, allowing a reducible one.
     *
     * Replaces the first balance equation by the normalization, which is the
     * equation it is redundant with (the columns of a generator sum to zero),
     * and solves the resulting square system. Unlike a null-space solve this
     * stays well posed when the chain is reducible with ONE closed class, which
     * the level-0 chain of a phase-expanded component routinely is: a phase-type
     * restarts in the support of alpha, so every service phase outside that
     * support is unreachable once the queue has emptied at least once.
     */
    private static double[] statVector(Matrix C) {
        int m = C.getNumRows();
        Matrix A = C.copy();
        for (int i = 0; i < m; i++) A.set(i, 0, 1.0);
        Matrix rhs = new Matrix(m, 1);
        rhs.set(0, 0, 1.0);
        Matrix v = new Matrix(m, 1);
        Matrix.solve(A.transpose(), rhs, v);
        double[] out = new double[m];
        for (int i = 0; i < m; i++) out[i] = v.get(i, 0);
        return out;
    }

    /**
     * Neuts' matrix-geometric solution of one open component with MPH phases per
     * level: R from logarithmic reduction, then the boundary equations of levels
     * 0 and 1,
     *   pi_0 B00 + pi_1 A2 = 0,   pi_0 A0 + pi_1 (A1 + R A2) = 0,
     * normalized by pi_0 e + pi_1 (I - R)^-1 e = 1. Returns null when R has no
     * sub-unit spectral radius, i.e. when the isolated component is unstable and
     * has no stationary tail to report.
     *
     * Logarithmic reduction rather than successive substitutions: this runs once
     * per component per fixed-point sweep, and the quadratic convergence is what
     * keeps that affordable.
     */
    private static QbdTail qbdMatrixTail(Matrix Qk, Matrix A0, Matrix A1, Matrix A2, int mph) {
        Matrix R;
        try {
            R = Qbd_R_logred.qbd_R_logred(A2, A1, A0, 500);
        } catch (RuntimeException e) {
            return null;
        }
        for (int i = 0; i < mph; i++) {
            for (int j = 0; j < mph; j++) {
                double v = R.get(i, j);
                // The minimal solution of a QBD is NON-NEGATIVE; anything else
                // is the iteration having failed rather than a rate matrix.
                if (Double.isNaN(v) || Double.isInfinite(v) || v < -1e-12) return null;
            }
        }

        Matrix B00 = levelBlock(Qk, 0, 0, mph);
        Matrix Sys = new Matrix(2 * mph, 2 * mph);
        Matrix lowerRight = A1.add(1.0, R.mult(A2));
        for (int i = 0; i < mph; i++) {
            for (int j = 0; j < mph; j++) {
                Sys.set(i, j, B00.get(i, j));
                Sys.set(i, mph + j, A0.get(i, j));
                Sys.set(mph + i, j, A2.get(i, j));
                Sys.set(mph + i, mph + j, lowerRight.get(i, j));
            }
        }
        Matrix IR = Matrix.eye(mph).add(-1.0, R);
        Matrix ones = new Matrix(mph, 1);
        for (int i = 0; i < mph; i++) ones.set(i, 0, 1.0);
        Matrix tailMass = new Matrix(mph, 1);
        try {
            Matrix.solve(IR, ones, tailMass);
        } catch (RuntimeException e) {
            return null;
        }
        // STABILITY WITHOUT AN EIGENSOLVER. (I-R)^-1 = I + R + R^2 + ...
        // converges exactly when the spectral radius is below one, and every row
        // of that series is e_i plus non-negative terms, so (I-R)^-1 e >= 1
        // entrywise. When the isolated component is unstable the series diverges
        // and the inverse picks up negative entries, so this is the spectral
        // condition without an eigensolve (the C++ twin has no LAPACK to call).
        for (int i = 0; i < mph; i++) {
            double w = tailMass.get(i, 0);
            if (Double.isNaN(w) || Double.isInfinite(w) || w < 1.0 - 1e-9) return null;
        }
        // Replace one column by the normalization pi_0 e + pi_1 (I-R)^-1 e = 1.
        for (int i = 0; i < mph; i++) {
            Sys.set(i, 0, 1.0);
            Sys.set(mph + i, 0, tailMass.get(i, 0));
        }
        Matrix rhs = new Matrix(2 * mph, 1);
        rhs.set(0, 0, 1.0);
        Matrix v = new Matrix(2 * mph, 1);
        try {
            Matrix.solve(Sys.transpose(), rhs, v);
        } catch (RuntimeException e) {
            return null;
        }
        Matrix pi0 = new Matrix(1, mph);
        Matrix pi1 = new Matrix(1, mph);
        for (int i = 0; i < mph; i++) {
            double a = v.get(i, 0);
            double b = v.get(mph + i, 0);
            if (Double.isNaN(a) || Double.isInfinite(a) || Double.isNaN(b) || Double.isInfinite(b)) {
                return null;
            }
            pi0.set(0, i, a);
            pi1.set(0, i, b);
        }
        Matrix IRinv = IR.inv();
        Matrix busy = pi1.mult(IRinv);              // sum_{n>=1} pi_n
        Matrix qlenVec = busy.mult(IRinv);          // E[N] = pi_1 (I-R)^-2 e
        double qlen = 0.0;
        for (int i = 0; i < mph; i++) qlen += qlenVec.get(0, i);
        return new QbdTail(R, pi0, pi1, busy, qlen);
    }

    /**
     * Materialize the matrix-geometric tail over NLEV levels, so the block norm
     * of the fixed point and the RCAT residual read one vector shape for every
     * component. The metrics use the closed forms in the tail instead.
     */
    private static Matrix qbdTailExpand(QbdTail g, int nlev, int mph) {
        Matrix pi = new Matrix(1, nlev * mph);
        for (int i = 0; i < mph; i++) pi.set(0, i, g.pi0.get(0, i));
        Matrix v = g.pi1;
        for (int n = 1; n < nlev; n++) {
            for (int i = 0; i < mph; i++) pi.set(0, n * mph + i, v.get(0, i));
            v = v.mult(g.R);
        }
        return pi;
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
