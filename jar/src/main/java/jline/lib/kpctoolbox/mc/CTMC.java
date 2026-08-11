/**
 * Continuous-Time Markov Chain (CTMC) analysis functions.
 *
 * Ported from MATLAB: matlab/lib/kpctoolbox/mc/ctmc_*.m
 */
package jline.lib.kpctoolbox.mc;

import java.util.ArrayDeque;
import java.util.ArrayList;
import java.util.Arrays;
import java.util.List;
import java.util.Random;

import jline.util.Pair;
import jline.util.matrix.Matrix;
import org.apache.commons.math3.linear.LUDecomposition;
import org.apache.commons.math3.linear.MatrixUtils;
import org.apache.commons.math3.linear.RealMatrix;
import org.apache.commons.math3.util.FastMath;

public final class CTMC {
    private CTMC() {}

    /**
     * Normalizes a matrix to be a valid infinitesimal generator.
     */
    public static Matrix ctmc_makeinfgen(Matrix Q) {
        int n = Q.getNumRows();
        Matrix result = new Matrix(n, n);

        for (int i = 0; i < n; i++) {
            double rowSum = 0.0;
            for (int j = 0; j < n; j++) {
                if (i != j) {
                    double value = Q.get(i, j);
                    result.set(i, j, value);
                    rowSum += value;
                }
            }
            result.set(i, i, -rowSum);
        }
        return result;
    }

    /**
     * Generates a random infinitesimal generator matrix.
     */
    public static Matrix ctmc_rand(int n) {
        Random random = new Random();
        Matrix Q = new Matrix(n, n);
        for (int i = 0; i < n; i++) {
            for (int j = 0; j < n; j++) {
                Q.set(i, j, random.nextDouble());
            }
        }
        return ctmc_makeinfgen(Q);
    }

    /**
     * Result of connected component analysis.
     */
    public static class ConnectedComponents {
        public final int numComponents;
        public final int[] componentAssignment;

        public ConnectedComponents(int numComponents, int[] componentAssignment) {
            this.numComponents = numComponents;
            this.componentAssignment = componentAssignment;
        }

        @Override
        public boolean equals(Object other) {
            if (this == other) return true;
            if (other == null || getClass() != other.getClass()) return false;
            ConnectedComponents o = (ConnectedComponents) other;
            return numComponents == o.numComponents && Arrays.equals(componentAssignment, o.componentAssignment);
        }

        @Override
        public int hashCode() {
            int result = numComponents;
            result = 31 * result + Arrays.hashCode(componentAssignment);
            return result;
        }
    }

    /**
     * Finds weakly connected components in a directed graph.
     */
    public static ConnectedComponents weaklyconncomp(Matrix G) {
        int n = G.getNumRows();

        Matrix adj = new Matrix(n, n);
        for (int i = 0; i < n; i++) {
            for (int j = 0; j < n; j++) {
                double gij = G.get(i, j);
                double gji = G.get(j, i);
                if (Double.isNaN(gij) || gij != 0.0 || Double.isNaN(gji) || gji != 0.0) {
                    adj.set(i, j, 1.0);
                    adj.set(j, i, 1.0);
                }
            }
        }

        boolean[] visited = new boolean[n];
        int[] component = new int[n];
        int numComponents = 0;

        for (int start = 0; start < n; start++) {
            if (!visited[start]) {
                numComponents++;
                ArrayDeque<Integer> queue = new ArrayDeque<Integer>();
                queue.add(start);
                visited[start] = true;
                component[start] = numComponents;

                while (!queue.isEmpty()) {
                    int current = queue.removeFirst();
                    for (int neighbor = 0; neighbor < n; neighbor++) {
                        if (!visited[neighbor] && adj.get(current, neighbor) != 0.0) {
                            visited[neighbor] = true;
                            component[neighbor] = numComponents;
                            queue.add(neighbor);
                        }
                    }
                }
            }
        }

        return new ConnectedComponents(numComponents, component);
    }

    /**
     * Result of CTMC solving.
     */
    public static class CTMCSolveResult {
        public final double[] equilibriumDistribution;
        public final Matrix generator;
        public final int numComponents;
        public final int[] componentAssignment;

        public CTMCSolveResult(double[] equilibriumDistribution, Matrix generator,
                               int numComponents, int[] componentAssignment) {
            this.equilibriumDistribution = equilibriumDistribution;
            this.generator = generator;
            this.numComponents = numComponents;
            this.componentAssignment = componentAssignment;
        }

        @Override
        public boolean equals(Object other) {
            if (this == other) return true;
            if (other == null || getClass() != other.getClass()) return false;
            CTMCSolveResult o = (CTMCSolveResult) other;
            return Arrays.equals(equilibriumDistribution, o.equilibriumDistribution)
                    && numComponents == o.numComponents
                    && Arrays.equals(componentAssignment, o.componentAssignment);
        }

        @Override
        public int hashCode() {
            int result = Arrays.hashCode(equilibriumDistribution);
            result = 31 * result + numComponents;
            result = 31 * result + Arrays.hashCode(componentAssignment);
            return result;
        }
    }

    /**
     * Computes the equilibrium distribution of a continuous-time Markov chain.
     */
    public static double[] ctmc_solve(Matrix Q) {
        return ctmc_solveFull(Q).equilibriumDistribution;
    }

    /**
     * Computes the equilibrium distribution with full details.
     */
    public static CTMCSolveResult ctmc_solveFull(Matrix Q) {
        int n = Q.getNumRows();

        if (n == 1) {
            return new CTMCSolveResult(new double[]{1.0}, Q, 1, new int[]{1});
        }

        Matrix normalizedQ = ctmc_makeinfgen(Q);

        Matrix symMatrix = new Matrix(n, n);
        for (int i = 0; i < n; i++) {
            for (int j = 0; j < n; j++) {
                if (FastMath.abs(normalizedQ.get(i, j)) + FastMath.abs(normalizedQ.get(j, i)) > 0) {
                    symMatrix.set(i, j, 1.0);
                }
            }
        }

        ConnectedComponents cc = weaklyconncomp(symMatrix);

        if (cc.numComponents > 1) {
            double[] p = new double[n];
            for (int c = 1; c <= cc.numComponents; c++) {
                List<Integer> indices = new ArrayList<Integer>();
                for (int i = 0; i < n; i++) {
                    if (cc.componentAssignment[i] == c) indices.add(i);
                }
                int m = indices.size();
                Matrix Qc = new Matrix(m, m);
                for (int ii = 0; ii < m; ii++) {
                    for (int jj = 0; jj < m; jj++) {
                        Qc.set(ii, jj, normalizedQ.get(indices.get(ii), indices.get(jj)));
                    }
                }
                double[] pc = ctmc_solve(ctmc_makeinfgen(Qc));
                for (int ii = 0; ii < m; ii++) {
                    p[indices.get(ii)] = pc[ii];
                }
            }

            double sum = 0.0;
            for (int i = 0; i < n; i++) sum += p[i];
            if (sum > 0) {
                for (int i = 0; i < n; i++) p[i] /= sum;
            }

            return new CTMCSolveResult(p, normalizedQ, cc.numComponents, cc.componentAssignment);
        }

        boolean allZero = true;
        for (int i = 0; i < n && allZero; i++) {
            for (int j = 0; j < n; j++) {
                if (normalizedQ.get(i, j) != 0.0) {
                    allZero = false;
                    break;
                }
            }
        }

        if (allZero) {
            double[] p = new double[n];
            int[] comp = new int[n];
            for (int i = 0; i < n; i++) { p[i] = 1.0 / n; comp[i] = 1; }
            return new CTMCSolveResult(p, normalizedQ, 1, comp);
        }

        int[] nnzel = new int[n];
        for (int i = 0; i < n; i++) nnzel[i] = i;
        Matrix Qnnz = normalizedQ;
        double[] bnnz = new double[n];
        Matrix Qnnz_prev = Qnnz;
        double[] bnnz_prev = bnnz.clone();

        boolean goon = true;
        while (goon) {
            int m = Qnnz.getNumRows();
            List<Integer> keep = new ArrayList<Integer>();
            for (int idx = 0; idx < m; idx++) {
                double colSum = 0.0, rowSum = 0.0;
                for (int k = 0; k < m; k++) {
                    colSum += FastMath.abs(Qnnz.get(k, idx));
                    rowSum += FastMath.abs(Qnnz.get(idx, k));
                }
                if (colSum != 0.0 && rowSum != 0.0) keep.add(idx);
            }

            int mKeep = keep.size();
            Matrix QnnzNew = new Matrix(mKeep, mKeep);
            double[] bnnzNew = new double[mKeep];
            for (int ii = 0; ii < mKeep; ii++) {
                bnnzNew[ii] = bnnz[keep.get(ii)];
                for (int jj = 0; jj < mKeep; jj++) {
                    QnnzNew.set(ii, jj, Qnnz.get(keep.get(ii), keep.get(jj)));
                }
            }

            Matrix QnnzNorm = ctmc_makeinfgen(QnnzNew);

            int[] newNnzel = new int[mKeep];
            for (int ii = 0; ii < mKeep; ii++) newNnzel[ii] = nnzel[keep.get(ii)];

            if (Qnnz_prev.getNumRows() == QnnzNorm.getNumRows() && bnnz_prev.length == bnnzNew.length) {
                goon = false;
            } else {
                Qnnz_prev = QnnzNorm;
                bnnz_prev = bnnzNew;
                nnzel = newNnzel;
            }

            Qnnz = QnnzNorm;
            bnnz = bnnzNew;
            nnzel = newNnzel;
        }

        if (Qnnz.getNumRows() == 0) {
            double[] p = new double[n];
            int[] comp = new int[n];
            for (int i = 0; i < n; i++) { p[i] = 1.0 / n; comp[i] = 1; }
            return new CTMCSolveResult(p, normalizedQ, 1, comp);
        }

        int mSolve = Qnnz.getNumRows();

        Matrix Qmod = new Matrix(mSolve, mSolve);
        for (int i = 0; i < mSolve; i++) {
            for (int j = 0; j < mSolve; j++) {
                Qmod.set(j, i, Qnnz.get(i, j));
            }
        }

        for (int i = 0; i < mSolve; i++) {
            Qmod.set(i, mSolve - 1, 1.0);
        }

        double[] b = new double[mSolve];
        b[mSolve - 1] = 1.0;

        RealMatrix realMatrix = MatrixUtils.createRealMatrix(mSolve, mSolve);
        for (int i = 0; i < mSolve; i++) {
            for (int j = 0; j < mSolve; j++) {
                realMatrix.setEntry(i, j, Qmod.get(i, j));
            }
        }

        double[] pSolve;
        try {
            pSolve = new LUDecomposition(realMatrix).getSolver().solve(MatrixUtils.createRealVector(b)).toArray();
        } catch (Exception e) {
            pSolve = new double[mSolve];
            for (int i = 0; i < mSolve; i++) pSolve[i] = Double.NaN;
        }

        double[] p = new double[n];
        for (int ii = 0; ii < nnzel.length; ii++) {
            p[nnzel[ii]] = pSolve[ii];
        }

        return new CTMCSolveResult(p, normalizedQ, cc.numComponents, cc.componentAssignment);
    }

    /**
     * Computes the time-reversed generator of a CTMC.
     */
    public static Matrix ctmc_timereverse(Matrix Q) {
        int n = Q.getNumRows();
        Matrix Qrev = new Matrix(n, n);
        double[] pie = ctmc_solve(Q);

        for (int i = 0; i < n; i++) {
            for (int j = 0; j < n; j++) {
                if (pie[j] != 0.0) {
                    Qrev.set(j, i, Q.get(i, j) * pie[i] / pie[j]);
                }
            }
        }
        return ctmc_makeinfgen(Qrev);
    }

    /**
     * Applies uniformization (randomization) to transform a CTMC into a DTMC.
     */
    public static Pair<Matrix, Double> ctmc_randomization(Matrix Q, Double q) {
        int n = Q.getNumRows();
        double rate = (q != null) ? q.doubleValue() : (maxAbsDiagonal(Q) + new Random().nextDouble());

        Matrix P = new Matrix(n, n);
        for (int i = 0; i < n; i++) {
            for (int j = 0; j < n; j++) {
                P.set(i, j, Q.get(i, j) / rate);
            }
            P.set(i, i, P.get(i, i) + 1.0);
        }
        return new Pair<Matrix, Double>(DTMC.dtmc_makestochastic(P), rate);
    }

    public static Pair<Matrix, Double> ctmc_randomization(Matrix Q) {
        return ctmc_randomization(Q, null);
    }

    /**
     * Computes transient probabilities using uniformization method.
     */
    public static Pair<double[], Integer> ctmc_uniformization(double[] pi0, Matrix Q, double t) {
        return ctmc_uniformization(pi0, Q, t, 1e-12, 100);
    }

    public static Pair<double[], Integer> ctmc_uniformization(double[] pi0, Matrix Q, double t, double tol, int maxiter) {
        int n = Q.getNumRows();
        double maxDiag = 0.0;
        for (int i = 0; i < n; i++) {
            maxDiag = Math.max(maxDiag, FastMath.abs(Q.get(i, i)));
        }
        double q = 1.1 * maxDiag;

        Matrix Qs = new Matrix(n, n);
        for (int i = 0; i < n; i++) {
            for (int j = 0; j < n; j++) {
                Qs.set(i, j, Q.get(i, j) / q);
            }
            Qs.set(i, i, Qs.get(i, i) + 1.0);
        }

        int k = 0;
        double s = 1.0;
        double r = 1.0;
        int kmax = 1;

        for (int iter = 0; iter < maxiter; iter++) {
            k++;
            r = r * (q * t) / k;
            s += r;
            if (1 - FastMath.exp(-q * t) * s <= tol) {
                kmax = k;
                break;
            }
        }

        double[] pi = new double[pi0.length];
        for (int i = 0; i < pi0.length; i++) pi[i] = pi0[i] * FastMath.exp(-q * t);
        double[] P = pi0.clone();
        double ri = FastMath.exp(-q * t);

        for (int j = 1; j <= kmax; j++) {
            double[] newP = new double[n];
            for (int i = 0; i < n; i++) {
                double sum = 0.0;
                for (int l = 0; l < n; l++) {
                    sum += P[l] * Qs.get(l, i);
                }
                newP[i] = sum;
            }
            P = newP;

            ri *= (q * t / j);
            for (int i = 0; i < n; i++) {
                pi[i] += ri * P[i];
            }
        }

        return new Pair<double[], Integer>(pi, kmax);
    }

    private static double maxAbsDiagonal(Matrix M) {
        int n = M.getNumRows();
        double maxVal = 0.0;
        for (int i = 0; i < n; i++) {
            maxVal = Math.max(maxVal, FastMath.abs(M.get(i, i)));
        }
        return maxVal;
    }

    /**
     * Computes the equilibrium distribution relative to a reference state.
     */
    public static double[] ctmc_relsolve(Matrix Q) {
        return ctmc_relsolve(Q, 0);
    }

    public static double[] ctmc_relsolve(Matrix Q, int refstate) {
        int n = Q.getNumRows();
        if (n == 1) {
            return new double[]{1.0};
        }
        Matrix normalizedQ = ctmc_makeinfgen(Q);

        Matrix Qmod = new Matrix(n, n);
        for (int i = 0; i < n; i++) {
            for (int j = 0; j < n; j++) {
                Qmod.set(j, i, normalizedQ.get(i, j));
            }
        }

        int lastRow = n - 1;
        for (int j = 0; j < n; j++) {
            Qmod.set(lastRow, j, 0.0);
        }
        Qmod.set(lastRow, refstate, 1.0);

        double[] b = new double[n];
        b[lastRow] = 1.0;

        RealMatrix realMatrix = MatrixUtils.createRealMatrix(n, n);
        for (int i = 0; i < n; i++) {
            for (int j = 0; j < n; j++) {
                realMatrix.setEntry(i, j, Qmod.get(i, j));
            }
        }

        try {
            return new LUDecomposition(realMatrix).getSolver().solve(MatrixUtils.createRealVector(b)).toArray();
        } catch (Exception e) {
            double[] result = new double[n];
            for (int i = 0; i < n; i++) result[i] = 1.0 / n;
            return result;
        }
    }
}
