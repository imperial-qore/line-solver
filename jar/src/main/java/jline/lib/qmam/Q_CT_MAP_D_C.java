/**
 * @file Q_CT_MAP_D_C - Continuous-Time MAP/D/c Queue Analyzer
 *
 * @since LINE 3.1.0
 */
package jline.lib.qmam;

import java.util.ArrayList;
import java.util.List;

import jline.lib.smc.NSFGHTOptions;
import jline.lib.smc.NSFPiOptions;
import jline.lib.smc.NSF_GHT;
import jline.lib.smc.Stat;
import jline.util.matrix.Matrix;

public final class Q_CT_MAP_D_C {
    private Q_CT_MAP_D_C() {}

    /**
     * Computes queue length and waiting time distribution for a MAP/D/c/FCFS queue.
     */
    public static MAPDcResult qCtMapDC(Matrix D0, Matrix D1, double s, int c, MAPDcOptions options) {
        int m = D0.getNumRows();

        if (D0.getNumCols() != m || D1.getNumRows() != m || D1.getNumCols() != m) {
            throw new IllegalArgumentException("D0 and D1 must be m x m matrices");
        }
        if (s <= 1e-14) throw new IllegalArgumentException("Service time s must be strictly positive");
        if (c < 1) throw new IllegalArgumentException("Number of servers c must be at least 1");

        double lambda = maxDiagonal(D0);
        Matrix P0 = D0.scale(1.0 / lambda).add(Matrix.eye(m));
        Matrix P1 = D1.scale(1.0 / lambda);

        Matrix thetaA = Stat.stat(P0.add(P1));
        double lambdaA = thetaA.mult(D1.sumRows()).get(0, 0);
        double epsilon = 1e-12;

        double load = lambdaA * s / c;
        if (load >= 1 - epsilon) {
            throw new IllegalArgumentException("The load " + load + " of the system exceeds one");
        }

        Matrix P0s = matrixExpm(D0.scale(s));
        Matrix Ptot = P0s.copy();

        double[] poissonTerms = computePoissonTerms(lambda * s, epsilon);
        int Pterms = poissonTerms.length;

        Matrix[] Kold = new Matrix[Pterms];
        for (int i = 0; i < Pterms; i++) {
            Kold[i] = (i == 0) ? Matrix.eye(m) : Matrix.zeros(m, m);
        }
        for (int j = 1; j < Pterms; j++) {
            Kold[j] = Kold[j - 1].mult(P0);
        }

        int k = 1;
        double probCum = minRowSum(Ptot);
        List<Matrix> Ps = new ArrayList<Matrix>();

        while (probCum < 1 - epsilon) {
            Matrix[] K = new Matrix[Pterms];
            for (int i = 0; i < Pterms; i++) K[i] = Matrix.zeros(m, m);

            K[0] = P1.mult(Kold[0]);
            double htemp = (k <= poissonTerms.length) ? poissonTerms[k - 1] : poissonPdf(k, lambda * s);
            Matrix Psk = K[0].scale(htemp);

            for (int j = 1; j < Pterms; j++) {
                K[j] = P0.mult(K[j - 1]).add(P1.mult(Kold[j]));
                htemp = htemp * lambda * s / (k - 1 + j + 1);
                Psk = Psk.add(K[j].scale(htemp));
            }

            Ps.add(Psk);
            Ptot = Ptot.add(Psk);
            Kold = K;
            probCum = minRowSum(Ptot);
            k++;
        }

        int numBlocks = 1 + Ps.size();
        Matrix A = new Matrix(m, m * numBlocks);

        for (int i = 0; i < m; i++) {
            for (int j = 0; j < m; j++) {
                A.set(i, j, P0s.get(i, j));
            }
        }

        for (int blk = 0; blk < Ps.size(); blk++) {
            Matrix Pblk = Ps.get(blk);
            for (int i = 0; i < m; i++) {
                for (int j = 0; j < m; j++) {
                    A.set(i, (blk + 1) * m + j, Pblk.get(i, j));
                }
            }
        }

        Matrix G = NSF_GHT.nsfGht(A, c, new NSFGHTOptions(options.getVerbose()));

        Matrix pi = NSF_GHT.nsfPi(null, A, G, new NSFPiOptions(options.getMaxNumComp(), options.getVerbose() > 0));

        int numLevels = pi.getNumCols() / m;
        Matrix ql = new Matrix(1, numLevels);
        for (int i = 0; i < numLevels; i++) {
            double levelSum = 0.0;
            for (int j = 0; j < m; j++) {
                levelSum += pi.get(0, i * m + j);
            }
            ql.set(0, i, levelSum);
        }

        List<Double> w = new ArrayList<Double>();
        double w0 = 0.0;
        int wEnd = Math.min(c * m, pi.getNumCols());
        for (int i = 0; i < wEnd; i++) {
            w0 += pi.get(0, i);
        }
        w.add(w0);

        double wtAccum = w0;
        int i = 2;
        while (wtAccum < 1 - 1e-10 && i * m * c < pi.getNumCols()) {
            double wLevel = w.get(i - 2);
            int startIdx = (i - 1) * m * c;
            int endIdx = Math.min(i * m * c, pi.getNumCols());
            for (int idx = startIdx; idx < endIdx; idx++) {
                wLevel += pi.get(0, idx);
            }
            w.add(wLevel);
            wtAccum = wLevel;
            i++;
        }

        if (options.getNumSteps() > 1) {
            List<Double> refinedW = computeRefinedWaitingTime(
                    pi, D0, P0, P1, lambda, s, c, m, thetaA, w, options.getNumSteps(), epsilon);
            Matrix wtMatrix = new Matrix(1, refinedW.size());
            for (int idx = 0; idx < refinedW.size(); idx++) {
                wtMatrix.set(0, idx, refinedW.get(idx));
            }
            return new MAPDcResult(ql, wtMatrix);
        }

        Matrix wtMatrix = new Matrix(1, w.size());
        for (int idx = 0; idx < w.size(); idx++) {
            wtMatrix.set(0, idx, w.get(idx));
        }
        return new MAPDcResult(ql, wtMatrix);
    }

    public static MAPDcResult qCtMapDC(Matrix D0, Matrix D1, double s, int c) {
        return qCtMapDC(D0, D1, s, c, new MAPDcOptions());
    }

    private static double[] computePoissonTerms(double lambdaS, double epsilon) {
        List<Double> terms = new ArrayList<Double>();
        double h = Math.exp(-lambdaS) * lambdaS;
        double sumH = h + Math.exp(-lambdaS);
        terms.add(h);
        while (sumH < 1 - epsilon) {
            double nextH = terms.get(terms.size() - 1) * lambdaS / (terms.size() + 1);
            terms.add(nextH);
            sumH += nextH;
        }
        double[] arr = new double[terms.size()];
        for (int i = 0; i < arr.length; i++) arr[i] = terms.get(i);
        return arr;
    }

    private static double poissonPdf(int k, double lambda) {
        double result = Math.exp(-lambda);
        for (int i = 1; i <= k; i++) {
            result *= lambda / i;
        }
        return result;
    }

    private static double maxDiagonal(Matrix M) {
        double maxVal = 0.0;
        for (int i = 0; i < M.getNumRows(); i++) {
            maxVal = Math.max(maxVal, Math.abs(M.get(i, i)));
        }
        return maxVal;
    }

    private static double minRowSum(Matrix M) {
        double minSum = Double.MAX_VALUE;
        for (int i = 0; i < M.getNumRows(); i++) {
            double rowSum = 0.0;
            for (int j = 0; j < M.getNumCols(); j++) {
                rowSum += M.get(i, j);
            }
            minSum = Math.min(minSum, rowSum);
        }
        return minSum;
    }

    private static Matrix matrixExpm(Matrix A) {
        return A.expm();
    }

    private static List<Double> computeRefinedWaitingTime(
            Matrix pi, Matrix D0, Matrix P0, Matrix P1, double lambda, double s, int c, int m,
            Matrix thetaA, List<Double> baseW, int numSteps, double epsilon) {
        // Simplified: returns base waiting time list (matches Kotlin behavior fallback).
        return baseW;
    }
}
