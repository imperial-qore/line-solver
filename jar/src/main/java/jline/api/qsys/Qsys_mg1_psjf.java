/**
 * @file M/G/1 queueing system analysis with PSJF scheduling
 *
 * @since LINE 3.0
 */
package jline.api.qsys;

import java.util.Arrays;
import java.util.Comparator;

import jline.io.Ret;
import jline.util.matrix.Matrix;

public final class Qsys_mg1_psjf {
    private Qsys_mg1_psjf() {}

    /**
     * Analyzes an M/G/1 queueing system with PSJF (Preemptive Shortest Job First) scheduling.
     */
    public static Ret.qsys_prio qsys_mg1_psjf(Matrix lambda, Matrix mu, Matrix cs) {
        double[] lambdaArr = lambda.toArray1D();
        double[] muArr = mu.toArray1D();
        double[] csArr = cs.toArray1D();

        if (lambdaArr.length != muArr.length || lambdaArr.length != csArr.length) {
            throw new IllegalArgumentException("lambda, mu, and cs must have the same length");
        }

        int K = lambdaArr.length;

        for (int i = 0; i < K; i++) {
            if (lambdaArr[i] <= 0.0) throw new IllegalArgumentException("lambda[" + i + "] must be positive");
            if (muArr[i] <= 0.0) throw new IllegalArgumentException("mu[" + i + "] must be positive");
            if (csArr[i] < 0.0) throw new IllegalArgumentException("cs[" + i + "] must be non-negative");
        }

        boolean allExp = true;
        for (int i = 0; i < K; i++) {
            if (Math.abs(csArr[i] - 1.0) >= 1e-6) {
                allExp = false;
                break;
            }
        }
        if (allExp) {
            // see _kb/03-api-layer.md for rationale
            double rhoCheck = 0.0;
            for (int i = 0; i < K; i++) {
                rhoCheck += lambdaArr[i] / muArr[i];
            }
            if (rhoCheck >= 1.0) {
                throw new IllegalStateException("System is unstable: utilization rho = " + rhoCheck + " >= 1");
            }
            double[] W = new double[K];
            for (int k = 0; k < K; k++) {
                double muK = muArr[k];
                double xMax = 20.0 / muK;
                int nIntervals = 4000; // even
                double h = xMax / nIntervals;
                double acc = 0.0;
                for (int s = 0; s <= nIntervals; s++) {
                    double x = s * h;
                    double t = psjfResponseExp(x, lambdaArr, muArr);
                    double fx = muK * Math.exp(-muK * x);
                    double wgt = (s == 0 || s == nIntervals) ? 1.0 : (s % 2 == 1 ? 4.0 : 2.0);
                    acc += wgt * t * fx;
                }
                W[k] = acc * h / 3.0;
            }
            double Qexp = 0.0;
            for (int i = 0; i < K; i++) {
                Qexp += lambdaArr[i] * W[i];
            }
            return new Ret.qsys_prio(new Matrix(W), Qexp / (1.0 + Qexp));
        }

        // Sort classes by mean service time (ascending)
        final double[] meanService = new double[K];
        for (int i = 0; i < K; i++) {
            meanService[i] = 1.0 / muArr[i];
        }
        Integer[] sortIdxObj = new Integer[K];
        for (int i = 0; i < K; i++) sortIdxObj[i] = i;
        Arrays.sort(sortIdxObj, new Comparator<Integer>() {
            @Override
            public int compare(Integer a, Integer b) {
                return Double.compare(meanService[a], meanService[b]);
            }
        });
        int[] sortIdx = new int[K];
        for (int i = 0; i < K; i++) sortIdx[i] = sortIdxObj[i];

        double[] lambdaSorted = new double[K];
        double[] muSorted = new double[K];
        double[] csSorted = new double[K];
        for (int i = 0; i < K; i++) {
            lambdaSorted[i] = lambdaArr[sortIdx[i]];
            muSorted[i] = muArr[sortIdx[i]];
            csSorted[i] = csArr[sortIdx[i]];
        }

        double[] rhoI = new double[K];
        double rhoTotal = 0.0;
        for (int i = 0; i < K; i++) {
            rhoI[i] = lambdaSorted[i] / muSorted[i];
            rhoTotal += rhoI[i];
        }
        if (rhoTotal >= 1.0) {
            throw new IllegalStateException("System is unstable: utilization rho = " + rhoTotal + " >= 1");
        }

        double[] WSorted = new double[K];
        for (int k = 0; k < K; k++) {
            double x = 1.0 / muSorted[k];
            double rhoX = 0.0;
            for (int i = 0; i <= k; i++) rhoX += rhoI[i];
            double m2X = 0.0;
            for (int i = 0; i <= k; i++) {
                double ES2i = (1.0 + csSorted[i] * csSorted[i]) / (muSorted[i] * muSorted[i]);
                m2X += lambdaSorted[i] * ES2i;
            }
            if (rhoX >= 1.0) {
                WSorted[k] = Double.POSITIVE_INFINITY;
            } else {
                double waitingTerm = m2X / (2.0 * (1.0 - rhoX) * (1.0 - rhoX));
                double serviceTerm = x / (1.0 - rhoX);
                WSorted[k] = waitingTerm + serviceTerm;
            }
        }

        int[] unsortIdx = new int[K];
        for (int i = 0; i < K; i++) {
            unsortIdx[sortIdx[i]] = i;
        }
        double[] WArr = new double[K];
        for (int i = 0; i < K; i++) {
            WArr[i] = WSorted[unsortIdx[i]];
        }

        double Q = 0.0;
        for (int i = 0; i < K; i++) {
            Q += lambdaArr[i] * WArr[i];
        }
        double rhohat = Q / (1.0 + Q);

        return new Ret.qsys_prio(new Matrix(WArr), rhohat);
    }

    /**
     * PSJF conditional mean response time T(x) for exponential job sizes,
     * matching MATLAB compute_psjf_response: only jobs of size at most x
     * contribute to the load and second moment seen by a size-x job.
     */
    private static double psjfResponseExp(double x, double[] lambdaArr, double[] muArr) {
        int K = lambdaArr.length;
        double rhoX = 0.0;
        double m2Scaled = 0.0;
        for (int i = 0; i < K; i++) {
            double muI = muArr[i];
            // int_0^x t f_i(t) dt
            double intT = 1.0 / muI - (1.0 / muI + x) * Math.exp(-muI * x);
            rhoX += lambdaArr[i] * intT;
            // int_0^x t^2 f_i(t) dt
            double intT2 = 2.0 / (muI * muI)
                    - (2.0 / (muI * muI) + 2.0 * x / muI + x * x) * Math.exp(-muI * x);
            m2Scaled += lambdaArr[i] * intT2;
        }
        if (rhoX >= 1.0) {
            return Double.POSITIVE_INFINITY;
        }
        return x / (1.0 - rhoX) + m2Scaled / (2.0 * (1.0 - rhoX) * (1.0 - rhoX));
    }
}
