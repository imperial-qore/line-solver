/**
 * Acyclic Phase-Type (APH) distribution functions.
 *
 * Ported from MATLAB: matlab/lib/kpctoolbox/aph/
 */
package jline.lib.kpctoolbox.aph;

import java.util.List;
import java.util.Random;

import jline.api.mam.Map_normalize;
import jline.api.mam.Map_renewal;
import jline.api.mam.Map_scale;
import jline.api.mam.Map_exponential;
import jline.lib.kpctoolbox.mc.DTMC;
import jline.util.Pair;
import jline.util.matrix.Matrix;
import jline.util.matrix.MatrixCell;
import org.apache.commons.math3.util.FastMath;

public final class APH {
    private APH() {}

    /**
     * Convolution patterns for APH simplification.
     */
    public enum ConvolutionPattern {
        SEQUENCE,
        PARALLEL,
        BRANCH
    }

    public static Pair<double[], Matrix> aph_simplify(double[] a1, Matrix T1, double[] a2, Matrix T2, ConvolutionPattern pattern) {
        return aph_simplify(a1, T1, a2, T2, 1.0, 1.0, pattern);
    }

    public static Pair<double[], Matrix> aph_simplify(double[] a1, Matrix T1, double[] a2, Matrix T2,
                                                       double p1, double p2, ConvolutionPattern pattern) {
        int order1 = a1.length;
        int order2 = a2.length;

        switch (pattern) {
            case SEQUENCE: {
                int n = order1 + order2;
                double a1Sum = 0.0;
                for (double v : a1) a1Sum += v;

                double[] alpha = new double[n];
                for (int i = 0; i < order1; i++) alpha[i] = a1[i];
                for (int i = 0; i < order2; i++) alpha[order1 + i] = (1.0 - a1Sum) * a2[i];

                Matrix T = new Matrix(n, n);
                for (int i = 0; i < order1; i++) {
                    for (int j = 0; j < order1; j++) T.set(i, j, T1.get(i, j));
                }
                double[] exitRates = new double[order1];
                for (int i = 0; i < order1; i++) {
                    double sum = 0.0;
                    for (int j = 0; j < order1; j++) sum += T1.get(i, j);
                    exitRates[i] = -sum;
                }
                for (int i = 0; i < order1; i++) {
                    for (int j = 0; j < order2; j++) {
                        T.set(i, order1 + j, exitRates[i] * a2[j]);
                    }
                }
                for (int i = 0; i < order2; i++) {
                    for (int j = 0; j < order2; j++) T.set(order1 + i, order1 + j, T2.get(i, j));
                }
                return new Pair<double[], Matrix>(alpha, T);
            }
            case PARALLEL: {
                int n = order1 * order2 + order1 + order2;
                double a1Sum = 0.0, a2Sum = 0.0;
                for (double v : a1) a1Sum += v;
                for (double v : a2) a2Sum += v;

                double[] alpha = new double[n];
                int idx = 0;
                for (int i = 0; i < order1; i++) {
                    for (int j = 0; j < order2; j++) {
                        alpha[idx++] = a1[i] * a2[j];
                    }
                }
                for (int i = 0; i < order1; i++) alpha[idx++] = (1.0 - a2Sum) * a1[i];
                for (int i = 0; i < order2; i++) alpha[idx++] = (1.0 - a1Sum) * a2[i];

                Matrix T = new Matrix(n, n);
                double[] exit1 = new double[order1];
                double[] exit2 = new double[order2];
                for (int i = 0; i < order1; i++) {
                    double s = 0.0;
                    for (int j = 0; j < order1; j++) s += T1.get(i, j);
                    exit1[i] = -s;
                }
                for (int i = 0; i < order2; i++) {
                    double s = 0.0;
                    for (int j = 0; j < order2; j++) s += T2.get(i, j);
                    exit2[i] = -s;
                }

                for (int i = 0; i < order1; i++) {
                    for (int j = 0; j < order2; j++) {
                        int rowIdx = i * order2 + j;
                        for (int k = 0; k < order1; k++) {
                            for (int l = 0; l < order2; l++) {
                                int colIdx = k * order2 + l;
                                double value = 0.0;
                                if (j == l) value += T1.get(i, k);
                                if (i == k) value += T2.get(j, l);
                                T.set(rowIdx, colIdx, value);
                            }
                        }
                    }
                }

                for (int i = 0; i < order1; i++) {
                    for (int j = 0; j < order2; j++) {
                        int rowIdx = i * order2 + j;
                        int colIdx = order1 * order2 + i;
                        T.set(rowIdx, colIdx, exit2[j]);
                    }
                }
                for (int i = 0; i < order1; i++) {
                    for (int j = 0; j < order2; j++) {
                        int rowIdx = i * order2 + j;
                        int colIdx = order1 * order2 + order1 + j;
                        T.set(rowIdx, colIdx, exit1[i]);
                    }
                }
                for (int i = 0; i < order1; i++) {
                    for (int j = 0; j < order1; j++) {
                        T.set(order1 * order2 + i, order1 * order2 + j, T1.get(i, j));
                    }
                }
                for (int i = 0; i < order2; i++) {
                    for (int j = 0; j < order2; j++) {
                        T.set(order1 * order2 + order1 + i, order1 * order2 + order1 + j, T2.get(i, j));
                    }
                }
                return new Pair<double[], Matrix>(alpha, T);
            }
            case BRANCH:
            default: {
                int n = order1 + order2;
                double[] alpha = new double[n];
                for (int i = 0; i < order1; i++) alpha[i] = p1 * a1[i];
                for (int i = 0; i < order2; i++) alpha[order1 + i] = p2 * a2[i];

                Matrix T = new Matrix(n, n);
                for (int i = 0; i < order1; i++) {
                    for (int j = 0; j < order1; j++) T.set(i, j, T1.get(i, j));
                }
                for (int i = 0; i < order2; i++) {
                    for (int j = 0; j < order2; j++) T.set(order1 + i, order1 + j, T2.get(i, j));
                }
                return new Pair<double[], Matrix>(alpha, T);
            }
        }
    }

    public static Pair<double[], Matrix> aph_convpara(List<Pair<double[], Matrix>> distributions) {
        if (distributions.isEmpty()) {
            throw new IllegalArgumentException("Need at least one distribution");
        }
        if (distributions.size() == 1) return distributions.get(0);

        Pair<double[], Matrix> r = aph_simplify(distributions.get(0).getLeft(), distributions.get(0).getRight(),
                distributions.get(1).getLeft(), distributions.get(1).getRight(), ConvolutionPattern.PARALLEL);
        double[] alpha = r.getLeft();
        Matrix T = r.getRight();
        for (int i = 2; i < distributions.size(); i++) {
            Pair<double[], Matrix> result = aph_simplify(alpha, T, distributions.get(i).getLeft(),
                    distributions.get(i).getRight(), ConvolutionPattern.SEQUENCE);
            alpha = result.getLeft();
            T = result.getRight();
        }
        return new Pair<double[], Matrix>(alpha, T);
    }

    public static Pair<double[], Matrix> aph_convseq(List<Pair<double[], Matrix>> distributions) {
        if (distributions.isEmpty()) {
            throw new IllegalArgumentException("Need at least one distribution");
        }
        if (distributions.size() == 1) return distributions.get(0);

        Pair<double[], Matrix> r = aph_simplify(distributions.get(0).getLeft(), distributions.get(0).getRight(),
                distributions.get(1).getLeft(), distributions.get(1).getRight(), ConvolutionPattern.SEQUENCE);
        double[] alpha = r.getLeft();
        Matrix T = r.getRight();
        for (int i = 2; i < distributions.size(); i++) {
            Pair<double[], Matrix> result = aph_simplify(alpha, T, distributions.get(i).getLeft(),
                    distributions.get(i).getRight(), ConvolutionPattern.SEQUENCE);
            alpha = result.getLeft();
            T = result.getRight();
        }
        return new Pair<double[], Matrix>(alpha, T);
    }

    public static MatrixCell aph_rand() {
        return aph_rand(2);
    }

    public static MatrixCell aph_rand(int K) {
        Random random = new Random();
        Matrix D0 = new Matrix(K, K);
        Matrix D1 = new Matrix(K, K);

        for (int i = 0; i < K; i++) {
            for (int j = i; j < K; j++) D0.set(i, j, random.nextDouble());
            for (int j = 0; j < i; j++) D0.set(i, j, 0.0);
        }
        for (int i = 0; i < K; i++) {
            for (int j = 0; j < K; j++) D1.set(i, j, random.nextDouble());
        }

        MatrixCell MAP = new MatrixCell(2);
        MAP.set(0, D0);
        MAP.set(1, D1);
        return Map_normalize.map_normalize(Map_renewal.map_renewal(MAP));
    }

    public static Pair<MatrixCell, Boolean> aph_fit(double e1, double e2, double e3) {
        return aph_fit(e1, e2, e3, 10);
    }

    public static Pair<MatrixCell, Boolean> aph_fit(double e1, double e2, double e3, int nmax) {
        boolean isExact = true;
        if (Double.isInfinite(e2) || Double.isInfinite(e3)) {
            return new Pair<MatrixCell, Boolean>(Map_exponential.map_exponential(e1), true);
        }

        double n2 = e2 / (e1 * e1);
        double n3 = e3 / (e1 * e2);

        boolean n2Feas = false, n3UbFeas = false, n3LbFeas = false;
        int n = 1;
        double un = 0.0;
        double unPrev = 0.0;

        while ((!n2Feas || !n3LbFeas || !n3UbFeas) && n < nmax) {
            n++;
            unPrev = un;

            double pn = ((n + 1) * (n2 - 2) / (3 * n2 * (n - 1))) *
                    (-2 * FastMath.sqrt(n + 1.0) / FastMath.sqrt(4.0 * (n + 1) - 3 * n * n2) - 1);
            double an = (n2 - 2) / (pn * (1 - n2) + FastMath.sqrt(pn * pn + pn * n * (n2 - 2) / (n - 1)));
            double ln = ((3 + an) * (n - 1) + 2 * an) / ((n - 1) * (1 + an * pn)) -
                    (2 * an * (n + 1)) / (2 * (n - 1) + an * pn * (n * an + 2 * n - 2));
            un = (1.0 / (n * n * n2)) * (2 * (n - 2) * (n * n2 - n - 1) *
                    FastMath.sqrt(1 + n * (n2 - 2) / (n - 1)) + (n + 2) * (3 * n * n2 - 2 * n - 2));

            if (n2 >= (n + 1.0) / n && n2 <= (n + 4.0) / (n + 1)) {
                n2Feas = true;
                if (n3 >= ln) n3LbFeas = true;
            } else if (n2 >= (n + 4.0) / (n + 1)) {
                n2Feas = true;
                if (n3 >= n2 * (n + 1) / n) n3LbFeas = true;
            }

            if (n2 >= (n + 1.0) / n && n2 <= (double) n / (n - 1)) {
                n2Feas = true;
                if (n3 <= un) n3UbFeas = true;
            } else if (n2 >= (double) n / (n - 1)) {
                n2Feas = true;
                n3UbFeas = true;
            }
        }

        double fitN2 = n2;
        double fitN3 = n3;

        if (!n2Feas || !n3LbFeas || !n3UbFeas || n >= nmax) {
            fitN2 = (n + 1.0) / n;
            fitN3 = 2 * fitN2 - 1;
            isExact = false;
        }

        MatrixCell MAP;
        if (fitN2 <= (double) n / (n - 1) || fitN3 <= 2 * fitN2 - 1) {
            double b = 2 * (4 - n * (3 * fitN2 - 4)) / (fitN2 * (4 + n - n * fitN3) +
                    FastMath.sqrt(n * fitN2) * FastMath.sqrt(12 * fitN2 * fitN2 * (n + 1) +
                            16 * fitN3 * (n + 1) + fitN2 * (n * (fitN3 - 15) * (fitN3 + 1) - 8 * (fitN3 + 3))));
            double a = (b * fitN2 - 2) * (n - 1) * b / ((b - 1) * n);
            double p = (b - 1) / a;
            double lambda = 1.0;
            double mu = lambda * (n - 1) / a;

            double[] alpha = new double[n];
            alpha[0] = p;
            alpha[n - 1] = 1 - p;

            Matrix T = new Matrix(n, n);
            for (int i = 0; i < n; i++) {
                T.set(i, i, -mu);
                if (i < n - 1) T.set(i, i + 1, mu);
            }
            T.set(n - 1, n - 1, -lambda);

            Matrix D0 = T.copy();
            Matrix D1 = new Matrix(n, n);
            for (int i = 0; i < n; i++) {
                double exitRate = 0.0;
                for (int j = 0; j < n; j++) exitRate += T.get(i, j);
                exitRate = -exitRate;
                for (int j = 0; j < n; j++) D1.set(i, j, exitRate * alpha[j]);
            }

            MAP = new MatrixCell(2);
            MAP.set(0, D0);
            MAP.set(1, D1);
        } else if (fitN2 > (double) n / (n - 1) && fitN3 > unPrev) {
            // Case 2 of 2 (long algebraic block); see Bobbio-Horvath-Telek 2005
            double K1 = n - 1.0;
            double K2 = n - 2.0;
            double K3 = 3 * fitN2 - 2 * fitN3;
            double K4 = fitN3 - 3;
            double K5 = n - fitN2;
            double K6 = 1 + fitN2 - fitN3;
            double K7 = n + fitN2 - n * fitN2;
            double K8 = 3 + 3 * fitN2 * fitN2 + fitN3 - 3 * fitN2 * fitN3;

            double innerSqrt = -16 * K1 * K1 * FastMath.pow(K7, 6.0) +
                    FastMath.pow(4 * K1 * K5 * K5 * K5 + K1 * K1 * K2 * K4 * K4 * n * fitN2 * fitN2 +
                            4 * K2 * n * fitN2 * (K4 * n * n - 3 * K6 * fitN2 + K8 * n), 2.0);
            double K9 = 108 * K1 * K1 * (4 * K2 * K2 * K3 * n * n * fitN2 +
                    K1 * K1 * K2 * K4 * K4 * n * fitN2 * fitN2 +
                    4 * K1 * K5 * (K5 * K5 - 3 * K2 * K6 * n * fitN2) +
                    FastMath.sqrt(innerSqrt));

            double K10 = K4 * K4 / (4 * K3 * K3) - K5 / (K1 * K3 * fitN2);
            double K9cbrt = (K9 >= 0) ? FastMath.pow(K9, 1.0 / 3.0) : -FastMath.pow(-K9, 1.0 / 3.0);
            double K11 = FastMath.pow(2.0, 1.0 / 3.0) * (3 * K5 * K5 + K2 * (K3 + 2 * K4) * n * fitN2) /
                    (K3 * K9cbrt * fitN2);
            double K12 = K9cbrt / (3 * FastMath.pow(2.0, 7.0 / 3.0) * K1 * K1 * K3 * fitN2);
            double K13 = FastMath.sqrt(K10 + K11 + K12);
            double K14 = (6 * K1 * K3 * K4 * K5 + 4 * K2 * K3 * K3 * n - K1 * K1 * K4 * K4 * K4 * fitN2) /
                    (4 * K1 * K1 * K3 * K3 * K3 * K13 * fitN2);
            double K15 = -K4 / (2 * K3);
            double K16 = FastMath.sqrt(2 * K10 - K11 - K12 - K14);
            double K17 = FastMath.sqrt(2 * K10 - K11 - K12 + K14);

            double innerSqrt18 = 81 * FastMath.pow(4 * K5 * K5 * K5 + 4 * K2 * K4 * K5 * n * fitN2 +
                    K1 * K2 * K4 * K4 * n * fitN2 * fitN2, 2.0) -
                    48 * FastMath.pow(3 * K5 * K5 + 2 * K2 * K4 * n * fitN2, 3.0);
            double K18 = 36 * K5 * K5 * K5 + 36 * K2 * K4 * K5 * n * fitN2 +
                    9 * K1 * K2 * K4 * K4 * n * fitN2 * fitN2 - FastMath.sqrt(innerSqrt18);
            double K18cbrt = (K18 >= 0) ? FastMath.pow(K18, 1.0 / 3.0) : -FastMath.pow(-K18, 1.0 / 3.0);
            double K19 = -K5 / (K1 * K4 * fitN2) -
                    FastMath.pow(2.0, 2.0 / 3.0) * (3 * K5 * K5 + 2 * K2 * K4 * n * fitN2) /
                            (FastMath.pow(3.0, 1.0 / 3.0) * K1 * K4 * fitN2 * K18cbrt) -
                    K18cbrt / (FastMath.pow(6.0, 2.0 / 3.0) * K1 * K4 * fitN2);
            double K20 = 6 * K1 * K3 * K4 * K5 + 4 * K2 * K3 * K3 * n - K1 * K1 * K4 * K4 * K4 * fitN2;
            double K21 = K11 + K12 + K5 / (2 * n * K1 * K3);
            double K22 = FastMath.sqrt(3 * K4 * K4 / (4 * K3 * K3) - 3 * K5 / (K1 * K3 * fitN2) +
                    FastMath.sqrt(4 * K21 * K21 - n * K2 / (fitN2 * K1 * K1 * K3)));

            double f;
            if (fitN3 > unPrev && fitN3 < 3 * fitN2 / 2) {
                f = K13 + K15 - K17;
            } else if (fitN3 == 2 * fitN2 / 2) {
                f = K19;
            } else if (fitN3 > 3 * fitN2 / 2 && K20 > 0) {
                f = -K13 + K15 + K16;
            } else if (K20 == 0.0) {
                f = K15 + K22;
            } else {
                f = K13 + K15 + K17;
            }

            double a = 2 * (f - 1) * (n - 1) / ((n - 1) * (fitN2 * f * f - 2 * f + 2) - n);
            double p = (f - 1) * a;
            double lambda = 1.0;
            double mu = lambda * (n - 1) / a;

            double[] alpha = new double[n];
            alpha[0] = p;
            alpha[1] = 1 - p;

            Matrix T = new Matrix(n, n);
            for (int i = 0; i < n; i++) {
                T.set(i, i, -mu);
                if (i < n - 1) T.set(i, i + 1, mu);
            }
            T.set(0, 0, -lambda);
            T.set(0, 1, lambda);

            Matrix D0 = T.copy();
            Matrix D1 = new Matrix(n, n);
            for (int i = 0; i < n; i++) {
                double exitRate = 0.0;
                for (int j = 0; j < n; j++) exitRate += T.get(i, j);
                exitRate = -exitRate;
                for (int j = 0; j < n; j++) D1.set(i, j, exitRate * alpha[j]);
            }

            MAP = new MatrixCell(2);
            MAP.set(0, D0);
            MAP.set(1, D1);
        } else {
            System.err.println("Warning: moment set cannot be matched with an APH distribution");
            isExact = false;
            double lambda = 1.0;
            double mu = lambda * n / e1;

            double[] alpha = new double[n];
            alpha[0] = 1.0;

            Matrix T = new Matrix(n, n);
            for (int i = 0; i < n; i++) {
                T.set(i, i, -mu);
                if (i < n - 1) T.set(i, i + 1, mu);
            }

            Matrix D0 = T.copy();
            Matrix D1 = new Matrix(n, n);
            for (int i = 0; i < n; i++) {
                double exitRate = 0.0;
                for (int j = 0; j < n; j++) exitRate += T.get(i, j);
                exitRate = -exitRate;
                for (int j = 0; j < n; j++) D1.set(i, j, exitRate * alpha[j]);
            }

            MAP = new MatrixCell(2);
            MAP.set(0, D0);
            MAP.set(1, D1);
        }

        return new Pair<MatrixCell, Boolean>(Map_scale.map_scale(Map_normalize.map_normalize(MAP), e1), isExact);
    }

    public static Pair<double[], double[]> ph2hyper(MatrixCell PH) {
        Matrix D0 = PH.get(0);
        int n = D0.getNumRows();
        for (int i = 0; i < n; i++) {
            for (int j = 0; j < n; j++) {
                if (i != j && FastMath.abs(D0.get(i, j)) > 1e-10) {
                    throw new IllegalArgumentException("The PH distribution is not hyper-exponential");
                }
            }
        }
        double[] lambda = new double[n];
        for (int i = 0; i < n; i++) lambda[i] = -D0.get(i, i);

        Matrix D1 = PH.get(1);
        Matrix P = new Matrix(n, n);
        for (int i = 0; i < n; i++) {
            for (int j = 0; j < n; j++) {
                P.set(i, j, D1.get(i, j) / lambda[i]);
            }
        }
        double[] prob = DTMC.dtmc_solve(P);
        return new Pair<double[], double[]>(lambda, prob);
    }

    public static double[] hyper_rand(double[] rates, double[] probs, int nSamples) {
        Random random = new Random();
        double[] samples = new double[nSamples];
        double[] cumProbs = new double[probs.length];
        cumProbs[0] = probs[0];
        for (int i = 1; i < probs.length; i++) cumProbs[i] = cumProbs[i - 1] + probs[i];

        for (int s = 0; s < nSamples; s++) {
            double u = random.nextDouble();
            int selected = 0;
            for (int i = 0; i < cumProbs.length; i++) {
                if (u <= cumProbs[i]) {
                    selected = i;
                    break;
                }
            }
            samples[s] = -FastMath.log(random.nextDouble()) / rates[selected];
        }
        return samples;
    }
}
