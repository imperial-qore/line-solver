/**
 * @file McKenna-Mitra integral with Gauss-Laguerre quadrature
 *
 * Implements Gauss-Laguerre quadrature integration for computing normalizing constants
 * in multi-class repairman models. Uses precomputed Gauss-Laguerre nodes and weights
 * for high-precision numerical integration of the McKenna-Mitra integral representation.
 *
 * The Gauss-Laguerre quadrature approximates integrals of the form:
 *   int_0^inf f(x) exp(-x) dx ~ sum_i w_i * f(x_i)
 *
 * @since LINE 3.0
 */
package jline.api.pfqn.nc;

import java.util.ArrayList;
import java.util.List;

import org.apache.commons.math3.util.FastMath;

import jline.io.Ret;
import jline.util.Maths;
import jline.util.matrix.Matrix;

/**
 * Compute the normalizing constant of a repairmen model using Gauss-Laguerre integration.
 */
public final class Pfqn_mmint2_gausslaguerre {
    private Pfqn_mmint2_gausslaguerre() {}

    public static Ret.pfqnNc pfqn_mmint2_gausslaguerre(Matrix L, Matrix N, Matrix Z) {
        return pfqn_mmint2_gausslaguerre(L, N, Z, 1);
    }

    public static Ret.pfqnNc pfqn_mmint2_gausslaguerre(Matrix L, Matrix N, Matrix Z, int m) {
        int nMax = 300;
        double[] nodes = getGaussLaguerreNodes(nMax);
        double[] weights = getGaussLaguerreWeights(nMax);

        List<Integer> nonzeroClasses = new ArrayList<Integer>();
        for (int i = 0; i < N.length(); i++) {
            if (N.get(i) > 0) {
                nonzeroClasses.add(i);
            }
        }

        int n = FastMath.min(nMax, 2 * (int) N.elementSum() + 1);

        // f(u) = N(nonzeroClasses) * log(Z(nonzeroClasses) + L(nonzeroClasses) * u)'
        double[] F = new double[n];
        for (int i = 0; i < n; i++) {
            double sum = 0.0;
            for (int j : nonzeroClasses) {
                sum += N.get(j) * FastMath.log(Z.get(j) + L.get(j) * nodes[i]);
            }
            F[i] = (m - 1) * FastMath.log(nodes[i]) + sum;
        }

        // g = log(w) + F - sum(factln(N)) - factln(m-1)
        Matrix g = new Matrix(1, n);
        double factlnSum = 0.0;
        for (int i = 0; i < N.length(); i++) {
            factlnSum += Maths.factln(N.get(i));
        }
        double factlnM1 = Maths.factln(m - 1);

        for (int i = 0; i < n; i++) {
            g.set(i, FastMath.log(weights[i]) + F[i] - factlnSum - factlnM1);
        }

        // lG = log(sum(exp(g)))
        double lG;
        double sumExpG = 0.0;
        for (int i = 0; i < n; i++) {
            sumExpG += FastMath.exp(g.get(i));
        }
        lG = FastMath.log(sumExpG);

        if (!Double.isFinite(lG)) {
            lG = Matrix.logsumexp(g);
        }

        double G = FastMath.exp(lG);
        return new Ret.pfqnNc(G, lG);
    }

    private static double[] getGaussLaguerreNodes(int n) {
        return computeGaussLaguerre(n)[0];
    }

    private static double[] getGaussLaguerreWeights(int n) {
        return computeGaussLaguerre(n)[1];
    }

    private static int cachedN = -1;
    private static double[] cachedNodes = null;
    private static double[] cachedWeights = null;

    private static synchronized double[][] computeGaussLaguerre(int n) {
        if (n == cachedN && cachedNodes != null && cachedWeights != null) {
            return new double[][]{cachedNodes, cachedWeights};
        }

        // Golub-Welsch: build symmetric tridiagonal matrix
        double[] diag = new double[n];
        double[] offdiag = new double[n - 1];

        for (int k = 0; k < n; k++) {
            diag[k] = (2.0 * k + 1.0);
        }
        for (int k = 0; k < n - 1; k++) {
            offdiag[k] = (k + 1);
        }

        Object[] eigenResult = symmetricTridiagonalEigen(diag, offdiag);
        double[] nodes = (double[]) eigenResult[0];
        double[][] eigenvectors = (double[][]) eigenResult[1];

        double[] weights = new double[n];
        for (int i = 0; i < n; i++) {
            weights[i] = eigenvectors[0][i] * eigenvectors[0][i];
        }

        Integer[] indices = new Integer[n];
        for (int i = 0; i < n; i++) indices[i] = i;
        final double[] nodesRef = nodes;
        java.util.Arrays.sort(indices, new java.util.Comparator<Integer>() {
            @Override
            public int compare(Integer a, Integer b) {
                return Double.compare(nodesRef[a], nodesRef[b]);
            }
        });
        double[] sortedNodes = new double[n];
        double[] sortedWeights = new double[n];
        for (int i = 0; i < n; i++) {
            sortedNodes[i] = nodes[indices[i]];
            sortedWeights[i] = weights[indices[i]];
        }

        cachedN = n;
        cachedNodes = sortedNodes;
        cachedWeights = sortedWeights;

        return new double[][]{sortedNodes, sortedWeights};
    }

    private static Object[] symmetricTridiagonalEigen(double[] diag, double[] offdiag) {
        int n = diag.length;
        double[] d = diag.clone();
        double[] e = new double[n];
        for (int i = 0; i < n - 1; i++) {
            e[i] = offdiag[i];
        }
        e[n - 1] = 0.0;

        // Initialize eigenvector matrix as identity
        double[][] z = new double[n][n];
        for (int i = 0; i < n; i++) z[i][i] = 1.0;

        int maxIter = 300;
        for (int l = 0; l < n; l++) {
            int iter = 0;
            while (true) {
                int m = l;
                while (m < n - 1) {
                    double dd = FastMath.abs(d[m]) + FastMath.abs(d[m + 1]);
                    if (FastMath.abs(e[m]) + dd == dd) break;
                    m++;
                }
                if (m == l) break;

                if (iter++ >= maxIter) break;

                double g = (d[l + 1] - d[l]) / (2.0 * e[l]);
                double r = FastMath.sqrt(g * g + 1.0);
                g = d[m] - d[l] + e[l] / (g + (g >= 0 ? FastMath.abs(r) : -FastMath.abs(r)));

                double s = 1.0;
                double c = 1.0;
                double p = 0.0;

                boolean rZero = false;
                for (int i = m - 1; i >= l; i--) {
                    double f = s * e[i];
                    double b = c * e[i];
                    r = FastMath.sqrt(f * f + g * g);
                    e[i + 1] = r;
                    if (r == 0.0) {
                        d[i + 1] -= p;
                        e[m] = 0.0;
                        rZero = true;
                        break;
                    }
                    s = f / r;
                    c = g / r;
                    g = d[i + 1] - p;
                    r = (d[i] - g) * s + 2.0 * c * b;
                    p = s * r;
                    d[i + 1] = g + p;
                    g = c * r - b;

                    for (int k = 0; k < n; k++) {
                        f = z[k][i + 1];
                        z[k][i + 1] = s * z[k][i] + c * f;
                        z[k][i] = c * z[k][i] - s * f;
                    }
                }

                if (rZero && m - 1 >= l) continue;

                d[l] -= p;
                e[l] = g;
                e[m] = 0.0;
            }
        }

        return new Object[]{d, z};
    }
}
