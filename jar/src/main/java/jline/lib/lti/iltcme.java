/*
 * Copyright (c) 2012-2026, QORE Lab, Imperial College London
 * All rights reserved.
 *
 * Unified inverse Laplace transform using Abate-Whitt framework.
 * Supports CME (Concentrated Matrix Exponential), Euler, and Gaver-Stehfest methods.
 * Based on: Horvath, Horvath, Almousa, Telek - "Numerical inverse Laplace
 * transformation using concentrated matrix-exponential distributions"
 */
package jline.lib.lti;

import java.io.BufferedReader;
import java.io.IOException;
import java.io.InputStream;
import java.io.InputStreamReader;
import java.lang.reflect.Type;
import java.nio.charset.StandardCharsets;
import java.util.List;
import java.util.function.Function;
import java.util.function.UnaryOperator;

import org.apache.commons.math3.complex.Complex;

import com.google.gson.Gson;
import com.google.gson.reflect.TypeToken;

public final class iltcme {

    private iltcme() {}

    // Lazy-loaded CME parameters from iltcme.json resource
    private static volatile List<CmeEntry> cmeParams;

    /**
     * One row of the pre-computed CME parameter table. Public because the same
     * rows define the CME distribution (see jline.lang.processes.CME), which
     * reads them as a matrix-exponential representation of order 2*n+1.
     */
    public static final class CmeEntry {
        public int n;
        public String optim;
        public List<Double> a;
        public List<Double> b;
        public double c;
        public double omega;
        public Object phi; // can be Double or List<Double>
        public double lognorm;
        public double mu1;
        public double mu2;
        public double cv2;
    }

    /**
     * Returns the pre-computed CME parameter table from iltcme.json.
     *
     * The table is loaded once and cached, and is shared by the inverse Laplace
     * transform below and by the CME distribution class. The returned list is
     * the cached object, so callers must not modify it.
     *
     * @return the CME parameter table
     */
    public static List<CmeEntry> parameters() {
        return loadParams();
    }

    private static List<CmeEntry> loadParams() {
        List<CmeEntry> params = cmeParams;
        if (params == null) {
            synchronized (iltcme.class) {
                params = cmeParams;
                if (params == null) {
                    InputStream stream = iltcme.class.getClassLoader().getResourceAsStream("iltcme.json");
                    if (stream == null) {
                        throw new RuntimeException("iltcme.json not found in classpath resources");
                    }
                    StringBuilder sb = new StringBuilder();
                    BufferedReader reader = new BufferedReader(new InputStreamReader(stream, StandardCharsets.UTF_8));
                    try {
                        String line;
                        while ((line = reader.readLine()) != null) {
                            sb.append(line).append('\n');
                        }
                    } catch (IOException e) {
                        throw new RuntimeException(e);
                    } finally {
                        try {
                            reader.close();
                        } catch (IOException ignored) {
                        }
                    }
                    String json = sb.toString();
                    Type type = new TypeToken<List<CmeEntry>>() {}.getType();
                    params = new Gson().fromJson(json, type);
                    cmeParams = params;
                }
            }
        }
        return params;
    }

    /**
     * Compute the inverse Laplace transform of fun at time points T.
     *
     * @param fun         Laplace transform function F(s), operating on Complex
     * @param T           array of time points (must be positive)
     * @param maxFnEvals  maximum number of function evaluations allowed
     * @param method      "cme" (default), "euler", or "gaver"
     * @return array of f(t) values
     */
    public static double[] ilt(UnaryOperator<Complex> fun, double[] T, int maxFnEvals, String method) {
        Complex[][] weights = abateWhittWeights(maxFnEvals, method);
        Complex[] eta = weights[0];
        Complex[] beta = weights[1];

        // Common Abate-Whitt evaluation: f(t) = (1/t) * sum(real(eta * F(beta/t)))
        double[] result = new double[T.length];
        for (int i = 0; i < T.length; i++) {
            double t = T[i];
            double sum = 0.0;
            for (int j = 0; j < eta.length; j++) {
                Complex s = beta[j].divide(t); // beta[j] / t
                Complex fVal = fun.apply(s);
                sum += eta[j].multiply(fVal).getReal();
            }
            result[i] = sum / t;
        }
        return result;
    }

    /**
     * Matrix-valued inverse Laplace transform. Evaluates the matrix-valued
     * transform {@code fun} once per Abate-Whitt node and accumulates the full
     * (complex) matrix, sharing the eta/beta weights with the scalar overload.
     *
     * @param fun         Laplace transform F(s) returning an nr x nc complex matrix
     * @param T           array of time points (must be positive)
     * @param maxFnEvals  maximum number of function evaluations allowed
     * @param method      "cme" (default), "euler", or "gaver"
     * @return f(t) as a double[T.length][nr][nc] array (real part)
     */
    public static double[][][] iltMatrix(Function<Complex, Complex[][]> fun, double[] T, int maxFnEvals, String method) {
        Complex[][] weights = abateWhittWeights(maxFnEvals, method);
        Complex[] eta = weights[0];
        Complex[] beta = weights[1];

        // Probe output dimensions with one cheap evaluation.
        Complex[][] probe = fun.apply(beta[0].divide(T[0]));
        int nr = probe.length;
        int nc = probe[0].length;

        double[][][] result = new double[T.length][nr][nc];
        for (int i = 0; i < T.length; i++) {
            double t = T[i];
            Complex[][] acc = new Complex[nr][nc];
            for (int r = 0; r < nr; r++) {
                for (int c = 0; c < nc; c++) {
                    acc[r][c] = Complex.ZERO;
                }
            }
            for (int j = 0; j < eta.length; j++) {
                Complex[][] fVal = fun.apply(beta[j].divide(t));
                for (int r = 0; r < nr; r++) {
                    for (int c = 0; c < nc; c++) {
                        acc[r][c] = acc[r][c].add(eta[j].multiply(fVal[r][c]));
                    }
                }
            }
            for (int r = 0; r < nr; r++) {
                for (int c = 0; c < nc; c++) {
                    result[i][r][c] = acc[r][c].getReal() / t;
                }
            }
        }
        return result;
    }

    public static double[][][] iltMatrix(Function<Complex, Complex[][]> fun, double[] T, int maxFnEvals) {
        return iltMatrix(fun, T, maxFnEvals, "cme");
    }

    /**
     * Compute the Abate-Whitt eta and beta weight vectors for the given method
     * and evaluation budget. Returned as {eta, beta}.
     */
    private static Complex[][] abateWhittWeights(int maxFnEvals, String method) {
        Complex[] eta;
        Complex[] beta;

        if ("cme".equals(method)) {
            List<CmeEntry> allParams = loadParams();
            // Find the most steep CME satisfying maxFnEvals
            CmeEntry best = allParams.get(0);
            for (int i = 1; i < allParams.size(); i++) {
                CmeEntry p = allParams.get(i);
                if (p.cv2 < best.cv2 && p.n + 1 <= maxFnEvals) {
                    best = p;
                }
            }
            int n = best.n;
            // eta = [c*mu1, (a + i*b)*mu1]
            eta = new Complex[n + 1];
            for (int idx = 0; idx <= n; idx++) {
                if (idx == 0) {
                    eta[idx] = new Complex(best.c * best.mu1);
                } else {
                    eta[idx] = new Complex(best.a.get(idx - 1) * best.mu1, best.b.get(idx - 1) * best.mu1);
                }
            }
            // beta = [1, 1 + i*(1:n)*omega] * mu1
            beta = new Complex[n + 1];
            for (int idx = 0; idx <= n; idx++) {
                if (idx == 0) {
                    beta[idx] = new Complex(best.mu1);
                } else {
                    beta[idx] = new Complex(best.mu1, idx * best.omega * best.mu1);
                }
            }
        } else if ("euler".equals(method)) {
            int nEuler = (int) Math.floor((maxFnEvals - 1) / 2.0);
            // Build eta array: [0.5, ones(n_euler), zeros(n_euler-1), 2^-n_euler]
            double[] etaRaw = new double[2 * nEuler + 1];
            etaRaw[0] = 0.5;
            for (int i = 1; i <= nEuler; i++) etaRaw[i] = 1.0;
            // zeros(n_euler-1) already 0
            etaRaw[2 * nEuler] = Math.pow(2.0, -nEuler);

            // Binomial accumulation using log-gamma for numerical stability
            for (int k = 1; k < nEuler; k++) {
                double logBinom = 0.0;
                for (int i = 1; i <= nEuler; i++) logBinom += Math.log(i);
                logBinom -= nEuler * Math.log(2.0);
                for (int i = 1; i <= k; i++) logBinom -= Math.log(i);
                for (int i = 1; i <= (nEuler - k); i++) logBinom -= Math.log(i);
                etaRaw[2 * nEuler - k] = etaRaw[2 * nEuler - k + 1] + Math.exp(logBinom);
            }

            // Apply scaling: eta = 10^(n_euler/3) * (-1)^k * eta
            double scale = Math.pow(10.0, nEuler / 3.0);
            eta = new Complex[2 * nEuler + 1];
            for (int k = 0; k <= 2 * nEuler; k++) {
                double sign = (k % 2 == 0) ? 1.0 : -1.0;
                eta[k] = new Complex(scale * sign * etaRaw[k]);
            }

            // beta = n_euler*log(10)/3 + i*pi*k
            beta = new Complex[2 * nEuler + 1];
            for (int k = 0; k <= 2 * nEuler; k++) {
                beta[k] = new Complex(nEuler * Math.log(10.0) / 3.0, Math.PI * k);
            }
        } else if ("gaver".equals(method)) {
            int mfe = maxFnEvals;
            if (mfe % 2 == 1) mfe--;
            int ndiv2 = mfe / 2;

            double[] etaArr = new double[mfe];
            double[] betaArr = new double[mfe];

            // Precompute log(1), log(1)+log(2), ...
            double[] logsum = new double[mfe + 1];
            logsum[0] = 0.0;
            for (int i = 1; i <= mfe; i++) logsum[i] = logsum[i - 1] + Math.log(i);

            for (int k = 1; k <= mfe; k++) {
                double insideSum = 0.0;
                int jStart = (int) Math.floor((k + 1) / 2.0);
                int jEnd = Math.min(k, ndiv2);
                for (int j = jStart; j <= jEnd; j++) {
                    insideSum += Math.exp(
                            (ndiv2 + 1) * Math.log(j)
                                    - logsum[ndiv2 - j]
                                    + logsum[2 * j]
                                    - 2 * logsum[j]
                                    - logsum[k - j]
                                    - logsum[2 * j - k]
                    );
                }
                double sign = ((k + ndiv2) % 2 == 0) ? 1.0 : -1.0;
                etaArr[k - 1] = Math.log(2.0) * sign * insideSum;
                betaArr[k - 1] = k * Math.log(2.0);
            }

            eta = new Complex[mfe];
            beta = new Complex[mfe];
            for (int it = 0; it < mfe; it++) {
                eta[it] = new Complex(etaArr[it]);
                beta[it] = new Complex(betaArr[it]);
            }
        } else {
            throw new IllegalArgumentException(
                    "Unknown inverse Laplace transform method: " + method + ". Supported: cme, euler, gaver"
            );
        }

        return new Complex[][]{eta, beta};
    }

    /**
     * ILT with the default method "cme".
     */
    public static double[] ilt(UnaryOperator<Complex> fun, double[] T, int maxFnEvals) {
        return ilt(fun, T, maxFnEvals, "cme");
    }
}
