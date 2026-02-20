/*
 * Copyright (c) 2012-2026, QORE Lab, Imperial College London
 * All rights reserved.
 */

package jline.util;

import static jline.GlobalConstants.Inf;
import static jline.GlobalConstants.NegInf;

import jline.util.matrix.ComplexMatrix;
import jline.util.matrix.Matrix;
import org.apache.commons.math3.analysis.polynomials.PolynomialFunction;
import org.apache.commons.math3.analysis.solvers.LaguerreSolver;
import org.apache.commons.math3.complex.Complex;
import org.apache.commons.math3.special.Gamma;
import org.apache.commons.math3.util.FastMath;
import org.ejml.data.DMatrixRMaj;
import org.ejml.dense.row.CommonOps_DDRM;
import org.ejml.dense.row.NormOps_DDRM;
import org.ejml.dense.row.factory.LinearSolverFactory_DDRM;
import org.ejml.interfaces.linsol.LinearSolverDense;

import java.util.*;
import java.util.function.Function;

import static jline.io.InputOutputKt.line_error;
import static jline.io.InputOutputKt.line_warning;
import static org.apache.commons.math3.primes.Primes.isPrime;

/**
 * Mathematical functions and utilities.
 */
public class Maths {

    // ================================================================================
    // RANDOM NUMBER GENERATION
    // Functions for generating random numbers and managing random seeds
    // ================================================================================
    
    // Random number generation now always uses MersenneTwister through RandomManager
    // This ensures reproducible results and consistent behavior across all solvers

    // ================================================================================
    // COMBINATORIAL FUNCTIONS
    // Functions for combinatorial mathematics including binomial coefficients,
    // factorials, combinations, permutations, and multichoose operations
    // ================================================================================

    /**
     * Computes the binomial coefficient "n choose k" by doing the work
     * in log-space then exponentiating with FastMath.exp().
     *
     * @param n total items
     * @param k items chosen
     * @return the binomial coefficient, or 0.0 if k<0 or k>n
     */
    public static double binomialCoeff(int n, int k) {
        double lg = logBinomial(n, k);
        if (Double.isInfinite(lg) && lg < 0) {
            return 0.0;
        }
        return FastMath.exp(lg);
    }

    // ================================================================================
    // MATRIX-RELATED UTILITIES
    // Functions for matrix operations including circulant matrices, matrix exponentials,
    // concatenation, and other specialized matrix computations
    // ================================================================================

    /**
     * Returns a circulant matrix of order c.
     * Creates a circulant matrix where each row is a cyclic shift of the previous row.
     * 
     * @param c the order of the circulant matrix
     * @return a c×c circulant matrix
     */
    public static Matrix circul(int c) {
        if (c == 1) {
            Matrix C = new Matrix(1, 1);
            C.set(0, 0, 1);
            return C;
        }

        Matrix R = new Matrix(c, c);
        R.set(0, c - 1, 1);
        for (int i = 1; i < c; i++) {
            R.set(i, i - 1, 1);
        }
        R = R.transpose();
        return R;
    }

    /**
     * Returns a circulant matrix from the given first column.
     * Creates a circulant matrix where the first column is given by the input matrix,
     * and each subsequent column is a cyclic shift of the previous one.
     * 
     * @param c the first column of the circulant matrix
     * @return a square circulant matrix
     */
    public static Matrix circul(Matrix c) {
        int n = c.length();
        Matrix R = new Matrix(n, n);
        R.set(0, n - 1, 1);
        for (int i = 1; i < n; i++) {
            R.set(i, i - 1, 1);
        }
        Matrix C = new Matrix(n, n);
        for (int t = 0; t < n; t++) {
            if (t > 1) {
                R = new Matrix(R.mult(R, null));
            }
            Matrix tmpC = new Matrix(0, 0);
            R.scaleEq(c.get(0, t), tmpC);
            C = C.add(1, tmpC);
        }
        return C;
    }

    // ================================================================================
    // BASIC ARITHMETIC OPERATIONS
    // Functions for basic mathematical operations including cumulative sums,
    // element-wise operations, min/max functions, and array manipulations
    // ================================================================================

    /**
     * Cumulative sum of an array, where the value at each index of the result is the
     * sum of all the previous values of the input.
     *
     * @param ar Array we want the cumulative sum of.
     * @return New array that is the cumulative sum of the input.
     */
    public static int[] cumSum(int[] ar) {
        int[] result = new int[ar.length];
        result[0] = ar[0];
        for (int i = 1; i < ar.length; i++) {
            result[i] = result[i - 1] + ar[i];
        }
        return result;
    }

    /**
     * Cumulative sum of an array, where the value at each index of the result is the
     * sum of all the previous values of the input.
     *
     * @param ar Array we want the cumulative sum of.
     * @return New array that is the cumulative sum of the input.
     */
    public static double[] cumSum(double[] ar) {
        double[] result = new double[ar.length];
        result[0] = ar[0];
        for (int i = 1; i < ar.length; i++) {
            result[i] = result[i - 1] + ar[i];
        }
        return result;
    }

    /**
     * Helper method that adds an integer to every element of an
     * integer array.
     *
     * @param ar Array to be added to.
     * @param i  Integer to add.
     */
    public static void elementAdd(int[] ar, int i) {
        for (int j = 0; j < ar.length; j++) {
            ar[j] += i;
        }
    }

    // ================================================================================
    // PROBABILITY DISTRIBUTIONS AND STATISTICAL FUNCTIONS
    // Functions for probability distributions, error functions, and statistical operations
    // ================================================================================

    /**
     * Computes the error function (erf) of a value.
     * The error function is the integral of the Gaussian distribution.
     * 
     * @param x the input value
     * @return the error function value erf(x)
     */
    public static double normalErf(double x) {
        // Constants
        double a1 = 0.254829592;
        double a2 = -0.284496736;
        double a3 = 1.421413741;
        double a4 = -1.453152027;
        double a5 = 1.061405429;
        double p = 0.3275911;

        // Save the sign of x
        int sign = (x < 0) ? -1 : 1;
        x = FastMath.abs(x);

        // A&S formula 7.1.26
        double t = 1.0 / (1.0 + p * x);
        double y = (((((a5 * t + a4) * t) + a3) * t + a2) * t + a1) * t;

        return sign * (1 - y * FastMath.exp(-x * x));
    }

    /**
     * Adapted from jblas and IHMC Original documentation:
     *
     * <p>Calculate matrix exponential of a square matrix.
     *
     * <p>A scaled Pade approximation algorithm from Golub and Van Loan "Matrix Computations",
     * algorithm 11.3.1 and 11.2.
     *
     * @param matrix
     * @return matrix exponential of the matrix
     */
    public static Matrix matrixExp(Matrix matrix) {
        double c0 = 1.0;
        double c1 = 0.5;
        double c2 = 0.12;
        double c3 = 0.01833333333333333;
        double c4 = 0.0019927536231884053;
        double c5 = 1.630434782608695E-4;
        double c6 = 1.0351966873706E-5;
        double c7 = 5.175983436853E-7;
        double c8 = 2.0431513566525E-8;
        double c9 = 6.306022705717593E-10;
        double c10 = 1.4837700484041396E-11;
        double c11 = 2.5291534915979653E-13;
        double c12 = 2.8101705462199615E-15;
        double c13 = 1.5440497506703084E-17;
        int size = matrix.getNumCols();
        DMatrixRMaj result = new DMatrixRMaj(size, size);
        DMatrixRMaj As;
        DMatrixRMaj As_2 = new DMatrixRMaj(size, size);
        DMatrixRMaj As_4 = new DMatrixRMaj(size, size);
        DMatrixRMaj As_6 = new DMatrixRMaj(size, size);
        DMatrixRMaj U;
        DMatrixRMaj V;
        DMatrixRMaj AV = new DMatrixRMaj(size, size);
        DMatrixRMaj N = new DMatrixRMaj(size, size), D = new DMatrixRMaj(size, size);
        DMatrixRMaj temp;
        LinearSolverDense<DMatrixRMaj> solver = LinearSolverFactory_DDRM.linear(size);

        As = new DMatrixRMaj(matrix.toArray2D());
        int j =
                FastMath.max(
                        0, 1 + (int) FastMath.floor(Math.log(NormOps_DDRM.normPInf(As)) / FastMath.log(2)));

        CommonOps_DDRM.scale(1.0 / FastMath.pow(2, j), As); // scaled version of A

        // calculate D and N using special Horner techniques
        CommonOps_DDRM.mult(As, As, As_2);
        CommonOps_DDRM.mult(As_2, As_2, As_4);
        CommonOps_DDRM.mult(As_4, As_2, As_6);

        // U = c0*I + c2*A^2 + c4*A^4 + (c6*I + c8*A^2 + c10*A^4 + c12*A^6)*A^6
        U = CommonOps_DDRM.identity(size);
        CommonOps_DDRM.scale(c0, U);
        CommonOps_DDRM.addEquals(U, c2, As_2);
        CommonOps_DDRM.addEquals(U, c4, As_4);

        temp = CommonOps_DDRM.identity(size);
        CommonOps_DDRM.scale(c6, temp);
        CommonOps_DDRM.addEquals(temp, c8, As_2);
        CommonOps_DDRM.addEquals(temp, c10, As_4);
        CommonOps_DDRM.addEquals(temp, c12, As_6);

        CommonOps_DDRM.multAdd(temp, As_6, U);

        // V = c1*I + c3*A^2 + c5*A^4 + (c7*I + c9*A^2 + c11*A^4 + c13*A^6)*A^6
        V = CommonOps_DDRM.identity(size);
        CommonOps_DDRM.scale(c1, V);
        CommonOps_DDRM.addEquals(V, c3, As_2);
        CommonOps_DDRM.addEquals(V, c5, As_4);

        temp = CommonOps_DDRM.identity(size);
        CommonOps_DDRM.scale(c7, temp);
        CommonOps_DDRM.addEquals(temp, c9, As_2);
        CommonOps_DDRM.addEquals(temp, c11, As_4);
        CommonOps_DDRM.addEquals(temp, c13, As_6);

        CommonOps_DDRM.multAdd(temp, As_6, V);

        CommonOps_DDRM.mult(As, V, AV);
        CommonOps_DDRM.add(U, AV, N);
        CommonOps_DDRM.subtract(U, AV, D);

        // solve DF = N for F
        solver.setA(D);
        solver.solve(N, result);

        // now square j times
        for (int k = 0; k < j; k++) {
            temp = new DMatrixRMaj(result);
            CommonOps_DDRM.mult(temp, temp, result);
        }

        Matrix array = new Matrix(result.numRows, result.numCols);
        for (int i = 0; i < result.numRows; i++) {
            for (int i1 = 0; i1 < result.numCols; i1++) {
                array.set(i, i1, result.get(i, i1));
            }
        }
        return array;
    }

    // ================================================================================
    // SPECIAL FUNCTIONS
    // Special mathematical functions including factorials, gamma functions,
    // and logarithmic operations
    // ================================================================================

    /**
     * Computes the factorial of an integer.
     * 
     * @param n the input value
     * @return n! (n factorial)
     */
    public static int fact(int n) {
        return (int) Gamma.gamma((double) 1 + n);
    }

    /**
     * Computes the factorial of a double value using the gamma function.
     * 
     * @param n the input value
     * @return n! = Γ(n+1)
     */
    public static double fact(double n) {
        return Gamma.gamma(1 + n);
    }

    /**
     * Computes the natural logarithm of the factorial.
     * 
     * @param n the input value
     * @return ln(n!)
     */
    public static double factln(int n) {
        return FastMath.log(Gamma.gamma(1 + n));
    }

    /**
     * Computes the natural logarithm of the factorial for a double value.
     * 
     * @param n the input value  
     * @return ln(n!) = ln(Γ(n+1))
     */
    public static double factln(double n) {
        return FastMath.log(Gamma.gamma(1 + n));
    }

    /**
     * Computes the gamma function via Lanczos approximation formula.
     * The gamma function is an extension of the factorial function to real and complex numbers.
     * 
     * @param x the input value
     * @return Γ(x)
     */
    public static double gammaFunction(double x) {
        double[] p = {0.99999999999980993, 676.5203681218851, -1259.1392167224028, 771.32342877765313,
                -176.61502916214059, 12.507343278686905, -0.13857109526572012, 9.9843695780195716e-6,
                1.5056327351493116e-7};
        int g = 7;
        if (x < 0.5) return FastMath.PI / (Math.sin(Math.PI * x) * gammaFunction(1 - x));
        x -= 1;
        double a = p[0];
        double t = x + g + 0.5;
        for (int i = 1; i < p.length; i++) {
            a += p[i] / (x + i);
        }
        return FastMath.sqrt(2 * FastMath.PI) * FastMath.pow(t, x + 0.5) * FastMath.exp(-t) * a;
    }

    // ================================================================================
    // NUMERICAL INTEGRATION AND OPTIMIZATION
    // Functions for numerical integration, quadrature, optimization algorithms,
    // and root-finding methods
    // ================================================================================

    /**
     * Grundmann-Moeller simplex integration
     */
    public static simplexQuadResult grnmol(Function<double[], Double> f, double[][] V, int s, double tol) {
        int n = V.length;
        int d = 0;
        double[] Q = new double[s + 1];
        double[] Qv = new double[s + 1];
        double Vol = 1.0 / Maths.fact(n);
        int nv = 0;

        while (true) {
            int m = n + 2 * d + 1;
            double[] al = new double[n];
            Arrays.fill(al, 1);
            int alz = 2 * d + 1;
            double Qs = 0;

            while (true) {
                double[] alzal = new double[n + 1];
                alzal[0] = alz;
                System.arraycopy(al, 0, alzal, 1, n);
                double[] Valzal = new double[n];
                for (int i = 0; i < n; i++) {
                    for (int j = 0; j < n + 1; j++) {
                        Valzal[i] += V[i][j] * alzal[j] / m;
                    }
                }
                Qs += f.apply(Valzal);
                nv++;
                for (int j = 0; j < n; j++) {
                    alz -= 2;
                    if (alz > 0) {
                        al[j] += 2;
                        break;
                    }
                    alz += al[j] + 1;
                    al[j] = 1;
                }
                if (alz == 2 * d + 1) {
                    break;
                }
            }

            d++;
            Qv[d - 1] = Vol * Qs;
            Q[d - 1] = 0;
            double den = 1;
            for (int i = n + 1; i <= m; i++) {
                den *= 2 * i;
            }
            double p = 2.0 / den;
            for (int i = 1; i <= d; i++) {
                Q[d - 1] += FastMath.pow(m + 2 - 2 * i, 2 * d - 1) * p * Qv[d + 1 - i - 1];
                p = -p * (m + 1 - i) / i;
            }

            if (d > s || (d > 1 && FastMath.abs(Q[d - 1] - Q[d - 1 - 1]) < tol * Q[d - 1 - 1])) {
                return new simplexQuadResult(0.0, Arrays.copyOfRange(Q, 0, d), nv);
            }
        }
    }

    // ================================================================================
    // STATISTICAL ANALYSIS AND DATA PROCESSING
    // Functions for histograms, data binning, and statistical analysis
    // ================================================================================

    /**
     * Implementation of MATLAB "hist": puts elements of v into k bins
     *
     * @param v      - vector of elements
     * @param minVal - lowest number in v
     * @param maxVal - highest number in v
     * @return - vector containing number of elements in each bin
     */
    public static Matrix binHist(Matrix v, int minVal, int maxVal) {
        Matrix result = new Matrix(1, maxVal - minVal + 1);
        for (int j = 0; j < v.getNumCols(); j++) {
            int value = (int) v.get(j);
            if (value >= minVal && value <= maxVal) {
                int index = value - minVal;
                result.set(0, index, result.get(0, index) + 1);
            }
        }
        return result;
    }

    // ================================================================================
    // LAPLACE APPROXIMATION METHODS
    // Functions for Laplace approximation and related numerical methods
    // ================================================================================

    public static laplaceApproxReturn laplaceapprox_h(Matrix x0, SerializableFunction<Matrix, Matrix> h) {
        int d = x0.getNumCols();
        double tol = 1e-5;
        double detnH = -1;
        Matrix H = new Matrix(1, 1);
        H.set(0, 0, -1);
        while (detnH < 0 && tol <= 1e-3) {
            H = num_hess_h(x0, tol, h);
            Matrix nH = H.copy();
            nH.scaleEq(-1);
            detnH = nH.det();
            tol *= 10;
        }
        if (detnH < 0) {
            line_warning("Maths", "laplaceapprox_h: det(-H)<0.");
        }
        Matrix infrad = h.apply(x0);
        double I = infrad.get(0) * (Math.pow(2 * FastMath.PI, d) / FastMath.sqrt(detnH));
        double logI = FastMath.log(infrad.get(0)) + ((double) d / 2 * FastMath.log(2 * FastMath.PI)) - FastMath.log(detnH);
        return new laplaceApproxReturn(H, I, logI);
    }

    public static laplaceApproxComplexReturn laplaceapprox_h_complex(Matrix x0, SerializableFunction<Matrix, ComplexMatrix> h) {
        int d = x0.getNumCols();
        double tol = 1e-5;
        Complex detnH = new Complex(-1);
        ComplexMatrix H = new ComplexMatrix(1, 1);
        H.set(0, 0, -1);
        while (detnH.getReal() < 0 && tol <= 1e-3) {
            H = num_hess_h_complex(x0, tol, h);
            ComplexMatrix nH = H.copy();
            nH.scale(-1);
            detnH = nH.det();
            tol *= 10;
        }
        if (detnH.getReal() < 0) {
            line_warning("Maths", "laplaceapprox_h_complex: det(-H)<0.");
        }
        ComplexMatrix infrad = h.apply(x0);
        Complex I = infrad.get(0).multiply(new Complex(Math.pow(2 * FastMath.PI, d)).divide(detnH).sqrt());
        Complex logI = infrad.get(0).log().add((new Complex((double) d / 2 * FastMath.log(2 * FastMath.PI)))).subtract(detnH.log());
        return new laplaceApproxComplexReturn(H, I, logI);
    }

    // ================================================================================
    // NUMERICAL UTILITIES
    // Functions for generating sequences, spaces, and other numerical utilities
    // ================================================================================

    public static double[] linSpace(double start, double end, int num) {
        double[] result = new double[num];
        if (num == 1) {
            result[0] = end;
            return result;
        }
        double step = (end - start) / (num - 1);
        for (int i = 0; i < num; i++) {
            result[i] = start + i * step;
        }
        return result;
    }

    /**
     * Computes log(n choose k) = log Γ(n+1) − log Γ(k+1) − log Γ(n−k+1).
     * Uses Commons-Math Gamma but FastMath for everything else.
     *
     * @param n total items
     * @param k items chosen
     * @return log of the binomial coefficient, or NegInf if k<0 or k>n
     */
    public static double logBinomial(int n, int k) {
        if (k < 0 || k > n) {
            return NegInf;
        }
        return Gamma.logGamma(n + 1.0)
                - Gamma.logGamma(k + 1.0)
                - Gamma.logGamma((n - k) + 1.0);
    }

    public static double logmeanexp(Matrix x) {
        return Matrix.logsumexp(x) - FastMath.log(x.getNumElements());
    }

    public static double[] logSpace(double min, double max, int n) {
        // Generate a logarithmically spaced vector
        double[] y = new double[n];
        double[] x = linSpace(min, max, n);
        for (int i = 0; i < n; i++) {
            y[i] = FastMath.pow(10, x[i]);
        }
        return y;
    }

    /**
     * Generates an integer list of n logarithmically spaced values.
     *
     * @param start First element of array.
     * @param stop  Last element of array.
     * @param n     Number of values.
     * @return Array of logarithmically spaced values.
     */
    public static int[] logSpacei(int start, int stop, int n) {
        // I know this sucks but im making it java 7 compatible on the last possible day, was using lambdas
        ArrayList<Integer> vector = new ArrayList<>(n);
        int i = 0;
        for (double log : logSpace(start, stop, n)) {
            vector.add((int) FastMath.round(log));
        }
        vector = new ArrayList<Integer>(new LinkedHashSet<Integer>(vector));
        int[] result = new int[n];
        int j = 0;
        for (Integer integer : vector) {
            result[j] = integer;
            j++;
        }
        return Arrays.copyOf(result, j);
    }

    /**
     * Returns the max of two numbers. If one of the two is NaN, it returns the other.
     *
     * @param x - the first number to be compared
     * @param y - the second number to be compared
     * @return - the max between the two
     */
    public static double max(double x, double y) {
        if (!Double.isNaN(x) && !Double.isNaN(y))
            return FastMath.max(x, y);
        else if (Double.isNaN(x))
            return y;
        return x;
    }

    /**
     * Returns the min of two numbers. If one of the two is NaN, it returns the other.
     *
     * @param x - the first number to be compared
     * @param y - the second number to be compared
     * @return - the min between the two
     */
    public static double min(double x, double y) {
        if (!Double.isNaN(x) && !Double.isNaN(y))
            return FastMath.min(x, y);
        else if (Double.isNaN(x))
            return y;
        return x;
    }

    /**
     * Generate all combinations of distributing n items into k bins
     */
    public static List<Matrix> multiChooseCombinations(int k, int n) {
        List<Matrix> results = new ArrayList<Matrix>();
        int[] current = new int[k];
        multiChooseCombinationsHelper(k, n, 0, current, results);
        return results;
    }

    // Helper method for multiChooseCombinations
    private static void multiChooseCombinationsHelper(int k, int n, int index, int[] current, List<Matrix> results) {
        if (index == k - 1) {
            current[index] = n;
            Matrix row = new Matrix(1, k);
            for (int i = 0; i < k; i++) {
                row.set(0, i, current[i]);
            }
            results.add(row);
            return;
        }

        for (int i = 0; i <= n; i++) {
            current[index] = i;
            multiChooseCombinationsHelper(k, n - i, index + 1, current, results);
        }
    }

    public static Matrix multiChooseCon(Matrix n, double S) {
        // n row vector
        // n[i] is number of elements at ith position we have to choose from
        Matrix v = new Matrix(0, 0);
        List<Integer> indices = new ArrayList<>();
        for (int i = 0; i < n.getNumElements(); i++) {
            if (n.get(i) != 0) {
                indices.add(i);
            }
        }
        if (S == 1) {
            for (int i : indices) {
                v.expandMatrix(v.getNumRows() + 1, n.getNumCols(),
                        (v.getNumRows() + 1) * n.getNumCols());
                int lastRow = v.getNumRows() - 1;
                for (int col = 0; col < n.getNumCols(); col++) {
                    v.set(lastRow, col, 0);
                }
                v.set(lastRow, i, 1);
            }
            return v;
        } else {
            for (int i : indices) {
                Matrix n1 = n.copy();
                n1.set(0, i, n1.get(i) - 1);
                Matrix T = multiChooseCon(n1, S - 1);
                Matrix y = new Matrix(T.getNumRows(), n.getNumCols());
                y.zero();
                for (int row = 0; row < y.getNumRows(); row++) {
                    y.set(row, i, 1);
                }
                if (v.isEmpty()) {
                    v = y.add(1, T);
                } else {
                    v = Matrix.concatRows(v, y.add(1, T), null);
                }
            }
            return v;
        }
    }

    public static Matrix multichoose(double n, double k) {

        Matrix v = new Matrix(1, (int) n);
        v.zero();

        if (n == 1) {
            v = new Matrix(1, 1);
            v.set(0, 0, k);
        } else if (k != 0) {


            List<Matrix> tmpSSRows = new ArrayList<>();
            for (int i = 0; i <= k; i++) {
                Matrix w = multichoose(n - 1, k - i);

                Matrix tmpSSRow = new Matrix(w.getNumRows(), w.getNumCols() + 1);
                for (int j = 0; j < w.getNumRows(); j++) {
                    tmpSSRow.set(j, 0, i);
                    for (int l = 1; l < w.getNumCols() + 1; l++) {
                        tmpSSRow.set(j, l, w.get(j, l - 1));
                    }
                }
                tmpSSRows.add(tmpSSRow);
            }
            int totalRows = 0;
            for (Matrix m : tmpSSRows) {
                totalRows += m.getNumRows();
            }
            v = new Matrix(totalRows, (int) n);
            int rowForV = 0;
            for (int i = 0; i < tmpSSRows.size(); i++) {
                int rowForTmpSSRows = 0;
                int initialRowForV = rowForV;
                for (int j = rowForV; j < initialRowForV + tmpSSRows.get(i).getNumRows(); j++) {
                    for (int l = 0; l < tmpSSRows.get(i).getNumCols(); l++) {
                        v.set(j, l, tmpSSRows.get(i).get(rowForTmpSSRows, l));
                    }
                    rowForV++;
                    rowForTmpSSRows++;
                }
            }
        }

        return v;
    }

    public static double multinomialln(Matrix n) {
        return factln(n.elementSum()) - Matrix.factln(n).elementSum();
    }

    public static double nCk(double n, double k) {
        return FastMath.exp(factln(n) - factln(k) - factln(n - k));
    }

    /**
     * Computes the combinations of the elements in v taken k at a time
     *
     * @param v - vector of elements
     * @param k - how many elements to pick at a time
     * @return - the combinations of the elements in v taken k at a time
     */
    public static Matrix nCk(Matrix v, int k) {
        int n = v.length();
        if (k < 0 || k > n) {
            return null;
        }
        int a = 1, b = 1, c = 1; // a == n!, b == k!, c == (n-k)!
        for (int i = 2; i <= n; i++) {
            a *= i;
            if (i <= k) {
                b *= i;
            }
            if (i <= n - k) {
                c *= i;
            }
        }
        Matrix res = new Matrix(a / (b * c), k);
        int row = 0;
        int[] indexes = new int[k];
        for (int i = 0; i < indexes.length; i++) {
            indexes[i] = i;
        }
        while (true) {
            for (int i = 0; i < k; i++) {
                res.set(row, i, v.get(indexes[i]));
            }
            row++;
            int last = k - 1;
            while (last >= 0 && indexes[last] == n - k + last) {
                last--;
            }
            if (last == -1) {
                break;
            }
            indexes[last]++;
            for (int i = last + 1; i < k; i++) {
                indexes[i] = indexes[i - 1] + 1;
            }
        }
        return res;
    }

    // ================================================================================
    // HELPER METHODS AND UTILITIES
    // Internal utility functions and helper methods
    // ================================================================================

    /**
     * Given an integer x input, find the next integer y greater than x
     * such that y is a power of 2.
     *
     * @param x Input to find next power
     * @return Closest integer greater than input that is a power of two.
     */
    protected static int nextPowerOfTwo(int x) {
        int highestOneBit = Integer.highestOneBit(x);
        return (x == highestOneBit) ? x : highestOneBit << 1;
    }

    // ================================================================================
    // NUMERICAL DIFFERENTIATION
    // Functions for computing numerical gradients and Hessians
    // ================================================================================

    public static Matrix num_grad_h(Matrix x0, double h, SerializableFunction<Matrix, Matrix> hfun) {
        Matrix df = new Matrix(x0.copy());
        df.zero();
        for (int i = 0; i < x0.getNumCols(); i++) {
            Matrix x1 = x0.copy();
            Matrix x2 = x0.copy();
            x1.set(i, x0.get(i) - h);
            x2.set(i, x0.get(i) + h);
            double y1;
            double y2;
            y1 = FastMath.log(hfun.apply(x1).get(0));
            y2 = FastMath.log(hfun.apply(x2).get(0));
            df.set(i, y2 - y1 / (2 * h));
        }
        return df;
    }

    public static ComplexMatrix num_grad_h_complex(Matrix x0, double h, SerializableFunction<Matrix, ComplexMatrix> hfun) {
        ComplexMatrix df = new ComplexMatrix(x0.copy(), x0.copy());
        df.zero();
        for (int i = 0; i < x0.getNumCols(); i++) {
            Matrix x1 = x0.copy();
            Matrix x2 = x0.copy();
            x1.set(i, x0.get(i) - h);
            x2.set(i, x0.get(i) + h);
            Complex y1;
            Complex y2;
            y1 = hfun.apply(x1).get(0).log();
            y2 = hfun.apply(x2).get(0).log();
            df.set(i, y2.subtract(y1).divide(2 * h));
        }
        return df;
    }

    public static Matrix num_hess_h(Matrix x0, double h, SerializableFunction<Matrix, Matrix> hfun) {
        Matrix H = new Matrix(0, x0.getNumElements());
        H.zero();
        for (int i = 0; i < x0.getNumElements(); i++) {
            Matrix x1 = x0.copy();
            x1.set(i, x1.get(i) - h);
            Matrix df1 = num_grad_h(x1, h, hfun);

            Matrix x2 = x0.copy();
            x2.set(i, x2.get(i) + h);
            Matrix df2 = num_grad_h(x2, h, hfun);
            Matrix d2f = new Matrix(1, x0.getNumElements());
            for (int j = 0; j < d2f.getNumElements(); j++) {
                d2f.set(j, (df2.get(j) - df1.get(j)) / (2 * h));
            }
            H = Matrix.concatRows(H, d2f, null);
        }
        return H;
    }

    public static ComplexMatrix num_hess_h_complex(Matrix x0, double h, SerializableFunction<Matrix, ComplexMatrix> hfun) {
        int d = x0.getNumElements();
        ComplexMatrix H = new ComplexMatrix(d, d);
        H.zero();

        // Function to compute log(real(h(x))) - take real part before log to match MATLAB
        java.util.function.Function<Matrix, Complex> logH = (Matrix x) -> {
            Complex hval = hfun.apply(x).get(0);
            // Take real part like MATLAB's y(i) = real(h(x(i,:)))
            return new Complex(FastMath.log(hval.getReal()), 0.0);
        };

        // Evaluate at center point
        Complex f0 = logH.apply(x0);

        for (int i = 0; i < d; i++) {
            for (int j = 0; j < d; j++) {
                if (i == j) {
                    // Diagonal: use direct second-order central difference
                    // H[i,i] = (f(x+h*e_i) - 2*f(x) + f(x-h*e_i)) / h^2
                    Matrix xPlus = x0.copy();
                    Matrix xMinus = x0.copy();
                    xPlus.set(i, x0.get(i) + h);
                    xMinus.set(i, x0.get(i) - h);

                    Complex fPlus = logH.apply(xPlus);
                    Complex fMinus = logH.apply(xMinus);

                    // H[i,i] = (fPlus - 2*f0 + fMinus) / h^2
                    Complex hii = fPlus.subtract(f0.multiply(2)).add(fMinus).divide(h * h);
                    H.set(i, i, hii);
                } else {
                    // Off-diagonal: use four-point stencil for mixed partials
                    // H[i,j] = (f(x+h*e_i+h*e_j) - f(x+h*e_i-h*e_j) - f(x-h*e_i+h*e_j) + f(x-h*e_i-h*e_j)) / (4*h^2)
                    Matrix xpp = x0.copy();
                    Matrix xpm = x0.copy();
                    Matrix xmp = x0.copy();
                    Matrix xmm = x0.copy();

                    xpp.set(i, x0.get(i) + h);
                    xpp.set(j, x0.get(j) + h);
                    xpm.set(i, x0.get(i) + h);
                    xpm.set(j, x0.get(j) - h);
                    xmp.set(i, x0.get(i) - h);
                    xmp.set(j, x0.get(j) + h);
                    xmm.set(i, x0.get(i) - h);
                    xmm.set(j, x0.get(j) - h);

                    Complex fpp = logH.apply(xpp);
                    Complex fpm = logH.apply(xpm);
                    Complex fmp = logH.apply(xmp);
                    Complex fmm = logH.apply(xmm);

                    // H[i,j] = (fpp - fpm - fmp + fmm) / (4*h^2)
                    Complex hij = fpp.subtract(fpm).subtract(fmp).add(fmm).divide(4 * h * h);
                    H.set(i, j, hij);
                }
            }
        }

        // Ensure symmetry
        for (int i = 0; i < d; i++) {
            for (int j = i + 1; j < d; j++) {
                Complex avg = H.get(i, j).add(H.get(j, i)).divide(2);
                H.set(i, j, avg);
                H.set(j, i, avg);
            }
        }

        return H;
    }

    public static Matrix permutations(Matrix vec) {
        // copy to avoid issues with needing vec in outer code as swap modifies inputted vector
        Matrix vecCopy = vec.copy();
        // use length to accept row or column vector
        Matrix indexes = new Matrix(1, vecCopy.length());
        indexes.zero();
        List<Matrix> lPermutations = new ArrayList<>();
        lPermutations.add(vecCopy.copy());

        int ind = 0;
        while (ind < vecCopy.length()) {
            if (indexes.get(ind) < ind) {
                swap(vecCopy, ind % 2 == 0 ? 0 : (int) indexes.get(ind), ind);
                lPermutations.add(vecCopy.copy());
                indexes.set(ind, indexes.get(ind) + 1);
                ind = 0;
            } else {
                indexes.set(ind, 0);
                ind++;
            }
        }
        lPermutations.sort((a, b) -> {
            for (int i = 0; i < a.length(); i++) {
                if (a.get(i) != b.get(i)) {
                    return Double.compare(b.get(i), a.get(i));
                }
            }
            return 0;
        });


        // Construct permuations matrix from perm vectors
        Matrix permutations = new Matrix(lPermutations.size(), vec.length());
        for (int i = 0; i < lPermutations.size(); i++) {
            Matrix permutation = lPermutations.get(i);
            for (int j = 0; j < permutations.getNumCols(); j++) {
                permutations.set(i, j, permutation.get(j));
            }
        }
        return permutations;
    }

    // ================================================================================
    // RANDOM NUMBER GENERATION FUNCTIONS
    // Core random number generation methods
    // ================================================================================

    // returns a random double in interval [0,1)
    public static double rand() {
        return RandomManager.nextDouble();
    }

    // returns a random Gaussian-distributed double
    public static double randn() {
        return RandomManager.nextGaussian();
    }

    // ================================================================================
    // PROBABILITY AND STOCHASTIC FUNCTIONS
    // Functions for probabilistic choices and stochastic operations
    // ================================================================================

    /**
     * Choose an element index according to probability vector.
     * Equivalent to MATLAB's probchoose function.
     *
     * @param p Probability vector (must sum to 1.0)
     * @return Index of chosen element (0-based)
     */
    public static int probchoose(List<Double> p) {
        if (p == null || p.isEmpty()) {
            throw new IllegalArgumentException("Probability vector cannot be null or empty");
        }
        
        // Compute cumulative sum
        List<Double> cumsum = new ArrayList<Double>(p.size());
        cumsum.add(p.get(0));
        for (int i = 1; i < p.size(); i++) {
            cumsum.add(cumsum.get(i - 1) + p.get(i));
        }
        
        // Generate random number
        double r = rand();
        
        // Find first position where r <= cumsum[i]
        for (int i = 0; i < cumsum.size(); i++) {
            if (r <= cumsum.get(i)) {
                return i;
            }
        }
        
        // If we get here due to numerical precision, return last element
        return p.size() - 1;
    }

    /**
     * Choose an element index according to probability vector.
     * Overloaded version for Matrix input.
     *
     * @param p Probability matrix (single row or column)
     * @return Index of chosen element (0-based)
     */
    public static int probchoose(Matrix p) {
        if (p == null || p.isEmpty()) {
            throw new IllegalArgumentException("Probability matrix cannot be null or empty");
        }
        
        List<Double> pList = new ArrayList<Double>();
        for (int i = 0; i < p.length(); i++) {
            pList.add(p.get(i));
        }
        return probchoose(pList);
    }

    // ================================================================================
    // NUMBER THEORY FUNCTIONS
    // Functions for prime numbers and other number-theoretic operations
    // ================================================================================

    public static List<Integer> primes(int N) {
        // Generate prime numbers up to N
        List<Integer> primes = new ArrayList<>();
        for (int i = 2; i <= N; i++) {
            if (isPrime(i)) {
                primes.add(i);
            }
        }
        return primes;
    }


    // ================================================================================
    // ARRAY POSITION AND INDEXING FUNCTIONS
    // Functions for finding positions of max/min values in arrays
    // ================================================================================

    /**
     * Returns the position of the maximum value in a vector.
     * Equivalent to MATLAB's maxpos function.
     *
     * @param v Input vector
     * @return Index of maximum value (0-based)
     */
    public static int maxpos(Matrix v) {
        if (v == null || v.isEmpty()) {
            throw new IllegalArgumentException("Input matrix cannot be null or empty");
        }
        
        double maxVal = NegInf;
        int maxIdx = 0;
        
        for (int i = 0; i < v.length(); i++) {
            if (v.get(i) > maxVal) {
                maxVal = v.get(i);
                maxIdx = i;
            }
        }
        
        return maxIdx;
    }

    /**
     * Returns the positions of the n largest values in a vector.
     * Equivalent to MATLAB's maxpos(v, n) function.
     *
     * @param v Input vector
     * @param n Number of positions to return
     * @return Array of indices of n largest values (0-based), sorted by value (descending)
     */
    public static int[] maxpos(Matrix v, int n) {
        if (v == null || v.isEmpty()) {
            throw new IllegalArgumentException("Input matrix cannot be null or empty");
        }
        if (n < 1 || n > v.length()) {
            throw new IllegalArgumentException("n must be between 1 and vector length");
        }

        int[] result = new int[n];
        boolean[] used = new boolean[v.length()];
        
        for (int count = 0; count < n; count++) {
            double maxVal = NegInf;
            int maxIdx = 0;
            
            for (int i = 0; i < v.length(); i++) {
                if (!used[i] && v.get(i) > maxVal) {
                    maxVal = v.get(i);
                    maxIdx = i;
                }
            }
            
            result[count] = maxIdx;
            used[maxIdx] = true;
        }
        
        return result;
    }

    /**
     * Returns the position of the minimum value in a vector.
     * Equivalent to MATLAB's minpos function.
     *
     * @param v Input vector
     * @return Index of minimum value (0-based)
     */
    public static int minpos(Matrix v) {
        if (v == null || v.isEmpty()) {
            throw new IllegalArgumentException("Input matrix cannot be null or empty");
        }
        
        double minVal = Inf;
        int minIdx = 0;
        
        for (int i = 0; i < v.length(); i++) {
            if (v.get(i) < minVal) {
                minVal = v.get(i);
                minIdx = i;
            }
        }
        
        return minIdx;
    }

    /**
     * Returns the positions of the n smallest values in a vector.
     * Equivalent to MATLAB's minpos(v, n) function.
     *
     * @param v Input vector
     * @param n Number of positions to return
     * @return Array of indices of n smallest values (0-based), sorted by value (ascending)
     */
    public static int[] minpos(Matrix v, int n) {
        if (v == null || v.isEmpty()) {
            throw new IllegalArgumentException("Input matrix cannot be null or empty");
        }
        if (n < 1 || n > v.length()) {
            throw new IllegalArgumentException("n must be between 1 and vector length");
        }

        int[] result = new int[n];
        boolean[] used = new boolean[v.length()];
        
        for (int count = 0; count < n; count++) {
            double minVal = Inf;
            int minIdx = 0;
            
            for (int i = 0; i < v.length(); i++) {
                if (!used[i] && v.get(i) < minVal) {
                    minVal = v.get(i);
                    minIdx = i;
                }
            }
            
            result[count] = minIdx;
            used[minIdx] = true;
        }
        
        return result;
    }

    public static boolean randomMatlabStyle() {
        return true; // Always true now - all random generation uses MersenneTwister through RandomManager
    }

    // ================================================================================
    // POLYNOMIAL AND ROOT-FINDING FUNCTIONS
    // Functions for polynomial operations and finding roots
    // ================================================================================

    public static Complex[] roots(Matrix coefficients) {
        return roots(new Matrix(coefficients));
    }

    public static Complex[] roots(double[] coefficients) {
        // This may be alternatively implemented using EJML as in here:
        // https://ejml.org/wiki/index.php?title=Example_Polynomial_Roots

        // Remove leading zeros from the coefficients array
        int firstNonZeroIndex = 0;
        while (firstNonZeroIndex < coefficients.length && coefficients[firstNonZeroIndex] == 0) {
            firstNonZeroIndex++;
        }

        if (firstNonZeroIndex == coefficients.length) {
            // All coefficients are zero
            throw new IllegalArgumentException("All coefficients are zero, no polynomial to solve.");
        }

        // Create a trimmed array of coefficients
        double[] trimmedCoefficients = new double[coefficients.length - firstNonZeroIndex];
        System.arraycopy(coefficients, firstNonZeroIndex, trimmedCoefficients, 0, trimmedCoefficients.length);

        PolynomialFunction polynomialFunction = new PolynomialFunction(trimmedCoefficients);
        LaguerreSolver solver = new LaguerreSolver();
        Complex[] complexRoots = solver.solveAllComplex(trimmedCoefficients, 0);

        return complexRoots;
    }

    // ================================================================================
    // RANDOM SEED MANAGEMENT
    // Functions for managing random number generator seeds
    // ================================================================================

    @SuppressWarnings("deprecation")
    public static void setMatlabRandomSeed(final long seed) {
        // For backward compatibility, use ThreadLocalRandom directly
        ThreadLocalRandom.setSeed((int) seed);
    }

    public static void setRandomNumbersMatlab(boolean matlab_style) {
        // No-op: always uses MersenneTwister through RandomManager now
        // Method retained for backward compatibility
    }

    /**
     * Computes a specialized function used in simplex-based optimization algorithms.
     * This function evaluates an exponential transformation of the input variables
     * combined with matrix operations involving L and N matrices.
     *
     * <p>The function performs the following computation:</p>
     * <ol>
     *   <li>Creates a vector v where v[i] = exp(x[i]) for i = 0..M-2, and v[M-1] = 1</li>
     *   <li>Computes tmp = log(v * L)^T (element-wise logarithm of matrix product)</li>
     *   <li>Returns exp(N * tmp + sum(x) - (sum(N) + M) * log(sum(v)))</li>
     * </ol>
     *
     * <p>This function is typically used in optimization contexts where the variables
     * are constrained to a simplex and require exponential parameterization.</p>
     *
     * @param x The input vector of optimization variables (length M-1)
     * @param L The transformation matrix for the exponential variables
     * @param N The weighting vector for the final computation
     * @return The computed function value as a double
     * @throws IllegalArgumentException if matrix dimensions are incompatible
     */
    public static double simplex_fun(double[] x, Matrix L, Matrix N) {
        int M = x.length + 1;
        Matrix v = new Matrix(1, M);
        for (int i = 0; i < M - 1; i++) {
            v.set(i, FastMath.exp(x[i]));
        }
        v.set(M - 1, 1.0);

        Matrix tmp = (v.mult(L)).transpose();
        for (int i = 0; i < tmp.length(); i++) {
            tmp.set(i, FastMath.log(tmp.get(i)));
        }

        double x_sum = 0.0;
        for (int i = 0; i < x.length; i++) {
            x_sum += x[i];
        }

        double f = FastMath.exp(N.mult(tmp).elementSum() + x_sum
                - (N.elementSum() + M) * FastMath.log(v.elementSum()));

        return f;
    }

    public static simplexQuadResult simplexquad(Function<double[], Double> f, int n, int order, double tol) {
        Matrix Iplus = Matrix.eye(n);
        Iplus = Iplus.concatCols(Matrix.zeros(n, 1));
        simplexQuadResult grnmolResult = grnmol(f, Iplus.toArray2D(), order, tol);
        double I = grnmolResult.Q[grnmolResult.Q.length - 1];
        return new simplexQuadResult(I, grnmolResult.Q, grnmolResult.nv);
    }

    /**
     * Softmin function.
     *
     * @param x     first term to compare
     * @param y     second term to compare
     * @param alpha softmin smoothing parameter
     * @return Softmin function value
     */
    public static double softmin(double x, double y, double alpha) {
        return -((-x) * FastMath.exp(-alpha * x) - y * FastMath.exp(-alpha * y))
                / (Math.exp(-alpha * x) + FastMath.exp(-alpha * y));
    }

    // ================================================================================
    // INDEX CONVERSION FUNCTIONS
    // Functions for converting between matrix subscripts and linear indices
    // ================================================================================

    public static List<Integer> sub2ind(Matrix dims, Matrix row, Matrix col) {
        List<Integer> result = new ArrayList<>();
        int numRows = (int) dims.get(0);
        int numCols = (int) dims.get(1);
        for (int i = 0; i < row.getNumElements(); i++) {

            int index = (int) (col.get(i) * numRows + row.get(i));
            result.add(index);
        }
        return result;
    }

    // Helper method that swaps elements in a matrix vector around
    private static void swap(Matrix vec, int a, int b) {
        double tmp = vec.get(a);
        vec.set(a, vec.get(b));
        vec.set(b, tmp);
    }

    // ================================================================================
    // ARRAY MANIPULATION UTILITIES
    // Functions for array operations, transposition, and data structure conversions
    // ================================================================================

    /**
     * Transposes a 2D array of values.
     *
     * @param data 2D array to be transposed.
     * @return New 2D array of transposed data.
     */
    public static double[][] transpose(double[][] data) {
        int n = data.length;
        int m = data[0].length;
        double[][] transpose = new double[m][n];
        for (int r = 0; r < m; r++) {
            for (int c = 0; c < n; c++) {
                transpose[r][c] = data[c][r];
            }
        }
        return transpose;
    }

    // ================================================================================
    // SORTING AND UNIQUENESS FUNCTIONS
    // Functions for sorting and finding unique elements
    // ================================================================================

    // "unique" function in matlab sorts in lexicographic, replicated here
    public static Matrix uniqueAndSort(Matrix space) {
        List<Matrix> uniqueRows = new ArrayList<>();
        for (int i = 0; i < space.getNumRows(); i++) {
            Matrix tmp = new Matrix(1, space.getNumCols());
            Matrix tmp2 = new Matrix(1, space.getNumCols());
            Matrix.extractRows(space, i, i + 1, tmp);
            boolean unique = true;
            for (int j = i + 1; j < space.getNumRows(); j++) {
                Matrix.extractRows(space, j, j + 1, tmp2);
                if (tmp.isEqualTo(tmp2)) {
                    unique = false;
                }
            }
            if (unique) {
                uniqueRows.add(tmp);
            }
        }

        Comparator<Matrix> lexico = (mat1, mat2) -> {
            for (int col = 0; col < mat1.getNumCols(); col++) {
                int result = Integer.compare((int) mat1.get(0, col), (int) mat2.get(0, col));
                if (result != 0) {
                    return result;
                }
            }
            return 0;
        };

        uniqueRows.sort(lexico);
        Matrix newSpace = new Matrix(uniqueRows.size(), uniqueRows.get(0).getNumCols());
        // So that states with jobs in phase 1 comes earlier
        for (int i = 0; i < uniqueRows.size(); i++) {
            for (int j = 0; j < uniqueRows.get(0).getNumCols(); j++) {
                newSpace.set(i, j, uniqueRows.get(i).get(0, j));
            }
        }

        return newSpace;
    }

    public static Matrix uniquePerms(Matrix vec) {

        // Vector is empty
        if (vec.isEmpty()) {
            return new Matrix(0, 0);
        }

        // Vector is not a single row
        if (vec.getNumRows() != 1) {
            throw new RuntimeException("Matrix passed to uniquePerms has more than one row. Unsupported.");
        }

        // Number of elements in the vector
        int n = vec.length();

        // Number of unique elements in the vector
        int nu = 0;
        LinkedList<Double> uniqueElements = new LinkedList<>();
        for (int i = 0; i < vec.length(); i++) {
            double element = vec.get(0, i);
            boolean unique = true;
            for (Double uniqueElement : uniqueElements) {
                if (element == uniqueElement) {
                    unique = false;
                    break;
                }
            }
            if (unique) {
                uniqueElements.add(element);
            }
        }
        nu = uniqueElements.size();

        // Only one unique element
        if (nu == 1) {
            return vec;
        }

        // Every element is unique
        if (n == nu) {
            return permutations(vec);
        }
        Matrix[] output = new Matrix[nu];
        for (int i = 0; i < nu; i++) {
            Matrix v = vec.copy();

            int ind = -1;
            double target = uniqueElements.get(i);
            for (int j = 0; j < v.length(); j++) {
                if (target == v.get(j)) {
                    ind = j;
                    break;
                }
            }
            Matrix newV = new Matrix(1, v.getNumCols() - 1);
            for (int j = 0, destCol = 0; j < v.getNumCols(); j++) {
                if (j != ind) {
                    newV.set(0, destCol++, v.get(j));
                }
            }
            v = newV.copy();
            Matrix temp = uniquePerms(v);
            Matrix repeated = new Matrix(temp.getNumRows(), 1);
            repeated.fill(uniqueElements.get(i));
            output[i] = repeated.concatCols(temp);
        }

        Matrix result = output[0];
        for (int i = 1; i < output.length; i++) {
            result = Matrix.concatRows(result, output[i], null).copy();
        }

        return result;
    }

    /**
     * Set difference operation for Matrix objects that are vectors.
     * Returns elements in A that are not in B.
     * 
     * @param A First vector (universe of values to check)
     * @param B Second vector (values to exclude from A)
     * @return Matrix containing elements in A that are not in B, or empty matrix if inputs are invalid
     */
    public static Matrix setdiff(Matrix A, Matrix B) {
        if (A == null || A.isEmpty()) {
            return new Matrix(0, 0);
        }
        
        // Check that both are vectors (either row or column vectors)
        boolean aIsVector = (A.getNumRows() == 1 || A.getNumCols() == 1);
        boolean bIsVector = (B.getNumRows() == 1 || B.getNumCols() == 1);
        
        if (!aIsVector || !bIsVector) {
            line_error("setdiff", "Both inputs must be vectors (either row or column vectors)");
            return new Matrix(0, 0);
        }
        
        // Check that vectors have identical sizes
        int aLength = A.length();
        int bLength = B.length();
        
        if (aLength != bLength) {
            line_error("setdiff", "Vectors must have identical sizes. A has " + aLength + " elements, B has " + bLength + " elements");
            return new Matrix(0, 0);
        }
        
        // Create list to store difference values
        List<Double> diffValues = new ArrayList<Double>();
        
        // Check each element in A
        for (int i = 0; i < aLength; i++) {
            double aValue = A.get(i);
            boolean found = false;
            
            // Check if this value exists in B
            for (int j = 0; j < bLength; j++) {
                if (aValue == B.get(j)) {
                    found = true;
                    break;
                }
            }
            
            // If not found in B, add to result
            if (!found) {
                diffValues.add(aValue);
            }
        }
        
        // Convert result to Matrix
        if (diffValues.isEmpty()) {
            return new Matrix(0, 0);
        }
        
        Matrix result = new Matrix(diffValues.size(), 1);
        for (int i = 0; i < diffValues.size(); i++) {
            result.set(i, 0, diffValues.get(i));
        }
        
        return result;
    }

    /**
     * Generates all combinations of n objects taken k at a time with repetition (multichoose).
     * Returns a Matrix for compatibility with existing code.
     *
     * @param n number of distinct objects
     * @param k total number of selections
     * @return Matrix where each row is a combination
     */
    public static Matrix multichoose(int n, int k) {
        List<int[]> combinations = multichooseList(n, k);
        if (combinations.isEmpty()) {
            return new Matrix(0, n);
        }
        
        Matrix result = new Matrix(combinations.size(), n);
        for (int i = 0; i < combinations.size(); i++) {
            int[] combo = combinations.get(i);
            for (int j = 0; j < n; j++) {
                result.set(i, j, combo[j]);
            }
        }
        return result;
    }
    
    /**
     * Generates all combinations of n objects taken k at a time with repetition (multichoose).
     * Each combination is an array where element i represents how many times object i is chosen.
     *
     * @param n number of distinct objects
     * @param k total number of selections
     * @return list of combinations as int arrays
     */
    public static List<int[]> multichooseList(int n, int k) {
        if (n == 1) {
            return Arrays.asList(new int[]{k});
        }
        if (k == 0) {
            int[] zeros = new int[n];
            return Arrays.asList(zeros);
        }

        List<int[]> result = new ArrayList<>();
        for (int i = 0; i <= k; i++) {
            List<int[]> subResults = multichooseList(n - 1, k - i);
            for (int[] sub : subResults) {
                int[] combination = new int[n];
                combination[0] = i;
                System.arraycopy(sub, 0, combination, 1, sub.length);
                result.add(combination);
            }
        }
        return result;
    }

    /**
     * Sorts combinations by number of non-zeros and their positions.
     *
     * @param combinations list of combinations to sort
     * @return sorted list of combinations
     */
    public static List<int[]> sortByNnzPos(List<int[]> combinations) {
        List<int[]> sorted = new ArrayList<>(combinations);
        sorted.sort(new NnzComparator());
        return sorted;
    }

    /**
     * Compares two arrays by number of zeros and their positions.
     * Returns 1 if a < b, -1 if a > b, 0 if equal.
     *
     * @param i1 first array
     * @param i2 second array
     * @return comparison result
     */
    public static int nnzcmp(int[] i1, int[] i2) {
        int nnz1 = countNonZeros(i1);
        int nnz2 = countNonZeros(i2);

        if (nnz1 > nnz2) {
            return 1; // i2 has more zeros and is thus greater
        } else if (nnz1 < nnz2) {
            return -1; // i1 has more zeros and is thus greater
        } else { // nnz1 == nnz2
            for (int j = 0; j < i1.length; j++) {
                if (i1[j] == 0 && i2[j] > 0) {
                    return 1; // i2 has the left-most zero
                } else if (i1[j] > 0 && i2[j] == 0) {
                    return -1;
                }
            }
            return 0;
        }
    }

    /**
     * Counts the number of non-zero elements in an array.
     *
     * @param array the array to count non-zeros in
     * @return number of non-zero elements
     */
    private static int countNonZeros(int[] array) {
        int count = 0;
        for (int val : array) {
            if (val != 0) count++;
        }
        return count;
    }

    /**
     * Finds the index of a matching row in a list.
     *
     * @param rows list of integer arrays
     * @param target array to find
     * @return index of matching row, or -1 if not found
     */
    public static int matchRow(List<int[]> rows, int[] target) {
        for (int i = 0; i < rows.size(); i++) {
            if (Arrays.equals(rows.get(i), target)) {
                return i;
            }
        }
        return -1;
    }

    /**
     * Computes binomial coefficient "n choose k".
     *
     * @param n total items
     * @param k items chosen
     * @return binomial coefficient
     */
    public static int nchoosek(int n, int k) {
        if (k > n) return 0;
        if (k == 0 || k == n) return 1;

        // Use integer arithmetic to avoid precision issues
        long result = 1;
        for (int i = 0; i < k; i++) {
            result = result * (n - i) / (i + 1);
        }
        return (int) result;
    }

    // ================================================================================
    // HELPER CLASSES AND DATA STRUCTURES
    // Supporting classes for returning complex results from mathematical functions
    // ================================================================================

    public static class laplaceApproxComplexReturn {
        public ComplexMatrix H;
        public Complex I;
        public Complex logI;

        public laplaceApproxComplexReturn(ComplexMatrix H, Complex I, Complex logI) {
            this.H = H;
            this.I = I;
            this.logI = logI;
        }
    }

    public static class laplaceApproxReturn {
        public Matrix H;
        public double I;
        public double logI;

        public laplaceApproxReturn(Matrix H, double I, double logI) {
            this.H = H;
            this.I = I;
            this.logI = logI;
        }
    }

    /**
     * Golden section search for finding the minimum of a unimodal function.
     * This is equivalent to Python's scipy.optimize.fminbound method.
     * 
     * @param function The function to minimize
     * @param a Lower bound of search interval
     * @param b Upper bound of search interval
     * @param tol Tolerance for convergence (default: 1e-5)
     * @return Array containing [optimal_x, optimal_value]
     */
    public static double[] goldenSectionSearch(Function<Double, Double> function, double a, double b, double tol) {
        final double gr = (Math.sqrt(5.0) + 1) / 2;  // Golden ratio
        
        double a_local = a;
        double b_local = b;
        double c = b_local - (b_local - a_local) / gr;
        double d = a_local + (b_local - a_local) / gr;
        
        while (Math.abs(b_local - a_local) > tol) {
            if (function.apply(c) < function.apply(d)) {
                b_local = d;
            } else {
                a_local = c;
            }
            
            c = b_local - (b_local - a_local) / gr;
            d = a_local + (b_local - a_local) / gr;
        }
        
        double p_opt = (b_local + a_local) / 2;
        return new double[]{p_opt, function.apply(p_opt)};
    }
    
    /**
     * Golden section search with default tolerance.
     * 
     * @param function The function to minimize
     * @param a Lower bound of search interval
     * @param b Upper bound of search interval
     * @return Array containing [optimal_x, optimal_value]
     */
    public static double[] goldenSectionSearch(Function<Double, Double> function, double a, double b) {
        return goldenSectionSearch(function, a, b, 1e-5);
    }

    public static class simplexQuadResult {
        public double I;
        public double[] Q;
        public int nv;

        public simplexQuadResult(double I, double[] Q, int nv) {
            this.I = I;
            this.Q = Q;
            this.nv = nv;
        }
    }

    /**
     * Comparator for sorting combinations by number of non-zeros and their positions.
     */
    private static class NnzComparator implements Comparator<int[]> {
        @Override
        public int compare(int[] a, int[] b) {
            return nnzcmp(a, b);
        }
    }
}