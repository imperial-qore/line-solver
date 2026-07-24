/**
 * @file M/G/1-type Stationary Distribution
 *
 * Computes the stationary distribution for M/G/1-type Markov chains.
 *
 * Based on the SMC Solver implementation by Benny Van Houdt.
 *
 * @since LINE 3.1.0
 */
package jline.lib.smc;

import jline.util.matrix.Matrix;

public final class MG1_pi {

    private MG1_pi() {}

    public static Matrix mg1_pi(Matrix B, Matrix A) {
        return mg1_pi(B, A, new MG1PiOptions());
    }

    public static Matrix mg1_pi(Matrix B, Matrix A, MG1PiOptions options) {
        int m = A.getNumRows();
        int dega = A.getNumCols() / m - 1;

        // Use boundary or default to A's first row structure
        Matrix boundary = (B != null) ? B : A;

        // First, compute G using appropriate solver
        Matrix G;
        if ("CR".equals(options.getSolver().toUpperCase())) {
            G = mg1_cr(A, new MG1CROptions(options.getMode(),
                    50, 2048, 1e-16,
                    options.getVerbose() ? 1 : 0,
                    "one"));
        } else {
            G = MG1_FI.mg1_fi(A, new MG1FIOptions(options.getVerbose() ? 1 : 0));
        }

        // Compute R matrix: R = sum(i=0 to dega) A_i * G^i
        Matrix R = A.extractCols(0, m).copy();
        Matrix Gpow = G.copy();
        for (int i = 1; i <= dega; i++) {
            R = R.add(A.extractCols(i * m, (i + 1) * m).mult(Gpow));
            Gpow = Gpow.mult(G);
        }

        // Compute stationary distribution
        int degb = boundary.getNumCols() / m - 1;

        Matrix Bhat = boundary.extractCols(0, m).copy();
        Gpow = G.copy();
        for (int i = 1; i <= degb; i++) {
            Bhat = Bhat.add(boundary.extractCols(i * m, (i + 1) * m).mult(Gpow));
            Gpow = Gpow.mult(G);
        }

        // Compute pi_0 as stationary distribution of Bhat
        Matrix pi0 = Stat.stat(Bhat);

        // Compute higher levels using pi_i = pi_0 * R^i
        java.util.List<Matrix> result = new java.util.ArrayList<Matrix>();
        result.add(pi0);

        Matrix Rpow = R.copy();
        for (int i = 1; i < options.getMaxNumComp(); i++) {
            Matrix pi_i = pi0.mult(Rpow);

            double mass = pi_i.elementSum();
            if (mass < 1e-15) {
                break;
            }

            result.add(pi_i);
            Rpow = Rpow.mult(R);
        }

        // Normalize the distribution
        double totalMass = 0.0;
        for (Matrix pi_i : result) {
            totalMass += pi_i.elementSum();
        }

        // Return as single row vector
        int totalSize = result.size() * m;
        Matrix piVec = new Matrix(1, totalSize);
        for (int i = 0; i < result.size(); i++) {
            for (int j = 0; j < m; j++) {
                piVec.set(0, i * m + j, result.get(i).get(0, j) / totalMass);
            }
        }

        return piVec;
    }

    /**
     * Computes the G matrix using Cyclic Reduction for M/G/1-type Markov Chains.
     */
    public static Matrix mg1_cr(Matrix A) {
        return mg1_cr(A, new MG1CROptions());
    }

    public static Matrix mg1_cr(Matrix A, MG1CROptions options) {
        int m = A.getNumRows();

        // Check whether G is known explicitly
        Matrix explicitG = MG1_EG.mg1_eg(A, options.getVerbose() > 0);
        if (explicitG != null) {
            return explicitG;
        }

        Matrix Dold = null;
        Matrix workA = A.copy();
        double drift = 0.0;
        double tau = 0.0;
        Matrix v = Matrix.zeros(m, 1);

        if (options.getMode().contains("ShiftPWCR")) {
            if (options.getVerbose() == 1) {
                Dold = workA.copy();
            }
            jline.lib.smc.Quadruple<Matrix, Double, Double, Matrix> shiftResult =
                    MG1_Shifts.mg1_shifts(workA, options.getShiftType());
            workA = shiftResult.getFirst();
            drift = shiftResult.getSecond();
            tau = shiftResult.getThird();
            v = shiftResult.getFourth();
        }

        // Start Cyclic Reduction
        int dega = workA.getNumCols() / m - 1;
        int paddedBlocks = (1 << (1 + floorLog2(dega))) + 1;
        Matrix D = new Matrix(paddedBlocks * m, m);
        // D = transpose of workA, stored as column blocks
        for (int i = 0; i <= dega; i++) {
            for (int r = 0; r < m; r++) {
                for (int c = 0; c < m; c++) {
                    D.set(i * m + r, c, workA.get(c, i * m + r));
                }
            }
        }

        int totalBlocks = paddedBlocks;
        Matrix Aeven = extractEvenBlocks(D, m, totalBlocks);
        Matrix Aodd = extractOddBlocks(D, m, totalBlocks);

        Matrix Ahatodd = new Matrix((Aeven.getNumRows() / m - 1 + 1) * m, m);
        for (int i = 1; i < Aeven.getNumRows() / m; i++) {
            for (int r = 0; r < m; r++) {
                for (int c = 0; c < m; c++) {
                    Ahatodd.set((i - 1) * m + r, c, Aeven.get(i * m + r, c));
                }
            }
        }
        for (int r = 0; r < m; r++) {
            for (int c = 0; c < m; c++) {
                Ahatodd.set((Aeven.getNumRows() / m - 1) * m + r, c, D.get((totalBlocks - 1) * m + r, c));
            }
        }
        Matrix Ahateven = Aodd.copy();

        // Compute Rj for stop criteria (PWCR mode)
        Matrix Rj = new Matrix(m, m);
        for (int i = 1; i < D.getNumRows() / m; i++) {
            Matrix block = D.extractRows(i * m, (i + 1) * m);
            Rj = Rj.add(block);
        }
        Rj = Matrix.eye(m).sub(Rj).inv().mult(D.extractRows(0, m));

        Matrix G = Matrix.zeros(m, m);
        int numit = 0;

        while (numit < options.getMaxNumIt()) {
            numit++;
            int nj = Aodd.getNumRows() / m - 1;

            Matrix Anew;
            Matrix Ahatnew;

            if (nj > 0) {
                int n = nj + 1;
                double omega = -2.0 * Math.PI / n;

                double[][][] ft1 = complexDftBlocks(Aodd, m, n, Aodd.getNumRows() / m);
                double[][][] ft2 = complexDftBlocks(Aeven, m, n, Aeven.getNumRows() / m);
                double[][][] ft3 = complexDftBlocks(Ahatodd, m, n, Ahatodd.getNumRows() / m);
                double[][][] ft4 = complexDftBlocks(Ahateven, m, n, Ahateven.getNumRows() / m);

                double[][][] ftAnew = new double[n][m][m * 2];
                double[][][] ftAhatnew = new double[n][m][m * 2];

                for (int cnt = 0; cnt < n; cnt++) {
                    double[][] Im = complexEye(m);
                    double[][] ImMinusOdd = complexMatSub(Im, ft1[cnt]);
                    double[][] invImMinusOdd = complexMatInv(ImMinusOdd, m);

                    double[][] prod1 = complexMatMul(ft2[cnt], invImMinusOdd, m);
                    double[][] prod2 = complexMatMul(prod1, ft3[cnt], m);
                    ftAhatnew[cnt] = complexMatAdd(ft4[cnt], prod2);

                    double[][] prod3 = complexMatMul(prod1, ft2[cnt], m);
                    double phase = omega * cnt;
                    double expReal = Math.cos(phase);
                    double expImag = Math.sin(phase);
                    double[][] scaledOdd = complexMatScale(ft1[cnt], expReal, expImag, m);
                    ftAnew[cnt] = complexMatAdd(scaledOdd, prod3);
                }

                Anew = complexIdftBlocks(ftAnew, m, n);
                Ahatnew = complexIdftBlocks(ftAhatnew, m, n);
            } else {
                Matrix ImOdd = Matrix.eye(m).sub(Aodd.extractRows(0, m));
                Matrix temp = Aeven.extractRows(0, m).mult(ImOdd.inv());
                Ahatnew = Ahateven.extractRows(0, m).add(temp.mult(Ahatodd.extractRows(0, m)));

                Matrix block1 = temp.mult(Aeven.extractRows(0, m));
                Anew = new Matrix(2 * m, m);
                for (int r = 0; r < m; r++) {
                    for (int c = 0; c < m; c++) {
                        Anew.set(r, c, block1.get(r, c));
                        Anew.set(m + r, c, Aodd.get(r, c));
                    }
                }
                Matrix temp2 = Ahatnew.copy();
                Ahatnew = new Matrix(m, m);
                for (int r = 0; r < m; r++) {
                    for (int c = 0; c < m; c++) {
                        Ahatnew.set(r, c, temp2.get(r, c));
                    }
                }
            }

            double nAnew = computeUpperHalfNorm(Anew, m);
            double nAhatnew = computeUpperHalfNorm(Ahatnew, m);
            int njCur = (nj > 0) ? nj : 0;

            while ((nAnew > (njCur + 1) * options.getEpsilonValue() ||
                    nAhatnew > (njCur + 1) * options.getEpsilonValue()) &&
                    njCur + 1 < options.getMaxNumRoot()) {
                njCur = 2 * (njCur + 1) - 1;
                int stopv = Math.min(njCur + 1, Aodd.getNumRows() / m);
                int n = njCur + 1;

                double omega2 = -2.0 * Math.PI / n;

                double[][][] ft1 = complexDftBlocks(Aodd, m, n, stopv);
                double[][][] ft2 = complexDftBlocks(Aeven, m, n, stopv);
                double[][][] ft3 = complexDftBlocks(Ahatodd, m, n, stopv);
                double[][][] ft4 = complexDftBlocks(Ahateven, m, n, stopv);

                double[][][] ftAnew2 = new double[n][m][m * 2];
                double[][][] ftAhatnew2 = new double[n][m][m * 2];

                for (int cnt = 0; cnt < n; cnt++) {
                    double[][] Im = complexEye(m);
                    double[][] ImMinusOdd = complexMatSub(Im, ft1[cnt]);
                    double[][] invImMinusOdd = complexMatInv(ImMinusOdd, m);

                    double[][] prod1 = complexMatMul(ft2[cnt], invImMinusOdd, m);
                    double[][] prod2 = complexMatMul(prod1, ft3[cnt], m);
                    ftAhatnew2[cnt] = complexMatAdd(ft4[cnt], prod2);

                    double[][] prod3 = complexMatMul(prod1, ft2[cnt], m);
                    double phase = omega2 * cnt;
                    double expReal = Math.cos(phase);
                    double expImag = Math.sin(phase);
                    double[][] scaledOdd = complexMatScale(ft1[cnt], expReal, expImag, m);
                    ftAnew2[cnt] = complexMatAdd(scaledOdd, prod3);
                }

                Anew = complexIdftBlocks(ftAnew2, m, n);
                Ahatnew = complexIdftBlocks(ftAhatnew2, m, n);

                nAnew = computeUpperHalfNorm(Anew, m);
                nAhatnew = computeUpperHalfNorm(Ahatnew, m);
            }

            if ((nAnew > (njCur + 1) * options.getEpsilonValue() ||
                    nAhatnew > (njCur + 1) * options.getEpsilonValue()) &&
                    njCur + 1 >= options.getMaxNumRoot()) {
                jline.io.InputOutput.line_warning(
                        "MG1_CR",
                        "Maximum number of '%d' reached, accuracy might be affected",
                        options.getMaxNumRoot());
            }

            if (njCur > 1) {
                int halfBlocks = (njCur + 1) / 2;
                Anew = Anew.extractRows(0, halfBlocks * m);
                Ahatnew = Ahatnew.extractRows(0, halfBlocks * m);
            }

            int numBlocksAnew = Anew.getNumRows() / m;
            Aeven = extractEvenBlocks(Anew, m, numBlocksAnew);
            Aodd = extractOddBlocks(Anew, m, numBlocksAnew);

            int numBlocksAhatnew = Ahatnew.getNumRows() / m;
            Ahateven = extractEvenBlocks(Ahatnew, m, numBlocksAhatnew);
            Ahatodd = extractOddBlocks(Ahatnew, m, numBlocksAhatnew);

            if (options.getVerbose() == 1) {
                String modeStr = "PWCR".equals(options.getMode()) ? "Point-wise" : "Shifted PWCR";
                System.out.println("The " + modeStr + " evaluation of Iteration " + numit + " required " + (njCur + 1) + " roots");
            }

            // Test stop criteria
            if ("PWCR".equals(options.getMode()) || "DCR".equals(options.getMode())) {
                Matrix Rnewj = Anew.extractRows(m, 2 * m);
                for (int i = 2; i < Anew.getNumRows() / m; i++) {
                    Rnewj = Rnewj.add(Anew.extractRows(i * m, (i + 1) * m));
                }
                Rnewj = Matrix.eye(m).sub(Rnewj).inv().mult(Anew.extractRows(0, m));

                if (Rj.sub(Rnewj).infinityNorm() < options.getEpsilonValue() ||
                        sumColsMaxElement(Matrix.eye(m).sub(Anew.extractRows(0, m).mult(Matrix.eye(m).sub(Anew.extractRows(m, 2 * m)).inv()))) < options.getEpsilonValue()) {
                    G = Ahatnew.extractRows(0, m);
                    for (int i = 1; i < Ahatnew.getNumRows() / m; i++) {
                        G = G.add(Rnewj.mult(Ahatnew.extractRows(i * m, (i + 1) * m)));
                    }
                    G = D.extractRows(0, m).mult(Matrix.eye(m).sub(G).inv());
                    break;
                }
                Rj = Rnewj;

                if (Anew.extractRows(0, m).infinityNorm() < options.getEpsilonValue() ||
                        elementAbsSum(Ahatnew.extractRows(m, Ahatnew.getNumRows())) < options.getEpsilonValue() ||
                        sumColsMaxElement(Matrix.eye(m).sub(D.extractRows(0, m).mult(Matrix.eye(m).sub(Ahatnew.extractRows(0, m)).inv()))) < options.getEpsilonValue()) {
                    G = D.extractRows(0, m).mult(Matrix.eye(m).sub(Ahatnew.extractRows(0, m)).inv());
                    break;
                }
            } else {
                // ShiftPWCR mode
                Matrix Gold = G.copy();
                G = D.extractRows(0, m).mult(Matrix.eye(m).sub(Ahatnew.extractRows(0, m)).inv());
                if (G.sub(Gold).infinityNorm() < options.getEpsilonValue() ||
                        (Ahatnew.getNumRows() > m && Ahatnew.extractRows(m, Ahatnew.getNumRows()).infinityNorm() < options.getEpsilonValue())) {
                    break;
                }
            }
        }

        if (numit == options.getMaxNumIt() && G.elementSum() == 0.0) {
            jline.io.InputOutput.line_warning("MG1_CR", "Maximum Number of Iterations %d reached", numit);
            Matrix arg;
            if (Ahateven.getNumRows() >= m) {
                arg = Ahateven.extractRows(0, m);
            } else {
                arg = Matrix.zeros(m, m);
            }
            G = D.extractRows(0, m).mult(Matrix.eye(m).sub(arg).inv());
        }

        // Transpose G
        G = G.transpose();

        // Apply shift correction
        if (options.getMode().contains("ShiftPWCR")) {
            String st = options.getShiftType();
            if ("one".equals(st)) {
                if (drift < 1) G = G.add(Matrix.ones(m, m).scale(1.0 / m));
            } else if ("tau".equals(st)) {
                if (drift > 1) G = G.add(v.mult(Matrix.ones(1, m)).scale(tau));
            } else if ("dbl".equals(st)) {
                if (drift < 1) G = G.add(Matrix.ones(m, m).scale(1.0 / m));
                if (drift > 1) G = G.add(v.mult(Matrix.ones(1, m)).scale(tau));
            }
        }

        if (options.getVerbose() == 1) {
            System.out.println("Final Residual Error for G: CR computation completed in " + numit + " iterations");
        }

        return G;
    }

    // Helper functions

    private static int floorLog2(int n) {
        if (n <= 1) return 0;
        int result = 0;
        int v = n;
        while (v > 1) {
            v = v >> 1;
            result++;
        }
        return result;
    }

    private static Matrix extractEvenBlocks(Matrix D, int m, int numBlocks) {
        int evenCount = (numBlocks + 1) / 2;
        Matrix result = new Matrix(evenCount * m, D.getNumCols());
        int idx = 0;
        for (int i = 0; i < numBlocks; i += 2) {
            for (int r = 0; r < m; r++) {
                for (int c = 0; c < D.getNumCols(); c++) {
                    result.set(idx * m + r, c, D.get(i * m + r, c));
                }
            }
            idx++;
        }
        return result;
    }

    private static Matrix extractOddBlocks(Matrix D, int m, int numBlocks) {
        int oddCount = numBlocks / 2;
        Matrix result = new Matrix(oddCount * m, D.getNumCols());
        int idx = 0;
        for (int i = 1; i < numBlocks; i += 2) {
            for (int r = 0; r < m; r++) {
                for (int c = 0; c < D.getNumCols(); c++) {
                    result.set(idx * m + r, c, D.get(i * m + r, c));
                }
            }
            idx++;
        }
        return result;
    }

    private static double[][][] complexDftBlocks(Matrix blocks, int m, int N, int maxBlocks) {
        int numBlocks = Math.min(maxBlocks, blocks.getNumRows() / m);
        double[][][] result = new double[N][m][m * 2];

        for (int k = 0; k < N; k++) {
            double angle = -2.0 * Math.PI * k / N;
            for (int n = 0; n < numBlocks; n++) {
                double theta = angle * n;
                double cosTheta = Math.cos(theta);
                double sinTheta = Math.sin(theta);
                for (int r = 0; r < m; r++) {
                    for (int c = 0; c < m; c++) {
                        double v = blocks.get(n * m + r, c);
                        result[k][r][c * 2] += v * cosTheta;
                        result[k][r][c * 2 + 1] += v * sinTheta;
                    }
                }
            }
        }
        return result;
    }

    private static Matrix complexIdftBlocks(double[][][] ft, int m, int N) {
        Matrix result = new Matrix(N * m, m);
        double invN = 1.0 / N;

        for (int n = 0; n < N; n++) {
            for (int k = 0; k < N; k++) {
                double angle = 2.0 * Math.PI * k * n / N;
                double cosAngle = Math.cos(angle);
                double sinAngle = Math.sin(angle);
                for (int r = 0; r < m; r++) {
                    for (int c = 0; c < m; c++) {
                        double re = ft[k][r][c * 2];
                        double im = ft[k][r][c * 2 + 1];
                        double cur = result.get(n * m + r, c);
                        result.set(n * m + r, c, cur + (re * cosAngle - im * sinAngle) * invN);
                    }
                }
            }
        }
        return result;
    }

    private static double[][] complexEye(int m) {
        double[][] result = new double[m][m * 2];
        for (int i = 0; i < m; i++) {
            result[i][i * 2] = 1.0;
        }
        return result;
    }

    private static double[][] complexMatSub(double[][] A, double[][] B) {
        int m = A.length;
        double[][] result = new double[m][];
        for (int r = 0; r < m; r++) {
            result[r] = new double[A[r].length];
            for (int c = 0; c < A[r].length; c++) {
                result[r][c] = A[r][c] - B[r][c];
            }
        }
        return result;
    }

    private static double[][] complexMatAdd(double[][] A, double[][] B) {
        int m = A.length;
        double[][] result = new double[m][];
        for (int r = 0; r < m; r++) {
            result[r] = new double[A[r].length];
            for (int c = 0; c < A[r].length; c++) {
                result[r][c] = A[r][c] + B[r][c];
            }
        }
        return result;
    }

    private static double[][] complexMatMul(double[][] A, double[][] B, int m) {
        double[][] result = new double[m][m * 2];
        for (int i = 0; i < m; i++) {
            for (int j = 0; j < m; j++) {
                double re = 0.0;
                double im = 0.0;
                for (int k = 0; k < m; k++) {
                    double aRe = A[i][k * 2];
                    double aIm = A[i][k * 2 + 1];
                    double bRe = B[k][j * 2];
                    double bIm = B[k][j * 2 + 1];
                    re += aRe * bRe - aIm * bIm;
                    im += aRe * bIm + aIm * bRe;
                }
                result[i][j * 2] = re;
                result[i][j * 2 + 1] = im;
            }
        }
        return result;
    }

    private static double[][] complexMatScale(double[][] A, double scaleRe, double scaleIm, int m) {
        double[][] result = new double[m][m * 2];
        for (int i = 0; i < m; i++) {
            for (int j = 0; j < m; j++) {
                double re = A[i][j * 2];
                double im = A[i][j * 2 + 1];
                result[i][j * 2] = re * scaleRe - im * scaleIm;
                result[i][j * 2 + 1] = re * scaleIm + im * scaleRe;
            }
        }
        return result;
    }

    private static double[][] complexMatInv(double[][] A, int m) {
        double[][] aug = new double[m][m * 4];
        for (int i = 0; i < m; i++) {
            for (int j = 0; j < m; j++) {
                aug[i][j * 2] = A[i][j * 2];
                aug[i][j * 2 + 1] = A[i][j * 2 + 1];
            }
            aug[i][(m + i) * 2] = 1.0;
        }

        for (int col = 0; col < m; col++) {
            double maxMag = 0.0;
            int pivotRow = col;
            for (int row = col; row < m; row++) {
                double re = aug[row][col * 2];
                double im = aug[row][col * 2 + 1];
                double mag = re * re + im * im;
                if (mag > maxMag) {
                    maxMag = mag;
                    pivotRow = row;
                }
            }

            if (pivotRow != col) {
                double[] temp = aug[col];
                aug[col] = aug[pivotRow];
                aug[pivotRow] = temp;
            }

            double pivRe = aug[col][col * 2];
            double pivIm = aug[col][col * 2 + 1];
            double pivMag2 = pivRe * pivRe + pivIm * pivIm;
            if (pivMag2 < 1e-30) continue;

            double invRe = pivRe / pivMag2;
            double invIm = -pivIm / pivMag2;
            for (int j = 0; j < 2 * m; j++) {
                double re = aug[col][j * 2];
                double im = aug[col][j * 2 + 1];
                aug[col][j * 2] = re * invRe - im * invIm;
                aug[col][j * 2 + 1] = re * invIm + im * invRe;
            }

            for (int row = 0; row < m; row++) {
                if (row == col) continue;
                double factRe = aug[row][col * 2];
                double factIm = aug[row][col * 2 + 1];
                for (int j = 0; j < 2 * m; j++) {
                    aug[row][j * 2] -= factRe * aug[col][j * 2] - factIm * aug[col][j * 2 + 1];
                    aug[row][j * 2 + 1] -= factRe * aug[col][j * 2 + 1] + factIm * aug[col][j * 2];
                }
            }
        }

        double[][] result = new double[m][m * 2];
        for (int i = 0; i < m; i++) {
            for (int j = 0; j < m; j++) {
                result[i][j * 2] = aug[i][(m + j) * 2];
                result[i][j * 2 + 1] = aug[i][(m + j) * 2 + 1];
            }
        }
        return result;
    }

    private static double computeUpperHalfNorm(Matrix A, int m) {
        int deg = A.getNumRows() / m;
        double maxNorm = 0.0;
        for (int i = deg / 2; i < deg; i++) {
            Matrix block = A.extractRows(i * m, (i + 1) * m);
            double norm = block.infinityNorm();
            if (norm > maxNorm) maxNorm = norm;
        }
        return maxNorm;
    }

    private static double sumColsMaxElement(Matrix m) {
        // sum cols and return max
        double max = Double.NEGATIVE_INFINITY;
        for (int j = 0; j < m.getNumCols(); j++) {
            double s = 0.0;
            for (int i = 0; i < m.getNumRows(); i++) {
                s += m.get(i, j);
            }
            if (s > max) max = s;
        }
        return max;
    }

    private static double elementAbsSum(Matrix m) {
        double sum = 0.0;
        for (int i = 0; i < m.getNumRows(); i++) {
            for (int j = 0; j < m.getNumCols(); j++) {
                sum += Math.abs(m.get(i, j));
            }
        }
        return sum;
    }
}
