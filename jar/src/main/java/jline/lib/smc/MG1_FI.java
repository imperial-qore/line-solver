/**
 * @file M/G/1-type Functional Iteration solver
 *
 * @since LINE 3.1.0
 */
package jline.lib.smc;

import jline.io.InputOutput;
import jline.util.matrix.Matrix;

public final class MG1_FI {
    private MG1_FI() {}

    /**
     * Functional Iterations for M/G/1-Type Markov Chains.
     */
    public static Matrix mg1_fi(Matrix A, MG1FIOptions options) {
        int m = A.getNumRows();
        int maxd = A.getNumCols() / m - 1;

        Matrix explicitG = MG1_EG.mg1_eg(A, options.getVerbose() > 0);
        if (explicitG != null) {
            return explicitG;
        }

        int numit = 0;
        double check = 1.0;
        Matrix G = (options.getStartValue() != null) ? options.getStartValue() : Matrix.zeros(m, m);

        Matrix workA = A.copy();
        double drift = 0.0;
        double tau = 0.0;
        Matrix v = null;

        boolean useShift = options.getMode().startsWith("Shift");
        if (useShift && options.getNonZeroBlocks() == null) {
            Quadruple<Matrix, Double, Double, Matrix> shiftResult = MG1_Shifts.mg1_shifts(workA, options.getShiftType());
            workA = shiftResult.getFirst();
            drift = shiftResult.getSecond();
            tau = shiftResult.getThird();
            v = shiftResult.getFourth();
        }

        String actualMode = useShift ? options.getMode().substring("Shift".length()) : options.getMode();

        if (options.getNonZeroBlocks() == null) {
            if (actualMode.contains("Natural")) {
                while (check > 1e-14 && numit < options.getMaxNumIt()) {
                    Matrix Gold = G.copy();
                    G = workA.extractCols(maxd * m, (maxd + 1) * m).copy();
                    for (int j = maxd - 1; j >= 0; j--) {
                        G = workA.extractCols(j * m, (j + 1) * m).add(G.mult(Gold));
                    }
                    check = G.sub(Gold).infinityNorm();
                    numit++;
                }
            } else if (actualMode.contains("Traditional")) {
                while (check > 1e-14 && numit < options.getMaxNumIt()) {
                    Matrix Gold = G.copy();
                    G = workA.extractCols(maxd * m, (maxd + 1) * m).copy();
                    for (int j = maxd - 1; j >= 2; j--) {
                        G = workA.extractCols(j * m, (j + 1) * m).add(G.mult(Gold));
                    }
                    G = workA.extractCols(0, m).add(G.mult(Gold).mult(Gold));
                    Matrix ImA1 = Matrix.eye(m).sub(workA.extractCols(m, 2 * m));
                    G = ImA1.inv().mult(G);
                    check = G.sub(Gold).infinityNorm();
                    numit++;
                }
            } else if (actualMode.contains("U-Based")) {
                while (check > 1e-14 && numit < options.getMaxNumIt()) {
                    Matrix Gold = G.copy();
                    G = workA.extractCols(maxd * m, (maxd + 1) * m).copy();
                    for (int j = maxd - 1; j >= 1; j--) {
                        G = workA.extractCols(j * m, (j + 1) * m).add(G.mult(Gold));
                    }
                    Matrix ImG = Matrix.eye(m).sub(G);
                    G = ImG.inv().mult(workA.extractCols(0, m));
                    check = G.sub(Gold).infinityNorm();
                    numit++;
                }
            }

            if (useShift && v != null) {
                if ("one".equals(options.getShiftType())) {
                    if (drift < 1) G = G.add(Matrix.ones(m, m).scale(1.0 / m));
                } else if ("tau".equals(options.getShiftType())) {
                    if (drift > 1) G = G.add(v.mult(Matrix.ones(1, m)).scale(tau));
                } else if ("dbl".equals(options.getShiftType())) {
                    if (drift < 1) G = G.add(Matrix.ones(m, m).scale(1.0 / m));
                    if (drift > 1) G = G.add(v.mult(Matrix.ones(1, m)).scale(tau));
                }
            }
        } else {
            int[] nzb = options.getNonZeroBlocks();
            int[] vec = new int[nzb.length + 1];
            vec[0] = 0;
            for (int i = 0; i < nzb.length; i++) vec[i + 1] = nzb[i];
            int[] vecDiff = new int[vec.length - 1];
            for (int i = 0; i < vecDiff.length; i++) vecDiff[i] = vec[i + 1] - vec[i];

            if (vecDiff[0] > 1 && "Traditional".equals(actualMode)) {
                MG1FIOptions newOpts = new MG1FIOptions("Natural", options.getMaxNumIt(),
                        options.getVerbose(), options.getShiftType(),
                        options.getStartValue(), options.getNonZeroBlocks());
                return mg1_fi(A, newOpts);
            }

            if (actualMode.contains("Natural")) {
                while (check > 1e-14 && numit < options.getMaxNumIt()) {
                    Matrix Gold = G.copy();
                    G = workA.extractCols(maxd * m, (maxd + 1) * m).copy();
                    for (int j = maxd - 1; j >= 0; j--) {
                        G = workA.extractCols(j * m, (j + 1) * m).add(G.mult(matrixPow(Gold, vecDiff[j + 1])));
                    }
                    check = G.sub(Gold).infinityNorm();
                    numit++;
                }
            } else if (actualMode.contains("Traditional")) {
                while (check > 1e-14 && numit < options.getMaxNumIt()) {
                    Matrix Gold = G.copy();
                    G = workA.extractCols(maxd * m, (maxd + 1) * m).copy();
                    for (int j = maxd - 1; j >= 2; j--) {
                        G = workA.extractCols(j * m, (j + 1) * m).add(G.mult(matrixPow(Gold, vecDiff[j + 1])));
                    }
                    G = workA.extractCols(0, m).add(G.mult(matrixPow(Gold, 1 + vecDiff[2])));
                    Matrix ImA1 = Matrix.eye(m).sub(workA.extractCols(m, 2 * m));
                    G = ImA1.inv().mult(G);
                    check = G.sub(Gold).infinityNorm();
                    numit++;
                }
            } else if (actualMode.contains("U-Based")) {
                while (check > 1e-14 && numit < options.getMaxNumIt()) {
                    Matrix Gold = G.copy();
                    G = workA.extractCols(maxd * m, (maxd + 1) * m).copy();
                    for (int j = maxd - 1; j >= 1; j--) {
                        G = workA.extractCols(j * m, (j + 1) * m).add(G.mult(matrixPow(Gold, vecDiff[j + 1])));
                    }
                    Matrix ImGpow = Matrix.eye(m).sub(G.mult(matrixPow(Gold, vecDiff[0] - 1)));
                    G = ImGpow.inv().mult(workA.extractCols(0, m));
                    check = G.sub(Gold).infinityNorm();
                    numit++;
                }
            }
        }

        if (numit == options.getMaxNumIt()) {
            InputOutput.line_warning("MG1_FI", "Maximum number of iterations %d reached", numit);
        }

        return G;
    }

    public static Matrix mg1_fi(Matrix A) {
        return mg1_fi(A, new MG1FIOptions());
    }

    private static Matrix matrixPow(Matrix m, int n) {
        if (n == 0) return Matrix.eye(m.getNumRows());
        if (n == 1) return m.copy();

        Matrix result = Matrix.eye(m.getNumRows());
        Matrix base = m.copy();
        int exp = n;

        while (exp > 0) {
            if (exp % 2 == 1) {
                result = result.mult(base);
            }
            base = base.mult(base);
            exp /= 2;
        }
        return result;
    }
}
