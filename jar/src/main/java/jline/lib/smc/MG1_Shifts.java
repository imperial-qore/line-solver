/**
 * @file M/G/1-type Shift Technique
 *
 * Implements the shift technique for M/G/1-type Markov chains to improve
 * convergence of iterative solvers.
 *
 * @since LINE 3.1.0
 */
package jline.lib.smc;

import jline.util.matrix.Matrix;

public final class MG1_Shifts {
    private MG1_Shifts() {}

    public static Quadruple<Matrix, Double, Double, Matrix> mg1_shifts(Matrix A) {
        return mg1_shifts(A, "one");
    }

    /**
     * Applies the shift technique to an M/G/1-type block matrix.
     */
    public static Quadruple<Matrix, Double, Double, Matrix> mg1_shifts(Matrix A, String shiftType) {
        int m = A.getNumRows();
        int maxd = A.getNumCols() / m - 1;

        // Compute sumA and beta for drift calculation
        Matrix sumA = A.extractCols(maxd * m, (maxd + 1) * m);
        Matrix beta = sumA.sumRows();
        for (int i = maxd - 1; i >= 1; i--) {
            sumA = sumA.add(A.extractCols(i * m, (i + 1) * m));
            beta = beta.add(sumA.sumRows());
        }
        sumA = sumA.add(A.extractCols(0, m));

        Matrix theta = Stat.stat(sumA);
        double drift = theta.mult(beta).get(0, 0);

        Matrix v = Matrix.zeros(m, 1);
        double tau = 1.0;

        Matrix workA = A.copy();
        Matrix hatA = Matrix.zeros(m, (maxd + 1) * m);

        if (drift < 1) {
            if ("tau".equals(shiftType) || "dbl".equals(shiftType)) {
                MG1DecayResult decayResult = MG1_Decay.mg1_decay(workA, true);
                tau = decayResult.getEta();
                Matrix uT = decayResult.getUT();
                if (uT == null) uT = Matrix.ones(1, m).scale(1.0 / m);

                for (int i = 0; i < m; i++) {
                    workA.set(i, m + i, workA.get(i, m + i) - 1.0);
                }

                double uTSum = uT.elementSum();
                if (uTSum != 0.0) {
                    uT = uT.scale(1.0 / uTSum);
                }

                Matrix rowhatA = new Matrix(1, (maxd + 1) * m);
                for (int c = 0; c < m; c++) {
                    double s = 0.0;
                    for (int k = 0; k < m; k++) {
                        s += uT.get(0, k) * workA.get(k, maxd * m + c);
                    }
                    rowhatA.set(0, maxd * m + c, s);
                }
                for (int i = maxd - 1; i >= 0; i--) {
                    for (int c = 0; c < m; c++) {
                        double s = 0.0;
                        for (int k = 0; k < m; k++) {
                            s += uT.get(0, k) * workA.get(k, i * m + c);
                        }
                        rowhatA.set(0, i * m + c, tau * rowhatA.get(0, (i + 1) * m + c) + s);
                    }
                }

                for (int i = 0; i <= maxd; i++) {
                    for (int r2 = 0; r2 < m; r2++) {
                        for (int c = 0; c < m; c++) {
                            hatA.set(r2, i * m + c, workA.get(r2, i * m + c) - rowhatA.get(0, i * m + c));
                        }
                    }
                }
            }

            if ("dbl".equals(shiftType)) {
                for (int i = 0; i <= maxd; i++) {
                    for (int r2 = 0; r2 < m; r2++) {
                        for (int c = 0; c < m; c++) {
                            workA.set(r2, i * m + c, hatA.get(r2, i * m + c));
                        }
                    }
                }
            }

            if ("one".equals(shiftType)) {
                for (int i = 0; i < m; i++) {
                    workA.set(i, m + i, workA.get(i, m + i) - 1.0);
                }
            }

            if ("one".equals(shiftType) || "dbl".equals(shiftType)) {
                Matrix colhatA = new Matrix(m, maxd + 1);
                for (int r2 = 0; r2 < m; r2++) {
                    double s = 0.0;
                    for (int c = 0; c < m; c++) {
                        s += workA.get(r2, c);
                    }
                    colhatA.set(r2, 0, s);
                }
                for (int i = 1; i <= maxd; i++) {
                    for (int r2 = 0; r2 < m; r2++) {
                        double s = 0.0;
                        for (int c = 0; c < m; c++) {
                            s += workA.get(r2, i * m + c);
                        }
                        colhatA.set(r2, i, colhatA.get(r2, i - 1) + s);
                    }
                }

                for (int i = 0; i <= maxd; i++) {
                    for (int r2 = 0; r2 < m; r2++) {
                        for (int c = 0; c < m; c++) {
                            hatA.set(r2, i * m + c, workA.get(r2, i * m + c) - colhatA.get(r2, i) / m);
                        }
                    }
                }
            }
        } else {
            if ("one".equals(shiftType) || "dbl".equals(shiftType)) {
                for (int i = 0; i < m; i++) {
                    workA.set(i, m + i, workA.get(i, m + i) - 1.0);
                }

                Matrix rowhatA = new Matrix(1, (maxd + 1) * m);
                for (int c = 0; c < m; c++) {
                    double s = 0.0;
                    for (int k = 0; k < m; k++) {
                        s += theta.get(0, k) * workA.get(k, maxd * m + c);
                    }
                    rowhatA.set(0, maxd * m + c, s);
                }
                for (int i = maxd - 1; i >= 0; i--) {
                    for (int c = 0; c < m; c++) {
                        double s = 0.0;
                        for (int k = 0; k < m; k++) {
                            s += theta.get(0, k) * workA.get(k, i * m + c);
                        }
                        rowhatA.set(0, i * m + c, rowhatA.get(0, (i + 1) * m + c) + s);
                    }
                }

                for (int i = 0; i <= maxd; i++) {
                    for (int r2 = 0; r2 < m; r2++) {
                        for (int c = 0; c < m; c++) {
                            hatA.set(r2, i * m + c, workA.get(r2, i * m + c) - rowhatA.get(0, i * m + c));
                        }
                    }
                }
            }

            if ("tau".equals(shiftType) || "dbl".equals(shiftType)) {
                Matrix inputA;
                if ("dbl".equals(shiftType)) {
                    Matrix tempA = hatA.copy();
                    for (int i = 0; i < m; i++) {
                        tempA.set(i, m + i, tempA.get(i, m + i) + 1.0);
                    }
                    inputA = tempA;
                } else {
                    inputA = workA;
                }

                GIM1CaudalResult caudalResult = GIM1_Caudal.gim1_caudal(inputA, false, true);
                tau = caudalResult.getEta();
                v = caudalResult.getV();
                if (v == null) v = Matrix.ones(m, 1).scale(1.0 / m);

                double vSum = v.elementSum();
                if (vSum > 0) {
                    v = v.scale(1.0 / vSum);
                }

                Matrix tauA = inputA.copy();
                for (int i = 0; i < m; i++) {
                    tauA.set(i, m + i, tauA.get(i, m + i) - 1.0);
                }

                Matrix colhatA2 = new Matrix(m, maxd + 1);
                for (int r2 = 0; r2 < m; r2++) {
                    double s = 0.0;
                    for (int c = 0; c < m; c++) {
                        s += tauA.get(r2, c) * v.get(c, 0);
                    }
                    colhatA2.set(r2, 0, s);
                }
                for (int i = 1; i <= maxd; i++) {
                    for (int r2 = 0; r2 < m; r2++) {
                        double s = 0.0;
                        for (int c = 0; c < m; c++) {
                            s += tauA.get(r2, i * m + c) * v.get(c, 0);
                        }
                        colhatA2.set(r2, i, colhatA2.get(r2, i - 1) / tau + s);
                    }
                }

                for (int i = 0; i <= maxd; i++) {
                    for (int r2 = 0; r2 < m; r2++) {
                        for (int c = 0; c < m; c++) {
                            hatA.set(r2, i * m + c, tauA.get(r2, i * m + c) - colhatA2.get(r2, i));
                        }
                    }
                }
            }
        }

        for (int i = 0; i < m; i++) {
            hatA.set(i, m + i, hatA.get(i, m + i) + 1.0);
        }

        return new Quadruple<Matrix, Double, Double, Matrix>(hatA, drift, tau, v);
    }
}
