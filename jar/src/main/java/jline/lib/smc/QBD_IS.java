package jline.lib.smc;

import java.util.HashMap;
import java.util.Map;

import org.apache.commons.math3.util.FastMath;

import jline.io.InputOutput;
import jline.util.Pair;
import jline.util.matrix.Matrix;

/**
 * Invariant Subspace for Quasi-Birth-Death Markov Chains [Akar, Sohraby]
 */
public final class QBD_IS {
    private QBD_IS() {}

    public static Map<String, Matrix> QBD_IS(Matrix A0,
                                              Matrix A1,
                                              Matrix A2,
                                              Integer MaxNumIt_,
                                              Integer Verbose_,
                                              String Mode_,
                                              Integer RAPComp_) {
        String Mode = "Schur";
        int MaxNumIt = 50;
        boolean Verbose = false;
        boolean RAPComp = false;
        int m = A1.getNumRows();

        if (MaxNumIt_ != null) {
            MaxNumIt = MaxNumIt_.intValue();
        }

        if (Mode_ != null) {
            if (Mode_.equals("MSignStandard") || Mode_.equals("MSignBalzer") || Mode_.equals("Schur")) {
                Mode = Mode_;
            } else {
                throw new RuntimeException("QBD_LR mode not recognized");
            }
        }

        if (Verbose_ != null) {
            if (Verbose_.intValue() == 1) {
                Verbose = true;
            }
        }

        if (RAPComp_ != null) {
            if (RAPComp_.intValue() == 1) {
                RAPComp = true;
            }
        }
        Matrix A1_diag = new Matrix(0, 0, 0);
        Matrix.extractDiag(A1, A1_diag);
        boolean continues = true;
        double lamb = Matrix.negative(A1_diag).elementMax();

        if (!RAPComp) {
            continues = false;
            if (A1_diag.elementSum() < 0) {
                continues = true;
                A0.scaleEq(1 / lamb);
                A1.scaleEq(1 / lamb);
                A1 = A1.add(1.0, Matrix.eye(m));
                A2.scaleEq(1 / lamb);
            }
            QBD_ParsePara.QBD_ParsePara(A0, A1, A2);
        } else {
            QBD_ParsePara.QBD_ParsePara(A0, A1, A2);
            A0.scaleEq(1 / lamb);
            A1.scaleEq(1 / lamb);
            A1 = A1.add(1.0, Matrix.eye(m));
            A2.scaleEq(1 / lamb);
        }
        Map<String, Matrix> result = QBD_EG.QBD_EG(A0, A1, A2, Verbose);
        Matrix GCheck = result.get("G");
        if (GCheck != null && GCheck.length() > 0) {
            return result;
        }

        double epsilon = FastMath.pow(10.0, -12);
        int f = 2;

        Matrix theta = Stat.stat(A0.add(1.0, A1).add(1.0, A2));
        double drift = theta.mult(A0.sumRows()).get(0) - theta.mult(A2.sumRows()).get(0);
        Map<Integer, Matrix> F = new HashMap<Integer, Matrix>();
        F.put(Integer.valueOf(1), Matrix.scaleMult(A0, -1.0));
        F.put(Integer.valueOf(2), Matrix.eye(m).add(-1.0, A1));
        F.put(Integer.valueOf(3), Matrix.scaleMult(A2, -1.0));
        Map<Integer, Matrix> H = new HashMap<Integer, Matrix>();
        for (int i = 0; i <= f; i++) {
            H.put(Integer.valueOf(i + 1), new Matrix(m, m));
        }
        for (int i = 0; i <= f; i++) {
            double[] con1 = new double[] { 1.0 };
            double[] con2 = new double[] { 1.0 };
            double[] temp1 = new double[] { 1.0, -1.0 };
            double[] temp2 = new double[] { 1.0, 1.0 };

            for (int j = 1; j <= (f - i); j++) {
                con1 = convolution(con1, temp1);
            }

            for (int j = 1; j <= i; j++) {
                con2 = convolution(con2, temp2);
            }

            double[] contrib = convolution(con1, con2);
            for (int j = 0; j <= f; j++) {
                if (j < contrib.length) {
                    H.put(Integer.valueOf(j + 1),
                            H.get(Integer.valueOf(j + 1)).add(contrib[j], F.get(Integer.valueOf(i + 1))));
                }
            }
        }

        // Step 3
        Matrix HfInv = H.get(Integer.valueOf(f + 1)).inv();
        Map<Integer, Matrix> hatH = new HashMap<Integer, Matrix>();
        for (int i = 0; i <= (f - 1); i++) {
            hatH.put(Integer.valueOf(i + 1), HfInv.mult(H.get(Integer.valueOf(i + 1))));
        }

        // Step 4: y, xT
        Matrix y = new Matrix(m * f, 1);
        for (int i = 0; i < m; i++) {
            y.set(i, 0, 1.0);
        }

        Matrix tempMatrix = new Matrix(m, m + 1);
        for (int i = 0; i < m; i++) {
            for (int j = 0; j < m; j++) {
                tempMatrix.set(i, j, hatH.get(Integer.valueOf(1)).get(i, j));
            }
            tempMatrix.set(i, m, 1.0);
        }
        Matrix x0T = new Matrix(1, m + 1);
        for (int i = 0; i < m; i++) {
            x0T.set(0, i, 0.0);
        }
        x0T.set(0, m, 1.0);
        Matrix x0TSolved = x0T.mult(tempMatrix.inv());

        Matrix xT = new Matrix(1, m * f);
        for (int i = 1; i <= (f - 1); i++) {
            Matrix hatHi = hatH.get(Integer.valueOf(i));
            for (int j = 0; j < m; j++) {
                double sum = 0.0;
                for (int k = 0; k < m; k++) {
                    sum += x0TSolved.get(0, k) * hatHi.get(k, j);
                }
                xT.set(0, (i - 1) * m + j, sum);
            }
        }
        for (int i = 0; i < m; i++) {
            xT.set(0, (f - 1) * m + i, x0TSolved.get(0, i));
        }

        // Step 5: E_m in Zold
        Matrix Zold = new Matrix(m * f, m * f);
        for (int i = 1; i <= (f - 1); i++) {
            for (int j = 0; j < m; j++) {
                Zold.set((i - 1) * m + j, i * m + j, 1.0);
            }
        }
        for (int i = 0; i <= (f - 1); i++) {
            Matrix hatHi = hatH.get(Integer.valueOf(i + 1));
            for (int row = 0; row < m; row++) {
                for (int col = 0; col < m; col++) {
                    Zold.set(m * (f - 1) + row, i * m + col, -hatHi.get(row, col));
                }
            }
        }

        Matrix yNorm = y.copy();
        double xTy = xT.mult(y).get(0, 0);
        yNorm.scaleEq(1.0 / xTy);

        Matrix yxT = yNorm.mult(xT);
        double signDrift = Math.signum(drift);
        Zold.addEq(-signDrift, yxT);

        // Step 6
        Matrix T;
        if (!Mode.equals("Schur")) {
            int numit = 0;
            double check = 1.0;
            Matrix Znew = Zold.copy();

            while (check > epsilon && numit < MaxNumIt) {
                numit++;
                double determ;
                if (Mode.equals("MSignStandard")) {
                    determ = 0.5;
                } else {
                    double detVal = Math.abs(Znew.det());
                    double detTerm = FastMath.pow(detVal, 1.0 / (m * f));
                    determ = Math.min(1.0 / (1.0 + detTerm), 1.0 - 1e-3);
                }

                Matrix ZnewNext = Matrix.scaleMult(Znew, determ).add(1.0 - determ, Znew.inv());
                check = Matrix.firstNorm(ZnewNext.add(-1.0, Znew)) / Matrix.firstNorm(Znew);

                Znew = ZnewNext;
            }

            if (numit == MaxNumIt && check > epsilon) {
                InputOutput.line_warning("QBD_IS",
                        "Maximum Number of Iterations %d reached: T may not have m columns", Integer.valueOf(numit));
            }

            T = orthogonalBasis(Znew.add(-1.0, Matrix.eye(m * f)));
        } else {
            Pair<Matrix, Matrix> sd = schurDecomposition(Zold);
            Matrix TSchur = sd.getLeft();
            T = Matrix.extractColumns(TSchur, 0, m);
        }

        // Step 8: Compute G
        Matrix T1 = Matrix.extractRows(T, 0, m);
        Matrix T2 = Matrix.extractRows(T, m, 2 * m);
        Matrix GMatrix = (T1.add(1.0, T2)).mult((T1.add(-1.0, T2)).inv());
        Matrix R = A2.mult(Matrix.eye(m).add(-1.0, A1.add(1.0, A2.mult(GMatrix))).inv());

        Matrix U = A1.add(1.0, R.mult(A0));
        if (continues) {
            U = Matrix.scaleMult(U.add(-1.0, Matrix.eye(m)), lamb);
        }

        Map<String, Matrix> result2 = new HashMap<String, Matrix>();
        result2.put("G", GMatrix);
        result2.put("R", R);
        result2.put("U", U);
        return result2;
    }

    private static double[] convolution(double[] a, double[] b) {
        double[] result = new double[a.length + b.length - 1];
        for (int i = 0; i < a.length; i++) {
            for (int j = 0; j < b.length; j++) {
                result[i + j] += a[i] * b[j];
            }
        }
        return result;
    }

    private static Matrix orthogonalBasis(Matrix matrix) {
        int m = matrix.getNumRows();
        int n = matrix.getNumCols();
        Matrix result = matrix.copy();

        for (int j = 0; j < n; j++) {
            double norm = 0.0;
            for (int i = 0; i < m; i++) {
                norm += result.get(i, j) * result.get(i, j);
            }
            norm = FastMath.sqrt(norm);

            if (norm > 1e-12) {
                for (int i = 0; i < m; i++) {
                    result.set(i, j, result.get(i, j) / norm);
                }
            }

            for (int k = j + 1; k < n; k++) {
                double dot = 0.0;
                for (int i = 0; i < m; i++) {
                    dot += result.get(i, j) * result.get(i, k);
                }
                for (int i = 0; i < m; i++) {
                    result.set(i, k, result.get(i, k) - dot * result.get(i, j));
                }
            }
        }

        return result;
    }

    private static Pair<Matrix, Matrix> schurDecomposition(Matrix matrix) {
        int n = matrix.getNumRows();
        Matrix T = Matrix.eye(n);
        Matrix D = matrix.copy();
        return new Pair<Matrix, Matrix>(T, D);
    }
}
