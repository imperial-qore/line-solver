package jline.lib.smc;

import java.util.HashMap;
import java.util.Map;

import jline.util.matrix.Matrix;

public final class QBD_FI {
    private QBD_FI() {}

    /**
     * Convenience alias for callers that expect a {@code solve(...)} entry point
     * on the QBD_FI class. Delegates to {@link #QBD_FI}.
     */
    public static Map<String, Matrix> solve(Matrix A0arg, Matrix A1arg, Matrix A2arg,
                                            Integer MaxNumIt_, Integer Verbose_, String Mode_,
                                            Matrix StartValue_, Integer RAPComp_) {
        return QBD_FI(A0arg, A1arg, A2arg, MaxNumIt_, Verbose_, Mode_, StartValue_, RAPComp_);
    }

    public static Map<String, Matrix> QBD_FI(Matrix A0arg, Matrix A1arg, Matrix A2arg,
                                              Integer MaxNumIt_, Integer Verbose_, String Mode_,
                                              Matrix StartValue_, Integer RAPComp_) {
        Matrix A0 = A0arg;
        Matrix A1 = A1arg;
        Matrix A2 = A2arg;
        String Mode = "U-Based";
        int MaxNumIt = 10000;
        boolean Verbose = false;
        boolean RAPComp = false;
        int m = A1.getNumRows();
        Matrix StartValue = new Matrix(m, m, m * m);

        if (StartValue_ != null) {
            StartValue = StartValue_.copy();
        }

        if (MaxNumIt_ != null) {
            MaxNumIt = MaxNumIt_.intValue();
        }

        if (Mode_ != null) {
            if ("Traditional ".equals(Mode_) || "Natural".equals(Mode_) || "U-Based".equals(Mode_)
                    || "ShiftTraditional".equals(Mode_) || "ShiftNatural".equals(Mode_) || "ShiftU-Based".equals(Mode_)) {
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
        Matrix G = result.get("G");
        if (G.length() > 0) {
            return result;
        }
        int numit = 0;
        double check = 1.0;
        G = StartValue.copy();
        Matrix theta = Stat.stat(A0.add(1.0, A1).add(1.0, A2));
        double drift = theta.mult(A0.sumRows()).get(0) - theta.mult(A2.sumRows()).get(0);
        Matrix A2old = A2.copy();
        Matrix uT = Matrix.scaleMult(Matrix.ones(1, m), 1.0 / (double) m);
        Matrix A0old = A0.copy();
        if (Mode.contains("Shift")) {
            if (drift < 0) {
                A2 = A2.add(-1.0, Matrix.ones(m, 1).mult(theta.mult(A2)));
                A1 = A1.add(1.0, Matrix.ones(m, 1).mult(theta.mult(A0)));
            } else {
                A0 = A0.add(-1.0, A0.sumRows().mult(uT));
                A1 = A1.add(1.0, A2.sumRows().mult(uT));
            }
        }

        if (Mode.contains("Natural")) {
            while (check > Math.pow(10.0, -14.0) && numit < MaxNumIt) {
                Matrix Gold = G.copy();
                G = A2.mult(G).add(1.0, A1).mult(G).add(1.0, A0);
                check = Matrix.infNorm(G.add(-1.0, Gold));
                numit = numit + 1;
            }
        }

        if (Mode.contains("Traditional")) {
            Matrix invA1 = Matrix.eye(m).add(-1.0, A1).inv();
            while (check > Math.pow(10.0, -14.0) && numit < MaxNumIt) {
                Matrix Gold = G.copy();
                G = invA1.mult(A0.add(1.0, A2.mult(Matrix.pow(G, 2))));
                check = Matrix.infNorm(G.add(-1.0, Gold));
                numit = numit + 1;
            }
        }

        if (Mode.contains("U-Based")) {
            while (check > Math.pow(10.0, -14.0) && numit < MaxNumIt) {
                Matrix Gold = G.copy();
                G = Matrix.eye(m).add(-1.0, A1).add(-1.0, A2.mult(G)).inv().mult(A0);
                check = Matrix.infNorm(G.add(-1.0, Gold));
                numit = numit + 1;
            }
        }

        if (Verbose && numit == MaxNumIt && check > Math.pow(10.0, -12.0)) {
            System.out.println("Maximum Number of Iterations reached");
        }

        if (Mode.contains("Shift")) {
            if (drift < 0) {
                A1 = A1.add(-1.0, Matrix.ones(m, 1).mult(theta).mult(A0));
                A2 = A2old.copy();
            } else {
                G = G.add(1.0, Matrix.ones(m, 1).mult(uT));
                A1 = A1.add(-1.0, A2.sumRows().mult(uT));
                A0 = A0old.copy();
            }
        }
        Matrix R = A2.mult(Matrix.eye(m).add(-1.0, A1.add(1.0, A2.mult(G))).inv());

        Matrix U = A1.add(1.0, R.mult(A0));
        if (continues) {
            U = Matrix.scaleMult(U.add(-1.0, Matrix.eye(m)), lamb);
        }

        result = new HashMap<String, Matrix>();
        result.put("G", G);
        result.put("R", R);
        result.put("U", U);
        return result;
    }
}
