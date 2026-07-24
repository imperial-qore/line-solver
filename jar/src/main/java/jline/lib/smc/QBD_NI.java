package jline.lib.smc;

import java.util.HashMap;
import java.util.Map;

import jline.util.matrix.Matrix;

public final class QBD_NI {
    private QBD_NI() {}

    public static Map<String, Matrix> QBD_NI(Matrix A0arg, Matrix A1arg, Matrix A2arg,
                                             Integer MaxNumIt_, Integer Verbose_,
                                             String Mode_, Integer RAPComp_) {
        // The uniformization below rescales the blocks in place, so work on
        // copies: MATLAB QBD_NI.m rebinds (A0=A0/lamb) and leaves the caller's
        // matrices untouched, as does QBD_CR here.
        Matrix A0 = A0arg.copy();
        Matrix A1 = A1arg.copy();
        Matrix A2 = A2arg.copy();
        String Mode = "Sylvest";
        int MaxNumIt = 50;
        boolean Verbose = false;
        boolean RAPComp = false;

        if (MaxNumIt_ != null) {
            MaxNumIt = MaxNumIt_;
        }

        if (Mode_ != null) {
            if ("Sylvest".equals(Mode_) || "Estimat".equals(Mode_) || "DirectSum".equals(Mode_)) {
                Mode = Mode_;
            } else {
                throw new RuntimeException("QBD_LR mode not recognized");
            }
        }

        if (Verbose_ != null) {
            if (Verbose_ == 1) {
                Verbose = true;
            }
        }

        if (RAPComp_ != null) {
            if (RAPComp_ == 1) {
                RAPComp = true;
            }
        }
        Matrix A1_diag = new Matrix(0, 0, 0);
        Matrix.extractDiag(A1, A1_diag);
        boolean continues = true;
        double lamb = Matrix.negative(A1_diag).elementMax();
        int m = A1.getNumRows();
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
        // Check whether G is known explicitly; if so QBD_EG has already solved
        // the QBD and its result stands. Only when G is unknown (empty) does the
        // Newton iteration below run. Mirrors MATLAB QBD_NI.m:
        // [G,R,U]=QBD_EG(...); if (~isempty(G)) return; end
        Map<String, Matrix> result = QBD_EG.QBD_EG(A0, A1, A2, Verbose);
        Matrix G = result.get("G");
        if (G.length() != 0) {
            return result;
        }
        Matrix R = new Matrix(m, m, m * m);
        double check = 1.0;
        int numit = 0;
        while (check > Math.pow(10.0, -12.0) && numit < MaxNumIt) {
            numit = numit + 1;
            Matrix YK;
            if (numit == 1) {
                YK = A2.mult(Matrix.eye(m).add(-1.0, A1).inv());
            } else {
                if ("Estimat".equals(Mode)) {
                    Matrix FRK = A2.add(1.0, R.mult(A1.add(-1.0, Matrix.eye(m)).add(1.0, R.mult(A0))));
                    Matrix ZK = FRK.mult(Matrix.eye(m).add(-1.0, A1)).inv();
                    YK = FRK.add(1.0, FRK.add(1.0, ZK.mult(A1)).add(1.0, R.mult(ZK).add(1.0, ZK.mult(R)).mult(A0)));
                } else {
                    Matrix D = Matrix.scaleMult(A2.add(1.0, R.mult(A1.add(-1.0, Matrix.eye(m)).add(1.0, R.mult(A0)))), -1.0);
                    Matrix C = A1.add(1.0, R.mult(A0)).add(-1.0, Matrix.eye(m));
                    if ("Sylvest".equals(Mode)) {
                        // Mirrors MATLAB QBD_NI.m: Yk = QBD_NI_Sylvest(A0',R',C',D')'.
                        // The solver returns the transposed unknown, so the result
                        // must be transposed back.
                        YK = QBD_NI_Sylvest.QBD_NI_Sylvest(A0.transpose(), R.transpose(),
                                C.transpose(), D.transpose()).transpose();
                    } else {
                        Matrix D_reshape = D.copy();
                        D_reshape.reshape(m * m, 1);
                        YK = A0.transpose().krons(R).add(1.0, C.transpose().krons(Matrix.eye(m)).inv().mult(D_reshape));
                        YK.reshape(m, m);
                    }
                }
            }
            R = R.add(1.0, YK);
            check = Matrix.infNorm(YK);
        }
        if (Verbose && numit == MaxNumIt && check > Math.pow(10.0, -12.0)) {
            System.out.println("Maximum Number of Iterations reached");
        }

        // G is derived from the R just computed by the Newton iteration, not the
        // reverse: mirrors MATLAB QBD_NI.m G=(eye(m)-(A1+R*A0))^(-1)*A0. R must
        // not be overwritten here -- it is the iteration's own fixed point.
        G = Matrix.eye(m).add(-1.0, A1.add(1.0, R.mult(A0))).inv().mult(A0);

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
