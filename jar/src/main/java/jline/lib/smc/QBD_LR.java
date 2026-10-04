package jline.lib.smc;

import java.util.HashMap;
import java.util.Map;

import org.apache.commons.math3.util.FastMath;

import jline.util.matrix.Matrix;

public final class QBD_LR {
    private QBD_LR() {}

    public static Map<String, Matrix> QBD_LR(Matrix A0arg, Matrix A1arg, Matrix A2arg,
                                             Integer MaxNumIt_, Integer Verbose_,
                                             String Mode_, Integer RAPComp_) {
        Matrix A0 = A0arg;
        Matrix A1 = A1arg;
        Matrix A2 = A2arg;
        String Mode = "Shift";
        int MaxNumIt = 50;
        boolean Verbose = false;
        boolean RAPComp = false;

        if (MaxNumIt_ != null) {
            MaxNumIt = MaxNumIt_;
        }
        if (Mode_ != null) {
            if (Mode_.equals("Shift") || Mode_.equals("Basic")) {
                Mode = Mode_;
            } else {
                throw new RuntimeException("QBD_LR mode not recognized");
            }
        }
        if (Verbose_ != null && Verbose_ == 1) {
            Verbose = true;
        }
        if (RAPComp_ != null && RAPComp_ == 1) {
            RAPComp = true;
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
        Map<String, Matrix> result = QBD_EG.QBD_EG(A0, A1, A2, Verbose);
        Matrix G = result.get("G");
        if (G != null && G.length() > 0) {
            return result;
        }
        Matrix theta = Stat.stat(A0.add(1.0, A1).add(1.0, A2));
        double drift = theta.mult(A0.sumRows()).get(0) - theta.mult(A2.sumRows()).get(0);
        Matrix A2old = A2.copy();
        Matrix uT = Matrix.scaleMult(Matrix.ones(1, m), 1.0 / (double) m);
        Matrix A0old = A0.copy();
        if (Mode.equals("Shift")) {
            if (drift < 0) {
                A2 = A2.add(-1.0, Matrix.ones(m, 1).mult(theta.mult(A2)));
                A1 = A1.add(1.0, Matrix.ones(m, 1).mult(theta.mult(A0)));
            } else {
                A0 = A0.add(-1.0, A0.sumRows().mult(uT));
                A1 = A1.add(1.0, A2.sumRows().mult(uT));
            }
        }

        Matrix B2 = Matrix.eye(m).add(-1.0, A1).inv();
        Matrix B0 = B2.mult(A2);
        B2 = B2.mult(A0);
        G = B2.copy();
        Matrix PI = B0.copy();
        double check = 1.0;
        int numit = 0;
        while (check > FastMath.pow(10.0, -14.0) && numit < MaxNumIt) {
            Matrix A1star = B2.mult(B0).add(1.0, B0.mult(B2));
            Matrix A0star = B0.mult(B0);
            Matrix A2star = B2.mult(B2);
            B0 = Matrix.eye(m).add(-1.0, A1star).inv();
            B2 = B0.mult(A2star);
            B0 = B0.mult(A0star);
            G = G.add(1.0, PI.mult(B2));
            PI = PI.mult(B0);
            check = FastMath.min(Matrix.infNorm(B0), Matrix.infNorm(B2));
            numit = numit + 1;
        }

        if (Verbose && numit == MaxNumIt && check > FastMath.pow(10.0, -14.0)) {
            System.out.println("Maximum Number of Iterations reached");
        }

        if (Mode.equals("Shift")) {
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

        Map<String, Matrix> retMap = new HashMap<String, Matrix>();
        retMap.put("G", G);
        retMap.put("R", R);
        retMap.put("U", U);
        return retMap;
    }
}
