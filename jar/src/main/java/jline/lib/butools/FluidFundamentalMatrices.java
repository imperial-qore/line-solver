package jline.lib.butools;

import java.util.HashMap;
import java.util.Map;

import org.apache.commons.math3.util.FastMath;

import jline.util.matrix.Matrix;

public final class FluidFundamentalMatrices {
    private FluidFundamentalMatrices() {}

    public static Map<String, Matrix> FluidFundamentalMatrices(
            Matrix Fpp, Matrix Fpm, Matrix Fmp, Matrix Fmm,
            Double precision_, Integer maxNumIt_, String method_) {
        Matrix Psi = Matrix.singleton(0.0);
        double precision = precision_ != null ? precision_ : 1e-14;
        int maxNumIt = maxNumIt_ != null ? maxNumIt_ : 150;
        String method = method_ != null ? method_ : "ADDA";
        int numit = 0;
        if (Fpp.getNumRows() == 0) {
            Psi = new Matrix(0, Fmm.getNumRows());
        } else if ("CR".equals(method)) {
            // unimplemented in source
        } else if ("ADDA".equals(method) || "SDA".equals(method)) {
            Matrix A = Fpp.copy();
            A.scaleEq(-1.0);
            Matrix B = Fpm.copy();
            Matrix C = Fmp.copy();
            Matrix D = Fmm.copy();
            D.scaleEq(-1.0);
            Matrix diag_A = Matrix.singleton(0.0);
            Matrix.extractDiag(A, diag_A);
            double gamma1 = diag_A.elementMax();
            Matrix diag_D = Matrix.singleton(0.0);
            Matrix.extractDiag(D, diag_D);
            double gamma2 = diag_D.elementMax();
            if ("SDA".equals(method)) {
                gamma1 = FastMath.max(gamma1, gamma2);
                gamma2 = gamma1;
            }
            int sA = A.getNumRows();
            int sD = D.getNumRows();
            Matrix IA = Matrix.eye(sA);
            Matrix ID = Matrix.eye(sD);
            Matrix gamma2IA = IA.copy();
            gamma2IA.scaleEq(gamma2);
            Matrix gamma1ID = ID.copy();
            gamma1ID.scaleEq(gamma1);
            A = A.add(1.0, gamma2IA);
            D = D.add(1.0, gamma1ID);
            Matrix Dginv = D.inv();
            Matrix Vginv = D.add(-1.0, C.mult(A.inv()).mult(B)).inv();
            Matrix Wginv = A.add(-1.0, B.mult(Dginv).mult(C)).inv();
            Matrix gammaVginv = Vginv.copy();
            gammaVginv.scaleEq(gamma1 + gamma2);
            Matrix Eg = ID.add(-1.0, gammaVginv);
            Matrix gammaWginV = Wginv.copy();
            gammaWginV.scaleEq(gamma1 + gamma2);
            Matrix Fg = IA.add(-1.0, gammaWginV);
            Matrix Gg = Dginv.mult(C).mult(Wginv);
            Gg.scaleEq(gamma1 + gamma2);
            Matrix Hg = Wginv.mult(B).mult(Dginv);
            Hg.scaleEq(gamma1 + gamma2);

            double diff = 1.0;
            while (diff > precision && numit < maxNumIt) {
                Vginv = Eg.mult(ID.add(-1.0, Gg.mult(Hg)).inv());
                Wginv = Fg.mult(IA.add(-1.0, Hg.mult(Gg)).inv());
                Gg = Gg.add(1.0, Vginv.mult(Gg).mult(Fg));
                Hg = Hg.add(1.0, Wginv.mult(Hg).mult(Eg));
                Eg = Vginv.mult(Eg);
                Fg = Wginv.mult(Fg);
                double neg = Matrix.firstNorm(Eg);
                double nfg = Matrix.firstNorm(Fg);
                if ("ADDA".equals(method)) {
                    double eta = FastMath.sqrt(nfg / neg);
                    Eg.scaleEq(eta);
                    Fg.scaleEq(1 / eta);
                    diff = neg * nfg;
                } else {
                    diff = FastMath.min(neg, nfg);
                }
                numit++;
            }
            Psi = Hg;
        }

        Map<String, Matrix> result = new HashMap<String, Matrix>();
        result.put("P", Psi);
        result.put("K", Fpp.add(1.0, Psi.mult(Fmp)));
        result.put("U", Fmm.add(1.0, Fmp.mult(Psi)));

        return result;
    }
}
