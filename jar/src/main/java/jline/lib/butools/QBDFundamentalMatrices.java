package jline.lib.butools;

import jline.lib.smc.QBD_CR;
import jline.lib.smc.QBD_FI;
import jline.lib.smc.QBD_IS;
import jline.lib.smc.QBD_LR;
import jline.lib.smc.QBD_NI;
import jline.util.matrix.Matrix;

import java.util.Map;

public final class QBDFundamentalMatrices {
    private QBDFundamentalMatrices() {}

    public static Map<String, Matrix> QBDFundamentalMatrices(Matrix B,
                                                             Matrix L,
                                                             Matrix F,
                                                             Double precision_,
                                                             Integer maxNumIt_,
                                                             String method_,
                                                             Integer Verbose_) {
        double precision = 1e-14;
        if (precision_ != null) {
            precision = precision_;
        }

        int maxNumIt = 50;
        if (maxNumIt_ != null) {
            maxNumIt = maxNumIt_;
        }

        String method = "CR";
        if (method_ != null) {
            method = method_;
        }

        int Verbose = 0;
        if (Verbose_ != null) {
            Verbose = Verbose_;
        }

        if ("LR".equals(method)) {
            return QBD_LR.QBD_LR(B, L, F, maxNumIt, Verbose, null, null);
        }

        if ("NI".equals(method)) {
            return QBD_NI.QBD_NI(B, L, F, maxNumIt, Verbose, null, null);
        }

        if ("IS".equals(method)) {
            return QBD_IS.QBD_IS(B, L, F, maxNumIt, Verbose, null, null);
        }

        if ("FI".equals(method)) {
            return QBD_FI.QBD_FI(B, L, F, maxNumIt, Verbose, null, null, null);
        }

        return QBD_CR.QBD_CR(B, L, F, maxNumIt, Verbose, null, null);
    }
}
