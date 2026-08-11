/*
 * Copyright (c) 2012-2026, QORE Lab, Imperial College London
 * All rights reserved.
 */
package jline.lib.smc;

import jline.io.Ret;
import jline.util.matrix.Matrix;

public final class QBD_CAUDAL {
    private QBD_CAUDAL() {}

    public static double QBD_CAUDAL(Matrix A0, Matrix A1, Matrix A2, boolean Dual) {
        if (Dual) {
            Matrix tmp = A2;
            A2 = A0;
            A0 = tmp;
        }

        double eta_min = 0.0;
        double eta_max = 1.0;
        double eta = 0.5;
        while (eta_max - eta_min > Math.pow(10.0, -15.0)) {
            Matrix A3 = A2.add(1.0, Matrix.scaleMult(A1, eta)).add(1.0, Matrix.scaleMult(A0, Math.pow(eta, 2.0)));
            Ret.Eigs eigs = A3.eigval();
            double new_eta = eigs.values.elementMax();
            if (new_eta > eta) {
                eta_min = eta;
            } else {
                eta_max = eta;
            }
            eta = (eta_max + eta_min) / 2;
        }
        return eta;
    }

    public static double QBD_CAUDAL(Matrix A0, Matrix A1, Matrix A2) {
        return QBD_CAUDAL(A0, A1, A2, false);
    }
}
