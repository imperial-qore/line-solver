/*
 * Copyright (c) 2012-2026, QORE Lab, Imperial College London
 * All rights reserved.
 * Ported from BUTools V2.0
 */
package jline.lib.butools.queues;

import jline.lib.butools.mam.GeneralFluidSolution;
import jline.lib.butools.mam.GeneralFluidSolve;
import jline.lib.butools.mc.CTMCSolve;
import jline.util.matrix.Matrix;

public final class FluFluQueue {
    private FluFluQueue() {}

    public static FluFluResult fluFluQueue(Matrix Qin, Matrix Rin, Matrix Qout, Matrix Rout, boolean srv0stop) {
        return fluFluQueue(Qin, Rin, Qout, Rout, srv0stop, 0, 0, 1e-14);
    }

    public static FluFluResult fluFluQueue(Matrix Qin, Matrix Rin, Matrix Qout, Matrix Rout,
                                           boolean srv0stop, int numFluidMoments, int numSojournMoments) {
        return fluFluQueue(Qin, Rin, Qout, Rout, srv0stop, numFluidMoments, numSojournMoments, 1e-14);
    }

    /**
     * Returns various performance measures of a fluid queue with independent
     * fluid arrival and service processes.
     */
    public static FluFluResult fluFluQueue(Matrix Qin, Matrix Rin, Matrix Qout, Matrix Rout,
                                           boolean srv0stop, int numFluidMoments, int numSojournMoments,
                                           double prec) {
        int Nin = Qin.getNumRows();
        int Nout = Qout.getNumRows();

        Matrix Iin = Matrix.eye(Nin);
        Matrix Iout = Matrix.eye(Nout);

        Matrix piIn = CTMCSolve.ctmcSolve(Qin, prec);
        Matrix piOut = CTMCSolve.ctmcSolve(Qout, prec);
        double lambda = piIn.mult(Rin).elementSum();
        double mu = piOut.mult(Rout).elementSum();

        GeneralFluidSolution fluidSolution = null;
        GeneralFluidSolution sojournSolution = null;
        double[] fluidMoments = null;
        double[] sojournMoments = null;

        if (numFluidMoments > 0) {
            Matrix Q = Qin.kron(Iout).add(1.0, Iin.kron(Qout));
            Matrix R = Rin.kron(Iout).add(-1.0, Iin.kron(Rout));

            Matrix Q0;
            if (srv0stop) {
                Matrix RoutPinv = Rout.pinv();
                Q0 = Qin.kron(Iout).add(1.0, Rin.kron(RoutPinv.mult(Qout)));
            } else {
                Q0 = null;
            }

            fluidSolution = GeneralFluidSolve.generalFluidSolve(Q, R, Q0, prec);

            Matrix ini = fluidSolution.getIni();
            Matrix K = fluidSolution.getK();
            Matrix clo = fluidSolution.getClo();

            Matrix negK = K.neg();
            Matrix invNegK = negK.inv();
            Matrix onesN = Matrix.ones(clo.getNumCols(), 1);

            fluidMoments = new double[numFluidMoments];
            Matrix invKPower = invNegK.copy();

            for (int m = 1; m <= numFluidMoments; m++) {
                invKPower = invKPower.mult(invNegK);
                double moment = factorial(m) * ini.mult(invKPower).mult(clo).mult(onesN).elementSum();
                fluidMoments[m - 1] = moment;
            }
        }

        if (numSojournMoments > 0) {
            Matrix Rh = Rin.kron(Iout).add(-1.0, Iin.kron(Rout));
            Matrix Qh = Qin.kron(Rout).add(1.0, Rin.kron(Qout));

            sojournSolution = GeneralFluidSolve.generalFluidSolve(Qh, Rh, null, prec);

            Matrix inih = sojournSolution.getIni();
            Matrix Kh = sojournSolution.getK();
            Matrix cloh = sojournSolution.getClo();

            Matrix kclo;
            if (srv0stop) {
                kclo = cloh.mult(Rin.kron(Rout)).scale(1.0 / (lambda * mu));
            } else {
                kclo = cloh.mult(Rin.kron(Iout)).scale(1.0 / lambda);
            }

            Matrix negKh = Kh.neg();
            Matrix invNegKh = negKh.inv();

            sojournMoments = new double[numSojournMoments];
            Matrix invKhPower = invNegKh.copy();

            for (int m = 1; m <= numSojournMoments; m++) {
                invKhPower = invKhPower.mult(invNegKh);
                double moment = factorial(m) * inih.mult(invKhPower).mult(kclo).elementSum();
                sojournMoments[m - 1] = moment;
            }
        }

        return new FluFluResult(fluidSolution, sojournSolution, fluidMoments, sojournMoments, lambda, mu);
    }

    private static double factorial(int n) {
        double result = 1.0;
        for (int i = 2; i <= n; i++) {
            result *= i;
        }
        return result;
    }
}
