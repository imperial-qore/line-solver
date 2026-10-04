/*
 * Copyright (c) 2012-2026, QORE Lab, Imperial College London
 * All rights reserved.
 * Ported from BUTools V2.0
 */
package jline.lib.butools.mam;

import java.util.ArrayList;
import java.util.List;
import java.util.Map;

import jline.lib.butools.FluidFundamentalMatrices;
import jline.util.matrix.Matrix;

public final class GeneralFluidSolve {
    private GeneralFluidSolve() {}

    public static GeneralFluidSolution generalFluidSolve(Matrix Q, Matrix R) {
        return generalFluidSolve(Q, R, null, 1e-14);
    }

    public static GeneralFluidSolution generalFluidSolve(Matrix Q, Matrix R, Matrix Q0) {
        return generalFluidSolve(Q, R, Q0, 1e-14);
    }

    /**
     * Returns the parameters of the matrix-exponentially distributed stationary
     * distribution of a general Markovian fluid model.
     */
    public static GeneralFluidSolution generalFluidSolve(Matrix Q, Matrix R, Matrix Q0, double prec) {
        int N = Q.getNumRows();

        List<Integer> ixz = new ArrayList<Integer>();
        List<Integer> ixp = new ArrayList<Integer>();
        List<Integer> ixn = new ArrayList<Integer>();

        for (int i = 0; i < N; i++) {
            double rate = R.get(i, i);
            if (Math.abs(rate) <= prec) ixz.add(i);
            else if (rate > prec) ixp.add(i);
            else if (rate < -prec) ixn.add(i);
        }

        int Nz = ixz.size();
        int Np = ixp.size();
        int Nn = ixn.size();

        Matrix P = Matrix.zeros(N, N);
        for (int i = 0; i < Nz; i++) P.set(i, ixz.get(i), 1.0);
        for (int i = 0; i < Np; i++) P.set(Nz + i, ixp.get(i), 1.0);
        for (int i = 0; i < Nn; i++) P.set(Nz + Np + i, ixn.get(i), 1.0);
        Matrix iP = P.inv();

        Matrix Qv = P.mult(Q).mult(iP);
        Matrix Rv = P.mult(R).mult(iP);

        Matrix Qv00 = (Nz > 0) ? Matrix.getSubMatrix(Qv, 0, Nz, 0, Nz) : Matrix.zeros(0, 0);
        Matrix Qv0pm = (Nz > 0) ? Matrix.getSubMatrix(Qv, 0, Nz, Nz, N) : Matrix.zeros(0, Np + Nn);
        Matrix Qvpm0 = (Nz > 0) ? Matrix.getSubMatrix(Qv, Nz, N, 0, Nz) : Matrix.zeros(Np + Nn, 0);
        Matrix Qvpmpm = Matrix.getSubMatrix(Qv, Nz, N, Nz, N);

        Matrix Qbar;
        if (Nz > 0) {
            Matrix negQv00 = Qv00.neg();
            Matrix iQv00 = negQv00.pinv();
            Qbar = Qvpmpm.add(1.0, Qvpm0.mult(iQv00).mult(Qv0pm));
        } else {
            Qbar = Qvpmpm.copy();
        }

        Matrix absRi = Matrix.zeros(Np + Nn, Np + Nn);
        for (int i = 0; i < Np + Nn; i++) {
            absRi.set(i, i, Math.abs(1.0 / Rv.get(Nz + i, Nz + i)));
        }

        Matrix Qz = absRi.mult(Qbar);

        Matrix Qzpp = (Np > 0) ? Matrix.getSubMatrix(Qz, 0, Np, 0, Np) : Matrix.zeros(0, 0);
        Matrix Qzpn = (Np > 0 && Nn > 0) ? Matrix.getSubMatrix(Qz, 0, Np, Np, Np + Nn) : Matrix.zeros(Np, Nn);
        Matrix Qznp = (Np > 0 && Nn > 0) ? Matrix.getSubMatrix(Qz, Np, Np + Nn, 0, Np) : Matrix.zeros(Nn, Np);
        Matrix Qznn = (Nn > 0) ? Matrix.getSubMatrix(Qz, Np, Np + Nn, Np, Np + Nn) : Matrix.zeros(0, 0);

        Map<String, Matrix> result = FluidFundamentalMatrices.FluidFundamentalMatrices(
                Qzpp, Qzpn, Qznp, Qznn, prec, null, null);
        Matrix Psi = result.get("P");
        Matrix K = result.get("K");
        Matrix U = result.get("U");

        Matrix Pm = Matrix.zeros(Np, Np + Nn);
        for (int i = 0; i < Np; i++) {
            Pm.set(i, i, 1.0);
            for (int j = 0; j < Nn; j++) {
                Pm.set(i, Np + j, Psi.get(i, j));
            }
        }

        Matrix iCn = (Nn > 0) ? Matrix.getSubMatrix(absRi, Np, Np + Nn, Np, Np + Nn) : Matrix.zeros(0, 0);
        Matrix iCp = (Np > 0) ? Matrix.getSubMatrix(absRi, 0, Np, 0, Np) : Matrix.zeros(0, 0);

        Matrix QvPz = (Nz > 0) ? Matrix.getSubMatrix(Qv, Nz, Nz + Np, 0, Nz) : Matrix.zeros(Np, 0);
        Matrix QvNz = (Nz > 0) ? Matrix.getSubMatrix(Qv, Nz + Np, N, 0, Nz) : Matrix.zeros(Nn, 0);

        Matrix clo;
        if (Nz > 0) {
            Matrix negQv00 = Qv00.neg();
            Matrix iQv00 = negQv00.pinv();
            Matrix leftPart = iCp.mult(QvPz).add(1.0, Psi.mult(iCn).mult(QvNz)).mult(iQv00);
            Matrix rightPart = Pm.mult(absRi);
            clo = Matrix.zeros(Np, N);
            for (int i = 0; i < Np; i++) {
                for (int j = 0; j < Nz; j++) {
                    clo.set(i, j, leftPart.get(i, j));
                }
                for (int j = 0; j < Np + Nn; j++) {
                    clo.set(i, Nz + j, rightPart.get(i, j));
                }
            }
        } else {
            clo = Pm.mult(absRi);
        }

        Matrix mass0;
        Matrix ini;

        if (Q0 == null) {
            clo = clo.mult(P);

            Matrix onesN = Matrix.ones(N, 1);
            Matrix onesNn = Matrix.ones(Nn, 1);
            Matrix onesNz = (Nz > 0) ? Matrix.ones(Nz, 1) : Matrix.zeros(0, 1);

            Matrix negK = K.neg();
            Matrix invNegK = negK.pinv();  // robust: -K is well-conditioned but high-dim det can be tiny

            Matrix Ua;
            if (Nz > 0) {
                Matrix negQv00 = Qv00.neg();
                Matrix iQv00 = negQv00.pinv();
                Matrix term1 = iCn.mult(QvNz).mult(iQv00).mult(onesNz);
                Matrix term2 = iCn.mult(onesNn);
                Matrix term3 = Qznp.mult(invNegK).mult(clo).mult(onesN);
                Ua = term1.add(1.0, term2).add(1.0, term3);
            } else {
                Matrix term2 = iCn.mult(onesNn);
                Matrix term3 = Qznp.mult(invNegK).mult(clo).mult(onesN);
                Ua = term2.add(1.0, term3);
            }

            Matrix UAugmented = Matrix.zeros(Nn, Nn + 1);
            for (int i = 0; i < Nn; i++) {
                for (int j = 0; j < Nn; j++) {
                    UAugmented.set(i, j, U.get(i, j));
                }
                UAugmented.set(i, Nn, Ua.get(i, 0));
            }
            Matrix rhs = Matrix.zeros(Nn + 1, 1);
            rhs.set(Nn, 0, 1.0);

            Matrix UAugT = UAugmented.transpose();
            Matrix pmT = UAugT.leftMatrixDivide(rhs);
            Matrix pm = pmT.transpose();

            Matrix mass0Unsorted = Matrix.zeros(1, N);
            if (Nz > 0) {
                Matrix negQv00 = Qv00.neg();
                Matrix iQv00 = negQv00.pinv();
                Matrix zeroTerm = pm.mult(iCn).mult(QvNz).mult(iQv00);
                for (int j = 0; j < Nz; j++) {
                    mass0Unsorted.set(0, j, zeroTerm.get(0, j));
                }
            }
            Matrix negTerm = pm.mult(iCn);
            for (int j = 0; j < Nn; j++) {
                mass0Unsorted.set(0, Nz + Np + j, negTerm.get(0, j));
            }
            mass0 = mass0Unsorted.mult(P);

            ini = pm.mult(Qznp);
        } else {
            Matrix Q0v = P.mult(Q0).mult(iP);

            Matrix M = Matrix.zeros(N, N);

            Matrix negCloRv = clo.mult(Rv).neg();
            for (int i = 0; i < Np; i++) {
                for (int j = 0; j < N; j++) {
                    M.set(i, j, negCloRv.get(i, j));
                }
            }

            for (int i = 0; i < Nn; i++) {
                for (int j = 0; j < N; j++) {
                    M.set(Np + i, j, Q0v.get(Nz + Np + i, j));
                }
            }

            for (int i = 0; i < Nz; i++) {
                for (int j = 0; j < N; j++) {
                    M.set(Np + Nn + i, j, Q0v.get(i, j));
                }
            }

            Matrix negK = K.neg();
            Matrix invNegK = negK.pinv();  // robust: -K is well-conditioned but high-dim det can be tiny
            Matrix invKClo = invNegK.mult(clo);
            Matrix Ma = Matrix.zeros(N, 1);
            for (int i = 0; i < Np; i++) {
                double rowSum = 0.0;
                for (int j = 0; j < N; j++) {
                    rowSum += invKClo.get(i, j);
                }
                Ma.set(i, 0, rowSum);
            }
            for (int i = Np; i < N; i++) {
                Ma.set(i, 0, 1.0);
            }

            Matrix MAugmented = Matrix.zeros(N, N + 1);
            for (int i = 0; i < N; i++) {
                for (int j = 0; j < N; j++) {
                    MAugmented.set(i, j, M.get(i, j));
                }
                MAugmented.set(i, N, Ma.get(i, 0));
            }
            Matrix rhs = Matrix.zeros(N + 1, 1);
            rhs.set(N, 0, 1.0);

            Matrix MAugT = MAugmented.transpose();
            Matrix solT = MAugT.leftMatrixDivide(rhs);
            Matrix sol = solT.transpose();

            ini = Matrix.zeros(1, Np);
            for (int j = 0; j < Np; j++) {
                ini.set(0, j, sol.get(0, j));
            }

            clo = clo.mult(P);

            Matrix mass0Unsorted = Matrix.zeros(1, N);
            for (int j = 0; j < Nz; j++) {
                mass0Unsorted.set(0, j, sol.get(0, Np + Nn + j));
            }
            for (int j = 0; j < Nn; j++) {
                mass0Unsorted.set(0, Nz + Np + j, sol.get(0, Np + j));
            }
            mass0 = mass0Unsorted.mult(P);
        }

        return new GeneralFluidSolution(mass0, ini, K, clo);
    }
}
