package jline.api.aoi;

import jline.lib.butools.mam.FluidTools;
import jline.util.matrix.Matrix;

/**
 * Single-buffer AoI solver using Markovian Fluid Queues.
 */
public final class Aoi_solve_singlebuffer {
    private Aoi_solve_singlebuffer() {}

    public static AoiMfqResult aoi_solve_singlebuffer(double lambda, Matrix sigma, Matrix S, double r) {
        int l = S.getNumCols();
        Matrix nu = S.mult(Matrix.ones(l, 1)).scale(-1.0);
        int z1 = l + 2;

        Matrix Q1 = Matrix.zeros(z1, z1);
        for (int i = 0; i < l; i++) {
            for (int j = 0; j < l; j++) Q1.set(i, j, S.get(i, j));
            Q1.set(i, l, nu.get(i, 0));
        }
        Q1.set(l, l, -lambda);
        Q1.set(l, l + 1, lambda);

        Matrix R1 = Matrix.eye(z1);
        R1.set(z1 - 2, z1 - 2, -1.0);
        R1.set(z1 - 1, z1 - 1, -1.0);

        Matrix Qtilde1 = Matrix.zeros(z1, z1);
        for (int j = 0; j < l; j++) Qtilde1.set(l, j, lambda * sigma.get(0, j));
        Qtilde1.set(l, l, -lambda);
        for (int j = 0; j < l; j++) Qtilde1.set(l + 1, j, sigma.get(0, j));
        Qtilde1.set(l + 1, l + 1, -1.0);

        Matrix e1 = Matrix.ones(z1, 1);
        Matrix QplusEE = Q1.add(1.0, e1.mult(e1.transpose()));
        Matrix pikT = new Matrix(z1, 1);
        Matrix.solveSafe(QplusEE.transpose(), e1, pikT);
        Matrix pik = pikT.transpose();

        Matrix QR1 = Q1.mult(Matrix.inv(R1));
        Matrix xR = R1.mult(e1);
        Matrix xL = pik;
        double xLxR = xL.mult(xR).get(0, 0);
        Matrix A1mat = QR1.add(1.0 / xLxR, xR.mult(xL));

        // ordschur(U,T,'rhp') as in MATLAB solveSingleBuffer
        Matrix P1 = FluidTools.orderedSchurRhp(A1mat)[0];

        int a1 = 2;
        int b1 = l;
        Matrix Ae1 = P1.transpose().mult(QR1).mult(P1);
        Matrix Amat1 = new Matrix(b1, b1);
        for (int i = 0; i < b1; i++) for (int j = 0; j < b1; j++) Amat1.set(i, j, Ae1.get(a1 + i, a1 + j));
        Matrix Hmat1 = new Matrix(b1, z1);
        for (int i = 0; i < b1; i++) for (int j = 0; j < z1; j++) Hmat1.set(i, j, P1.get(j, a1 + i));
        Matrix Qts1 = new Matrix(a1, z1);
        for (int i = 0; i < a1; i++) for (int j = 0; j < z1; j++) Qts1.set(i, j, Qtilde1.get(b1 + i, j));

        Matrix Ainv1 = Matrix.inv(Amat1);
        Matrix eqn1 = new Matrix(z1, z1 + 1);
        Matrix HR1 = Hmat1.mult(R1);
        Matrix AinvHones1 = Ainv1.mult(Hmat1.mult(Matrix.ones(z1, 1))).scale(-1.0);
        for (int i = 0; i < b1; i++) {
            for (int j = 0; j < z1; j++) eqn1.set(i, j, HR1.get(i, j));
            eqn1.set(i, z1, AinvHones1.get(i, 0));
        }
        for (int i = 0; i < a1; i++) {
            for (int j = 0; j < z1; j++) eqn1.set(b1 + i, j, -Qts1.get(i, j));
            eqn1.set(b1 + i, z1, 1.0);
        }
        Matrix rhs1 = new Matrix(z1 + 1, 1);
        rhs1.set(z1, 0, 1.0);
        // Overdetermined (z1+1) x z1 system: least-squares solve as in MATLAB mldivide
        Matrix sol1 = eqn1.transpose().leftMatrixDivide(rhs1);

        Matrix wait_g = new Matrix(1, b1);
        for (int i = 0; i < b1; i++) wait_g.set(0, i, sol1.get(i, 0));
        Matrix wait_d = new Matrix(a1, 1);
        for (int i = 0; i < a1; i++) wait_d.set(i, 0, sol1.get(b1 + i, 0));
        double c_0 = wait_d.get(0, 0);

        Matrix wait_A = Amat1.add(-r * lambda, Matrix.eye(l));
        Matrix selVec = new Matrix(z1, 1);
        selVec.set(l, 0, 1.0);
        selVec.set(l + 1, 0, r);
        Matrix wait_H = Hmat1.mult(selVec);

        Matrix wait_Ainv = Matrix.inv(wait_A);
        double n_1 = 1.0 / (-wait_g.mult(wait_Ainv).mult(wait_H).get(0, 0) + c_0);
        Matrix wait_g_scaled = wait_g.scale(n_1);

        Matrix Mdiag = wait_Ainv.mult(wait_H).scale(-1.0);
        Matrix Mmat = Matrix.zeros(l, l);
        for (int i = 0; i < l; i++) Mmat.set(i, i, Mdiag.get(i, 0));
        Matrix Minv = Matrix.inv(Mmat);
        Matrix B = Minv.mult(wait_A).mult(Mmat);
        Matrix beta = wait_g_scaled.mult(Mmat);
        double betaSum = 0.0;
        for (int j = 0; j < l; j++) betaSum += beta.get(0, j);
        double beta_0 = 1.0 - betaSum;
        Matrix psi = B.mult(Matrix.ones(l, 1)).scale(-1.0);

        int z2 = 4 * l + 2;
        int a2 = 1;
        int b2 = z2 - 1;

        Matrix Q2 = Matrix.zeros(z2, z2);
        for (int i = 0; i < l; i++) for (int j = 0; j < l; j++) Q2.set(i, j, B.get(i, j));
        Matrix psiSigma = psi.kron(sigma);
        for (int i = 0; i < l; i++) for (int j = 0; j < l; j++) Q2.set(i, l + j, psiSigma.get(i, j));
        Matrix SlI = S.add(-lambda, Matrix.eye(l));
        Matrix lamI = Matrix.eye(l).scale(lambda);
        for (int i = 0; i < l; i++) {
            for (int j = 0; j < l; j++) Q2.set(l + i, l + j, SlI.get(i, j));
            for (int j = 0; j < l; j++) Q2.set(l + i, 2 * l + j, lamI.get(i, j));
            Q2.set(l + i, 3 * l, nu.get(i, 0));
        }
        Matrix nuSigma = nu.kron(sigma);
        for (int i = 0; i < l; i++) {
            for (int j = 0; j < l; j++) Q2.set(2 * l + i, 2 * l + j, S.get(i, j));
            for (int j = 0; j < l; j++) Q2.set(2 * l + i, 3 * l + 1 + j, nuSigma.get(i, j));
        }
        Q2.set(3 * l, 3 * l, -lambda);
        for (int j = 0; j < l; j++) Q2.set(3 * l, 3 * l + 1 + j, lambda * sigma.get(0, j));
        for (int i = 0; i < l; i++) {
            for (int j = 0; j < l; j++) Q2.set(3 * l + 1 + i, 3 * l + 1 + j, S.get(i, j));
            Q2.set(3 * l + 1 + i, 4 * l + 1, nu.get(i, 0));
        }

        Matrix R2 = Matrix.eye(z2);
        R2.set(z2 - 1, z2 - 1, -1.0);
        Matrix Qtilde2 = Q2.copy();
        for (int j = 0; j < l; j++) Qtilde2.set(z2 - 1, j, beta.get(0, j));
        for (int j = 0; j < l; j++) Qtilde2.set(z2 - 1, l + j, beta_0 * sigma.get(0, j));
        Qtilde2.set(z2 - 1, z2 - 1, -1.0);

        Matrix u1 = new Matrix(z2, 1);
        for (int i = 0; i < z2 - 1; i++) u1.set(i, 0, 1.0);
        u1.set(z2 - 1, 0, -1.0);
        double norm_u1 = Math.sqrt(u1.transpose().mult(u1).get(0, 0));
        Matrix u = u1.copy();
        u.set(0, 0, u.get(0, 0) - norm_u1);
        double utu = u.transpose().mult(u).get(0, 0);
        Matrix P2 = Matrix.eye(z2).add(-2.0 / utu, u.mult(u.transpose()));

        Matrix QR2 = Q2.mult(Matrix.inv(R2));
        Matrix Ae2 = P2.transpose().mult(QR2).mult(P2);
        Matrix Amat2 = new Matrix(b2, b2);
        for (int i = 0; i < b2; i++) for (int j = 0; j < b2; j++) Amat2.set(i, j, Ae2.get(a2 + i, a2 + j));
        Matrix Hmat2 = new Matrix(b2, z2);
        for (int i = 0; i < b2; i++) for (int j = 0; j < z2; j++) Hmat2.set(i, j, P2.get(j, a2 + i));
        Matrix Qts2 = new Matrix(a2, z2);
        for (int i = 0; i < a2; i++) for (int j = 0; j < z2; j++) Qts2.set(i, j, Qtilde2.get(b2 + i, j));

        Matrix Ainv2 = Matrix.inv(Amat2);
        Matrix HR2 = Hmat2.mult(R2);
        Matrix AinvHones2 = Ainv2.mult(Hmat2.mult(Matrix.ones(z2, 1))).scale(-1.0);
        Matrix eqn2 = new Matrix(z2, z2 + 1);
        for (int i = 0; i < b2; i++) {
            for (int j = 0; j < z2; j++) eqn2.set(i, j, HR2.get(i, j));
            eqn2.set(i, z2, AinvHones2.get(i, 0));
        }
        for (int i = 0; i < a2; i++) {
            for (int j = 0; j < z2; j++) eqn2.set(b2 + i, j, -Qts2.get(i, j));
            eqn2.set(b2 + i, z2, 1.0);
        }
        Matrix rhs2 = new Matrix(z2 + 1, 1);
        rhs2.set(z2, 0, 1.0);
        // Overdetermined (z2+1) x z2 system: least-squares solve as in MATLAB mldivide
        Matrix sol2 = eqn2.transpose().leftMatrixDivide(rhs2);

        Matrix g = new Matrix(1, b2);
        for (int i = 0; i < b2; i++) g.set(0, i, sol2.get(i, 0));

        Matrix selectedAoI = new Matrix(z2, 1);
        for (int i = 3 * l; i < 4 * l + 1; i++) selectedAoI.set(i, 0, 1.0);

        Matrix aoiH = Hmat2.mult(selectedAoI);
        double normConst = g.mult(Ainv2).mult(aoiH).scale(-1.0).get(0, 0);
        Matrix aoiG = g.scale(1.0 / normConst);
        Matrix aoiA = Amat2;
        Matrix Ainv2_2 = Ainv2.mult(Ainv2);
        Matrix Ainv2_3 = Ainv2_2.mult(Ainv2);
        double aoiMean = aoiG.mult(Ainv2_2).mult(aoiH).get(0, 0);
        double aoiVar = -2.0 * aoiG.mult(Ainv2_3).mult(aoiH).get(0, 0) - aoiMean * aoiMean;

        Matrix selectedPAoI = new Matrix(z2, 1);
        for (int i = 0; i < l; i++) selectedPAoI.set(3 * l + 1 + i, 0, nu.get(i, 0));
        Matrix paoiH = Hmat2.mult(selectedPAoI);
        double normConstP = g.mult(Ainv2).mult(paoiH).scale(-1.0).get(0, 0);
        Matrix paoiG = g.scale(1.0 / normConstP);
        Matrix paoiA = Amat2;
        double paoiMean = paoiG.mult(Ainv2_2).mult(paoiH).get(0, 0);
        double paoiVar = -2.0 * paoiG.mult(Ainv2_3).mult(paoiH).get(0, 0) - paoiMean * paoiMean;

        return new AoiMfqResult(
                aoiG, aoiA, aoiH,
                aoiMean, Math.max(0.0, aoiVar),
                paoiG, paoiA, paoiH,
                paoiMean, Math.max(0.0, paoiVar),
                "singlebuffer", r);
    }
}
