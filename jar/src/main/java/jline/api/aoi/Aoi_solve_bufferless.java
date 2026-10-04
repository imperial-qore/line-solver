/**
 * @file Bufferless AoI solver using Markovian Fluid Queues
 *
 * Computes matrix exponential (ME) representations for Age of Information
 * and Peak AoI distributions in PH/PH/1/1 (or PH/PH/1/1*) bufferless systems.
 *
 * Based on: aoi-fluid toolbox by Ozancan Dogan, Nail Akar, Eray Unsal Atay
 * BSD 2-Clause License, 2020
 *
 * Reference: Eqns (8), (13), and Algorithm 1 from the aoi-fluid paper.
 *
 * @since LINE 3.0
 */
package jline.api.aoi;

import jline.util.matrix.Matrix;

public final class Aoi_solve_bufferless {
    private Aoi_solve_bufferless() {}

    /**
     * Solve bufferless (PH/PH/1/1 or PH/PH/1/1*) AoI system using MFQ.
     *
     * Computes matrix exponential parameters for AoI and Peak AoI distributions
     * along with their first two moments.
     *
     * @param tau Initial probability vector of the arrival process (1 x k)
     * @param T Sub-generator matrix of the arrival process (k x k)
     * @param sigma Initial probability vector of the service process (1 x l)
     * @param S Sub-generator matrix of the service process (l x l)
     * @param p Packet preemption probability (0 = FCFS, 1 = preemptive)
     * @return AoiMfqResult with ME parameters and moments
     */
    public static AoiMfqResult aoi_solve_bufferless(Matrix tau, Matrix T, Matrix sigma, Matrix S, double p) {
        int k = T.getNumCols();
        int l = S.getNumCols();
        Matrix kappa = T.mult(Matrix.ones(k, 1)).scale(-1.0);
        Matrix nu = S.mult(Matrix.ones(l, 1)).scale(-1.0);
        int z = 2 * k * l + k + 1;
        int a = 1;
        int b = z - 1;

        // Q11 = kron(eye(k),S) + kron(T,eye(l)) + (1-p)*kron(kron(kappa,tau),eye(l))
        Matrix Q11 = Matrix.eye(k).kron(S)
                .add(1.0, T.kron(Matrix.eye(l)))
                .add(1.0 - p, kappa.kron(tau).kron(Matrix.eye(l)));

        // Q33 = Q11 + p*kron(kappa, kron(ones(l,1), kron(tau, sigma)))
        Matrix Q33 = Q11.copy().add(p, kappa.kron(Matrix.ones(l, 1).kron(tau.kron(sigma))));

        // Construct full Q matrix (z x z)
        Matrix Q = Matrix.zeros(z, z);

        Matrix block12 = Matrix.eye(k).kron(nu);
        Matrix block14 = kappa.kron(Matrix.ones(l, 1));
        for (int i = 0; i < k * l; i++) {
            for (int j = 0; j < k * l; j++) Q.set(i, j, Q11.get(i, j));
            for (int j = 0; j < k; j++) Q.set(i, k * l + j, block12.get(i, j));
            Q.set(i, 2 * k * l + k, p * block14.get(i, 0));
        }

        Matrix block23 = kappa.kron(tau.kron(sigma));
        for (int i = 0; i < k; i++) {
            for (int j = 0; j < k; j++) Q.set(k * l + i, k * l + j, T.get(i, j));
            for (int j = 0; j < k * l; j++) Q.set(k * l + i, k * l + k + j, block23.get(i, j));
        }

        Matrix block34 = Matrix.ones(k, 1).kron(nu);
        for (int i = 0; i < k * l; i++) {
            for (int j = 0; j < k * l; j++) Q.set(k * l + k + i, k * l + k + j, Q33.get(i, j));
            Q.set(k * l + k + i, 2 * k * l + k, block34.get(i, 0));
        }

        Matrix R = Matrix.eye(z);
        R.set(z - 1, z - 1, -1.0);

        Matrix Qtilde = Q.copy();
        Matrix tauSigma = tau.kron(sigma);
        for (int j = 0; j < k * l; j++) {
            Qtilde.set(z - 1, j, tauSigma.get(0, j));
        }
        Qtilde.set(z - 1, z - 1, -1.0);

        // Construction of orthogonal matrix P (Householder reflection)
        Matrix u1 = new Matrix(z, 1);
        for (int i = 0; i < z - 1; i++) u1.set(i, 0, 1.0);
        u1.set(z - 1, 0, -1.0);
        double norm_u1 = Math.sqrt(u1.transpose().mult(u1).get(0, 0));

        Matrix u = u1.copy();
        u.set(0, 0, u.get(0, 0) - norm_u1);

        double utu = u.transpose().mult(u).get(0, 0);
        Matrix uut = u.mult(u.transpose());
        Matrix P = Matrix.eye(z).add(-2.0 / utu, uut);

        Matrix QR = Q.mult(Matrix.inv(R));
        Matrix Ae = P.transpose().mult(QR).mult(P);

        Matrix Amat = new Matrix(b, b);
        for (int i = 0; i < b; i++) {
            for (int j = 0; j < b; j++) {
                Amat.set(i, j, Ae.get(a + i, a + j));
            }
        }

        Matrix Hmat = new Matrix(b, z);
        for (int i = 0; i < b; i++) {
            for (int j = 0; j < z; j++) {
                Hmat.set(i, j, P.get(j, a + i));
            }
        }

        Matrix Qtildestar = new Matrix(a, z);
        for (int i = 0; i < a; i++) {
            for (int j = 0; j < z; j++) {
                Qtildestar.set(i, j, Qtilde.get(b + i, j));
            }
        }

        Matrix HR = Hmat.mult(R);
        Matrix Ainv = Matrix.inv(Amat);
        Matrix AinvHones = Ainv.mult(Hmat.mult(Matrix.ones(z, 1))).scale(-1.0);

        Matrix eqnMatrix2 = new Matrix(z, z + 1);
        for (int i = 0; i < b; i++) {
            for (int j = 0; j < z; j++) eqnMatrix2.set(i, j, HR.get(i, j));
            eqnMatrix2.set(i, z, AinvHones.get(i, 0));
        }
        for (int i = 0; i < a; i++) {
            for (int j = 0; j < z; j++) eqnMatrix2.set(b + i, j, -Qtildestar.get(i, j));
            eqnMatrix2.set(b + i, z, 1.0);
        }

        Matrix rhs = new Matrix(z + 1, 1);
        rhs.set(z, 0, 1.0);

        // Overdetermined (z+1) x z system: least-squares solve as in MATLAB mldivide
        Matrix solution = eqnMatrix2.transpose().leftMatrixDivide(rhs);

        Matrix g = new Matrix(1, b);
        for (int i = 0; i < b; i++) g.set(0, i, solution.get(i, 0));

        // AoI: restrict to states in Phases 2 and 3 only
        Matrix selectedAoI = new Matrix(z, 1);
        for (int i = k * l; i < 2 * k * l + k; i++) selectedAoI.set(i, 0, 1.0);

        Matrix aoiH = Hmat.mult(selectedAoI);
        double normConst = g.mult(Ainv).mult(aoiH).scale(-1.0).get(0, 0);
        Matrix aoiG = g.scale(1.0 / normConst);
        Matrix aoiA = Amat;

        Matrix Ainv2 = Ainv.mult(Ainv);
        Matrix Ainv3 = Ainv2.mult(Ainv);
        double aoiMean = aoiG.mult(Ainv2).mult(aoiH).get(0, 0);
        double aoiVar = -2.0 * aoiG.mult(Ainv3).mult(aoiH).get(0, 0) - aoiMean * aoiMean;

        // PAoI: restrict to states in Phase 3 scaled by nu
        Matrix selectedPAoI = new Matrix(z, 1);
        Matrix kronOnesNu = Matrix.ones(k, 1).kron(nu);
        for (int i = 0; i < k * l; i++) selectedPAoI.set(k * l + k + i, 0, kronOnesNu.get(i, 0));

        Matrix paoiH = Hmat.mult(selectedPAoI);
        double normConstP = g.mult(Ainv).mult(paoiH).scale(-1.0).get(0, 0);
        Matrix paoiG = g.scale(1.0 / normConstP);
        Matrix paoiA = Amat;

        double paoiMean = paoiG.mult(Ainv2).mult(paoiH).get(0, 0);
        double paoiVar = -2.0 * paoiG.mult(Ainv3).mult(paoiH).get(0, 0) - paoiMean * paoiMean;

        return new AoiMfqResult(
                aoiG, aoiA, aoiH,
                aoiMean, Math.max(0.0, aoiVar),
                paoiG, paoiA, paoiH,
                paoiMean, Math.max(0.0, paoiVar),
                "bufferless", p);
    }
}
