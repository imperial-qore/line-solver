/**
 * @file Q_CT_PH_PH_1 - Continuous-Time PH/PH/1 Queue Analyzer
 *
 * Computes queue length and waiting time distribution for a continuous-time
 * PH/PH/1/FCFS queue.
 *
 * Based on the Q-MAM library by Benny Van Houdt.
 *
 * @since LINE 3.1.0
 */
package jline.lib.qmam;

import java.util.ArrayList;
import java.util.List;

import jline.io.InputOutput;
import jline.util.matrix.Matrix;

public final class Q_CT_PH_PH_1 {
    private Q_CT_PH_PH_1() {}

    /**
     * Computes queue length and waiting time distribution for a PH/PH/1/FCFS queue.
     */
    public static PHPH1Result qCtPhPh1(Matrix alpha, Matrix T, Matrix beta, Matrix S, PHPH1Options options) {
        // Validate dimensions
        int ma = alpha.getNumCols();
        int ms = beta.getNumCols();

        if (T.getNumRows() != ma || T.getNumCols() != ma) {
            throw new IllegalArgumentException("T matrix dimensions must match alpha");
        }
        if (S.getNumRows() != ms || S.getNumCols() != ms) {
            throw new IllegalArgumentException("S matrix dimensions must match beta");
        }

        // Arrival process
        Matrix t = T.mult(Matrix.ones(ma, 1)).scale(-1.0);
        Matrix invNegT = T.scale(-1.0).inv();
        double avgT = alpha.mult(invNegT).mult(Matrix.ones(ma, 1)).get(0, 0);

        // Service process
        Matrix s = S.mult(Matrix.ones(ms, 1)).scale(-1.0);
        Matrix invNegS = S.scale(-1.0).inv();
        double avgS = beta.mult(invNegS).mult(Matrix.ones(ms, 1)).get(0, 0);

        int mtot = ms * ma;
        double rho = avgS / avgT;

        if (rho >= 1) {
            throw new IllegalArgumentException("The load " + rho + " of the system exceeds one");
        }

        // Compute classic QBD blocks A0, A1, A2
        Matrix eyeMa = Matrix.eye(ma);
        Matrix eyeMs = Matrix.eye(ms);

        Matrix A0 = t.mult(alpha).kron(eyeMs);
        Matrix A1 = T.kron(eyeMs).add(eyeMa.kron(S));

        // Compute QBD blocks in Latouche-Ramaswami approach
        Matrix invmA1 = A1.scale(-1.0).inv();
        Matrix tEye = t.kron(eyeMs);
        Matrix alphaEye = alpha.kron(eyeMs);
        Matrix eyeBeta = eyeMa.kron(beta);
        Matrix eyeS = eyeMa.kron(s);

        Matrix A0pp = alphaEye.mult(invmA1).mult(tEye);
        Matrix A0mp = eyeBeta.mult(invmA1).mult(tEye);
        Matrix A2pm = alphaEye.mult(invmA1).mult(eyeS);
        Matrix A2mm = eyeBeta.mult(invmA1).mult(eyeS);

        // Construct block matrices
        Matrix A0n = new Matrix(ms + ma, ms + ma);
        Matrix A2n = new Matrix(ms + ma, ms + ma);

        // A0n = [A0pp zeros(ms,ma); A0mp zeros(ma)]
        for (int i = 0; i < ms; i++) {
            for (int j = 0; j < ms; j++) {
                A0n.set(i, j, A0pp.get(i, j));
            }
        }
        for (int i = 0; i < ma; i++) {
            for (int j = 0; j < ms; j++) {
                A0n.set(ms + i, j, A0mp.get(i, j));
            }
        }

        // A2n = [zeros(ms) A2pm; zeros(ma,ms) A2mm]
        for (int i = 0; i < ms; i++) {
            for (int j = 0; j < ma; j++) {
                A2n.set(i, ms + j, A2pm.get(i, j));
            }
        }
        for (int i = 0; i < ma; i++) {
            for (int j = 0; j < ma; j++) {
                A2n.set(ms + i, ms + j, A2mm.get(i, j));
            }
        }

        // Compute matrix Gamma: NE corner of matrix G using cyclic reduction
        Matrix itB0 = A0n.copy();
        Matrix itB2 = A2n.copy();
        Matrix Gamma = new Matrix(ms, ma);

        // Initialize Gamma from itB2
        for (int i = 0; i < ms; i++) {
            for (int j = 0; j < ma; j++) {
                Gamma.set(i, j, itB2.get(i, ms + j));
            }
        }

        Matrix itT = itB0.copy();
        double check = 1.0;
        int numit = 1;
        Matrix eyeMaMsTotal = Matrix.eye(ma + ms);

        while (check > 1e-13) {
            Matrix itA1 = itB0.mult(itB2).add(itB2.mult(itB0));
            Matrix invFactor = eyeMaMsTotal.sub(itA1).inv();
            itB0 = invFactor.mult(itB0.mult(itB0));
            itB2 = invFactor.mult(itB2.mult(itB2));

            Matrix tmp = itT.mult(itB2);
            for (int i = 0; i < ms; i++) {
                for (int j = 0; j < ma; j++) {
                    Gamma.set(i, j, Gamma.get(i, j) + tmp.get(i, ms + j));
                }
            }
            itT = itT.mult(itB0);

            check = Matrix.ones(1, ms).mult(Gamma).mult(Matrix.ones(ma, 1)).sub(Matrix.singleton(1.0)).infinityNorm();
            numit++;
        }

        Matrix Gm = Matrix.eye(ma).sub(A0mp.mult(Gamma)).inv().mult(A2mm);
        Matrix RGam = A0pp.mult(Matrix.eye(ms).sub(Gamma.mult(A0mp)).inv());

        // Compute queue length distribution
        Matrix Gstar = invmA1.mult(eyeS.mult(eyeBeta))
                .add(invmA1.mult(tEye).mult(Gamma).mult(Gm).mult(eyeBeta));
        Matrix Rstar = A0.mult(A1.add(A0.mult(Gstar)).scale(-1.0).inv());

        // Compute pi_0
        Matrix betaGamma = beta.mult(Gamma);
        Matrix betaGammaInvNegT = betaGamma.mult(invNegT);
        double normFactor = betaGammaInvNegT.mult(Matrix.ones(ma, 1)).get(0, 0);
        Matrix pi0Unnorm = betaGammaInvNegT.scale((1 - rho) / normFactor);
        Matrix pi0 = pi0Unnorm.kron(beta);

        // Compute pi_1, pi_2, ...
        List<Matrix> piLevels = new ArrayList<Matrix>();
        piLevels.add(pi0);

        double sumpi = pi0.elementSum();
        numit = 1;

        while (sumpi < 1 - 1e-10 && numit < 1 + options.maxNumComp) {
            Matrix piNext = piLevels.get(piLevels.size() - 1).mult(Rstar);
            piLevels.add(piNext);
            numit++;
            sumpi += piNext.elementSum();

            if (options.verbose > 0 && numit % options.verbose == 0) {
                System.out.println("Accumulated mass after " + numit + " iterations: " + sumpi);
            }
        }

        // Compute queue length distribution
        Matrix ql = new Matrix(1, piLevels.size());
        for (int i = 0; i < piLevels.size(); i++) {
            ql.set(0, i, piLevels.get(i).elementSum());
        }

        if (numit == 1 + options.maxNumComp) {
            InputOutput.line_warning("Q_CT_PH_PH_1", "Maximum Number of Components " + (numit - 1) + " reached");
        }

        // Compute waiting time PH representation
        Matrix sigtilde = beta.mult(invNegS).scale(1.0 / beta.mult(invNegS).mult(Matrix.ones(ms, 1)).get(0, 0));
        Matrix Delta = Matrix.diag(sigtilde.getRow(0).toArray1D());

        // wait_T = inv(Delta) * (S + RGam * s * beta)' * Delta
        Matrix innerMatrix = S.add(RGam.mult(s).mult(beta));
        Matrix DeltaInv = Delta.inv();
        Matrix waitT = DeltaInv.mult(innerMatrix.transpose()).mult(Delta);

        // Compute wait_alpha
        Matrix sigrho = sigtilde.scale(rho);
        Matrix theta = s.transpose().mult(Delta).scale(-1.0 / beta.mult(invNegS).mult(Matrix.ones(ms, 1)).get(0, 0));
        Matrix D = DeltaInv.mult(RGam.transpose()).mult(Delta);
        Matrix thetaD = theta.mult(D);
        double thetaDOnes = thetaD.mult(Matrix.ones(ms, 1)).get(0, 0);

        double prob0 = beta.mult(Matrix.eye(ms).sub(RGam).inv()).mult(Matrix.ones(ms, 1)).get(0, 0);
        Matrix waitAlpha = thetaD.scale((1 - 1.0 / prob0) / thetaDOnes);

        return new PHPH1Result(ql, waitAlpha, waitT);
    }

    public static PHPH1Result qCtPhPh1(Matrix alpha, Matrix T, Matrix beta, Matrix S) {
        return qCtPhPh1(alpha, T, beta, S, new PHPH1Options());
    }
}
