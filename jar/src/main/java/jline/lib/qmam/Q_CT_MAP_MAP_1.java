/**
 * @file Q_CT_MAP_MAP_1 - Continuous-Time MAP/MAP/1 Queue Analyzer
 *
 * Computes queue length, sojourn time and waiting time distribution for a
 * continuous-time MAP/MAP/1/FCFS queue.
 *
 * Based on the Q-MAM library by Benny Van Houdt.
 *
 * @since LINE 3.1.0
 */
package jline.lib.qmam;

import java.util.ArrayList;
import java.util.List;
import java.util.Map;

import jline.api.mc.Ctmc_solve;
import jline.lib.smc.QBD_CR;
import jline.lib.smc.QBD_FI;
import jline.lib.smc.QBD_pi;
import jline.lib.smc.Stat;
import jline.util.matrix.ComplexMatrix;
import jline.util.matrix.Matrix;

public final class Q_CT_MAP_MAP_1 {
    private Q_CT_MAP_MAP_1() {}

    /**
     * Computes queue length and time distributions for a MAP/MAP/1/FCFS queue.
     */
    public static MAPMAP1Result qCtMapMap1(Matrix C0, Matrix C1, Matrix D0, Matrix D1, MAPMAP1Options options) {
        int ma = C0.getNumRows();
        int ms = D0.getNumRows();
        int mtot = ma * ms;

        // Validate dimensions
        if (C0.getNumCols() != ma || C1.getNumRows() != ma || C1.getNumCols() != ma) {
            throw new IllegalArgumentException("Arrival process matrices must be ma x ma");
        }
        if (D0.getNumCols() != ms || D1.getNumRows() != ms || D1.getNumCols() != ms) {
            throw new IllegalArgumentException("Service process matrices must be ms x ms");
        }

        // Test the load of the queue
        Matrix invC0 = C0.scale(-1.0).inv();
        Matrix piA = Stat.stat(C1.mult(invC0));
        double lambda = piA.mult(C1).elementSum();

        Matrix invD0 = D0.scale(-1.0).inv();
        Matrix piS = Stat.stat(D1.mult(invD0));
        double mu = piS.mult(D1).elementSum();

        double load = lambda / mu;
        if (load >= 1) {
            throw new IllegalArgumentException("The load " + load + " of the system exceeds one");
        }

        // Compute classic QBD blocks A0, A1, A2
        Matrix eyeMa = Matrix.eye(ma);
        Matrix eyeMs = Matrix.eye(ms);

        Matrix Am1 = eyeMa.kron(D1);
        Matrix A0 = eyeMa.kron(D0).add(C0.kron(eyeMs));
        Matrix A1 = C1.kron(eyeMs);
        Matrix B0 = C0.kron(eyeMs);

        // Compute G and R using appropriate solver
        Map<String, Matrix> qbdResult;
        if (options.mode.contains("FI")) {
            qbdResult = QBD_FI.QBD_FI(Am1, A0, A1, null, options.verbose > 0 ? Integer.valueOf(1) : null, null, null, null);
        } else {
            qbdResult = QBD_CR.QBD_CR(Am1, A0, A1, null, options.verbose > 0 ? Integer.valueOf(1) : null, null, null);
        }

        Matrix R = qbdResult.get("R");

        // Stationary distribution. The boundary vector pi0 is the stationary
        // vector of the level-0 censored generator (B0 + R*Am1); B0 = kron(C0, I)
        // is the empty-queue local block and Am1 = kron(I, D1) a service
        // completion returning to level 0. NB: QBD_pi computes stat(B1 + R*B0)
        // after a CT->DT conversion, which does not yield this boundary vector
        // for continuous-time blocks (it leaves a large residual and over-weights
        // the tail), so the boundary is solved directly here.
        Matrix Mbound = B0.add(R.mult(Am1));
        Matrix pi0 = Ctmc_solve.ctmc_solve(Mbound); // (1 x mtot), sums to 1
        Matrix ImRinv = Matrix.eye(mtot).sub(R).inv();
        double norm = pi0.mult(ImRinv).mult(Matrix.ones(mtot, 1)).get(0, 0);
        pi0 = pi0.scale(1.0 / norm);

        // Queue length distribution: pi_n = pi0 * R^n, aggregated over phases.
        java.util.List<Double> qlList = new java.util.ArrayList<Double>();
        Matrix pk = pi0.copy();
        double mass = pk.elementSum();
        qlList.add(mass);
        int guard = 1;
        while (mass < 1.0 - 1e-12 && guard < options.maxNumComp) {
            pk = pk.mult(R);
            double s = pk.elementSum();
            qlList.add(s);
            mass += s;
            guard++;
        }
        Matrix ql = new Matrix(1, qlList.size());
        for (int i = 0; i < qlList.size(); i++) {
            ql.set(0, i, qlList.get(i));
        }

        // Compute Sojourn and Waiting time PH representation
        Matrix LM = C1.kron(D1);
        Matrix eyeD0 = eyeMa.kron(D0);
        Matrix C0eye = C0.kron(eyeMs);

        // Compute T iteratively
        Matrix Told = Matrix.zeros(mtot, mtot);
        Matrix Tnew = eyeD0.copy();
        // Pre-loop default so L is definitely assigned; overwritten on the first
        // (always-executed, since Tnew=eyeD0 != Told=0) iteration of either branch.
        Matrix L = Matrix.zeros(mtot, mtot);

        if (options.mode.contains("Direct")) {
            // Direct method
            Matrix eyeC0eye = Matrix.eye(mtot).kron(C0eye);
            while (Matrix.infNorm(Told.sub(Tnew)) > 1e-10) {
                Told = Tnew.copy();

                // Solve: vec(L) = -vec(I) * (kron(T',I) + kron(I,C0eye))^{-1}
                Matrix kronTI = Tnew.transpose().kron(Matrix.eye(mtot));
                Matrix system = kronTI.add(eyeC0eye);
                Matrix negEyeVec = new Matrix(mtot * mtot, 1);
                for (int i = 0; i < mtot; i++) {
                    negEyeVec.set(i * mtot + i, 0, -1.0);
                }
                Matrix Lvec = new Matrix(mtot * mtot, 1);
                Matrix.solve(system, negEyeVec, Lvec);

                L = new Matrix(mtot, mtot);
                for (int i = 0; i < mtot; i++) {
                    for (int j = 0; j < mtot; j++) {
                        L.set(i, j, Lvec.get(j * mtot + i, 0));
                    }
                }

                Tnew = eyeD0.add(L.mult(LM));
            }
            // L already holds reshape(Lvec,mtot,mtot)' from the converged iteration
            // (per MATLAB Q_CT_MAP_MAP_1.m Direct branch: L=reshape(L,mtot,mtot)').
        } else {
            // Sylvester method using the complex Schur decomposition, mirroring
            // MATLAB Q_CT_MAP_MAP_1.m: [U,Tr]=schur(C0,'complex').
            java.util.Map<String, ComplexMatrix> schur = Q_Sylvest.schurDecomposition(C0);
            ComplexMatrix U = schur.get("U");
            ComplexMatrix Tr = schur.get("T");
            ComplexMatrix Ukron = Q_Sylvest.kronEye(U, ms);
            ComplexMatrix Trkron = Q_Sylvest.kronEye(Tr, ms);

            while (Matrix.infNorm(Told.sub(Tnew)) > 1e-10) {
                Told = Tnew.copy();
                L = Q_Sylvest.qSylvest(Ukron, Trkron, Tnew);
                Tnew = eyeD0.add(L.mult(LM));
            }
        }

        // Compute Smat
        Matrix thetaTot = piA.mult(C1).kron(piS).scale(1.0 / (mu * load));

        // Find non-zero entries
        List<Integer> nonzIndices = new ArrayList<Integer>();
        for (int i = 0; i < mtot; i++) {
            if (thetaTot.get(0, i) > 0) {
                nonzIndices.add(i);
            }
        }

        if (nonzIndices.isEmpty()) {
            return new MAPMAP1Result(ql, null, null, null);
        }

        int nz = nonzIndices.size();
        Matrix thetaTotRed = new Matrix(1, nz);
        for (int i = 0; i < nz; i++) {
            thetaTotRed.set(0, i, thetaTot.get(0, nonzIndices.get(i)));
        }

        Matrix TnewReduced = new Matrix(nz, nz);
        for (int i = 0; i < nz; i++) {
            for (int j = 0; j < nz; j++) {
                TnewReduced.set(i, j, Tnew.get(nonzIndices.get(i), nonzIndices.get(j)));
            }
        }

        Matrix diagThetaRed = Matrix.diag(thetaTotRed.getRow(0).toArray1D());
        Matrix diagThetaRedInv = diagThetaRed.inv();
        Matrix Smat = diagThetaRedInv.mult(TnewReduced.transpose()).mult(diagThetaRed);

        // Alpha vector of PH representation of Sojourn time
        Matrix D1sumCol = D1.sumRows();
        Matrix sojAlphaFull = Matrix.diag(thetaTot.getRow(0).toArray1D())
                .mult(Matrix.ones(ma, 1).kron(D1sumCol))
                .scale(load / lambda);
        Matrix sojAlpha = new Matrix(1, nz);
        for (int i = 0; i < nz; i++) {
            sojAlpha.set(0, i, sojAlphaFull.get(nonzIndices.get(i), 0));
        }

        // Alpha vector of PH representation of Waiting time
        Matrix C1sumCol = C1.sumRows();
        Matrix waitAlphaFull = Matrix.diag(thetaTot.getRow(0).toArray1D())
                .mult(L)
                .mult(C1sumCol.kron(D1sumCol))
                .scale(load / lambda);
        Matrix waitAlpha = new Matrix(1, nz);
        for (int i = 0; i < nz; i++) {
            waitAlpha.set(0, i, waitAlphaFull.get(nonzIndices.get(i), 0));
        }

        return new MAPMAP1Result(ql, sojAlpha.transpose(), waitAlpha.transpose(), Smat);
    }

    public static MAPMAP1Result qCtMapMap1(Matrix C0, Matrix C1, Matrix D0, Matrix D1) {
        return qCtMapMap1(C0, C1, D0, D1, new MAPMAP1Options());
    }
}
