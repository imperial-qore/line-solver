package jline.lib.qmam;

import java.util.ArrayList;
import java.util.Collections;
import java.util.List;

import jline.lib.smc.Stat;
import jline.util.matrix.ComplexMatrix;
import jline.util.matrix.Matrix;

/**
 * Q_CT_MMAPK_PHK_1 - Continuous-Time MMAP[K]/PH[K]/1 Queue Analyzer.
 * Based on the Q-MAM library by Benny Van Houdt.
 */
public final class Q_CT_MMAPK_PHK_1 {
    private Q_CT_MMAPK_PHK_1() {}

    public static MMAPKPHK1Result qCtMmapkPhk1(Matrix D0, List<Matrix> D, List<Matrix> alpha, List<Matrix> S) {
        return qCtMmapkPhk1(D0, D, alpha, S, new MMAPKPHK1Options());
    }

    public static MMAPKPHK1Result qCtMmapkPhk1(Matrix D0, List<Matrix> D, List<Matrix> alpha, List<Matrix> S, MMAPKPHK1Options options) {
        int K = alpha.size();
        int m = D0.getNumRows();
        if (D.size() != K) throw new IllegalArgumentException("Number of D matrices must equal number of service types");
        if (S.size() != K) throw new IllegalArgumentException("Number of S matrices must equal number of service types");

        int[] smk = new int[K + 1];
        smk[0] = 0;
        for (int i = 0; i < K; i++) smk[i + 1] = smk[i] + alpha.get(i).getNumCols();
        int mser = smk[K];
        int mtot = mser * m;

        Matrix Dsum = D0.copy();
        for (int i = 0; i < K; i++) Dsum = Dsum.add(D.get(i));

        Matrix diagDsum = new Matrix(m, 1);
        Matrix.extractDiag(Dsum, diagDsum);
        double maxDiagDsum = diagDsum.scale(-1.0).elementMax();
        Matrix piMat = Stat.stat(Matrix.eye(m).sub(Dsum.scale(1.0 / maxDiagDsum)));

        double[] lambdas = new double[K];
        double[] mus = new double[K];
        List<Matrix> beta = new ArrayList<Matrix>();
        for (int i = 0; i < K; i++) {
            lambdas[i] = piMat.mult(D.get(i)).elementSum();
            Matrix sSum = S.get(i).sumRows();
            Matrix temp = S.get(i).sub(sSum.mult(alpha.get(i)));
            if (S.get(i).getNumRows() > 1) {
                Matrix diagTemp = new Matrix(S.get(i).getNumRows(), 1);
                Matrix.extractDiag(temp, diagTemp);
                double maxDiagTemp = diagTemp.scale(-1.0).elementMax();
                beta.add(Stat.stat(Matrix.eye(S.get(i).getNumRows()).sub(temp.scale(1.0 / maxDiagTemp))));
            } else {
                beta.add(Matrix.singleton(1.0));
            }
            mus[i] = -beta.get(i).mult(S.get(i).sumRows()).get(0, 0);
        }

        double lambda = piMat.mult(Dsum.sub(D0)).elementSum();
        double load = 0.0;
        for (int i = 0; i < K; i++) load += lambdas[i] / mus[i];
        if (load >= 1) throw new IllegalArgumentException("The load " + load + " of the system exceeds one");

        Matrix Tser = Matrix.zeros(mser, mser);
        Matrix tser = Matrix.zeros(mser, 1);
        for (int i = 0; i < K; i++) {
            for (int r = 0; r < S.get(i).getNumRows(); r++) {
                for (int c = 0; c < S.get(i).getNumCols(); c++) {
                    Tser.set(smk[i] + r, smk[i] + c, S.get(i).get(r, c));
                }
            }
            Matrix sSum = S.get(i).sumRows().scale(-1.0);
            for (int r = 0; r < S.get(i).getNumRows(); r++) {
                tser.set(smk[i] + r, 0, sSum.get(r, 0));
            }
        }

        Matrix LM = Matrix.zeros(mtot, mtot);
        for (int i = 0; i < K; i++) {
            Matrix selector = new Matrix(1, mser);
            for (int j = 0; j < alpha.get(i).getNumCols(); j++) {
                selector.set(0, smk[i] + j, alpha.get(i).get(0, j));
            }
            Matrix tserSelector = tser.mult(selector);
            Matrix kronPart = D.get(i).kron(tserSelector);
            LM.addEq(1.0, kronPart);
        }

        Matrix eyeTser = Matrix.eye(m).kron(Tser);
        Matrix D0eye = D0.kron(Matrix.eye(mser));
        Matrix Told = Matrix.zeros(mtot, mtot);
        Matrix Tnew = eyeTser.copy();
        Matrix L;

        if ("Direct".equals(options.mode)) {
            Matrix eyeD0eye = Matrix.eye(mtot).kron(D0eye);
            while (Matrix.infNorm(Told.sub(Tnew)) > 1e-10) {
                Told = Tnew.copy();
                Matrix kronTI = Tnew.transpose().kron(Matrix.eye(mtot));
                Matrix system = kronTI.add(eyeD0eye);
                Matrix negEyeVec = new Matrix(mtot * mtot, 1);
                for (int i = 0; i < mtot; i++) negEyeVec.set(i * mtot + i, 0, -1.0);
                Matrix Lvec = new Matrix(mtot * mtot, 1);
                Matrix.solve(system, negEyeVec, Lvec);
                L = new Matrix(mtot, mtot);
                for (int i = 0; i < mtot; i++) {
                    for (int j = 0; j < mtot; j++) L.set(i, j, Lvec.get(j * mtot + i, 0));
                }
                Tnew = eyeTser.add(L.mult(LM));
            }
        } else {
            java.util.Map<String, ComplexMatrix> schur = Q_Sylvest.schurDecomposition(D0);
            ComplexMatrix Ukron = Q_Sylvest.kronEye(schur.get("U"), mser);
            ComplexMatrix Trkron = Q_Sylvest.kronEye(schur.get("T"), mser);
            while (Matrix.infNorm(Told.sub(Tnew)) > 1e-10) {
                Told = Tnew.copy();
                L = Q_Sylvest.qSylvest(Ukron, Trkron, Tnew);
                Tnew = eyeTser.add(L.mult(LM));
            }
        }

        java.util.Map<String, ComplexMatrix> schur = Q_Sylvest.schurDecomposition(D0);
        L = Q_Sylvest.qSylvest(Q_Sylvest.kronEye(schur.get("U"), mser),
                Q_Sylvest.kronEye(schur.get("T"), mser), Tnew);

        Matrix thetaTot = new Matrix(1, mtot);
        for (int i = 0; i < K; i++) {
            Matrix betaVec = new Matrix(1, mser);
            for (int j = 0; j < beta.get(i).getNumCols(); j++) {
                betaVec.set(0, smk[i] + j, beta.get(i).get(0, j) / mus[i]);
            }
            Matrix piDi = piMat.mult(D.get(i));
            Matrix contrib = piDi.kron(betaVec);
            thetaTot.addEq(1.0, contrib);
        }
        thetaTot.scaleEq(1.0 / load);

        List<Integer> nonzIndices = new ArrayList<Integer>();
        for (int i = 0; i < mtot; i++) {
            if (thetaTot.get(0, i) > 0) nonzIndices.add(i);
        }
        if (nonzIndices.isEmpty()) {
            return new MMAPKPHK1Result(
                    Collections.<Matrix>emptyList(),
                    new Matrix(1, 1),
                    Collections.<Matrix>emptyList(),
                    Collections.<Matrix>emptyList(),
                    null);
        }

        int nz = nonzIndices.size();
        Matrix thetaTotRed = new Matrix(1, nz);
        for (int i = 0; i < nz; i++) thetaTotRed.set(0, i, thetaTot.get(0, nonzIndices.get(i)));
        Matrix TnewReduced = new Matrix(nz, nz);
        for (int i = 0; i < nz; i++) {
            for (int j = 0; j < nz; j++) TnewReduced.set(i, j, Tnew.get(nonzIndices.get(i), nonzIndices.get(j)));
        }

        Matrix diagThetaRed = Matrix.diag(thetaTotRed.getRow(0).toArray1D());
        Matrix diagThetaRedInv = diagThetaRed.inv();
        Matrix Smat = diagThetaRedInv.mult(TnewReduced.transpose()).mult(diagThetaRed);

        List<Matrix> sojAlpha = new ArrayList<Matrix>();
        Matrix diagThetaTot = Matrix.diag(thetaTot.getRow(0).toArray1D());
        Matrix onesM = Matrix.ones(m, 1);

        for (int i = 0; i < K; i++) {
            Matrix tempTser = Matrix.zeros(mser, 1);
            for (int r = 0; r < S.get(i).getNumRows(); r++) {
                tempTser.set(smk[i] + r, 0, tser.get(smk[i] + r, 0));
            }
            Matrix sojFull = diagThetaTot.mult(onesM.kron(tempTser)).scale(load / lambdas[i]);
            Matrix sojReduced = new Matrix(1, nz);
            for (int j = 0; j < nz; j++) sojReduced.set(0, j, sojFull.get(nonzIndices.get(j), 0));
            sojAlpha.add(sojReduced);
        }
        Matrix sojOverallFull = diagThetaTot.mult(onesM.kron(tser)).scale(load / lambda);
        Matrix sojOverall = new Matrix(1, nz);
        for (int j = 0; j < nz; j++) sojOverall.set(0, j, sojOverallFull.get(nonzIndices.get(j), 0));
        sojAlpha.add(sojOverall);

        List<Matrix> waitAlpha = new ArrayList<Matrix>();
        for (int i = 0; i < K; i++) {
            Matrix DiSumCol = D.get(i).sumRows();
            Matrix waitFull = diagThetaTot.mult(L).mult(DiSumCol.kron(tser)).scale(load / lambdas[i]);
            Matrix waitReduced = new Matrix(1, nz);
            for (int j = 0; j < nz; j++) waitReduced.set(0, j, waitFull.get(nonzIndices.get(j), 0));
            waitAlpha.add(waitReduced);
        }
        Matrix DsumMinusD0SumCol = Dsum.sub(D0).sumRows();
        Matrix waitOverallFull = diagThetaTot.mult(L).mult(DsumMinusD0SumCol.kron(tser)).scale(load / lambda);
        Matrix waitOverall = new Matrix(1, nz);
        for (int j = 0; j < nz; j++) waitOverall.set(0, j, waitOverallFull.get(nonzIndices.get(j), 0));
        waitAlpha.add(waitOverall);

        List<Double> qlTList = new ArrayList<Double>();
        qlTList.add(1 - load);

        Matrix Cn = Tnew.copy();
        Matrix Ln = Matrix.zeros(mtot, mtot);
        Matrix DsumEye = Dsum.sub(D0).kron(Matrix.eye(mser));
        java.util.Map<String, ComplexMatrix> schurTotal = Q_Sylvest.schurDecomposition(D0);
        ComplexMatrix Utot = Q_Sylvest.kronEye(schurTotal.get("U"), mser);
        ComplexMatrix Trtot = Q_Sylvest.kronEye(schurTotal.get("T"), mser);

        int n = 1;
        double sumQl = 1 - load;
        while (sumQl < 1 - 1e-10 && n < 1 + options.maxNumComp) {
            Ln = Q_Sylvest.qSylvest(Utot, Trtot, Cn);
            Cn = Ln.mult(DsumEye);
            double qlVal = -load * thetaTot.mult(Ln.sumRows()).get(0, 0);
            qlTList.add(qlVal);
            sumQl += qlVal;
            n++;
        }

        Matrix qlTotal = new Matrix(1, qlTList.size());
        for (int i = 0; i < qlTList.size(); i++) qlTotal.set(0, i, qlTList.get(i));

        List<Matrix> qlPerType = new ArrayList<Matrix>();
        for (int t = 0; t < K; t++) {
            double scaleFactor = lambdas[t] / lambda;
            Matrix qlType = new Matrix(1, qlTotal.getNumCols());
            for (int i = 0; i < qlTotal.getNumCols(); i++) qlType.set(0, i, qlTotal.get(0, i) * scaleFactor);
            qlPerType.add(qlType);
        }

        return new MMAPKPHK1Result(qlPerType, qlTotal, sojAlpha, waitAlpha, Smat);
    }
}
