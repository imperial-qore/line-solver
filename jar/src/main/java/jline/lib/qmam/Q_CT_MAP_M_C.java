package jline.lib.qmam;

import java.util.ArrayList;
import java.util.List;
import java.util.Map;

import jline.lib.smc.QBD_CR;
import jline.lib.smc.QBD_FI;
import jline.lib.smc.Stat;
import jline.util.matrix.ComplexMatrix;
import jline.util.matrix.Matrix;

/**
 * Q_CT_MAP_M_C - Continuous-Time MAP/M/c Queue Analyzer.
 * Computes queue length and waiting time distribution for a continuous-time MAP/M/c/FCFS queue.
 */
public final class Q_CT_MAP_M_C {
    private Q_CT_MAP_M_C() {}

    public static MAPMcResult qCtMapMC(Matrix D0, Matrix D1, double mu, int c) {
        return qCtMapMC(D0, D1, mu, c, new MAPMcOptions());
    }

    public static MAPMcResult qCtMapMC(Matrix D0, Matrix D1, double mu, int c, MAPMcOptions options) {
        int m = D0.getNumRows();
        if (D0.getNumCols() != m || D1.getNumRows() != m || D1.getNumCols() != m) {
            throw new IllegalArgumentException("D0 and D1 must be " + m + "x" + m + " matrices");
        }
        if (mu <= 1e-14) throw new IllegalArgumentException("Service rate mu must be strictly positive");
        if (c < 1) throw new IllegalArgumentException("Number of servers c must be at least 1");

        Matrix invD0 = D0.scale(-1.0).inv();
        Matrix piA = Stat.stat(D1.mult(invD0));
        double lambda = piA.mult(D1).elementSum();
        double load = lambda / (mu * c);
        if (load >= 1.0) throw new IllegalArgumentException("The load " + load + " of the system exceeds one");

        Matrix eyeM = Matrix.eye(m);
        Matrix A0 = eyeM.scale(c * mu);
        Matrix A1 = D0.sub(eyeM.scale(c * mu));
        Matrix A2 = D1.copy();

        Map<String, Matrix> qbdResult;
        if (options.mode.contains("FI")) {
            qbdResult = QBD_FI.solve(A0, A1, A2, null, (options.verbose > 0) ? Integer.valueOf(1) : null, null, null, null);
        } else {
            qbdResult = QBD_CR.solve(A0, A1, A2, null, (options.verbose > 0) ? Integer.valueOf(1) : null, null, null);
        }
        Matrix R = qbdResult.get("R");

        Matrix piGJL = new Matrix(1, c * m);
        if (c > 1) {
            List<Matrix> invC = new ArrayList<Matrix>();
            invC.add(D0.scale(-1.0).inv());
            for (int i = 2; i < c; i++) {
                Matrix innerMatrix = D0.sub(eyeM.scale((i - 1) * mu))
                        .add(invC.get(i - 2).mult(D1).scale((i - 1) * mu))
                        .scale(-1.0);
                invC.add(innerMatrix.inv());
            }
            Matrix boundaryMatrix = D0.sub(eyeM.scale((c - 1) * mu))
                    .add(R.mult(A0))
                    .add(invC.get(c - 2).mult(D1).scale((c - 1) * mu))
                    .add(eyeM);
            Matrix piCm1 = Stat.stat(boundaryMatrix);
            for (int j = 0; j < m; j++) piGJL.set(0, (c - 1) * m + j, piCm1.get(0, j));

            for (int i = c - 1; i >= 1; i--) {
                Matrix piNext = new Matrix(1, m);
                for (int j = 0; j < m; j++) piNext.set(0, j, piGJL.get(0, i * m + j));
                Matrix piPrev = piNext.mult(invC.get(i - 1)).scale(i * mu);
                for (int j = 0; j < m; j++) piGJL.set(0, (i - 1) * m + j, piPrev.get(0, j));
            }
        } else {
            Matrix boundaryMatrix = D0.add(R.mult(A0)).add(eyeM);
            Matrix pi0 = Stat.stat(boundaryMatrix);
            for (int j = 0; j < m; j++) piGJL.set(0, j, pi0.get(0, j));
        }

        Matrix ImR = eyeM.sub(R);
        Matrix ImRinv = ImR.inv();

        double K = 0.0;
        for (int i = 0; i < (c - 1) * m; i++) K += piGJL.get(0, i);
        Matrix piCm1 = new Matrix(1, m);
        for (int j = 0; j < m; j++) piCm1.set(0, j, piGJL.get(0, (c - 1) * m + j));
        K += piCm1.mult(ImRinv).mult(Matrix.ones(m, 1)).get(0, 0);
        piGJL.scaleEq(1.0 / K);

        Matrix piC1 = new Matrix(1, m);
        for (int j = 0; j < m; j++) piC1.set(0, j, piGJL.get(0, (c - 1) * m + j));
        List<Matrix> piLevels = new ArrayList<Matrix>();
        piLevels.add(piC1.copy());
        double sumpi = piGJL.elementSum();
        int numit = 1;
        while (sumpi < 1 - 1e-10 && numit < 1 + options.maxNumComp - c) {
            Matrix piNext = piLevels.get(piLevels.size() - 1).mult(R);
            piLevels.add(piNext);
            numit++;
            sumpi += piNext.elementSum();
            if (options.verbose > 0 && numit % options.verbose == 0) {
                System.out.println("Accumulated mass after " + numit + " iterations: " + sumpi);
            }
        }

        Matrix ql = new Matrix(1, c - 1 + piLevels.size());
        for (int i = 0; i < c - 1; i++) {
            double levelSum = 0.0;
            for (int j = 0; j < m; j++) levelSum += piGJL.get(0, i * m + j);
            ql.set(0, i, levelSum);
        }
        for (int i = 0; i < piLevels.size(); i++) ql.set(0, c - 1 + i, piLevels.get(i).elementSum());

        double probZero = 0.0;
        for (int i = 0; i < c * m; i++) probZero += piGJL.get(0, i);
        Matrix D1sumCol = D1.sumRows();

        double numerator = 0.0;
        for (int i = 0; i < c; i++) {
            Matrix piLevel = new Matrix(1, m);
            for (int j = 0; j < m; j++) piLevel.set(0, j, piGJL.get(0, i * m + j));
            numerator += piLevel.mult(D1sumCol).get(0, 0);
        }

        int piTsize = (c - 1) * m + piLevels.size() * m;
        Matrix piT = new Matrix(1, piTsize);
        for (int i = 0; i < (c - 1) * m; i++) piT.set(0, i, piGJL.get(0, i));
        for (int i = 0; i < piLevels.size(); i++) {
            for (int j = 0; j < m; j++) piT.set(0, (c - 1) * m + i * m + j, piLevels.get(i).get(0, j));
        }
        double denominator = 0.0;
        int numPiTLevels = piTsize / m;
        for (int i = 0; i < numPiTLevels; i++) {
            Matrix piLevel = new Matrix(1, m);
            for (int j = 0; j < m; j++) piLevel.set(0, j, piT.get(0, i * m + j));
            denominator += piLevel.mult(D1sumCol).get(0, 0);
        }
        probZero = numerator / denominator;

        Matrix piCm1ForAlpha = new Matrix(1, m);
        for (int j = 0; j < m; j++) piCm1ForAlpha.set(0, j, piGJL.get(0, (c - 1) * m + j));
        Matrix temp = piCm1ForAlpha.mult(ImRinv).mult(D1);
        Matrix alphaVec = temp.scale(1.0 / temp.elementSum());

        Matrix Told = Matrix.zeros(m, m);
        Matrix Tnew = A0.scale(-1.0).copy();
        Matrix L;
        int TIterCount = 0;
        int maxTIterations = 1000;

        if (options.mode.contains("Direct")) {
            Matrix eyeD0eye = eyeM.kron(D0);
            double tNorm = Matrix.infNorm(Told.sub(Tnew));
            while (tNorm > 1e-10 && TIterCount < maxTIterations) {
                Told = Tnew.copy();
                Matrix kronTI = Tnew.transpose().kron(eyeM);
                Matrix system = kronTI.add(eyeD0eye);
                Matrix negEyeVec = new Matrix(m * m, 1);
                for (int i = 0; i < m; i++) negEyeVec.set(i * m + i, 0, -1.0);
                Matrix Lvec = new Matrix(m * m, 1);
                Matrix.solve(system, negEyeVec, Lvec);
                L = new Matrix(m, m);
                for (int i = 0; i < m; i++) {
                    for (int j = 0; j < m; j++) L.set(i, j, Lvec.get(j * m + i, 0));
                }
                Tnew = A0.scale(-1.0).add(L.mult(D1).scale(mu * c));
                TIterCount++;
                tNorm = Matrix.infNorm(Told.sub(Tnew));
            }
        } else {
            java.util.Map<String, ComplexMatrix> schur = Q_Sylvest.schurDecomposition(D0);
            ComplexMatrix U = schur.get("U");
            ComplexMatrix Tr = schur.get("T");
            double tNorm = Matrix.infNorm(Told.sub(Tnew));
            while (tNorm > 1e-10 && TIterCount < maxTIterations) {
                Told = Tnew.copy();
                L = Q_Sylvest.qSylvest(U, Tr, Tnew);
                Tnew = A0.scale(-1.0).add(L.mult(D1).scale(mu * c));
                TIterCount++;
                tNorm = Matrix.infNorm(Told.sub(Tnew));
            }
            if (TIterCount >= maxTIterations && tNorm > 1e-10) {
                if (m == 1) {
                    double a = -A0.get(0, 0);
                    double b = D1.get(0, 0) * (mu * c);
                    double d = D0.get(0, 0);
                    double p = d - a;
                    double q = -(a * d - b);
                    double discriminant = p * p + 4 * q;
                    if (discriminant >= 0) {
                        double sqrtD = Math.sqrt(discriminant);
                        double x1 = (-p + sqrtD) / 2;
                        double x2 = (-p - sqrtD) / 2;
                        Tnew.set(0, 0, (x2 < 0) ? x2 : x1);
                    }
                }
            }
        }

        List<Integer> nonzIndices = new ArrayList<Integer>();
        for (int i = 0; i < m; i++) {
            if (alphaVec.get(0, i) > 0) nonzIndices.add(i);
        }
        if (nonzIndices.isEmpty()) return new MAPMcResult(ql, null, null);

        int nz = nonzIndices.size();
        Matrix thetaTotRed = new Matrix(1, nz);
        for (int i = 0; i < nz; i++) thetaTotRed.set(0, i, alphaVec.get(0, nonzIndices.get(i)));
        Matrix TnewReduced = new Matrix(nz, nz);
        for (int i = 0; i < nz; i++) {
            for (int j = 0; j < nz; j++) TnewReduced.set(i, j, Tnew.get(nonzIndices.get(i), nonzIndices.get(j)));
        }

        double[] thetaArr = thetaTotRed.getRow(0).toArray1D();
        Matrix diagThetaRed = Matrix.diag(thetaArr);
        Matrix diagThetaRedInv = diagThetaRed.inv();
        Matrix Smat = diagThetaRedInv.mult(TnewReduced.transpose()).mult(diagThetaRed);

        Matrix rhoVec = Tnew.add(A0).sumRows();
        double alphaRhoSum = 0.0;
        for (int i = 0; i < m; i++) alphaRhoSum += alphaVec.get(0, i) * rhoVec.get(i, 0);

        Matrix waitAlphaFull = new Matrix(1, m);
        for (int i = 0; i < m; i++) {
            waitAlphaFull.set(0, i, (1 - probZero) * alphaVec.get(0, i) * rhoVec.get(i, 0) / alphaRhoSum);
        }
        Matrix waitAlpha = new Matrix(1, nz);
        for (int i = 0; i < nz; i++) waitAlpha.set(0, i, waitAlphaFull.get(0, nonzIndices.get(i)));

        return new MAPMcResult(ql, waitAlpha, Smat);
    }
}
