/**
 * @file Q_RAP_RAP_1 - RAP/RAP/1 Queue Analyzer
 *
 * Computes queue length distribution for a RAP/RAP/1/FCFS queue.
 *
 * Based on the Q-MAM library by Benny Van Houdt.
 *
 * @since LINE 3.1.0
 */
package jline.lib.qmam;

import java.util.Map;
import java.util.Objects;

import jline.lib.smc.QBD_CR;
import jline.lib.smc.QBD_FI;
import jline.lib.smc.QBD_pi;
import jline.lib.smc.Stat;
import jline.util.matrix.Matrix;

public final class Q_RAP_RAP_1 {
    private Q_RAP_RAP_1() {}

    /**
     * Result of RAP/RAP/1 queue analysis
     */
    public static final class RAPRAP1Result {
        public final Matrix queueLength;

        public RAPRAP1Result(Matrix queueLength) {
            this.queueLength = queueLength;
        }

        @Override
        public boolean equals(Object o) {
            if (this == o) return true;
            if (!(o instanceof RAPRAP1Result)) return false;
            RAPRAP1Result that = (RAPRAP1Result) o;
            return Objects.equals(queueLength, that.queueLength);
        }

        @Override
        public int hashCode() {
            return Objects.hash(queueLength);
        }
    }

    /**
     * Options for RAP/RAP/1 queue analysis
     */
    public static final class RAPRAP1Options {
        public final String mode;
        public final int maxNumComp;
        public final int verbose;

        public RAPRAP1Options() {
            this("CR", 50, 0);
        }

        public RAPRAP1Options(String mode, int maxNumComp, int verbose) {
            this.mode = mode;
            this.maxNumComp = maxNumComp;
            this.verbose = verbose;
        }
    }

    public static RAPRAP1Result qRapRap1(Matrix C0, Matrix C1, Matrix D0, Matrix D1) {
        return qRapRap1(C0, C1, D0, D1, new RAPRAP1Options());
    }

    /**
     * Computes queue length distribution for a RAP/RAP/1/FCFS queue.
     */
    public static RAPRAP1Result qRapRap1(Matrix C0, Matrix C1, Matrix D0, Matrix D1, RAPRAP1Options options) {
        int mA = C0.getNumRows();
        int mS = D0.getNumRows();
        int mtot = mA * mS;

        // Validate dimensions
        if (!(C0.getNumCols() == mA && C1.getNumRows() == mA && C1.getNumCols() == mA)) {
            throw new IllegalArgumentException("Arrival process matrices must be mA x mA");
        }
        if (!(D0.getNumCols() == mS && D1.getNumRows() == mS && D1.getNumCols() == mS)) {
            throw new IllegalArgumentException("Service process matrices must be mS x mS");
        }

        // Test the load of the queue
        Matrix eyeA = Matrix.eye(mA);
        Matrix eyeS = Matrix.eye(mS);

        Matrix sumA = C0.add(C1).add(eyeA);
        Matrix sumS = D0.add(D1).add(eyeS);

        Matrix piA = Stat.stat(sumA);
        Matrix piS = Stat.stat(sumS);

        double lambda = piA.mult(C1).elementSum();
        double mu = piS.mult(D1).elementSum();

        double load = lambda / mu;
        if (!(load < 1)) {
            throw new IllegalArgumentException("The load " + load + " of the system exceeds one");
        }

        // Compute QBD blocks A0, A1, A2
        Matrix A0 = eyeA.kron(D1);
        Matrix A1 = C0.kron(eyeS).add(eyeA.kron(D0));
        Matrix A2 = C1.kron(eyeS);

        Matrix B0 = A0.copy();
        Matrix B1 = C0.kron(eyeS);

        // Compute G and R using appropriate solver with RAPComp=1
        Map<String, Matrix> qbdResult;
        Integer verboseArg = options.verbose > 0 ? Integer.valueOf(1) : null;
        if (options.mode.contains("FI")) {
            qbdResult = QBD_FI.QBD_FI(A0, A1, A2, null, verboseArg, null, null, 1);
        } else {
            qbdResult = QBD_CR.QBD_CR(A0, A1, A2, null, verboseArg, null, 1);
        }

        Matrix R = qbdResult.get("R");

        // Compute stationary distribution
        Matrix pi = QBD_pi.QBD_pi(B0, B1, R, options.maxNumComp, options.verbose, null, 1);

        // Compute queue length distribution
        int numLevels = pi.getNumCols() / mtot;
        Matrix ql = new Matrix(1, numLevels);
        for (int i = 0; i < numLevels; i++) {
            double sum = 0.0;
            for (int j = 0; j < mtot; j++) {
                sum += pi.get(0, i * mtot + j);
            }
            ql.set(0, i, sum);
        }

        return new RAPRAP1Result(ql);
    }
}
