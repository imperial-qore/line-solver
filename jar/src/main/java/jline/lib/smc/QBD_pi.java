package jline.lib.smc;

import java.util.ArrayList;
import java.util.List;

import jline.util.matrix.Matrix;

public final class QBD_pi {
    private QBD_pi() {}

    /**
     * QBD_pi: Stationary vector of a Quasi-Birth-Death Markov Chains [Neuts]
     */
    public static Matrix QBD_pi(Matrix B0, Matrix B1, Matrix R) {
        return QBD_pi(B0, B1, R, 500, 0, null, 0);
    }

    public static Matrix QBD_pi(Matrix B0, Matrix B1, Matrix R, int MaxNumComp, int Verbose, Matrix Boundary, int RAPComp) {
        Matrix B1_work = B1.copy();
        Matrix B0_work = B0.copy();
        Matrix Boundary_work = Boundary == null ? null : Boundary.copy();

        int m = R.getNumRows();

        // Convert to discrete time problem, if needed
        Matrix B1_diag = new Matrix(m, 1, 0);
        Matrix.extractDiag(B1_work, B1_diag);

        if (B1_diag.elementMin() < 0 || RAPComp == 1) { // continuous time
            double lamb = -B1_diag.elementMin();

            // Simplified normalization for continuous time case
            B1_work = B1_work.scale(1.0 / lamb).add(1.0, Matrix.eye(m));

            B0_work.scaleEq(1.0 / lamb);
        }

        Matrix temp = Matrix.eye(m).add(-1.0, R).inv();
        // Positive recurrence is sp(R) < 1. For a Markovian QBD, R is
        // nonnegative and sp(R) < 1 is equivalent to (I-R)^-1 >= 0, which is
        // the test QBD_pi.m performs; with RAP components R may have negative
        // entries, (I-R)^-1 then legitimately has negative entries too, and
        // reading a single negative entry as non-recurrence rejects stable
        // RAP/RAP/1 queues. The eigenvalue test is the condition itself and
        // agrees with the entrywise one whenever R is nonnegative.
        double spR = 0.0;
        for (org.apache.commons.math3.complex.Complex ev : R.eig()) {
            spR = Math.max(spR, ev.abs());
        }
        if (spR >= 1.0) {
            throw new RuntimeException("The spectral radius of R is not below 1: QBD is not pos. recurrent");
        }

        // Simplified implementation - in practice this would be more complex
        Matrix pi0 = Stat.stat(B1_work.add(1.0, R.mult(B0_work)));
        double normalizer = pi0.mult(temp).mult(Matrix.ones(m, 1)).get(0);
        pi0.scaleEq(1.0 / normalizer);

        List<Matrix> pi_components = new ArrayList<Matrix>();
        pi_components.add(pi0.copy());

        double sumpi = pi0.elementSum();
        int numit = 1;

        while (sumpi < 1 - 1e-10 && numit < MaxNumComp) {
            Matrix pi_next = pi_components.get(pi_components.size() - 1).mult(R);
            pi_components.add(pi_next);
            numit++;
            sumpi += pi_next.elementSum();

            if (Verbose > 0 && numit % Verbose == 0) {
                System.out.println("Accumulated mass after " + numit + " iterations: " + sumpi);
            }
        }

        if (numit == MaxNumComp) {
            System.out.println("Maximum Number of Components " + numit + " reached");
        }

        // Concatenate all components horizontally
        Matrix result = pi_components.get(0).copy();
        for (int i = 1; i < pi_components.size(); i++) {
            result = Matrix.concatColumns(result, pi_components.get(i), null);
        }

        return result;
    }
}
