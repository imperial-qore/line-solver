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

            if (Boundary_work != null) {
                // The boundary carries [B2; A1+R*A0]: its lower block is a
                // generator too, so the uniformization constant has to cover
                // it as well, and only that lower block takes the +I shift.
                int mbCt = B1_work.getNumRows();
                Matrix lower = Matrix.extract(Boundary_work, mbCt, Boundary_work.getNumRows(),
                        0, Boundary_work.getNumCols());
                Matrix lowerDiag = new Matrix(m, 1, 0);
                Matrix.extractDiag(lower, lowerDiag);
                lamb = Math.max(lamb, -lowerDiag.elementMin());
                Boundary_work = Boundary_work.scale(1.0 / lamb);
                for (int i = 0; i < m; i++) {
                    Boundary_work.set(mbCt + i, i, Boundary_work.get(mbCt + i, i) + 1.0);
                }
                B1_work = B1_work.scale(1.0 / lamb).add(1.0, Matrix.eye(mbCt));
            } else {
                B1_work = B1_work.scale(1.0 / lamb).add(1.0, Matrix.eye(m));
            }

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

        List<Matrix> pi_components = new ArrayList<Matrix>();
        double sumpi;
        int numit = 1;
        Matrix level0 = null;

        if (Boundary_work == null) {
            Matrix pi0 = Stat.stat(B1_work.add(1.0, R.mult(B0_work)));
            double normalizer = pi0.mult(temp).mult(Matrix.ones(m, 1)).get(0);
            pi0.scaleEq(1.0 / normalizer);
            pi_components.add(pi0.copy());
            sumpi = pi0.elementSum();
        } else {
            // General boundary, QBD_pi.m else-branch: level 0 and level 1 are
            // solved together on [[B1; B0] Boundary], because B0 and B2 need
            // not be square and the level-0 block then differs in size from
            // the repeating one. Before this branch existed the Boundary
            // argument was accepted and DISCARDED, so every general-boundary
            // caller silently got the default-boundary answer.
            int mb = B1_work.getNumRows();
            Matrix leftCol = Matrix.concatRows(B1_work, B0_work, null);
            Matrix joint = Matrix.concatColumns(leftCol, Boundary_work, null);
            Matrix pi01 = Stat.stat(joint);

            Matrix pi0 = Matrix.extract(pi01, 0, 1, 0, mb);
            Matrix pi1 = Matrix.extract(pi01, 0, 1, mb, mb + m);
            double normalizer = pi0.elementSum()
                    + pi1.mult(temp).mult(Matrix.ones(m, 1)).get(0);
            pi0.scaleEq(1.0 / normalizer);
            pi1.scaleEq(1.0 / normalizer);

            level0 = pi0;
            pi_components.add(pi1.copy());
            sumpi = pi0.elementSum() + pi1.elementSum();
        }

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

        // Concatenate all components horizontally, level 0 first
        Matrix result = level0 == null ? pi_components.get(0).copy() : level0.copy();
        int first = level0 == null ? 1 : 0;
        for (int i = first; i < pi_components.size(); i++) {
            result = Matrix.concatColumns(result, pi_components.get(i), null);
        }

        return result;
    }
}
