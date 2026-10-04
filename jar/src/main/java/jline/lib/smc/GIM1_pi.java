/**
 * @file GI/M/1-type Stationary Distribution
 *
 * Computes the stationary distribution for GI/M/1-type Markov chains.
 *
 * Based on the SMC Solver implementation by Benny Van Houdt.
 *
 * @since LINE 3.1.0
 */
package jline.lib.smc;

import java.util.ArrayList;
import java.util.List;

import jline.io.InputOutput;
import jline.util.matrix.Matrix;

public final class GIM1_pi {
    private GIM1_pi() {}

    /**
     * Computes the stationary distribution of a GI/M/1-type Markov chain.
     */
    public static Matrix gim1_pi(Matrix B, Matrix R, GIM1PiOptions options) {
        int m = R.getNumRows();

        // Check spectral radius of R
        Matrix ImR = Matrix.eye(m).sub(R);
        Matrix ImRinv;
        try {
            ImRinv = ImR.inv();
        } catch (Exception e) {
            throw new IllegalArgumentException("The spectral radius of R is not below 1: QBD is not positive recurrent");
        }

        // Verify (I-R)^{-1} is nonnegative
        for (int i = 0; i < m; i++) {
            for (int j = 0; j < m; j++) {
                if (ImRinv.get(i, j) < -100 * 1e-15) {
                    throw new IllegalArgumentException("The spectral radius of R is not below 1: QBD is not positive recurrent");
                }
            }
        }

        if (options.boundary == null) {
            // Standard boundary case
            return computePiStandard(B, R, ImRinv, m, options);
        } else {
            // General boundary case
            return computePiGeneralBoundary(B, R, ImRinv, m, options);
        }
    }

    public static Matrix gim1_pi(Matrix B, Matrix R) {
        return gim1_pi(B, R, new GIM1PiOptions());
    }

    /**
     * Compute pi for standard boundary case.
     */
    private static Matrix computePiStandard(Matrix B, Matrix R, Matrix ImRinv, int m, GIM1PiOptions options) {
        int maxb = B.getNumRows() / m;

        // Compute BR = B1 + R*B2 + R^2*B3 + ... + R^(maxb-1)*B_maxb
        // Using Horner's method: BR = R*(...R*(R*B_maxb + B_{maxb-1}) + ... + B_2) + B_1
        Matrix BR = B.extractRows((maxb - 1) * m, maxb * m);
        for (int i = maxb - 1; i >= 1; i--) {
            BR = R.mult(BR).add(B.extractRows((i - 1) * m, i * m));
        }

        // Compute pi_0 as stationary distribution of BR
        Matrix pi0 = Stat.stat(BR);

        // Normalize: pi_0 = pi_0 / (pi_0 * (I-R)^{-1} * e)
        Matrix e = Matrix.ones(m, 1);
        double normFactor = pi0.mult(ImRinv).mult(e).get(0, 0);
        pi0 = pi0.scale(1.0 / normFactor);

        // Compute higher levels using pi_i = pi_{i-1} * R
        List<Matrix> piLevels = new ArrayList<Matrix>();
        piLevels.add(pi0);

        double sumPi = pi0.elementSum();
        int numit = 1;

        while (sumPi < 1 - 1e-10 && numit < 1 + options.maxNumComp) {
            Matrix piNext = piLevels.get(piLevels.size() - 1).mult(R);
            piLevels.add(piNext);
            numit++;
            sumPi += piNext.elementSum();

            if (options.verbose > 0 && numit % options.verbose == 0) {
                System.out.println("Accumulated mass after " + numit + " iterations: " + sumPi);
            }
        }

        if (numit == 1 + options.maxNumComp) {
            InputOutput.line_warning("GIM1_pi", "Maximum Number of Components " + (numit - 1) + " reached");
        }

        // Reshape to row vector
        int totalCols = piLevels.size() * m;
        Matrix result = new Matrix(1, totalCols);
        for (int i = 0; i < piLevels.size(); i++) {
            for (int j = 0; j < m; j++) {
                result.set(0, i * m + j, piLevels.get(i).get(0, j));
            }
        }

        return result;
    }

    /**
     * Compute pi for general boundary case.
     */
    private static Matrix computePiGeneralBoundary(Matrix B, Matrix R, Matrix ImRinv, int m, GIM1PiOptions options) {
        Matrix boundary = options.boundary;
        int mb = B.getNumCols(); // number of states of boundary level
        int maxbm1 = (B.getNumRows() - mb) / m; // maxb - 1

        // Compute BR1 from B blocks (excluding first block)
        Matrix BR1 = B.extractRows(mb + (maxbm1 - 1) * m, mb + maxbm1 * m);
        for (int i = maxbm1 - 1; i >= 1; i--) {
            BR1 = R.mult(BR1).add(B.extractRows(mb + (i - 1) * m, mb + i * m));
        }

        // Compute BR0 from boundary blocks
        int maxa = (boundary.getNumRows() - mb) / m;
        Matrix BR0 = boundary.extractRows(mb + (maxa - 1) * m, mb + maxa * m);
        for (int i = maxa - 1; i >= 1; i--) {
            BR0 = R.mult(BR0).add(boundary.extractRows(mb + (i - 1) * m, mb + i * m));
        }

        // Construct combined matrix and compute pi0, pi1
        Matrix B0 = B.extractRows(0, mb);
        Matrix BoundaryTop = boundary.extractRows(0, mb);
        Matrix combinedTop = Matrix.concatColumns(B0, BoundaryTop, null);
        Matrix combinedBottom = Matrix.concatColumns(BR0, BR1, null);
        Matrix combined = Matrix.concatRows(combinedTop, combinedBottom, null);

        Matrix pi01 = Stat.stat(combined);

        // Normalize: pi0 * e + pi1 * (I-R)^{-1} * e = 1
        Matrix e = Matrix.ones(m, 1);
        Matrix emb = Matrix.ones(mb, 1);
        Matrix pi0part = pi01.extractCols(0, mb);
        Matrix pi1part = pi01.extractCols(mb, mb + m);
        double normFactor = pi0part.mult(emb).get(0, 0) + pi1part.mult(ImRinv).mult(e).get(0, 0);
        pi01 = pi01.scale(1.0 / normFactor);

        // Extract normalized pi0 and pi1
        Matrix pi0 = pi01.extractCols(0, mb);
        Matrix piCurrent = pi01.extractCols(mb, mb + m);

        // Compute higher levels
        List<Matrix> piLevels = new ArrayList<Matrix>();
        piLevels.add(piCurrent);

        double sumPi = pi0.elementSum() + piCurrent.elementSum();
        int numit = 1;

        while (sumPi < 1 - 1e-10 && numit < options.maxNumComp) {
            Matrix piNext = piLevels.get(piLevels.size() - 1).mult(R);
            piLevels.add(piNext);
            numit++;
            sumPi += piNext.elementSum();

            if (options.verbose > 0 && numit % options.verbose == 0) {
                System.out.println("Accumulated mass after " + numit + " iterations: " + sumPi);
            }
        }

        if (numit == options.maxNumComp) {
            InputOutput.line_warning("GIM1_pi", "Maximum Number of Components " + numit + " reached");
        }

        // Reshape to row vector: [pi0 pi1 pi2 ...]
        int totalCols = mb + piLevels.size() * m;
        Matrix result = new Matrix(1, totalCols);

        // Copy pi0
        for (int j = 0; j < mb; j++) {
            result.set(0, j, pi0.get(0, j));
        }

        // Copy pi levels
        for (int i = 0; i < piLevels.size(); i++) {
            for (int j = 0; j < m; j++) {
                result.set(0, mb + i * m + j, piLevels.get(i).get(0, j));
            }
        }

        return result;
    }
}
