package jline.api.mc;

import java.util.ArrayList;
import java.util.List;

import jline.solvers.ctmc.SolverCTMC;
import jline.util.matrix.Matrix;

public final class Ctmc_pseudostochcomp {
    private Ctmc_pseudostochcomp() {}

    public static SolverCTMC.StochCompResult ctmc_pseudostochcomp(Matrix Q, List<Double> I_list) {
        Matrix I = new Matrix(I_list);
        if (I_list.isEmpty()) {
            // Default: use first half of states
            int halfSize = (int) Math.ceil(Q.getNumCols() / 2.0);
            List<Double> defaultIndices = new ArrayList<Double>();
            for (int i = 0; i < halfSize; i++) {
                defaultIndices.add((double) i);
            }
            I = new Matrix(defaultIndices);
        }

        // Compute Ic as complement of I
        List<Double> diff_values = new ArrayList<Double>();
        for (int checkValue = 0; checkValue < Q.getNumCols(); checkValue++) {
            boolean check = false;
            for (int I_row = 0; I_row < I.getNumRows(); I_row++) {
                if (I.get(I_row) == (double) checkValue) {
                    check = true;
                    break;
                }
            }
            if (!check) {
                diff_values.add((double) checkValue);
            }
        }
        Matrix Ic = new Matrix(diff_values);

        // Extract submatrices
        Matrix Q11 = Q.getSubMatrix(I, I);
        Matrix Q12 = Q.getSubMatrix(I, Ic);
        Matrix Q21 = Q.getSubMatrix(Ic, I);
        Matrix Q22 = Q.getSubMatrix(Ic, Ic);

        // Solve for stationary distribution
        Matrix abar = Ctmc_solve.ctmc_solve(Q);

        // Extract abar(Ic) - stationary probabilities for complement states
        List<Double> zeroList = new ArrayList<Double>();
        zeroList.add(0.0);
        Matrix abar_Ic = abar.getSubMatrix(Ic, new Matrix(zeroList));

        // Y = abar(Ic) * Q21
        Matrix Y = abar_Ic.transpose().mult(Q21);

        // Create ones vector of size Q12.columns
        Matrix ones = Matrix.ones(Q12.getNumCols(), 1);

        // S = Q11 + Q12 * ones * Y / sum(Y)
        double sumY = Y.sumCols().get(0);
        Matrix scaledY = Y.scale(1.0 / sumY);
        Matrix Q12_ones = Q12.mult(ones);
        Matrix addTerm = Q12_ones.mult(scaledY);
        Matrix S = Q11.add(1.0, addTerm);

        // For compatibility with existing StochCompResult, we can set T = Q12 * ones * Y / sum(Y)
        Matrix T = addTerm;

        return new SolverCTMC.StochCompResult(S, Q11, Q12, Q21, Q22, T);
    }
}
