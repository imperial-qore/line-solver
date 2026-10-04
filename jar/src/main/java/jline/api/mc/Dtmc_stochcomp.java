package jline.api.mc;

import java.util.ArrayList;
import java.util.Arrays;
import java.util.List;
import java.util.stream.Collectors;
import java.util.stream.IntStream;

import org.apache.commons.math3.util.FastMath;

import jline.util.matrix.Matrix;

public final class Dtmc_stochcomp {
    private Dtmc_stochcomp() {}

    /**
     * Returns the stochastic complement of a DTMC.
     *
     * @param P Transition matrix of the DTMC
     * @param I Indexes of states to be kept in the stochastic complement
     * @return Transition matrix of the stochastic complement
     * @throws RuntimeException If the transition matrix is null
     */
    public static Matrix dtmc_stochcomp(Matrix P, List<Integer> I) {
        // Note that in this function, List is used instead of Matrix for performance consideration

        int lengthP = FastMath.max(P.getNumCols(), P.getNumRows());
        if (I == null || I.size() == 0) {
            I = new ArrayList<Integer>();
            for (int i = 0; i < (int) FastMath.ceil(lengthP / 2.0); i++) {
                I.add(i);
            }
        }

        List<Integer> Ic = IntStream.rangeClosed(0, lengthP - 1).boxed().collect(Collectors.toList());
        Ic.removeAll(I);

        Matrix P11 = new Matrix(I.size(), I.size());
        Matrix P12 = new Matrix(I.size(), Ic.size());
        Matrix P21 = new Matrix(Ic.size(), I.size());
        Matrix P22 = new Matrix(Ic.size(), Ic.size());

        for (int colIdx = 0; colIdx < P.getNumCols(); colIdx++) {
            for (int rowIdx = 0; rowIdx < P.getNumRows(); rowIdx++) {
                double value = P.get(rowIdx, colIdx);
                if (value > 0.0 || Double.isNaN(value)) {
                    if (I.contains(colIdx)) {
                        if (I.contains(rowIdx)) {
                            P11.set(I.indexOf(rowIdx), I.indexOf(colIdx), value);
                        } else {
                            P21.set(Ic.indexOf(rowIdx), I.indexOf(colIdx), value);
                        }
                    } else {
                        if (I.contains(rowIdx)) {
                            P12.set(I.indexOf(rowIdx), Ic.indexOf(colIdx), value);
                        } else {
                            P22.set(Ic.indexOf(rowIdx), Ic.indexOf(colIdx), value);
                        }
                    }
                }
            }
        }

        double[] values = new double[Ic.size()];
        Arrays.fill(values, 1.0);
        Matrix S2 = Matrix.diagMatrix(null, values, 0, values.length).sub(P22);

        // S=P11+P12*(S2 \ P21);
        Matrix s2_p21 = new Matrix(S2.getNumRows(), P21.getNumCols());

        // Check if S2 is empty
        if (S2.getNumRows() == 0 || S2.getNumCols() == 0) {
            return P11;
        }

        // see _kb/03-api-layer.md for rationale
        if (S2.getNumRows() > Ctmc_solve.GMRES_MIN_STATES) {
            Matrix it = Ctmc_gmres.ctmc_gmres(S2, P21, 0.0, 0, 0);
            if (it == null) {
                // Short-recurrence retry before the direct factorization, as in Ctmc_solve.
                it = Ctmc_bicgstab.ctmc_bicgstab(S2, P21, 0.0, 0);
            }
            if (it != null) {
                return P11.add(1.0, P12.mult(it, null));
            }
        }

        // Use solveSafe to handle singular I-P22 (matches MATLAB's mldivide behavior)
        if (!Matrix.solveSafe(S2, P21, s2_p21)) {
            // S2 is singular — regularize by adding small epsilon to diagonal
            double eps = 1e-10;
            Matrix S2reg = S2.copy();
            for (int i = 0; i < S2reg.getNumRows(); i++) {
                S2reg.set(i, i, S2reg.get(i, i) + eps);
            }
            Matrix.solve(S2reg, P21, s2_p21);
        }
        Matrix S = P11.add(1.0, P12.mult(s2_p21, null));
        return S;
    }

    /** Stochastic complement together with the blocks it was built from. */
    public static final class DtmcStochCompResult {
        public final Matrix S;
        public final Matrix P11;
        public final Matrix P12;
        public final Matrix P21;
        public final Matrix P22;

        public DtmcStochCompResult(Matrix S, Matrix P11, Matrix P12, Matrix P21, Matrix P22) {
            this.S = S;
            this.P11 = P11;
            this.P12 = P12;
            this.P21 = P21;
            this.P22 = P22;
        }
    }

    /**
     * Stochastic complement of a DTMC together with the four blocks of the
     * transition matrix partitioned by the kept states and their complement.
     * Twin of the MATLAB [S,P11,P12,P21,P22] = dtmc_stochcomp(P,I) and of the
     * Python dtmc_stochcomp_full.
     *
     * @param P Transition matrix of the DTMC
     * @param I Indexes of states to be kept in the stochastic complement
     * @return the complement and the blocks it was built from
     */
    public static DtmcStochCompResult dtmc_stochcomp_full(Matrix P, List<Integer> I) {
        int lengthP = FastMath.max(P.getNumCols(), P.getNumRows());
        if (I == null || I.size() == 0) {
            I = new ArrayList<Integer>();
            for (int i = 0; i < (int) FastMath.ceil(lengthP / 2.0); i++) {
                I.add(i);
            }
        }
        List<Integer> Ic = IntStream.rangeClosed(0, lengthP - 1).boxed().collect(Collectors.toList());
        Ic.removeAll(I);

        Matrix P11 = new Matrix(I.size(), I.size());
        Matrix P12 = new Matrix(I.size(), Ic.size());
        Matrix P21 = new Matrix(Ic.size(), I.size());
        Matrix P22 = new Matrix(Ic.size(), Ic.size());
        for (int rowIdx = 0; rowIdx < I.size(); rowIdx++) {
            for (int colIdx = 0; colIdx < I.size(); colIdx++) {
                P11.set(rowIdx, colIdx, P.get(I.get(rowIdx), I.get(colIdx)));
            }
            for (int colIdx = 0; colIdx < Ic.size(); colIdx++) {
                P12.set(rowIdx, colIdx, P.get(I.get(rowIdx), Ic.get(colIdx)));
            }
        }
        for (int rowIdx = 0; rowIdx < Ic.size(); rowIdx++) {
            for (int colIdx = 0; colIdx < I.size(); colIdx++) {
                P21.set(rowIdx, colIdx, P.get(Ic.get(rowIdx), I.get(colIdx)));
            }
            for (int colIdx = 0; colIdx < Ic.size(); colIdx++) {
                P22.set(rowIdx, colIdx, P.get(Ic.get(rowIdx), Ic.get(colIdx)));
            }
        }
        return new DtmcStochCompResult(dtmc_stochcomp(P, I), P11, P12, P21, P22);
    }
}
