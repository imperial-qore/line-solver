/**
 * @file The conserved quantities of a stochastic Petri net's fluid limit.
 *
 * Port of {@code matlab/src/solvers/FLD/fluid_petri_conservation.m}.
 *
 * @since LINE 3.0
 */
package jline.solvers.fluid.petri;

import java.util.ArrayList;
import java.util.List;

import jline.util.matrix.Matrix;

/**
 * The conserved quantities, as equations: {@code u'D = 0} implies {@code u'x} is constant.
 *
 * <p>On the marking coordinates those u are the net's P-invariants; on a mode's phase block the
 * all-ones vector is one of them, which is the statement that the phase coordinates are a
 * distribution. Both come out of the SAME null space, so the phase normalisation needs no separate
 * row.
 *
 * <p>AN OPEN NET LOSES THE ROWS ITS ARRIVALS BREAK, automatically: the arrival columns are part of
 * D, so a u an arrival moves is not in the null space and never appears here.
 */
public final class PetriConservation {

    /** (nrows x nstate) conservation rows. */
    public Matrix C;
    /** The value each row holds, evaluated at the initial marking. */
    public double[] N;
    /** {@code max|C D|}: how far the rows fail to be invariants, and zero when they are. */
    public double leak;
    /** One human-readable statement per row, for the report. */
    public List<String> label = new ArrayList<String>();

    private PetriConservation() {}

    /**
     * Build the conservation rows of a term set.
     *
     * @param t the assembled terms
     * @return the rows, their values and their residual leak
     */
    public static PetriConservation conservation(PetriTerms t) {
        PetriConservation cons = new PetriConservation();
        if (t.D == null || t.D.isEmpty() || t.D.getNumCols() == 0) {
            cons.C = new Matrix(0, t.nstate);
            cons.N = new double[0];
            cons.leak = 0.0;
            return cons;
        }

        // null(D', 'r'), transposed: one row per invariant.
        Matrix Z = nullRational(t.D.transpose());
        Matrix C = Z.transpose();
        // A rational basis has exact zeros; anything below this is a rounding
        // artefact of the elimination, not a coefficient.
        for (int i = 0; i < C.getNumRows(); i++) {
            for (int j = 0; j < C.getNumCols(); j++) {
                if (Math.abs(C.get(i, j)) < 1e-12) {
                    C.set(i, j, 0.0);
                }
            }
        }
        // Scale each row by its smallest nonzero magnitude, so a row reads as a
        // statement about named places with small integer weights.
        List<Integer> keep = new ArrayList<Integer>();
        for (int i = 0; i < C.getNumRows(); i++) {
            double smallest = Double.POSITIVE_INFINITY;
            boolean any = false;
            for (int j = 0; j < C.getNumCols(); j++) {
                double v = Math.abs(C.get(i, j));
                if (v != 0.0) {
                    any = true;
                    smallest = Math.min(smallest, v);
                }
            }
            if (!any) {
                continue;
            }
            for (int j = 0; j < C.getNumCols(); j++) {
                C.set(i, j, C.get(i, j) / smallest);
            }
            keep.add(i);
        }

        Matrix Ck = new Matrix(keep.size(), t.nstate);
        for (int r = 0; r < keep.size(); r++) {
            for (int j = 0; j < t.nstate; j++) {
                Ck.set(r, j, C.get(keep.get(r), j));
            }
        }
        cons.C = Ck;
        cons.N = new double[Ck.getNumRows()];
        for (int r = 0; r < Ck.getNumRows(); r++) {
            double v = 0.0;
            for (int j = 0; j < t.nstate; j++) {
                v += Ck.get(r, j) * t.x0[j];
            }
            cons.N[r] = v;
        }
        cons.leak = 0.0;
        for (int r = 0; r < Ck.getNumRows(); r++) {
            for (int c = 0; c < t.D.getNumCols(); c++) {
                double v = 0.0;
                for (int j = 0; j < t.nstate; j++) {
                    v += Ck.get(r, j) * t.D.get(j, c);
                }
                cons.leak = Math.max(cons.leak, Math.abs(v));
            }
        }
        for (int r = 0; r < Ck.getNumRows(); r++) {
            cons.label.add(rowLabel(t, Ck, r));
        }
        return cons;
    }

    /**
     * A rational basis of the null space, as MATLAB's {@code null(A,'r')} returns.
     *
     * <p>D is INTEGRAL -- arc multiplicities and unit phase moves -- so the reduced row echelon
     * form gives exact rational rows, which keeps each conservation row readable as a statement
     * about named places instead of an arbitrary orthogonal mixture. An SVD basis would span the
     * same space and say nothing.
     */
    private static Matrix nullRational(Matrix A) {
        int rows = A.getNumRows();
        int cols = A.getNumCols();
        if (rows == 0 || cols == 0) {
            return new Matrix(cols, 0);
        }
        Matrix R = A.copy();
        List<Integer> piv = new ArrayList<Integer>();
        int r = 0;
        final double tol = 1e-12;
        for (int c = 0; c < cols && r < rows; c++) {
            int k = r;
            double best = Math.abs(R.get(r, c));
            for (int i = r + 1; i < rows; i++) {
                if (Math.abs(R.get(i, c)) > best) {
                    best = Math.abs(R.get(i, c));
                    k = i;
                }
            }
            if (best <= tol) {
                for (int i = r; i < rows; i++) {
                    R.set(i, c, 0.0);
                }
                continue;
            }
            if (k != r) {
                for (int j = 0; j < cols; j++) {
                    double tmp = R.get(r, j);
                    R.set(r, j, R.get(k, j));
                    R.set(k, j, tmp);
                }
            }
            double p = R.get(r, c);
            for (int j = 0; j < cols; j++) {
                R.set(r, j, R.get(r, j) / p);
            }
            for (int i = 0; i < rows; i++) {
                if (i == r) {
                    continue;
                }
                double f = R.get(i, c);
                if (f == 0.0) {
                    continue;
                }
                for (int j = 0; j < cols; j++) {
                    R.set(i, j, R.get(i, j) - f * R.get(r, j));
                }
            }
            piv.add(c);
            r++;
        }
        List<Integer> free = new ArrayList<Integer>();
        for (int c = 0; c < cols; c++) {
            if (!piv.contains(c)) {
                free.add(c);
            }
        }
        Matrix Z = new Matrix(cols, free.size());
        for (int i = 0; i < free.size(); i++) {
            int f = free.get(i);
            Z.set(f, i, 1.0);
            for (int rr = 0; rr < piv.size(); rr++) {
                Z.set(piv.get(rr), i, -R.get(rr, f));
            }
        }
        return Z;
    }

    /** One conservation row, written out over the coordinates it touches. */
    private static String rowLabel(PetriTerms t, Matrix C, int r) {
        StringBuilder sb = new StringBuilder();
        for (int s = 0; s < C.getNumCols(); s++) {
            double w = C.get(r, s);
            if (w == 0.0) {
                continue;
            }
            String nm;
            if (s < t.nm) {
                nm = t.namesNode.get(t.coordNode[s]) + "(class " + (t.coordClass[s] + 1) + ")";
            } else {
                nm = "phase";
                for (PetriMode md : t.modes) {
                    if (md.zblk == null) {
                        continue;
                    }
                    for (int q = 0; q < md.zblk.length; q++) {
                        if (md.zblk[q] == s) {
                            nm = md.label + " phase " + (q + 1);
                            break;
                        }
                    }
                }
            }
            if (sb.length() > 0) {
                sb.append(" + ");
            }
            sb.append(w == 1.0 ? nm : (trim(w) + "*" + nm));
        }
        return sb.toString();
    }

    /** %g formatting, so an integral weight prints without a decimal tail. */
    private static String trim(double v) {
        if (v == Math.rint(v) && Math.abs(v) < 1e15) {
            return Long.toString((long) v);
        }
        return Double.toString(v);
    }
}
