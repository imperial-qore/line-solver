/**
 * @file Finite place capacities of a stochastic Petri net, as linear rows.
 *
 * Port of {@code matlab/src/solvers/FLD/fluid_petri_constraints.m}.
 *
 * @since LINE 3.0
 */
package jline.solvers.fluid.petri;

import java.util.ArrayList;
import java.util.Arrays;
import java.util.List;

import jline.lang.NetworkStruct;
import jline.util.matrix.Matrix;

/**
 * Every finite place capacity as a linear row {@code A x <= b}.
 *
 * <p>THE GATE IS A LOSS ON THE DEPOSIT: LINE loses the tokens a firing would push past a place's
 * capacity, so the fluid analogue scales the DEPOSIT leg of every event adding mass to the capped
 * place and leaves the removal leg alone. The rows here name which coordinates are capped; the
 * scaling itself lives in the tangent clamp.
 */
public final class PetriConstraints {

    /** (nrows x nstate) capacity rows. */
    public Matrix A;
    /** The bound of each row. */
    public double[] b;
    /** One human-readable statement per row, for the report. */
    public List<String> label = new ArrayList<String>();
    /** {@code cover[r][s]} is true where row r caps coordinate s. */
    public boolean[][] cover;

    private PetriConstraints() {}

    /**
     * Build the capacity rows of a term set.
     *
     * @param sn the network the places belong to
     * @param t  the assembled terms
     * @return the rows and their bounds
     */
    public static PetriConstraints constraints(NetworkStruct sn, PetriTerms t) {
        List<double[]> rows = new ArrayList<double[]>();
        List<Double> bs = new ArrayList<Double>();
        List<String> labels = new ArrayList<String>();

        for (int ind : t.places) {
            int ist = (int) sn.nodeToStation.get(ind);
            List<Integer> slots = new ArrayList<Integer>();
            for (int k = 0; k < t.K; k++) {
                if (t.pidx[ind][k] >= 0) {
                    slots.add(t.pidx[ind][k]);
                }
            }
            if (slots.isEmpty()) {
                continue;
            }
            double cap = (sn.cap != null && ist >= 0 && ist < sn.cap.length())
                    ? sn.cap.get(ist) : Double.POSITIVE_INFINITY;
            if (Double.isFinite(cap)) {
                double[] row = new double[t.nstate];
                for (int s : slots) {
                    row[s] = 1.0;
                }
                rows.add(row);
                bs.add(cap);
                labels.add("capacity " + fmt(cap) + " of place " + t.namesNode.get(ind));
            }
            for (int k = 0; k < t.K; k++) {
                int s = t.pidx[ind][k];
                if (s < 0 || sn.classcap == null || ist < 0 || ist >= sn.classcap.getNumRows()) {
                    continue;
                }
                double ccap = sn.classcap.get(ist, k);
                if (!Double.isFinite(ccap)) {
                    continue;
                }
                double[] row = new double[t.nstate];
                row[s] = 1.0;
                rows.add(row);
                bs.add(ccap);
                labels.add("class-" + (k + 1) + " capacity " + fmt(ccap) + " of place "
                        + t.namesNode.get(ind));
            }
        }

        // TWO ROWS THAT SAY THE SAME THING ARE A SINGULAR NEWTON SYSTEM, not a
        // redundancy the least squares absorbs, so an exact duplicate is pruned
        // and the tighter bound survives.
        boolean[] keep = new boolean[rows.size()];
        Arrays.fill(keep, true);
        for (int c = 0; c < rows.size(); c++) {
            if (!keep[c]) {
                continue;
            }
            for (int d = c + 1; d < rows.size(); d++) {
                if (keep[d] && Arrays.equals(rows.get(c), rows.get(d))) {
                    if (bs.get(d) < bs.get(c)) {
                        keep[c] = false;
                        break;
                    }
                    keep[d] = false;
                }
            }
        }

        List<Integer> idx = new ArrayList<Integer>();
        for (int i = 0; i < keep.length; i++) {
            if (keep[i]) {
                idx.add(i);
            }
        }
        PetriConstraints con = new PetriConstraints();
        con.A = new Matrix(idx.size(), t.nstate);
        con.b = new double[idx.size()];
        con.cover = new boolean[idx.size()][t.nstate];
        for (int r = 0; r < idx.size(); r++) {
            double[] row = rows.get(idx.get(r));
            for (int s = 0; s < t.nstate; s++) {
                con.A.set(r, s, row[s]);
                con.cover[r][s] = row[s] > 0.0;
            }
            con.b[r] = bs.get(idx.get(r));
            con.label.add(labels.get(idx.get(r)));
        }
        return con;
    }

    /** %g formatting, so an integral capacity prints without a decimal tail. */
    private static String fmt(double v) {
        if (v == Math.rint(v) && Math.abs(v) < 1e15) {
            return Long.toString((long) v);
        }
        return Double.toString(v);
    }
}
