package jline.api.mam;

import java.util.ArrayList;
import java.util.Arrays;
import java.util.Comparator;
import java.util.List;
import java.util.Map;

import jline.util.matrix.Matrix;
import jline.util.matrix.MatrixCell;
import jline.GlobalConstants;

public final class Mmap_super_safe {
    private Mmap_super_safe() {}

    /**
     * Safely combines multiple MMAPs into a single superposed MMAP while considering order constraints.
     *
     * @param MMAPS    a map of MMAPs to be combined
     * @param maxorder the maximum allowed order for the resulting superposed MMAP
     * @param method   the method for combining MMAPs; "default" or "match"
     * @return a MatrixCell representing the combined superposed MMAP
     */
    public static MatrixCell mmap_super_safe(Map<Integer, MatrixCell> MMAPS, int maxorder, String method) {
        // Handle empty hashmap case
        if (MMAPS.isEmpty()) {
            Matrix lambda = new Matrix(1, 1);
            lambda.set(0, 0, 1.0);
            return Mmap_exponential.mmap_exponential(lambda);
        }

        // A component with an all-zero aggregate arrival matrix (D1 == 0) has zero
        // arrival rate: it contributes nothing to the superposition. When such a
        // component carries more than one phase (e.g. the transient slow phase of a
        // high-SCV APH fit that reaches a station the flow never visits), its phase
        // generator D0+D1 is absorbing, so map_scv/map_pie/map_prob -> ctmc_solve
        // fail with "no recurrent state". Canonicalize it to the equivalent order-1
        // null (matching the marking count), whose moments are well defined and
        // whose superposition is an identity. Detection uses the arrival-matrix norm,
        // not mmap_lambda, because mmap_lambda itself calls map_prob -> ctmc_solve.
        // Mirrors matlab/lib/m3a/m3a/mmap/mmap_super_safe.m.
        for (int i = 0; i < MMAPS.size(); i++) {
            MatrixCell flow = MMAPS.get(i);
            if (flow.get(0).getNumRows() > 1 && flow.get(1).norm() < 1e-13) {
                int Kmarks = flow.size() - 2;
                MMAPS.put(i, Mmap_exponential.mmap_exponential(new Matrix(1, Kmarks), 1));
            }
        }

        MatrixCell sup = new MatrixCell();
        List<Double> scv_unmarked = new ArrayList<Double>();
        for (int i = 0; i < MMAPS.size(); i++) {
            scv_unmarked.add(Map_scv.map_scv(MMAPS.get(i).get(0), MMAPS.get(i).get(1)));
        }

        Integer[] indices = new Integer[scv_unmarked.size()];
        for (int i = 0; i < indices.length; i++) {
            indices[i] = i;
        }
        final List<Double> scvFinal = scv_unmarked;
        // Sorting flows by SCV is a numerical heuristic (superpose the smoothest flow
        // first). It must NOT decide the mark order of the result: mmap_super
        // concatenates the marks of its operands in superposition order and callers
        // read back mark k as class k, so a fold order other than caller order renames
        // the classes. That is repaired unconditionally below, by mark provenance.
        // The tolerance here remains so that meaningless floating-point noise does not
        // reshuffle the Kronecker factors either: two Poisson streams both have SCV 1,
        // yet an aggregate accumulated through super/scale carries a few ulp of error
        // (e.g. 1.0000000000000004) and would sort after an exactly-1.0 flow. MATLAB's
        // sort() is likewise stable and, getting an exact 1.0, keeps the same order.
        // See matlab/lib/m3a/m3a/mmap/mmap_super_safe.m.
        Arrays.sort(indices, new Comparator<Integer>() {
            @Override
            public int compare(Integer a, Integer b) {
                double sa = scvFinal.get(a);
                double sb = scvFinal.get(b);
                double scale = Math.max(1.0, Math.max(Math.abs(sa), Math.abs(sb)));
                if (Math.abs(sa - sb) <= GlobalConstants.FineTol * scale) {
                    return 0; // indistinguishable: keep caller order (stable sort)
                }
                return Double.compare(sa, sb);
            }
        });
        int[] sortedIndices = new int[indices.length];
        for (int i = 0; i < indices.length; i++) {
            sortedIndices[i] = indices[i].intValue();
        }

        // Mark provenance: a zero-rate component has SCV NaN and sorts last, so a
        // chain that never visits the station used to push its marks ahead of one
        // that does, renaming both chains' classes.
        int[] markbase = new int[MMAPS.size() + 1];
        for (int i = 0; i < MMAPS.size(); i++) {
            markbase[i + 1] = markbase[i] + (MMAPS.get(i).size() - 2);
        }
        List<Integer> outorder = new ArrayList<Integer>();

        for (int i = 0; i < sortedIndices.length; i++) {
            int smallest_value = sortedIndices[i];
            MatrixCell flow = MMAPS.get(smallest_value);
            for (int j = markbase[smallest_value]; j < markbase[smallest_value + 1]; j++) {
                outorder.add(Integer.valueOf(j));
            }
            // Bound the order of each individual flow to maxorder. A single flow
            // whose order already exceeds maxorder (e.g. a high-order Erlang from a
            // near-deterministic APH fit) would otherwise pass through uncapped as
            // the superposition base and blow up downstream matrix-analytic solves.
            if (flow.get(0).length() > maxorder) {
                if (maxorder >= 2) {
                    flow = Mamap2m_fit_gamma_fb_mmap.mamap2m_fit_gamma_fb_mmap(flow);
                } else {
                    flow = Mmap_exponential.mmap_exponential(Mmap_lambda.mmap_lambda(flow));
                }
            }
            if (sup.isEmpty()) {
                sup = flow;
                if (maxorder == 1) {
                    sup = Mmap_exponential.mmap_exponential(Mmap_lambda.mmap_lambda(flow));
                }
            } else {
                if (sup.get(0).length() * flow.get(0).length() > maxorder) {
                    if (sup.get(0).length() * 2 < maxorder) {
                        sup = Mmap_super.mmap_super(sup, Mamap2m_fit_gamma_fb_mmap.mamap2m_fit_gamma_fb_mmap(flow), method);
                    } else {
                        sup = Mmap_super.mmap_super(sup, Mmap_exponential.mmap_exponential(Mmap_lambda.mmap_lambda(flow)), method);
                    }
                } else {
                    sup = Mmap_super.mmap_super(sup, flow, method);
                }
            }
        }

        // Restore the caller's mark order; "match" keeps one mark per class and
        // fails the count test, so it is left alone.
        if (sup.size() - 2 == outorder.size()) {
            boolean sorted = true;
            for (int j = 1; j < outorder.size(); j++) {
                if (outorder.get(j).intValue() < outorder.get(j - 1).intValue()) {
                    sorted = false;
                    break;
                }
            }
            if (!sorted) {
                MatrixCell out = new MatrixCell();
                out.set(0, sup.get(0));
                out.set(1, sup.get(1));
                for (int j = 0; j < outorder.size(); j++) {
                    out.set(2 + outorder.get(j).intValue(), sup.get(2 + j));
                }
                sup = out;
            }
        }

        return sup;
    }

    public static MatrixCell mmap_super_safe(Map<Integer, MatrixCell> MMAPS, int maxorder) {
        return mmap_super_safe(MMAPS, maxorder, "default");
    }
}
