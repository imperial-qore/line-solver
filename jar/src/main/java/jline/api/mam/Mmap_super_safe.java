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
        // first), but it also fixes the MARK order of the result: mmap_super
        // concatenates the marks of its operands in superposition order, and callers
        // read back mark k as class k. Ordering by a raw Double.compare therefore lets
        // meaningless floating-point noise permute the class-to-mark mapping: two
        // Poisson streams both have SCV 1, yet an aggregate accumulated through
        // super/scale carries a few ulp of error (e.g. 1.0000000000000004), which
        // sorts it AFTER an exactly-1.0 flow and silently swaps two classes' results.
        // Treat SCVs that agree to within tolerance as equal so the (stable) sort keeps
        // them in caller order, which is class order. MATLAB's sort() is likewise
        // stable and, getting an exact 1.0, keeps the same order; this makes the
        // agreement robust rather than accidental.
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

        for (int i = 0; i < sortedIndices.length; i++) {
            int smallest_value = sortedIndices[i];
            MatrixCell flow = MMAPS.get(smallest_value);
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

        return sup;
    }

    public static MatrixCell mmap_super_safe(Map<Integer, MatrixCell> MMAPS, int maxorder) {
        return mmap_super_safe(MMAPS, maxorder, "default");
    }
}
