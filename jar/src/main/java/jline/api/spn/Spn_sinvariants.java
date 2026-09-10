package jline.api.spn;

import java.util.ArrayList;
import java.util.List;

import jline.lang.NetworkStruct;
import jline.lang.NodeParam;
import jline.lang.constant.NodeType;
import jline.lang.nodeparam.TransitionNodeParam;
import jline.lang.nodes.Node;
import jline.util.matrix.Matrix;

/**
 * Minimal-support S-invariants (P-invariants) of a stochastic Petri net, and the
 * load vector V = S m0.
 *
 * <p>An S-invariant is a non-negative left null vector of the incidence matrix,
 * U' C = 0, so U' m is conserved by every firing. The minimal-support ones form
 * a basis of all of them (S. Balsamo, A. Marin, I. Stojic, FGCS 111 (2020),
 * Sec. 3.1) and are what the convolution algorithm {@code Spn_conv} decomposes
 * the reachability set along; {@code Spn_mdd} uses a single positive invariant
 * for a much weaker purpose, to bound each place a priori.</p>
 *
 * <p>FARKAS' ALGORITHM, on [C | I]: for each transition column in turn, keep the
 * rows that already annihilate it and add, for every pair of rows of opposite
 * sign in it, the positive combination that cancels it; then drop every row
 * whose support strictly contains another's, which is what leaves the minimal
 * supports. Rows are kept in integer arithmetic and divided by their gcd, so a
 * multiplicity is never lost to rounding and two invariants that differ only by
 * a positive scale are the same row.</p>
 *
 * <p>ARC MULTIPLICITIES MUST BE INTEGRAL. A fractional arc has no Petri-net
 * meaning and would make the gcd normalisation and the ILP-free convolution both
 * wrong, so it is refused rather than rounded.</p>
 *
 * <p>Levels are the (place, class) pairs of {@code Spn_mdd}, place-major, so the
 * invariants come out in the coordinates the decision diagram and
 * {@code Spn_conv} both use.</p>
 *
 * <p>MATLAB twin: {@code spn_sinvariants.m}. Python twin:
 * {@code api/spn/sinvariants.py}.</p>
 */
public class Spn_sinvariants {

    private Spn_sinvariants() {}

    /** The invariant basis of a net, in place-level coordinates. */
    public static class SpnInvariants {
        /** Node indices of the places, in level order. */
        public int[] places;
        /** S[i][p]: weight of place level p in minimal-support invariant i. */
        public long[][] S;
        /** V = S m0, the load vector. */
        public long[] V;
        /** The initial marking the load vector was taken against. */
        public long[] m0;
    }

    private static long gcd(long a, long b) {
        a = a < 0 ? -a : a;
        b = b < 0 ? -b : b;
        while (b != 0) {
            long t = a % b;
            a = b;
            b = t;
        }
        return a;
    }

    /** An arc multiplicity, refused unless integral. */
    private static long asInteger(double x, String what) {
        double r = Math.floor(x + 0.5);
        if (Math.abs(x - r) > 1e-9) {
            throw new RuntimeException("spn_sinvariants: " + what + " is not integral; a fractional "
                    + "arc multiplicity has no Petri-net meaning and no invariant basis over the "
                    + "integers");
        }
        return (long) r;
    }

    /** True when the support of a is contained in the support of b. */
    private static boolean supportSubset(long[] a, long[] b, int off, int n) {
        for (int k = 0; k < n; k++) {
            if (a[off + k] != 0 && b[off + k] == 0) {
                return false;
            }
        }
        return true;
    }

    private static boolean supportEqual(long[] a, long[] b, int off, int n) {
        for (int k = 0; k < n; k++) {
            if ((a[off + k] != 0) != (b[off + k] != 0)) {
                return false;
            }
        }
        return true;
    }

    /** Minimal-support S-invariants with the marking taken from the model. */
    public static SpnInvariants spn_sinvariants(NetworkStruct sn) {
        return spn_sinvariants(sn, null);
    }

    /**
     * Minimal-support S-invariants and the load vector of a net.
     *
     * @param sn a NetworkStruct holding Places and Transitions
     * @param init initial tokens per place level, place-major; null takes them
     *        from the reference station of each closed class, as
     *        {@code Spn_mdd} does
     * @return the invariant basis and the load vector
     */
    public static SpnInvariants spn_sinvariants(NetworkStruct sn, double[] init) {
        List<Integer> places = new ArrayList<Integer>();
        List<Integer> transitions = new ArrayList<Integer>();
        for (int i = 0; i < sn.nodetype.size(); i++) {
            if (sn.nodetype.get(i) == NodeType.Place) {
                places.add(Integer.valueOf(i));
            } else if (sn.nodetype.get(i) == NodeType.Transition) {
                transitions.add(Integer.valueOf(i));
            }
        }
        if (places.isEmpty() || transitions.isEmpty()) {
            throw new RuntimeException(
                    "spn_sinvariants: the model holds no Place or no Transition node");
        }
        int R = sn.nclasses;
        int P = places.size();
        int n = P * R;

        // ---- incidence C[level][mode] = post - pre, one column per (transition, mode)
        List<long[]> cols = new ArrayList<long[]>();
        for (int t = 0; t < transitions.size(); t++) {
            int ind = transitions.get(t).intValue();
            Node node = sn.nodes.get(ind);
            NodeParam param = sn.nodeparam.get(node);
            if (!(param instanceof TransitionNodeParam)) {
                continue;
            }
            TransitionNodeParam tp = (TransitionNodeParam) param;
            for (int m = 0; m < tp.nmodes; m++) {
                long[] col = new long[n];
                for (int pp = 0; pp < P; pp++) {
                    for (int k = 0; k < R; k++) {
                        double pre = arc(tp.enabling, m, places.get(pp).intValue(), k);
                        double post = arc(tp.firing, m, places.get(pp).intValue(), k);
                        col[pp * R + k] = asInteger(post, "a firing arc")
                                - asInteger(pre, "an enabling arc");
                    }
                }
                cols.add(col);
            }
        }
        int ncols = cols.size();

        // ---- Farkas on [C | I]: row p starts as (C[p], e_p)
        List<long[]> rows = new ArrayList<long[]>();
        for (int p = 0; p < n; p++) {
            long[] row = new long[ncols + n];
            for (int c = 0; c < ncols; c++) {
                row[c] = cols.get(c)[p];
            }
            row[ncols + p] = 1;
            rows.add(row);
        }
        for (int c = 0; c < ncols; c++) {
            List<long[]> next = new ArrayList<long[]>();
            for (int r = 0; r < rows.size(); r++) {
                if (rows.get(r)[c] == 0) {
                    next.add(rows.get(r));
                }
            }
            for (int a = 0; a < rows.size(); a++) {
                if (rows.get(a)[c] <= 0) {
                    continue;
                }
                for (int b = 0; b < rows.size(); b++) {
                    if (rows.get(b)[c] >= 0) {
                        continue;
                    }
                    long pa = rows.get(a)[c];
                    long nb = -rows.get(b)[c];
                    long d = gcd(pa, nb);
                    long fa = nb / d;
                    long fb = pa / d;
                    long[] combo = new long[ncols + n];
                    long gall = 0;
                    for (int k = 0; k < combo.length; k++) {
                        combo[k] = fa * rows.get(a)[k] + fb * rows.get(b)[k];
                        gall = gcd(gall, combo[k]);
                    }
                    if (gall > 1) {
                        for (int k = 0; k < combo.length; k++) {
                            combo[k] /= gall;
                        }
                    }
                    boolean nonzero = false;
                    for (int k = 0; k < n; k++) {
                        nonzero = nonzero || combo[ncols + k] != 0;
                    }
                    if (nonzero) {
                        next.add(combo);
                    }
                }
            }
            // support-minimality filter, applied at every step so the row set
            // cannot grow combinatorially on the way to the answer
            List<long[]> keep = new ArrayList<long[]>();
            for (int r = 0; r < next.size(); r++) {
                boolean dominated = false;
                for (int s = 0; s < next.size() && !dominated; s++) {
                    if (s == r) {
                        continue;
                    }
                    if (!supportSubset(next.get(s), next.get(r), ncols, n)) {
                        continue;
                    }
                    boolean same = supportEqual(next.get(s), next.get(r), ncols, n);
                    if (!same || s < r) {          // keep the first of equal supports
                        dominated = true;
                    }
                }
                if (!dominated) {
                    keep.add(next.get(r));
                }
            }
            rows = keep;
        }

        SpnInvariants out = new SpnInvariants();
        out.places = new int[P];
        for (int pp = 0; pp < P; pp++) {
            out.places[pp] = places.get(pp).intValue();
        }
        List<long[]> basis = new ArrayList<long[]>();
        for (int r = 0; r < rows.size(); r++) {
            boolean nonneg = true;
            for (int k = 0; k < n; k++) {
                nonneg = nonneg && rows.get(r)[ncols + k] >= 0;
            }
            if (!nonneg) {                          // an S-invariant is non-negative by definition
                continue;
            }
            long[] y = new long[n];
            for (int k = 0; k < n; k++) {
                y[k] = rows.get(r)[ncols + k];
            }
            basis.add(y);
        }
        out.S = new long[basis.size()][];
        for (int i = 0; i < basis.size(); i++) {
            out.S[i] = basis.get(i);
        }

        // ---- initial marking and the load vector V = S m0
        out.m0 = new long[n];
        if (init != null) {
            if (init.length != n) {
                throw new RuntimeException(
                        "spn_sinvariants: init must hold one token count per place level");
            }
            for (int i = 0; i < n; i++) {
                out.m0[i] = asInteger(init[i], "an initial marking");
            }
        } else {
            for (int k = 0; k < R; k++) {
                double njobs = sn.njobs.get(k);
                if (Double.isInfinite(njobs) || Double.isNaN(njobs)) {
                    throw new RuntimeException("spn_sinvariants: class " + (k + 1) + " is open, so "
                            + "the net has no finite load vector");
                }
                int ref = (int) sn.refstat.get(k);
                for (int pp = 0; pp < P; pp++) {
                    if ((int) sn.nodeToStation.get(places.get(pp).intValue()) == ref) {
                        out.m0[pp * R + k] += asInteger(njobs, "a class population");
                    }
                }
            }
        }
        out.V = new long[out.S.length];
        for (int i = 0; i < out.S.length; i++) {
            for (int p = 0; p < n; p++) {
                out.V[i] += out.S[i][p] * out.m0[p];
            }
        }
        return out;
    }

    private static double arc(List<Matrix> mats, int m, int node, int cls) {
        if (mats == null || mats.size() <= m || mats.get(m) == null) {
            return 0.0;
        }
        return Math.max(0.0, mats.get(m).get(node, cls));
    }
}
