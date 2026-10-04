/**
 * @file The active set of the immediate transitions and the equations pinning their flows.
 *
 * Port of {@code matlab/src/solvers/FLD/fluid_petri_immediate.m}.
 *
 * @since LINE 3.0
 */
package jline.solvers.fluid.petri;

import java.util.ArrayList;
import java.util.Arrays;
import java.util.List;
import java.util.TreeSet;

/**
 * The active set of the IMMEDIATE transitions, and the equations that pin their flows.
 *
 * <p>An immediate transition has no rate: its fluid limit is a FLOW, an algebraic unknown pinned by
 * the constraint that its binding input place holds no mass,
 *
 * <pre>    phi_j &gt;= 0,   x_b = 0 for the coordinate b that binds mode j</pre>
 *
 * with the GSPN conflict rule supplying the extra equation when two modes drain one place:
 * {@code phi_j*weight_l = phi_l*weight_j} among the enabled modes of highest firing priority, and
 * {@code phi = 0} below it.
 *
 * <p>THE COUNT IS SQUARE BY CONSTRUCTION: V pins plus (F-V) ratio rows is F equations for F flows.
 */
public final class PetriImmediate {

    /** Kinds of pinning equation. */
    public static final int PIN = 0;
    public static final int RATIO = 1;
    public static final int ZERO = 2;

    /** One pinning equation. */
    public static final class Row {
        public final int kind;
        /** PIN: the pinned coordinate. RATIO/ZERO: an immediate-mode index. */
        public final int a;
        /** RATIO: the second mode. Unused otherwise. */
        public final int b;
        public final double wa;
        public final double wb;

        public Row(int kind, int a, int b, double wa, double wb) {
            this.kind = kind;
            this.a = a;
            this.b = b;
            this.wa = wa;
            this.wb = wb;
        }
    }

    /** Number of immediate modes. */
    public int n;
    /** Which of them are currently firing. */
    public boolean[] active;
    /** The input coordinate each active mode is bound by; -1 when inactive. */
    public int[] bind;
    /** The distinct bound coordinates, ascending. */
    public int[] pins;
    /** The equations. */
    public List<Row> rows = new ArrayList<Row>();

    /**
     * True when an inhibitor arc of this mode has reached its threshold.
     *
     * <p>A HARD TEST ON THE MEAN, AND A KNOWN WRONG ANSWER WHEN THE MEAN SITS ON THE THRESHOLD -- a
     * timed mode closes the same indicator smoothly, as {@code Phi((thr-m)/sd)} in
     * {@link PetriSystem#theta}; this path does not. See {@code _kb/06-solver-catalog.md} for the
     * measurement and for what a real fix costs.
     */
    private static boolean inhibited(PetriMode md, double[] x) {
        for (int b = 0; b < md.inhSlot.size(); b++) {
            if (x[md.inhSlot.get(b)] >= md.inhThr.get(b)) {
                return true;
            }
        }
        return false;
    }

    /**
     * Assign bindings and build the pinning equations at the marking {@code x}.
     *
     * @param t   the assembled terms
     * @param x   the current state
     * @param imm the previous immediate state, or null to start one
     * @return the immediate state, mutated in place when one was given
     */
    public static PetriImmediate immediate(PetriTerms t, double[] x, PetriImmediate imm) {
        final int n = t.immIdx.size();
        if (imm == null) {
            imm = new PetriImmediate();
            imm.n = n;
            imm.active = new boolean[n];
            Arrays.fill(imm.active, true);
            imm.bind = new int[n];
            Arrays.fill(imm.bind, -1);
            imm.pins = new int[0];
            // An inhibited mode never fires, so it neither carries a flow nor
            // empties a place. An EMPTY INPUT PLACE IS NOT A REASON TO
            // DEACTIVATE -- that is the normal state of an enabled immediate
            // mode, and the whole content of its pin.
            for (int k = 0; k < n; k++) {
                PetriMode md = t.modes.get(t.immIdx.get(k));
                if (md.arcSlot.isEmpty()) {
                    throw new RuntimeException("Immediate mode " + md.label
                            + " has no enabling arc, so nothing bounds its firing flow and the net "
                            + "has no fluid limit. Give it an input place, or make it timed.");
                }
                if (inhibited(md, x)) {
                    imm.active[k] = false;
                }
            }
        }

        // ---- the assignment: each active mode binds the input arc it is shortest of
        for (int k = 0; k < n; k++) {
            if (!imm.active[k]) {
                imm.bind[k] = -1;
                continue;
            }
            PetriMode md = t.modes.get(t.immIdx.get(k));
            if (imm.bind[k] >= 0 && md.arcSlot.contains(Integer.valueOf(imm.bind[k]))) {
                continue;  // a binding the caller set explicitly is kept
            }
            int best = md.arcSlot.get(0);
            double bestLev = x[md.arcSlot.get(0)] / md.arcW.get(0);
            for (int a = 1; a < md.arcSlot.size(); a++) {
                double lev = x[md.arcSlot.get(a)] / md.arcW.get(a);
                if (lev < bestLev) {
                    bestLev = lev;
                    best = md.arcSlot.get(a);
                }
            }
            imm.bind[k] = best;
        }

        // ---- the equations
        TreeSet<Integer> pinset = new TreeSet<Integer>();
        for (int k = 0; k < n; k++) {
            if (imm.active[k] && imm.bind[k] >= 0) {
                pinset.add(imm.bind[k]);
            }
        }
        imm.pins = new int[pinset.size()];
        int at = 0;
        for (int p : pinset) {
            imm.pins[at++] = p;
        }
        List<Row> rows = new ArrayList<Row>();
        for (int p : imm.pins) {
            rows.add(new Row(PIN, p, 0, 0.0, 0.0));
            List<Integer> grp = new ArrayList<Integer>();
            for (int k = 0; k < n; k++) {
                if (imm.active[k] && imm.bind[k] == p) {
                    grp.add(k);
                }
            }
            if (grp.size() <= 1) {
                continue;
            }
            int topPrio = Integer.MIN_VALUE;
            for (int k : grp) {
                topPrio = Math.max(topPrio, t.modes.get(t.immIdx.get(k)).prio);
            }
            List<Integer> top = new ArrayList<Integer>();
            List<Integer> low = new ArrayList<Integer>();
            for (int k : grp) {
                if (t.modes.get(t.immIdx.get(k)).prio == topPrio) {
                    top.add(k);
                } else {
                    low.add(k);
                }
            }
            double w0 = t.modes.get(t.immIdx.get(top.get(0))).weight;
            for (int i = 1; i < top.size(); i++) {
                rows.add(new Row(RATIO, top.get(0), top.get(i), w0,
                        t.modes.get(t.immIdx.get(top.get(i))).weight));
            }
            for (int k : low) {
                rows.add(new Row(ZERO, k, 0, 0.0, 0.0));
            }
        }
        for (int k = 0; k < n; k++) {
            if (!imm.active[k]) {
                rows.add(new Row(ZERO, k, 0, 0.0, 0.0));
            }
        }
        imm.rows = rows;
        return imm;
    }
}
