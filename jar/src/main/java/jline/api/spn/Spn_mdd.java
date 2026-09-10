package jline.api.spn;

import java.util.ArrayList;
import java.util.Arrays;
import java.util.LinkedHashSet;
import java.util.List;
import java.util.Map;
import java.util.Set;

import jline.api.mdd.MDD;
import jline.api.mdd.MddDescriptor;
import jline.api.mdd.MddEvent;
import jline.api.mdd.MddLocalMatrix;
import jline.api.mdd.MddStruct;
import jline.lang.JobClass;
import jline.lang.Network;
import jline.lang.NetworkStruct;
import jline.lang.NodeParam;
import jline.lang.constant.NodeType;
import jline.lang.constant.TimingStrategy;
import jline.lang.nodeparam.TransitionNodeParam;
import jline.lang.nodes.Node;
import jline.util.matrix.Matrix;
import jline.util.matrix.MatrixCell;

/**
 * Decision-diagram reachable set and Kronecker rate descriptor of a stochastic
 * Petri net, so that {@link jline.api.mdd.Mdd_mcd} can analyse it.
 *
 * <p>Levels are of two kinds. <b>Place levels</b>: one per (Place, class) pair,
 * holding a token count. A multiclass net therefore has P*R of them, ordered
 * place-major, so level p*R+k is class k in place p. <b>Phase levels</b>: one
 * per mode whose firing time has more than one phase, holding the phase the
 * running server occupies.</p>
 *
 * <p>The rate structure factorises exactly under single-server firing
 * semantics: a mode fires at a constant rate whenever every input level holds
 * its enabling multiplicity and no inhibitor level has reached its threshold, so
 * W_l^e[i, i + fire(l) - enab(l)] = 1 for enab(l) &lt;= i &lt; inhib(l) at every
 * place level. A phase-type mode contributes two event families on its phase
 * level, the internal phase changes D0 (marking unchanged) and the firings D1
 * (marking moved), each gated by the same per-level enabling indicators. Both
 * are products of per-level terms, which is what Eq. 1 of the paper
 * requires.</p>
 *
 * <p><b>Phase-type firing and the memory policy.</b> LINE discards a running
 * server's phase when its mode becomes disabled, i.e. preemptive repeat.
 * Resetting a mode's phase is then triggered by a JOINT condition on the place
 * levels, which is not a product of per-level terms and has no Kronecker form.
 * What this descriptor encodes is preemptive resume: a disabled mode's phase
 * freezes and continues when the mode is re-enabled. The two policies coincide
 * exactly when a mode is never disabled while running, so reachability records,
 * for free, whether any phase-type mode was ever found disabled. phmemory
 * "exact" (the default) errors when one was; "resume" proceeds deliberately
 * with the resume semantics.</p>
 *
 * <p><b>Other restrictions</b>, each an error and never a silent approximation:
 * no immediate transitions (they make vanishing states, which must be
 * eliminated before a Kronecker rate descriptor exists) and no
 * marking-dependent firing rates. A multi-server mode is accepted only when its
 * enabling touches ONE level, because the enabling degree
 * min_l floor(m(l)/enab(l)) is otherwise not a product of per-level terms.</p>
 *
 * <p>MATLAB twin: {@code spn_mdd.m}. Python twin: {@code api/spn/mdd.py}.</p>
 */
public class Spn_mdd {

    private Spn_mdd() {}

    /** One (transition, mode) pair of the net, in level coordinates. */
    public static class SpnMode {
        /** Node index of the transition. */
        public int trans;
        /** Mode index within the transition, 0-based. */
        public int mode;
        /** Enabling multiplicity per place level. */
        public double[] enab;
        /** Inhibition threshold per place level; infinite when absent. */
        public double[] inhib;
        /** Firing outcome per place level. */
        public double[] fire;
        public double[][] D0;
        public double[][] D1;
        public double[] pie;
        public int nph;
        public double srv;
        /** Marking-dependent firing-rate multiplier; null for the unit one. */
        public jline.util.SerializableFunction<Matrix, Double> dep;
    }

    /** Everything the caller needs alongside the descriptor. */
    public static class SpnInfo {
        /** Node indices of the places. */
        public int[] places;
        public String[] placenames;
        public String[] classnames;
        /** 1 for a place level, 2 for a phase level. */
        public int[] levelkind;
        public String[] levelname;
        public List<SpnMode> modes;
        public int[] init;
        public MDD mdd;
        /** 1-based phase level of each mode, 0 when the mode has one phase. */
        public int[] phaseof;
        /** Whether each mode was ever found disabled in a reachable marking. */
        public boolean[] everDisabled;
        public int nplacelevels;
        /** Node count of the model, so an arc-shaped marking matrix can be rebuilt. */
        public int nnodes;
        public int nclasses;
        /** Whether the Kronecker descriptor was built. */
        public boolean descriptor;
    }

    /** Descriptor, diagram and metadata returned together. */
    public static class SpnResult {
        public MddStruct mdds;
        public MddDescriptor desc;
        public SpnInfo info;
    }

    /** Options of the translation. */
    public static class SpnOptions {
        /** Per-place-level token bound; null infers it from a place invariant. */
        public double[] bound;
        /** "exact" (default) or "resume". */
        public String phmemory = "exact";
        /**
         * Build the Kronecker rate descriptor (default true). Pass false for the
         * MDD-rec route, which reads only the reachable set: the restrictions that
         * exist purely because a Kronecker form must factorise per level
         * (marking-dependent firing rates, multi-server modes drawing from several
         * places) are then lifted, in exchange for the firing times having to be
         * exponential.
         */
        public boolean descriptor = true;
        public boolean verbose = false;
    }

    /** Translate the net with the default options. */
    public static SpnResult spn_mdd(Network model) {
        return spn_mdd(model, new SpnOptions());
    }

    /**
     * Build the reachable set and Kronecker descriptor of a stochastic Petri net.
     *
     * @param model a Network holding Places and Transitions
     * @param options translation options; null takes the defaults
     */
    public static SpnResult spn_mdd(Network model, SpnOptions options) {
        if (options == null) {
            options = new SpnOptions();
        }
        NetworkStruct sn = model.getStruct();
        int R = sn.nclasses;
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
            throw new RuntimeException("spn_mdd: the model holds no Place or no Transition node");
        }
        int P = places.size();
        int L = P * R;                                 // place levels, place-major

        // ---- collect the (transition, mode) pairs
        List<SpnMode> md = new ArrayList<SpnMode>();
        for (int t = 0; t < transitions.size(); t++) {
            int ind = transitions.get(t).intValue();
            Node node = sn.nodes.get(ind);
            NodeParam param = sn.nodeparam.get(node);
            if (!(param instanceof TransitionNodeParam)) {
                continue;
            }
            TransitionNodeParam tp = (TransitionNodeParam) param;
            for (int m = 0; m < tp.nmodes; m++) {
                if (tp.timing != null && tp.timing.size() > m
                        && tp.timing.get(m) == TimingStrategy.IMMEDIATE) {
                    throw new RuntimeException("spn_mdd: mode " + (m + 1) + " of node " + (ind + 1)
                            + " is IMMEDIATE; vanishing states must be eliminated before the net "
                            + "has a Kronecker rate descriptor");
                }
                jline.util.SerializableFunction<Matrix, Double> dep =
                        (tp.firingdep != null && tp.firingdep.size() > m)
                                ? tp.firingdep.get(m) : null;
                if (dep != null && options.descriptor) {
                    throw new RuntimeException("spn_mdd: mode " + (m + 1) + " of node " + (ind + 1)
                            + " has a marking-dependent firing rate; g(marking) is not a product "
                            + "of per-level terms");
                }
                MatrixCell proc = firingProc(tp, m);
                if (proc == null || proc.size() < 2) {
                    throw new RuntimeException("spn_mdd: mode " + (m + 1) + " of node " + (ind + 1)
                            + " has no Markovian firing process; a general distribution has no "
                            + "finite phase level");
                }
                SpnMode e = new SpnMode();
                e.trans = ind;
                e.mode = m;
                e.enab = arcvec(tp.enabling, m, places, R, 0.0);
                e.inhib = arcvec(tp.inhibiting, m, places, R, Double.POSITIVE_INFINITY);
                e.fire = arcvec(tp.firing, m, places, R, 0.0);
                e.D0 = toArray(proc.get(0));
                e.D1 = toArray(proc.get(1));
                e.nph = e.D0.length;
                e.pie = firingPie(tp, m, e.nph);
                e.srv = (tp.nmodeservers != null && tp.nmodeservers.getNumElements() > m)
                        ? tp.nmodeservers.get(m) : 1.0;
                e.dep = dep;
                if (!options.descriptor && e.nph > 1) {
                    throw new RuntimeException("spn_mdd: mode " + (m + 1) + " of node " + (ind + 1)
                            + " has a phase-type firing time; the reachable-set-only mode carries "
                            + "no phase level, and a product-form marking process must be "
                            + "memoryless in the marking alone");
                }
                md.add(e);
            }
        }
        int E = md.size();
        for (int e = 0; options.descriptor && e < E; e++) {
            int nz = 0;
            for (int l = 0; l < L; l++) {
                if (md.get(e).enab[l] != 0) {
                    nz++;
                }
            }
            if (md.get(e).srv != 1 && nz > 1) {
                throw new RuntimeException("spn_mdd: mode " + (md.get(e).mode + 1) + " of node "
                        + (md.get(e).trans + 1) + " has " + md.get(e).srv + " servers and draws "
                        + "from " + nz + " levels; the enabling degree min_l floor(m(l)/enab(l)) "
                        + "is then not a product of per-level terms and admits no Kronecker form");
            }
            if (md.get(e).srv != 1 && nz == 0) {
                throw new RuntimeException("spn_mdd: mode " + (md.get(e).mode + 1) + " of node "
                        + (md.get(e).trans + 1) + " has " + md.get(e).srv + " servers but consumes "
                        + "from no place, so its enabling degree is unbounded and its firing rate "
                        + "undefined");
            }
        }

        // ---- phase levels for the multi-phase modes
        int[] phaseof = new int[E];                    // 0 = none; else 1-based level
        int Q = 0;
        for (int e = 0; options.descriptor && e < E; e++) {
            if (md.get(e).nph > 1) {
                Q++;
                phaseof[e] = L + Q;
            }
        }
        int K = L + Q;

        double[][] netm = new double[E][L];
        for (int e = 0; e < E; e++) {
            for (int l = 0; l < L; l++) {
                netm[e][l] = md.get(e).fire[l] - md.get(e).enab[l];
            }
        }

        // ---- initial marking and per-level bounds
        double[] init0 = initMarking(model, sn, places, R);
        double[] winv = placeInvariantWeights(netm, L);
        double vinv = 0;
        if (winv != null) {
            for (int l = 0; l < L; l++) {
                vinv += winv[l] * init0[l];
            }
        }
        double[] bound;
        if (options.bound != null && options.bound.length > 0) {
            bound = new double[L];
            for (int l = 0; l < L; l++) {
                bound[l] = options.bound.length == 1 ? options.bound[0] : options.bound[l];
            }
        } else if (winv != null) {
            bound = new double[L];
            for (int l = 0; l < L; l++) {
                bound[l] = winv[l] > 0 ? Math.floor(vinv / winv[l]) : vinv;
            }
        } else {
            throw new RuntimeException("spn_mdd: the net has no place invariant with positive "
                    + "weights, so the marking is not bounded a priori; pass options.bound");
        }

        int[] domain = new int[K];
        for (int l = 0; l < L; l++) {
            domain[l] = (int) bound[l] + 1;
        }
        for (int e = 0; e < E; e++) {
            if (phaseof[e] > 0) {
                domain[phaseof[e] - 1] = md.get(e).nph;
            }
        }

        int[] init = new int[K];
        for (int l = 0; l < L; l++) {
            init[l] = (int) init0[l];
        }
        for (int e = 0; e < E; e++) {
            if (phaseof[e] > 0) {
                for (int a = 0; a < md.get(e).nph; a++) {
                    if (md.get(e).pie[a] > 0) {
                        init[phaseof[e] - 1] = a;
                        break;
                    }
                }
            }
        }

        // ---- reachable set; the closure records which modes were ever disabled,
        // so the phase-memory question is answered without a second pass over |S|
        boolean[] everDisabled = new boolean[E];
        MDD mdd = new MDD(domain);
        mdd.insert(init);
        List<int[]> frontier = new ArrayList<int[]>();
        frontier.add(init);
        int head = 0;
        while (head < frontier.size()) {
            int[] s = frontier.get(head);
            head++;
            List<int[]> succ = new ArrayList<int[]>();
            successors(s, md, netm, phaseof, domain, L, succ, everDisabled);
            for (int r = 0; r < succ.size(); r++) {
                int[] t = succ.get(r);
                if (!mdd.member(t)) {
                    mdd.insert(t);
                    frontier.add(t);
                }
            }
            if (head > 1024 && 2 * head > frontier.size()) {
                frontier = new ArrayList<int[]>(frontier.subList(head, frontier.size()));
                head = 0;
            }
        }
        mdd.compact();

        for (int e = 0; e < E; e++) {
            if (phaseof[e] > 0 && options.descriptor && everDisabled[e]
                    && !"resume".equalsIgnoreCase(options.phmemory)) {
                throw new RuntimeException("spn_mdd: mode " + (md.get(e).mode + 1) + " of node "
                        + (md.get(e).trans + 1) + " has a phase-type firing time AND is disabled "
                        + "in some reachable marking. LINE discards the phase on disabling "
                        + "(preemptive repeat) but that reset is a joint condition on the place "
                        + "levels and has no Kronecker form, so this descriptor would encode "
                        + "preemptive resume instead and disagree with SolverCTMC. Pass "
                        + "phmemory=\"resume\" to accept the resume semantics.");
            }
        }

        // ---- Kronecker event matrices
        List<MddEvent> events = new ArrayList<MddEvent>();
        for (int e = 0; options.descriptor && e < E; e++) {
            SpnMode mde = md.get(e);
            Set<Integer> gate = new LinkedHashSet<Integer>();
            Set<Integer> touched = new LinkedHashSet<Integer>();
            for (int l = 0; l < L; l++) {
                if (mde.enab[l] > 0 || !Double.isInfinite(mde.inhib[l])) {
                    gate.add(Integer.valueOf(l));
                    touched.add(Integer.valueOf(l));
                }
            }
            for (int l = 0; l < L; l++) {
                if (netm[e][l] != 0) {
                    touched.add(Integer.valueOf(l));
                }
            }
            List<Integer> touchedSorted = new ArrayList<Integer>(touched);
            java.util.Collections.sort(touchedSorted);
            int degl = -1;
            if (mde.srv != 1) {
                for (int l = 0; l < L; l++) {
                    if (mde.enab[l] > 0) {
                        degl = l;                       // level carrying the degree
                        break;
                    }
                }
            }
            if (phaseof[e] == 0) {
                if (touchedSorted.isEmpty()) {
                    continue;
                }
                int[] lev = new int[touchedSorted.size()];
                MddLocalMatrix[] W = new MddLocalMatrix[touchedSorted.size()];
                for (int t = 0; t < touchedSorted.size(); t++) {
                    int l = touchedSorted.get(t).intValue();
                    lev[t] = l;
                    W[t] = placemat(l, mde, netm[e][l], domain[l], l == degl, 1.0);
                }
                // the scalar firing rate rides on the first touched level
                int l0 = touchedSorted.get(0).intValue();
                W[0] = placemat(l0, mde, netm[e][l0], domain[l0], l0 == degl, mde.D1[0][0]);
                events.add(new MddEvent(mde.trans, mde.mode, lev, W));
            } else {
                int q = phaseof[e] - 1;
                // (1) internal phase changes: marking unchanged, gated by enabling
                int nnzOff = 0;
                for (int a = 0; a < mde.nph; a++) {
                    for (int b = 0; b < mde.nph; b++) {
                        if (a != b && mde.D0[a][b] != 0) {
                            nnzOff++;
                        }
                    }
                }
                if (nnzOff > 0) {
                    List<Integer> g = new ArrayList<Integer>(gate);
                    java.util.Collections.sort(g);
                    int[] lev = new int[g.size() + 1];
                    MddLocalMatrix[] W = new MddLocalMatrix[g.size() + 1];
                    for (int t = 0; t < g.size(); t++) {
                        int l = g.get(t).intValue();
                        lev[t] = l;
                        W[t] = placemat(l, mde, 0.0, domain[l], false, 1.0);
                    }
                    MddLocalMatrix.Builder bld = new MddLocalMatrix.Builder(mde.nph);
                    for (int a = 0; a < mde.nph; a++) {
                        for (int b = 0; b < mde.nph; b++) {
                            if (a != b) {
                                bld.add(a, b, mde.D0[a][b]);
                            }
                        }
                    }
                    lev[g.size()] = q;
                    W[g.size()] = bld.build();
                    events.add(new MddEvent(mde.trans, mde.mode, lev, W));
                }
                // (2) firings: marking moved, phase redrawn through D1
                int[] lev = new int[touchedSorted.size() + 1];
                MddLocalMatrix[] W = new MddLocalMatrix[touchedSorted.size() + 1];
                for (int t = 0; t < touchedSorted.size(); t++) {
                    int l = touchedSorted.get(t).intValue();
                    lev[t] = l;
                    W[t] = placemat(l, mde, netm[e][l], domain[l], l == degl, 1.0);
                }
                MddLocalMatrix.Builder bld = new MddLocalMatrix.Builder(mde.nph);
                for (int a = 0; a < mde.nph; a++) {
                    for (int b = 0; b < mde.nph; b++) {
                        bld.add(a, b, mde.D1[a][b]);
                    }
                }
                lev[touchedSorted.size()] = q;
                W[touchedSorted.size()] = bld.build();
                events.add(new MddEvent(mde.trans, mde.mode, lev, W));
            }
        }

        MddDescriptor desc = new MddDescriptor();
        desc.K = K;
        desc.N = 0;                                    // the invariant below replaces it
        desc.domain = domain;
        desc.events = events;
        if (winv != null) {
            double[] w = new double[K];
            System.arraycopy(winv, 0, w, 0, L);
            desc.invariantWeights = w;
            desc.invariantValue = vinv;
        }

        // ---- descriptive information
        SpnInfo info = new SpnInfo();
        info.places = new int[P];
        info.placenames = new String[P];
        for (int pp = 0; pp < P; pp++) {
            info.places[pp] = places.get(pp).intValue();
            info.placenames[pp] = sn.nodenames.get(info.places[pp]);
        }
        info.classnames = new String[R];
        for (int k = 0; k < R; k++) {
            info.classnames[k] = sn.classnames.get(k);
        }
        info.levelkind = new int[K];
        for (int l = 0; l < L; l++) {
            info.levelkind[l] = 1;
        }
        for (int l = L; l < K; l++) {
            info.levelkind[l] = 2;
        }
        info.levelname = new String[K];
        for (int pp = 0; pp < P; pp++) {
            for (int k = 0; k < R; k++) {
                info.levelname[pp * R + k] = info.placenames[pp] + "." + info.classnames[k];
            }
        }
        for (int e = 0; e < E; e++) {
            if (phaseof[e] > 0) {
                info.levelname[phaseof[e] - 1] = "phase("
                        + sn.nodenames.get(md.get(e).trans) + ".m" + (md.get(e).mode + 1) + ")";
            }
        }
        info.modes = md;
        info.init = init;
        info.mdd = mdd;
        info.phaseof = phaseof;
        info.everDisabled = everDisabled;
        info.nplacelevels = L;
        info.nnodes = sn.nnodes;
        info.nclasses = R;
        info.descriptor = options.descriptor;

        if (options.verbose) {
            System.out.format("%nSPN -> MDD: %d places x %d classes = %d place levels, "
                    + "%d phase levels%n", P, R, L, Q);
            System.out.format("  modes = %d, |S| = %d%n", E, mdd.cardinality());
        }

        SpnResult res = new SpnResult();
        res.mdds = mdd.toStruct();
        res.desc = desc;
        res.info = info;
        return res;
    }

    // -----------------------------------------------------------------------
    private static MatrixCell firingProc(TransitionNodeParam tp, int m) {
        if (tp.firingproc == null) {
            return null;
        }
        int i = 0;
        for (Map.Entry<jline.lang.Mode, MatrixCell> en : tp.firingproc.entrySet()) {
            if (i == m) {
                return en.getValue();
            }
            i++;
        }
        return null;
    }

    private static double[] firingPie(TransitionNodeParam tp, int m, int nph) {
        double[] pie = new double[nph];
        Matrix pv = null;
        if (tp.firingpie != null) {
            int i = 0;
            for (Map.Entry<jline.lang.Mode, Matrix> en : tp.firingpie.entrySet()) {
                if (i == m) {
                    pv = en.getValue();
                    break;
                }
                i++;
            }
        }
        if (pv == null || pv.getNumElements() == 0) {
            pie[0] = 1.0;
            return pie;
        }
        double s = 0;
        for (int a = 0; a < nph && a < pv.getNumElements(); a++) {
            pie[a] = pv.get(a);
            s += pie[a];
        }
        if (s <= 0) {
            Arrays.fill(pie, 0.0);
            pie[0] = 1.0;
            return pie;
        }
        for (int a = 0; a < nph; a++) {
            pie[a] /= s;
        }
        return pie;
    }

    private static double[][] toArray(Matrix m) {
        int r = m.getNumRows();
        int c = m.getNumCols();
        double[][] a = new double[r][c];
        for (int i = 0; i < r; i++) {
            for (int j = 0; j < c; j++) {
                a[i][j] = m.get(i, j);
            }
        }
        return a;
    }

    /** (nnodes x nclasses) arc matrix -> place-major length-(P*R) level vector. */
    private static double[] arcvec(List<Matrix> mats, int m, List<Integer> places, int R,
                                   double fillval) {
        int P = places.size();
        double[] v = new double[P * R];
        Arrays.fill(v, fillval);
        if (mats == null || mats.size() <= m || mats.get(m) == null) {
            return v;
        }
        Matrix mat = mats.get(m);
        for (int pp = 0; pp < P; pp++) {
            for (int k = 0; k < R; k++) {
                double x = mat.get(places.get(pp).intValue(), k);
                if (!Double.isInfinite(fillval) && fillval == 0) {
                    x = Math.max(0.0, x);
                }
                v[pp * R + k] = x;
            }
        }
        return v;
    }

    /**
     * Local matrix of one mode at place level l.
     *
     * <p>Move the count by net, but only from local states that satisfy this
     * level's enabling and inhibition. applyDegree scales each row by the
     * enabling degree min(floor(m/enab), srv), i.e. the number of concurrently
     * firing servers; it is carried by the single enabling level of a
     * multi-server mode. scale carries the scalar firing rate on whichever level
     * the caller chose to put it.</p>
     */
    private static MddLocalMatrix placemat(int l, SpnMode mde, double net, int d,
                                           boolean applyDegree, double scale) {
        MddLocalMatrix.Builder bld = new MddLocalMatrix.Builder(d);
        for (int i = 0; i < d; i++) {
            if (!(i >= mde.enab[l] && i < mde.inhib[l])) {
                continue;
            }
            int j = i + (int) net;
            if (j < 0 || j > d - 1) {
                continue;
            }
            double val = 1.0;
            if (applyDegree && mde.enab[l] > 0) {
                val = Math.min(Math.floor(i / mde.enab[l]), mde.srv);
            }
            bld.add(i, j, scale * val);
        }
        return bld.build();
    }

    /**
     * Tokens per (place, class) at time zero.
     *
     * <p>Taken from the model state when set, and from the reference station
     * otherwise.</p>
     */
    private static double[] initMarking(Network model, NetworkStruct sn, List<Integer> places,
                                        int R) {
        int P = places.size();
        double[] init = new double[P * R];
        boolean anyset = false;
        for (int pp = 0; pp < P; pp++) {
            Node node = sn.nodes.get(places.get(pp).intValue());
            Matrix st = null;
            if (node instanceof jline.lang.nodes.StatefulNode) {
                st = ((jline.lang.nodes.StatefulNode) node).getState();
            }
            if (st != null && st.getNumElements() > 0) {
                anyset = true;
                for (int k = 0; k < R && k < st.getNumElements(); k++) {
                    init[pp * R + k] = st.get(k);
                }
            }
        }
        boolean allZero = true;
        for (int i = 0; i < init.length; i++) {
            if (init[i] != 0) {
                allZero = false;
                break;
            }
        }
        if (!anyset || allZero) {
            for (int k = 0; k < R; k++) {
                int ref = (int) sn.refstat.get(k);
                for (int pp = 0; pp < P; pp++) {
                    if ((int) sn.nodeToStation.get(places.get(pp).intValue()) == ref) {
                        init[pp * R + k] = sn.njobs.get(k);
                    }
                }
            }
        }
        return init;
    }

    /**
     * A place invariant w'*m = const, as a positive left null vector of the
     * incidence matrix over the place levels.
     */
    /**
     * A strictly positive place invariant w, i.e. w &gt;= 0 with netm w = 0.
     *
     * <p>Scanning a null-space BASIS for a positive vector is not enough, and
     * the fork-join net is the counterexample: its two minimal-support
     * invariants are (1,1,0,1) and (1,0,1,1), neither positive, while their sum
     * (2,1,1,2) is. Which basis a codebase's null space routine returns then
     * decides whether the net is accepted, which is how MATLAB and the JAR came
     * to disagree on it. So the non-negative generators are computed directly,
     * by Farkas' algorithm on [netm' | I] -- the same construction
     * {@link Spn_sinvariants} uses on the incidence matrix -- and summed.</p>
     */
    private static double[] placeInvariantWeights(double[][] netm, int L) {
        int E = netm.length;
        boolean conservative = true;
        for (int e = 0; e < E; e++) {
            double s = 0;
            for (int l = 0; l < L; l++) {
                s += netm[e][l];
            }
            if (Math.abs(s) > 1e-12) {
                conservative = false;
                break;
            }
        }
        if (conservative) {
            double[] w = new double[L];
            Arrays.fill(w, 1.0);
            return w;
        }
        if (E == 0) {
            return null;
        }
        List<double[]> M = new ArrayList<double[]>();   // L x E, level l's column
        List<double[]> B = new ArrayList<double[]>();   // the combination it stands for
        for (int l = 0; l < L; l++) {
            double[] mrow = new double[E];
            for (int e = 0; e < E; e++) {
                mrow[e] = netm[e][l];
            }
            double[] brow = new double[L];
            brow[l] = 1.0;
            M.add(mrow);
            B.add(brow);
        }
        for (int e = 0; e < E; e++) {
            List<double[]> Mn = new ArrayList<double[]>();
            List<double[]> Bn = new ArrayList<double[]>();
            for (int i = 0; i < M.size(); i++) {
                if (Math.abs(M.get(i)[e]) < 1e-12) {
                    Mn.add(M.get(i));
                    Bn.add(B.get(i));
                }
            }
            for (int a = 0; a < M.size(); a++) {
                if (!(M.get(a)[e] > 1e-12)) {
                    continue;
                }
                for (int b = 0; b < M.size(); b++) {
                    if (!(M.get(b)[e] < -1e-12)) {
                        continue;
                    }
                    double ca = -M.get(b)[e];
                    double cb = M.get(a)[e];
                    double[] c = new double[E];
                    double[] d = new double[L];
                    for (int k = 0; k < E; k++) {
                        c[k] = ca * M.get(a)[k] + cb * M.get(b)[k];
                    }
                    for (int k = 0; k < L; k++) {
                        d[k] = ca * B.get(a)[k] + cb * B.get(b)[k];
                    }
                    long g = gcdVec(c, d);
                    if (g > 0) {
                        for (int k = 0; k < E; k++) {
                            c[k] /= g;
                        }
                        for (int k = 0; k < L; k++) {
                            d[k] /= g;
                        }
                    }
                    Mn.add(c);
                    Bn.add(d);
                }
            }
            minimalSupport(Mn, Bn, L);
            M = Mn;
            B = Bn;
        }
        if (B.isEmpty()) {
            return null;
        }
        double[] s = new double[L];
        for (int i = 0; i < B.size(); i++) {
            for (int l = 0; l < L; l++) {
                s[l] += B.get(i)[l];
            }
        }
        double mn = Double.MAX_VALUE;
        for (int l = 0; l < L; l++) {
            if (!(s[l] > 1e-9)) {
                return null;
            }
            mn = Math.min(mn, s[l]);
        }
        double[] w = new double[L];
        for (int l = 0; l < L; l++) {
            w[l] = s[l] / mn;
        }
        return w;
    }

    /** gcd of the integral entries of both vectors, 0 when any is not integral. */
    private static long gcdVec(double[] a, double[] b) {
        long g = 0;
        double[][] all = {a, b};
        for (int i = 0; i < 2; i++) {
            for (int k = 0; k < all[i].length; k++) {
                double x = all[i][k];
                if (Math.abs(x - Math.rint(x)) > 1e-9) {
                    return 0;
                }
                long y = Math.abs((long) Math.rint(x));
                while (y != 0) {
                    long t = g % y;
                    g = y;
                    y = t;
                }
            }
        }
        return g;
    }

    /**
     * Drop every row whose support strictly contains another's, which is what
     * leaves the minimal supports and stops the pair expansion from blowing up.
     */
    private static void minimalSupport(List<double[]> M, List<double[]> B, int L) {
        int n = B.size();
        boolean[] drop = new boolean[n];
        for (int i = 0; i < n; i++) {
            if (drop[i]) {
                continue;
            }
            for (int j = 0; j < n; j++) {
                if (i == j || drop[j]) {
                    continue;
                }
                boolean contained = true;
                boolean strict = false;
                for (int l = 0; l < L && contained; l++) {
                    boolean si = Math.abs(B.get(i)[l]) > 1e-12;
                    boolean sj = Math.abs(B.get(j)[l]) > 1e-12;
                    if (sj && !si) {
                        contained = false;
                    }
                    if (si && !sj) {
                        strict = true;
                    }
                }
                if (contained && strict) {
                    drop[i] = true;
                    break;
                }
            }
        }
        for (int i = n - 1; i >= 0; i--) {
            if (drop[i]) {
                M.remove(i);
                B.remove(i);
            }
        }
    }

    /** Successors of s, recording which modes were found disabled here. */
    private static void successors(int[] s, List<SpnMode> md, double[][] netm, int[] phaseof,
                                   int[] domain, int L, List<int[]> out, boolean[] everDisabled) {
        int E = md.size();
        for (int e = 0; e < E; e++) {
            SpnMode mde = md.get(e);
            boolean enabled = true;
            for (int l = 0; l < L; l++) {
                if (s[l] < mde.enab[l] || s[l] >= mde.inhib[l]) {
                    enabled = false;
                    break;
                }
            }
            if (!enabled) {
                everDisabled[e] = true;
                continue;
            }
            if (phaseof[e] == 0) {
                int[] t = s.clone();
                boolean ok = true;
                for (int l = 0; l < L; l++) {
                    int nv = s[l] + (int) netm[e][l];
                    if (nv < 0 || nv > domain[l] - 1) {
                        ok = false;
                        break;
                    }
                    t[l] = nv;
                }
                if (ok) {
                    out.add(t);
                }
            } else {
                int q = phaseof[e] - 1;
                int ph = s[q];
                for (int j = 0; j < mde.nph; j++) {
                    if (j != ph && mde.D0[ph][j] != 0) {
                        int[] t = s.clone();
                        t[q] = j;
                        out.add(t);
                    }
                }
                int[] moved = s.clone();
                boolean ok = true;
                for (int l = 0; l < L; l++) {
                    int nv = s[l] + (int) netm[e][l];
                    if (nv < 0 || nv > domain[l] - 1) {
                        ok = false;
                        break;
                    }
                    moved[l] = nv;
                }
                if (ok) {
                    for (int j = 0; j < mde.nph; j++) {
                        if (mde.D1[ph][j] != 0) {
                            int[] t = moved.clone();
                            t[q] = j;
                            out.add(t);
                        }
                    }
                }
            }
        }
    }
}
