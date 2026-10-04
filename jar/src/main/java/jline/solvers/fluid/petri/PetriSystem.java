/*
 * Copyright (c) 2012-2026, QORE Lab, Imperial College London
 * All rights reserved.
 */

package jline.solvers.fluid.petri;

import jline.lang.constant.TimingStrategy;
import jline.util.matrix.Matrix;

import java.util.ArrayList;
import java.util.List;
import java.util.Map;
import java.util.TreeMap;

/**
 * The closed enabling term of every transition mode, the rate of every event
 * column, and the drift Jacobian. Java twin of the MATLAB
 * {@code fluid_petri_theta}, {@code fluid_petri_rates} and
 * {@code fluid_petri_jacobian}.
 */
public final class PetriSystem {

    private PetriSystem() {
    }

    /** The enabling terms and their derivatives, for one iterate. */
    public static final class Theta {
        public double[] theta;
        public int[][] dslot;
        public double[][] dval;
        public double[] dep;
        public int[][] depslot;
        public double[][] depval;
    }

    /** One entry of the closure covariance, zero where the drift never reads it. */
    private static double sig(PetriTerms t, double[] s2, int a, int b) {
        if (a >= t.pairIndex.length || b >= t.pairIndex.length) {
            return 0.0;
        }
        int i = t.pairIndex[a][b];
        return i < 0 ? 0.0 : s2[i];
    }

    /**
     * The closed enabling term of every mode, and its derivative.
     *
     * <pre>  theta_j = ( prod_b Phi((thr_b - m_b)/sd_b) ) * E[ min_a(m_a/w_a), c_j ]</pre>
     *
     * <p>Both factors collapse to their first-order form at zero variance -- Phi
     * becomes the hard indicator and the min closure becomes min() -- so the
     * mean-field limit is one code path, not two.
     *
     * <p>THE INHIBITOR GATE IS WHY A PETRI NET NEEDS A SMOOTHED CLOSURE AT ALL,
     * quite apart from accuracy: the indicator is a step, and a Newton solver
     * has no derivative to descend on a step.
     *
     * <p>THE VARIANCES ARE UNKNOWNS, NOT FUNCTIONS OF X: s2 is pinned by its own
     * consistency row in the DAE, so the derivative is with respect to the MEANS
     * only.
     */
    public static Theta theta(PetriTerms t, double[] x, double[] s2) {
        int nmod = t.modes.size();
        Theta th = new Theta();
        th.theta = new double[nmod];
        th.dslot = new int[nmod][];
        th.dval = new double[nmod][];
        th.dep = new double[nmod];
        th.depslot = new int[nmod][];
        th.depval = new double[nmod][];
        if (s2 == null || s2.length == 0) {
            s2 = new double[Math.max(t.npair, 1)];
        }

        for (int j = 0; j < nmod; j++) {
            PetriMode md = t.modes.get(j);
            th.dep[j] = 1.0;
            th.depslot[j] = new int[0];
            th.depval[j] = new double[0];

            int A = md.arcSlot.size();
            double[] mu = new double[A];
            for (int a = 0; a < A; a++) {
                mu[a] = x[md.arcSlot.get(a)] / md.arcW.get(a);
            }
            Matrix Sarg = null;
            if (md.closable && A > 0) {
                Sarg = new Matrix(A, A);
                Sarg.fill(0.0);
                for (int a = 0; a < A; a++) {
                    for (int b = 0; b < A; b++) {
                        double v = sig(t, s2, md.arcSlot.get(a), md.arcSlot.get(b))
                                / (md.arcW.get(a) * md.arcW.get(b));
                        Sarg.set(a, b, v);
                    }
                }
            }
            PetriClosures.MinMulti mm = PetriClosures.minMultiClosure(mu, Sarg, md.c);

            int nb = md.inhSlot.size();
            double[] gate = new double[nb];
            double[] dgate = new double[nb];
            double ginh = 1.0;
            for (int b = 0; b < nb; b++) {
                double mb = x[md.inhSlot.get(b)];
                double thr = md.inhThr.get(b);
                double vb = sig(t, s2, md.inhSlot.get(b), md.inhSlot.get(b));
                if (vb > 0) {
                    double sd = Math.sqrt(vb);
                    double zb = (thr - mb) / sd;
                    gate[b] = PetriClosures.normCdf(zb);
                    dgate[b] = -PetriClosures.normPdf(zb) / sd;
                } else {
                    gate[b] = (mb < thr - PetriClosures.FINE_TOL * Math.max(1.0, thr)) ? 1.0 : 0.0;
                    dgate[b] = 0.0;
                }
                ginh *= gate[b];
            }

            List<Integer> slots = new ArrayList<Integer>();
            List<Double> vals = new ArrayList<Double>();
            for (int a = 0; a < A; a++) {
                slots.add(md.arcSlot.get(a));
                vals.add(ginh * mm.g[a] / md.arcW.get(a));
            }
            for (int b = 0; b < nb; b++) {
                double others;
                if (gate[b] != 0) {
                    others = ginh / gate[b];
                } else {
                    others = 1.0;
                    for (int q = 0; q < nb; q++) {
                        if (q != b) {
                            others *= gate[q];
                        }
                    }
                }
                slots.add(md.inhSlot.get(b));
                vals.add(mm.h * others * dgate[b]);
            }
            // an arc and an inhibitor arc may share a coordinate, so accumulate
            TreeMap<Integer, Double> acc = new TreeMap<Integer, Double>();
            for (int q = 0; q < slots.size(); q++) {
                Integer s = slots.get(q);
                Double prev = acc.get(s);
                acc.put(s, (prev == null ? 0.0 : prev.doubleValue()) + vals.get(q).doubleValue());
            }
            int[] us = new int[acc.size()];
            double[] uv = new double[acc.size()];
            int p = 0;
            for (Map.Entry<Integer, Double> e : acc.entrySet()) {
                us[p] = e.getKey().intValue();
                uv[p] = e.getValue().doubleValue();
                p++;
            }
            th.theta[j] = ginh * mm.h;
            th.dslot[j] = us;
            th.dval[j] = uv;

            if (md.dep != null) {
                depOf(t, md, x, th, j);
            }
        }
        return th;
    }

    /**
     * The marking-dependent firing multiplier and its gradient.
     *
     * <p>g is a user function of the (nnodes x nclasses) marking, so it is
     * evaluated at the MEAN marking -- a first-order closure of g, the same
     * order at which SolverCTMC evaluates it per state and the only one
     * available without the distribution of the marking. Its gradient has no
     * analytic form, so it is taken by central differences.
     */
    private static void depOf(PetriTerms t, PetriMode md, double[] x, Theta th, int j) {
        Matrix mm = new Matrix(t.I, t.K);
        mm.fill(0.0);
        for (int s = 0; s < t.nm; s++) {
            mm.set(t.coordNode[s], t.coordClass[s], x[s]);
        }
        th.dep[j] = md.dep.apply(mm).doubleValue();
        List<Integer> slots = new ArrayList<Integer>();
        List<Double> vals = new ArrayList<Double>();
        for (int s = 0; s < t.nm; s++) {
            double h = Math.max(1e-6 * Math.abs(x[s]), 1e-6);
            int i = t.coordNode[s];
            int k = t.coordClass[s];
            Matrix mp = mm.copy();
            mp.set(i, k, mm.get(i, k) + h);
            Matrix mn = mm.copy();
            mn.set(i, k, Math.max(0.0, mm.get(i, k) - h));
            double hh = mp.get(i, k) - mn.get(i, k);
            if (hh <= 0) {
                continue;
            }
            double d = (md.dep.apply(mp).doubleValue() - md.dep.apply(mn).doubleValue()) / hh;
            if (d != 0) {
                slots.add(s);
                vals.add(d);
            }
        }
        int[] sl = new int[slots.size()];
        double[] vl = new double[vals.size()];
        for (int q = 0; q < sl.length; q++) {
            sl[q] = slots.get(q).intValue();
            vl[q] = vals.get(q).doubleValue();
        }
        th.depslot[j] = sl;
        th.depval[j] = vl;
    }

    /**
     * The rate of every event column.
     *
     * <pre>
     *   kind 1  firing of mode j    single phase: rateBase*theta*dep
     *                               multi  phase: rateBase*y(j,h)
     *   kind 2  internal phase change             rateBase*y(j,h)
     *   kind 3  exogenous arrival                 a constant
     *   kind 4  firing of an IMMEDIATE mode       phi_j, an algebraic unknown
     *   kind 5  the server latch                  mu_j, a free-sign unknown
     * </pre>
     */
    public static double[] rates(PetriTerms t, double[] x, double[] phi, double[] mu, Theta th) {
        int nmod = t.modes.size();
        int[] immPos = new int[nmod];
        int[] latchPos = new int[nmod];
        for (int j = 0; j < nmod; j++) {
            immPos[j] = -1;
            latchPos[j] = -1;
        }
        for (int q = 0; q < t.immIdx.size(); q++) {
            immPos[t.immIdx.get(q)] = q;
        }
        for (int q = 0; q < t.latchMode.size(); q++) {
            latchPos[t.latchMode.get(q)] = q;
        }
        double[] r = new double[t.nev];
        for (int e = 0; e < t.nev; e++) {
            int k = t.evKind[e];
            if (k == 3) {
                r[e] = t.rateBase[e];
            } else if (k == 4) {
                r[e] = (phi == null) ? 0.0 : phi[immPos[t.evMode[e]]];
            } else if (k == 5) {
                r[e] = (mu == null) ? 0.0 : mu[latchPos[t.evMode[e]]];
            } else {
                int j = t.evMode[e];
                PetriMode md = t.modes.get(j);
                if (md.nph == 1) {
                    r[e] = t.rateBase[e] * th.theta[j] * th.dep[j];
                } else {
                    r[e] = t.rateBase[e] * x[md.zblk[t.evPhase[e]]];
                }
            }
        }
        return r;
    }

    /**
     * Drift Jacobian A = D * dR/dX.
     *
     * <p>This is what the Lyapunov equation of the linear noise approximation is
     * written about, so it has to be the derivative of the SAME rate vector
     * {@link #rates} returns: a covariance solved about an inconsistent Jacobian
     * is not the covariance of anything. THE VARIANCES ARE HELD.
     */
    public static Matrix jacobian(PetriTerms t, Theta th) {
        Matrix Jr = new Matrix(Math.max(t.nev, 1), t.nstate);
        Jr.fill(0.0);
        for (int e = 0; e < t.nev; e++) {
            int k = t.evKind[e];
            if (k == 3 || k == 4 || k == 5) {
                continue;
            }
            int j = t.evMode[e];
            PetriMode md = t.modes.get(j);
            double base = t.rateBase[e];
            if (md.nph == 1) {
                for (int q = 0; q < th.dslot[j].length; q++) {
                    int s = th.dslot[j][q];
                    Jr.set(e, s, Jr.get(e, s) + base * th.dep[j] * th.dval[j][q]);
                }
                for (int q = 0; q < th.depslot[j].length; q++) {
                    int s = th.depslot[j][q];
                    Jr.set(e, s, Jr.get(e, s) + base * th.theta[j] * th.depval[j][q]);
                }
            } else {
                int zc = md.zblk[t.evPhase[e]];
                Jr.set(e, zc, Jr.get(e, zc) + base);
            }
        }
        return t.D.mult(Jr);
    }

    /** Whether a mode is timed, spelled once. */
    public static boolean isTimed(PetriMode md) {
        return md.timing == TimingStrategy.TIMED;
    }
}
