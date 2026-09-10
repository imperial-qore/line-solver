package jline.api.spn;

import jline.api.mdd.Mdd_rec;
import jline.api.mdd.MddStruct;
import jline.api.spn.Spn_mdd.SpnInfo;
import jline.api.spn.Spn_mdd.SpnMode;
import jline.api.spn.Spn_rec_enabled.SpnEnabling;

/**
 * Stationary measures of a product-form stochastic Petri net from the MDD-rec
 * masses.
 *
 * <p>S. Balsamo, A. Marin, I. Stojic, FGCS 111 (2020) 475-490, Sec. 3.1 for the
 * definitions and Sec. 5.3 for the recursions they are read off.</p>
 *
 * <pre>
 *   n(P_j) = sum_k k P(m_j = k)                    mean tokens
 *   u(P_j) = 1 - P(m_j = 0)                        place utilization
 *   u(T_j) = P(e_j &gt;= 1)                           transition utilization
 *   x(T_j) = sum_k min(k, c_j) W(T_j) P(e_j = k)   throughput
 *   x(P_j) = sum_T I_j(T) x(T)                     tokens removed per unit time
 * </pre>
 *
 * <p>ONE DEVIATION FROM THE PAPER'S x(T_j), AND IT IS A GENERALISATION. The
 * paper writes x(T_j) = sum_k k W(T_j) P(e_j = k), which is INFINITE-SERVER
 * firing semantics -- every enabling set fires in parallel. LINE's own rate law
 * is min(enabling degree, nmodeservers) * W(T), so c_j above is the mode's
 * server count: c_j = 1 recovers single-server semantics,
 * x = W(T) P(e &gt;= 1), and c_j = infinity recovers the paper's formula exactly.
 * Using the paper's form for a single-server mode would report a throughput that
 * grows with the token population of a net whose transition can only fire one
 * set at a time.</p>
 *
 * <p>The measures come out of ONE reachable set and ONE set of g_l, so they are
 * mutually consistent by construction: no per-measure fixed point, no
 * iteration.</p>
 *
 * <p>MATLAB twin: {@code spn_metrics.m}. Python twin:
 * {@code api/spn/metrics.py}.</p>
 */
public class Spn_metrics {

    private Spn_metrics() {}

    /** The stationary measures of Sec. 3.1, per place level and per mode. */
    public static class SpnMetricsResult {
        /** The normalising constant the measures are taken against. */
        public double G;
        /** Mean tokens per place level. */
        public double[] tokens;
        /** Place utilization, P(m_j &gt; 0). */
        public double[] placeUtil;
        /** Place throughput, tokens removed per unit time. */
        public double[] placeTput;
        /** Transition (mode) utilization, P(e_j &gt;= 1). */
        public double[] modeUtil;
        /** Transition (mode) throughput. */
        public double[] modeTput;
        /** marginal[l][k] = P(m_l = k). */
        public double[][] marginal;
    }

    /**
     * Every measure of Sec. 3.1 from one diagram and one product form.
     *
     * @param mdds the reachable set built by {@code Spn_mdd}
     * @param g per-level product-form factors g_l(v)
     * @param info the metadata {@code Spn_mdd} returned alongside the diagram
     * @return the stationary measures
     */
    public static SpnMetricsResult spn_metrics(MddStruct mdds, double[][] g, SpnInfo info) {
        int L = info.nplacelevels;
        SpnMetricsResult out = new SpnMetricsResult();
        out.G = Mdd_rec.mdd_rec(mdds, g);
        if (!(out.G > 0)) {
            throw new RuntimeException("spn_metrics: the normalising constant is not positive; the "
                    + "g_l passed do not describe a product form over this reachable set");
        }

        out.marginal = new double[L][];
        out.tokens = new double[L];
        out.placeUtil = new double[L];
        out.placeTput = new double[L];
        for (int l = 0; l < L; l++) {
            double[] mass = Mdd_rec.mdd_rec_marginal(mdds, g, l);
            out.marginal[l] = new double[mass.length];
            for (int k = 0; k < mass.length; k++) {
                out.marginal[l][k] = mass[k] / out.G;
                out.tokens[l] += k * out.marginal[l][k];
            }
            out.placeUtil[l] = 1.0 - out.marginal[l][0];
        }

        int E = info.modes.size();
        out.modeUtil = new double[E];
        out.modeTput = new double[E];
        for (int e = 0; e < E; e++) {
            SpnMode mde = info.modes.get(e);
            SpnEnabling en = Spn_rec_enabled.spn_rec_enabled(mdds, g, mde, L);
            out.modeUtil[e] = en.ge[1] / out.G;
            // W(T) is the scalar firing rate of the mode; a phase-type firing
            // time has no single rate, so its throughput is left to the
            // phase-level marginal rather than reported through this formula.
            if (mde.nph > 1) {
                throw new RuntimeException("spn_metrics: mode " + (mde.mode + 1) + " of node "
                        + mde.trans + " has a phase-type firing time, whose throughput is not W(T) "
                        + "times an enabling probability; read it from the phase-level marginal "
                        + "instead");
            }
            // The formula below is W(T)*E[min(enabling degree, servers)], which
            // is the rate law only when no marking-dependent multiplier is in
            // play. With one, the firing rate is not a function of the enabling
            // degree at all, so the enabling-degree law is the wrong summary to
            // take it from. The MARGINALS above are unaffected -- they come from
            // the product form, not the rates.
            if (mde.dep != null) {
                throw new RuntimeException("spn_metrics: mode " + (mde.mode + 1) + " of node "
                        + mde.trans + " has a marking-dependent firing rate, so its throughput is "
                        + "not W(T) times a function of the enabling degree and cannot be read "
                        + "from the enabling-degree law. The token marginals are still exact");
            }
            double rate = mde.D1[0][0];
            double x = 0.0;
            for (int k = 1; k < en.eq.length; k++) {
                double served = Double.isInfinite(mde.srv) ? k : Math.min(k, mde.srv);
                x += served * rate * (en.eq[k] / out.G);
            }
            out.modeTput[e] = x;
            for (int l = 0; l < L; l++) {
                if (mde.enab[l] > 0) {
                    out.placeTput[l] += mde.enab[l] * x;
                }
            }
        }
        return out;
    }
}
