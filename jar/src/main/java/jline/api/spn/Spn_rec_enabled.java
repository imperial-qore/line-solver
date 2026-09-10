package jline.api.spn;

import jline.api.mdd.Mdd_rec;
import jline.api.mdd.MddStruct;
import jline.api.spn.Spn_mdd.SpnMode;

/**
 * Enabling-degree distribution of one mode of a product-form stochastic Petri
 * net, by the masked MDD-rec recursion.
 *
 * <p>S. Balsamo, A. Marin, I. Stojic, FGCS 111 (2020) 475-490, Sec. 5.3.</p>
 *
 * <p>The enabling degree of a mode in marking m is
 * e(m) = min_{l : I_l &gt; 0} floor(m_l / I_l), zero when any inhibitor threshold
 * is met. P(e &gt;= k) is therefore the mass of the marking subset in which EVERY
 * input level holds at least k*I_l tokens and no inhibitor fires, which is a
 * per-level restriction and so exactly what {@code Mdd_rec.mdd_rec_masked}
 * computes: the paper's second modified recurrence is the same walk under a
 * different mask, not a second algorithm.</p>
 *
 * <p>The masses returned are UNNORMALISED, as in the paper; divide by G from
 * {@code Mdd_rec.mdd_rec} for probabilities. {@code Spn_metrics} does that and
 * turns them into the transition measures.</p>
 *
 * <p>MATLAB twin: {@code spn_rec_enabled.m}. Python twin:
 * {@code api/spn/rec_enabled.py}.</p>
 */
public class Spn_rec_enabled {

    private Spn_rec_enabled() {}

    /** Unnormalised enabling-degree masses of one mode. */
    public static class SpnEnabling {
        /** ge[k] is the mass of {e &gt;= k}; ge[0] is the whole reachable set. */
        public double[] ge;
        /** eq[k] is the mass of {e == k}, i.e. ge[k] - ge[k+1]. */
        public double[] eq;
        /** The largest enabling degree the place bounds permit, E_j in the paper. */
        public int maxDegree;
    }

    /**
     * Enabling-degree masses of one mode over the reachable set in {@code mdds}.
     *
     * @param mdds the reachable set built by {@code Spn_mdd}
     * @param g per-level product-form factors, one row per level
     * @param mde the mode, as returned in {@code SpnInfo.modes}
     * @param nplacelevels how many leading levels are place levels
     * @return the unnormalised enabling-degree masses
     */
    public static SpnEnabling spn_rec_enabled(MddStruct mdds, double[][] g, SpnMode mde,
                                              int nplacelevels) {
        if (nplacelevels > mdds.K) {
            throw new RuntimeException("spn_rec_enabled: more place levels than diagram levels");
        }
        // E_j: the enabling degree cannot exceed what the tightest input place
        // bound allows. A mode with no input place has no bound and is refused
        // rather than silently truncated, matching Spn_mdd's own refusal.
        int emax = 0;
        boolean hasInput = false;
        for (int l = 0; l < nplacelevels; l++) {
            if (!(mde.enab[l] > 0)) {
                continue;
            }
            int cap = (int) Math.floor((mdds.domain[l] - 1) / mde.enab[l]);
            emax = hasInput ? Math.min(emax, cap) : cap;
            hasInput = true;
        }
        if (!hasInput) {
            throw new RuntimeException("spn_rec_enabled: the mode consumes from no place, so its "
                    + "enabling degree is unbounded");
        }

        SpnEnabling out = new SpnEnabling();
        out.maxDegree = emax;
        out.ge = new double[emax + 2];
        out.eq = new double[emax + 2];
        for (int k = 0; k <= emax; k++) {
            boolean[][] mask = new boolean[mdds.K][];
            for (int j = 0; j < mdds.K; j++) {
                mask[j] = new boolean[mdds.domain[j]];
                for (int v = 0; v < mdds.domain[j]; v++) {
                    mask[j][v] = true;
                }
            }
            for (int l = 0; l < nplacelevels; l++) {
                double need = mde.enab[l] * k;
                for (int v = 0; v < mdds.domain[l]; v++) {
                    boolean shortOfTokens = v < need;
                    boolean inhibited = v >= mde.inhib[l];
                    // k = 0 asks only that the marking exist, so the inhibitor
                    // test belongs to k >= 1: e = 0 covers the inhibited
                    // markings too.
                    if (shortOfTokens || (k > 0 && inhibited)) {
                        mask[l][v] = false;
                    }
                }
            }
            out.ge[k] = Mdd_rec.mdd_rec_masked(mdds, g, mask);
        }
        for (int k = 0; k <= emax; k++) {
            out.eq[k] = out.ge[k] - out.ge[k + 1];
        }
        return out;
    }
}
