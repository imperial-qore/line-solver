/**
 * @file Characteristic times for h-LRU / LRU(m) cache lists
 *
 * Solves the TTL (characteristic-time) fixed point of the list-based h-LRU
 * (LRU(m)) replacement policy (Gast and Van Houdt, SIGMETRICS 2015): h LRU
 * lists, a miss inserts at the head of list 1, a hit in list l exchanges the
 * item with the tail of list l+1. Under the approximation the level process
 * of an item with request rate lam is a birth-death chain with up-probability
 * 1-e(l) and down-probability e(l), e(l)=exp(-lam*T(l)), so the stationary
 * list probabilities are pi_l proportional to prod_{s<=l} (1-e_s)/e_s.
 *
 * @since LINE 3.0
 */
package jline.api.cache;

import org.apache.commons.math3.util.FastMath;

import jline.util.matrix.Matrix;

public final class Cache_t_hlru {
    private Cache_t_hlru() {}

    /**
     * Characteristic time of each list of an h-LRU cache, solved by per-list
     * bisection with Gauss-Seidel sweeps.
     *
     * @param gamma per-item request rates; an (n x h) matrix is accepted for
     *              backward compatibility (first column used)
     * @param m     (1 x h) list capacities
     * @return (1 x h) characteristic time of each list
     */
    public static Matrix cache_t_hlru(final Matrix gamma, final Matrix m) {
        int n = gamma.getNumRows();
        int h = m.length();
        double[] lam = new double[n];
        double lamMean = 0.0;
        for (int k = 0; k < n; k++) {
            lam[k] = gamma.get(k, 0);
            lamMean += lam[k];
        }
        lamMean = FastMath.max(lamMean / n, 1e-14);

        double[] t = new double[h];
        for (int l = 0; l < h; l++) {
            t[l] = 1.0 / lamMean;
        }
        for (int sweep = 0; sweep < 200; sweep++) {
            double maxRel = 0.0;
            for (int l = 0; l < h; l++) {
                double told = t[l];
                double lo = 0.0;
                double hi = FastMath.max(t[l], 1.0 / lamMean);
                while (occupancy(lam, t, l, hi, h) < m.get(l) && hi < 1e12) {
                    hi = 2 * hi;
                }
                for (int it = 0; it < 100; it++) {
                    double mid = 0.5 * (lo + hi);
                    if (occupancy(lam, t, l, mid, h) < m.get(l)) {
                        lo = mid;
                    } else {
                        hi = mid;
                    }
                }
                t[l] = 0.5 * (lo + hi);
                maxRel = FastMath.max(maxRel, FastMath.abs(t[l] - told) / FastMath.max(told, 1e-14));
            }
            if (maxRel < 1e-8) {
                break;
            }
        }
        Matrix res = new Matrix(1, h);
        for (int l = 0; l < h; l++) {
            res.set(0, l, t[l]);
        }
        return res;
    }

    /** Total occupancy of list l when its characteristic time is tl. */
    static double occupancy(double[] lam, double[] t, int l, double tl, int h) {
        double[] tt = t.clone();
        tt[l] = tl;
        double occ = 0.0;
        for (int k = 0; k < lam.length; k++) {
            double[] w = levelWeights(lam[k], tt, h);
            double sum = 0.0;
            for (int s = 0; s <= h; s++) {
                sum += w[s];
            }
            occ += w[1 + l] / sum;
        }
        return occ;
    }

    /** Unnormalized level weights: w_l = prod_{s<=l} (1-e_s)/e_s, w_0 = 1. */
    static double[] levelWeights(double lam, double[] t, int h) {
        double[] w = new double[h + 1];
        w[0] = 1.0;
        for (int l = 0; l < h; l++) {
            double e = FastMath.exp(-lam * t[l]);
            w[1 + l] = w[l] * (1 - e) / FastMath.max(e, 1e-300);
        }
        return w;
    }
}
