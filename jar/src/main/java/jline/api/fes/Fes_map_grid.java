/**
 * @file Population grid of a MAP flow-equivalent server
 *
 * @since LINE 3.0
 */
package jline.api.fes;

/**
 * Returns the population levels at which the inter-departure MAP is evaluated.
 *
 * Fitting one MAP per population is wasteful because the processes of neighbouring
 * populations are similar. Section 5.2.2 of Casale, Mi, Cherkasova and Smirni, IEEE Trans.
 * Soft. Eng. 37(5), 2011, evaluates the first ten populations and ten further equispaced
 * points, which is what this function returns.
 *
 * Copyright (c) 2012-2026, Imperial College London
 * All rights reserved.
 */
public final class Fes_map_grid {
    private Fes_map_grid() {}

    /** Leading populations kept in full. */
    public static final int NHEAD = 10;
    /** Equispaced points above them. */
    public static final int NTAIL = 10;

    /**
     * Returns the default grid for a largest population n.
     *
     * @param n largest population
     * @return sorted populations to evaluate
     */
    public static int[] fes_map_grid(int n) {
        return fes_map_grid(n, NHEAD, NTAIL);
    }

    /**
     * Returns the grid for a largest population n.
     *
     * @param n     largest population
     * @param nhead leading populations kept in full
     * @param ntail equispaced points above them
     * @return sorted populations to evaluate
     */
    public static int[] fes_map_grid(int n, int nhead, int ntail) {
        if (n <= nhead + ntail) {
            int[] all = new int[n];
            for (int i = 0; i < n; i++) {
                all[i] = i + 1;
            }
            return all;
        }
        boolean[] taken = new boolean[n + 1];
        for (int i = 1; i <= nhead; i++) {
            taken[i] = true;
        }
        for (int j = 0; j < ntail; j++) {
            double t = (ntail == 1) ? n : (nhead + 1) + (double) j * (n - nhead - 1) / (ntail - 1);
            taken[(int) Math.round(t)] = true;
        }
        int count = 0;
        for (int i = 1; i <= n; i++) {
            if (taken[i]) {
                count++;
            }
        }
        int[] grid = new int[count];
        int p = 0;
        for (int i = 1; i <= n; i++) {
            if (taken[i]) {
                grid[p++] = i;
            }
        }
        return grid;
    }
}
