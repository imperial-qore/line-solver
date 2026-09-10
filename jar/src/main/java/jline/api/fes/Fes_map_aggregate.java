/**
 * @file Recursive MAP flow-equivalent server for a station subset
 *
 * @since LINE 3.0
 */
package jline.api.fes;

import java.util.ArrayList;
import java.util.List;

import jline.api.mam.Map2_fit_idc;
import jline.api.mam.Map_idc;
import jline.api.mam.Map_moment;
import jline.io.Ret;
import jline.util.matrix.MatrixCell;

/**
 * Aggregates a subnetwork into a load-dependent MAP flow-equivalent server that reproduces
 * mean, variability and burstiness of its output.
 *
 * Implements the recursion of Section 5.2.1 of Casale, Mi, Cherkasova and Smirni, IEEE
 * Trans. Soft. Eng. 37(5), 2011. The first station seeds the flow-equivalent server; every
 * further station is folded against the running server by building the inter-departure MAP
 * of the resulting pair at each population level and fitting a MAP(2) to its first three
 * moments and index of dispersion. Unlike the classic flow-equivalent server, which keeps
 * only the mean throughput, this one also carries the burstiness of the departure stream,
 * so a bottleneck switch across the aggregated resources stays visible to the rest of the
 * model.
 *
 * Levels are evaluated on a grid and the four descriptors are interpolated between grid
 * points. The MAP is refitted at every level from the interpolated descriptors, never
 * interpolated entrywise.
 *
 * Copyright (c) 2012-2026, Imperial College London
 * All rights reserved.
 */
public final class Fes_map_aggregate {
    private Fes_map_aggregate() {}

    /**
     * Aggregates a subnetwork on the default population grid.
     *
     * @param maps    service process of each station, already scaled by its visit ratio
     * @param servers number of servers of each station, infinite for a delay
     * @param n       largest population the flow-equivalent server must serve
     * @return the load-dependent MAP and the descriptors it was fitted from
     */
    public static FesMapAggregateResult fes_map_aggregate(List<MatrixCell> maps, double[] servers, int n) {
        return fes_map_aggregate(maps, servers, n, Fes_map_grid.fes_map_grid(n), "ssolve");
    }

    /**
     * Aggregates a subnetwork.
     *
     * @param maps    service process of each station, already scaled by its visit ratio
     * @param servers number of servers of each station, infinite for a delay
     * @param n       largest population the flow-equivalent server must serve
     * @param grid    populations at which the inter-departure MAP is evaluated
     * @param method  moment evaluation method, "ssolve" or "euler"
     * @return the load-dependent MAP and the descriptors it was fitted from
     */
    public static FesMapAggregateResult fes_map_aggregate(List<MatrixCell> maps, double[] servers, int n,
                                                          int[] grid, String method) {
        int M = maps.size();
        if (M < 1) {
            throw new IllegalArgumentException("At least one station is required.");
        }
        if (servers.length != M) {
            throw new IllegalArgumentException("One server count per station is required.");
        }

        List<MatrixCell> fes = Fes_map_levels.fes_map_levels(maps.get(0), n, servers[0]);
        double[][] moments = new double[4][n];
        int[] status = new int[n];
        for (int k = 1; k <= n; k++) {
            MatrixCell f = fes.get(k - 1);
            moments[0][k - 1] = Map_moment.map_moment(f, 1);
            moments[1][k - 1] = Map_moment.map_moment(f, 2);
            moments[2][k - 1] = Map_moment.map_moment(f, 3);
            moments[3][k - 1] = Map_idc.map_idc(f);
        }

        for (int i = 1; i < M; i++) {
            double[][] gmom = new double[4][grid.length];
            List<MatrixCell> stationLev = Fes_map_levels.fes_map_levels(maps.get(i), n, servers[i]);
            for (int g = 0; g < grid.length; g++) {
                int k = grid[g];
                MatrixCell T = Fes_map_interdeparture.fes_map_interdeparture(stationLev, fes, k);
                FesMapMomentsResult mom = Fes_map_moments.fes_map_moments(T.get(0), T.get(1), method);
                gmom[0][g] = mom.e1;
                gmom[1][g] = mom.e2;
                gmom[2][g] = mom.e3;
                gmom[3][g] = mom.idc;
            }

            if (grid.length < n) {
                double[] xs = new double[grid.length];
                for (int g = 0; g < grid.length; g++) {
                    xs[g] = grid[g];
                }
                double[] xq = new double[n];
                for (int k = 0; k < n; k++) {
                    xq[k] = k + 1;
                }
                for (int r = 0; r < 4; r++) {
                    moments[r] = Fes_map_interp.fes_map_interp(xs, gmom[r], xq);
                }
            } else {
                moments = gmom;
            }

            List<MatrixCell> newFes = new ArrayList<MatrixCell>(n);
            for (int k = 0; k < n; k++) {
                Ret.mamMAPFitIdcReturn fit = Map2_fit_idc.map2_fit_idc(
                        moments[0][k], moments[1][k], moments[2][k], moments[3][k]);
                newFes.add(fit.MAP);
                status[k] = fit.status;
            }
            fes = newFes;
        }

        double[] throughput = new double[n];
        for (int k = 0; k < n; k++) {
            throughput[k] = 1.0 / moments[0][k];
        }
        return new FesMapAggregateResult(fes, throughput, moments, status, grid);
    }
}
