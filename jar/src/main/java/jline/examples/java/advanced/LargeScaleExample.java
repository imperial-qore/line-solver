/*
 * Copyright (c) 2012-2026, QORE Lab, Imperial College London
 * All rights reserved.
 */
package jline.examples.java.advanced;

import jline.lang.ClosedClass;
import jline.lang.Network;
import jline.lang.RoutingMatrix;
import jline.lang.nodes.Station;
import jline.lang.constant.SchedStrategy;
import jline.lang.nodes.Delay;
import jline.lang.nodes.Queue;
import jline.lang.processes.Erlang;
import jline.solvers.NetworkAvgTable;
import jline.solvers.SolverOptions;
import jline.solvers.fluid.FLD;
import jline.solvers.fluid.SolverFluid;

import java.util.ArrayList;
import java.util.List;
import java.util.Random;

/**
 * Large-scale transient fluid analysis with trajectory-based iteration (TBI).
 *
 * <p>This example builds a vehicle-sharing-style closed network with M queueing
 * stations and one Erlang transit delay per ordered station pair, giving O(M^2)
 * stations and, with 16 Erlang phases per delay, a fluid ODE with about 1500
 * state variables. At this size the monolithic stiff fluid solution (method
 * {@code closing}) takes several minutes, dominated by the Jacobian
 * factorizations, while trajectory-based iteration solves each station cell
 * separately against frozen inbound trajectories and completes in seconds.</p>
 *
 * <p>M. Sheldon, D. Tuncer, G. Casale, "TBI: Transient Hierarchical Modeling of
 * Large-Scale Vehicle Sharing Systems", IEEE Transactions on Intelligent
 * Transportation Systems.</p>
 *
 * <p>Each cell holds one queueing station together with its outbound transit
 * delays, mirroring the spatial submodels of the paper.</p>
 */
public class LargeScaleExample {

    /**
     * Transient fluid solution of an O(M^2)-station closed network by TBI
     * (largescale_tbi.m).
     */
    public static void largescale_tbi() {
        final int M = 10;    // queueing stations; the model has M + M*(M-1) in total
        final int N = 250;   // closed population (vehicles)
        final int kph = 16;  // Erlang phases per transit delay
        final double Tend = 20;

        Network model = new Network("tbi_largescale");
        // the reference seeds MATLAB's twister; the demands here are drawn from
        // the same ranges rather than from the same stream
        Random rng = new Random(1);
        Queue[] Q = new Queue[M];
        for (int i = 0; i < M; i++) {
            Q[i] = new Queue(model, "Q" + (i + 1), SchedStrategy.PS);
        }
        Delay[][] D = new Delay[M][M];
        for (int i = 0; i < M; i++) {
            for (int j = 0; j < M; j++) {
                if (i != j) {
                    D[i][j] = new Delay(model, "D" + (i + 1) + "_" + (j + 1));
                }
            }
        }
        ClosedClass job = new ClosedClass(model, "C1", N, Q[0]);
        double[][] lambda = new double[M][M];
        for (int i = 0; i < M; i++) {
            Q[i].setService(job, Erlang.fitMeanAndOrder(1.0 / (1 + 3 * rng.nextDouble()), 4));
            for (int j = 0; j < M; j++) {
                if (i != j) {
                    D[i][j].setService(job,
                            Erlang.fitMeanAndOrder(0.2 + 2 * rng.nextDouble(), kph));
                    lambda[i][j] = rng.nextDouble();
                }
            }
        }
        RoutingMatrix P = model.initRoutingMatrix();
        for (int i = 0; i < M; i++) {
            double tot = 0;
            for (int j = 0; j < M; j++) {
                tot += lambda[i][j];
            }
            for (int j = 0; j < M; j++) {
                if (i != j) {
                    P.set(job, job, Q[i], D[i][j], lambda[i][j] / tot);
                    P.set(job, job, D[i][j], Q[j], 1.0);
                }
            }
        }
        model.link(P);
        model.initDefault();

        // one cell per queueing station plus its outbound transit delays
        List<Station> stations = model.getStations();
        List<String> names = new ArrayList<String>();
        for (Station s : stations) {
            names.add(s.getName());
        }
        List<int[]> cells = new ArrayList<int[]>();
        for (int i = 0; i < M; i++) {
            int[] idx = new int[M];
            idx[0] = names.indexOf("Q" + (i + 1));
            int k = 1;
            for (int j = 0; j < M; j++) {
                if (i != j) {
                    idx[k++] = names.indexOf("D" + (i + 1) + "_" + (j + 1));
                }
            }
            cells.add(idx);
        }

        SolverOptions options = FLD.defaultOptions();
        options.method = "tbi";
        options.timespan = new double[]{0, Tend};
        options.stiff = true;
        options.config.tbi_cells = cells;
        SolverFluid solver = new SolverFluid(model, options);

        long t0 = System.nanoTime();
        NetworkAvgTable avgTable = solver.getAvgTable();
        double elapsed = (System.nanoTime() - t0) / 1e9;
        System.out.printf("TBI solved %d stations (%d ODE variables) in %.1f seconds.%n",
                model.getNumberOfStations(), M * (M - 1) * kph + M * 4, elapsed);

        // queue-length summary at the queueing stations
        List<Double> qlen = avgTable.getQLen();
        List<String> rows = avgTable.getStationNames();
        System.out.printf("%n%-10s %12s%n", "Station", "QLen");
        for (int i = 0; i < rows.size(); i++) {
            if (rows.get(i).startsWith("Q")) {
                System.out.printf("%-10s %12.5f%n", rows.get(i), qlen.get(i));
            }
        }

        // For comparison, the undecomposed solution of the same model,
        //   options.method = "closing";
        // takes several minutes on the same machine (about 80 seconds already at
        // M=8, and beyond 10 minutes at M=12).
    }

    public static void main(String[] args) {
        largescale_tbi();
    }
}
