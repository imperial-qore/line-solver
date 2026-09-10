/*
 * Copyright (c) 2012-2026, QORE Lab, Imperial College London
 * All rights reserved.
 */

package jline.api.qsys;

/**
 * Tandem network Lindley recursion on a sample path.
 *
 * <p>Propagates the waiting times of a series of K single-server FCFS stations in
 * tandem, driven by the primitives of the sample path: the interarrival times at
 * the first station and the per-station service times.
 *
 * <p>At the first station this is Lindley's recursion,
 * {@code W(n+1,1) = max(W(n,1) + S(n,1) - A(n), 0)}. Downstream the interarrival
 * time is not a primitive: the arrival epoch of customer n at station k is its
 * departure epoch from station k-1, so the interarrival time at station k is the
 * interdeparture time upstream. Writing {@code G(n,k)} for the interarrival time at
 * station k between customers n and n+1, with {@code G(n,1) = A(n)}, the exact
 * interdeparture identity is
 *
 * <pre>
 *   G(n,k+1) = G(n,k) + W(n+1,k) - W(n,k) + S(n+1,k) - S(n,k),
 * </pre>
 *
 * <p>equivalently and more transparently
 *
 * <pre>
 *   G(n,k+1) = max(G(n,k) - W(n,k) - S(n,k), 0) + S(n+1,k),
 * </pre>
 *
 * <p>an idle period at station k followed by the next customer's service there. The
 * recursion at station k is then
 * {@code W(n+1,k) = max(W(n,k) + S(n,k) - G(n,k), 0)}.
 *
 * <p>Nothing here is distributional, so the recursion is exact for arbitrary
 * interarrival and service times, dependent or not, and is the reference a
 * simulated tandem sample path can be checked against directly. It reproduces a
 * direct event-driven tandem simulation to 1e-12 over four stations.
 *
 * <p>Note that proposition 1 of the reference states this identity without the
 * {@code S(n+1,k) - S(n,k)} term, which makes it wrong as a sample-path identity:
 * the omitted difference has mean zero, so the mean interdeparture time survives,
 * but individual waiting times do not. Implementing it as published gives
 * station-1 waiting times that are correct and downstream ones that are not, by up
 * to several mean service times. The form above is used instead.
 *
 * <p>Port of MATLAB qsys_tandem_lindley.m.
 *
 * <p>Reference: S. Palomo, J. Pender, "Learning the Tandem Network Lindley
 * Recursion", Proc. Winter Simulation Conference, 2021, equations 2 and 3 and
 * proposition 1, the last corrected as described above; D. V. Lindley, "The Theory
 * of Queues with a Single Server", Proc. Camb. Phil. Soc. 48, 1952.
 *
 * @since LINE 3.1.0
 */
public final class Qsys_tandem_lindley {
    private Qsys_tandem_lindley() {}

    /**
     * Propagates the waiting times from an empty network.
     *
     * @param A interarrival times at the first station, {@code A[n]} separating
     *          customers n and n+1
     * @param S service times, {@code S[n][k]} for customer n at station k
     * @return the waiting times and derived epochs
     */
    public static QsysTandemPathResult qsys_tandem_lindley(double[] A, double[][] S) {
        return qsys_tandem_lindley(A, S, null);
    }

    /**
     * Propagates the waiting times from a given initial state.
     *
     * @param A  interarrival times at the first station
     * @param S  service times, {@code S[n][k]} for customer n at station k
     * @param W0 waiting times of customer 1 at each station, null for an empty
     *           network
     * @return the waiting times and derived epochs
     */
    public static QsysTandemPathResult qsys_tandem_lindley(double[] A, double[][] S,
                                                           double[] W0) {
        if (A == null || A.length < 1) {
            throw new IllegalArgumentException("A must hold at least one interarrival time");
        }
        int N = A.length;
        if (S == null || S.length != N) {
            throw new IllegalArgumentException("S must have one row per customer, "
                    + (S == null ? "got null" : "got " + S.length + " rows for " + N
                    + " interarrivals"));
        }
        int K = S[0] == null ? 0 : S[0].length;
        if (K < 1) {
            throw new IllegalArgumentException("S must have at least one column, one per station");
        }
        for (int n = 0; n < N; n++) {
            if (!Double.isFinite(A[n]) || A[n] < 0.0) {
                throw new IllegalArgumentException(
                        "A must hold finite nonnegative interarrival times");
            }
            if (S[n] == null || S[n].length != K) {
                throw new IllegalArgumentException("Every row of S must have " + K + " columns");
            }
            for (int k = 0; k < K; k++) {
                if (!Double.isFinite(S[n][k]) || S[n][k] < 0.0) {
                    throw new IllegalArgumentException(
                            "S must hold finite nonnegative service times");
                }
            }
        }

        double[] start = new double[K];
        if (W0 != null) {
            if (W0.length != K) {
                throw new IllegalArgumentException(
                        "W0 must hold one waiting time per station, " + K + " of them");
            }
            for (int k = 0; k < K; k++) {
                if (!Double.isFinite(W0[k]) || W0[k] < 0.0) {
                    throw new IllegalArgumentException(
                            "W0 must hold finite nonnegative waiting times");
                }
                start[k] = W0[k];
            }
        }

        double[][] W = new double[N][K];
        double[][] G = new double[N][K];
        System.arraycopy(start, 0, W[0], 0, K);
        for (int k = 0; k < K; k++) {
            G[N - 1][k] = Double.NaN;
        }

        for (int n = 0; n < N - 1; n++) {
            double gap = A[n];
            for (int k = 0; k < K; k++) {
                G[n][k] = gap;
                W[n + 1][k] = Math.max(W[n][k] + S[n][k] - gap, 0.0);
                // the interdeparture time here is the interarrival time one station on
                gap = Math.max(gap - W[n][k] - S[n][k], 0.0) + S[n + 1][k];
            }
        }

        double[][] T = new double[N][K];
        double[][] departure = new double[N][K];
        double epoch = 0.0;
        for (int n = 0; n < N; n++) {
            for (int k = 0; k < K; k++) {
                T[n][k] = W[n][k] + S[n][k];
            }
            departure[n][0] = epoch + T[n][0];
            for (int k = 1; k < K; k++) {
                departure[n][k] = departure[n][k - 1] + T[n][k];
            }
            if (n < N - 1) {
                epoch += A[n];
            }
        }

        return new QsysTandemPathResult(W, G, T, departure);
    }
}
