/**
 * @file Per-station metrics behind a MAP flow-equivalent server
 *
 * @since LINE 3.0
 */
package jline.api.fes;

import jline.GlobalConstants;
import jline.api.pfqn.mva.Pfqn_mvams;
import jline.io.Ret;
import jline.util.matrix.Matrix;

/**
 * Recovers the per-station metrics of an aggregated subnetwork by conditioning on the
 * population held by the flow-equivalent server.
 *
 * fes_map_solve returns the distribution pk of the jobs held by the aggregate. The metrics
 * of the stations behind it follow by conditioning, E[Y_i] = sum_k pk(k) Y_i(k), with
 * Y_i(k) the metric of station i when the isolated subnetwork holds k jobs. This is the
 * decomposition step of the hierarchical analysis of Chandy, Herzog and Woo, IBM J. Res.
 * Dev. 19(1), 1975, and it is exact for a product-form subnetwork. It is an approximation
 * whenever the burstiness that the MAP flow-equivalent server carries also matters inside
 * the subnetwork, because the conditional solve is the product-form one; the aggregate
 * metrics returned by fes_map_solve do not rely on it.
 *
 * Copyright (c) 2012-2026, Imperial College London
 * All rights reserved.
 */
public final class Fes_map_deaggregate {
    private Fes_map_deaggregate() {}

    /**
     * @param pk      distribution of the jobs held by the aggregate, index k = P(k jobs)
     * @param L       service demands of the isolated subnetwork
     * @param mi      servers per station, infinite for a delay
     * @param isDelay true where the station is a pure delay
     * @return per-station queue length, utilization, throughput and residence time
     */
    public static FesMapDeaggregateResult fes_map_deaggregate(double[] pk, double[] L, double[] mi,
                                                              boolean[] isDelay) {
        int M = L.length;
        int n = pk.length - 1;
        double[] QN = new double[M];
        double[] UN = new double[M];
        double[] XN = new double[M];
        double[] RN = new double[M];

        int nq = 0;
        double Z = 0;
        for (int i = 0; i < M; i++) {
            if (isDelay[i]) {
                Z += L[i];
            } else {
                nq++;
            }
        }
        int[] queueIdx = new int[nq];
        Matrix Lq = new Matrix(nq, 1);
        Matrix miq = new Matrix(1, nq);
        int p = 0;
        for (int i = 0; i < M; i++) {
            if (!isDelay[i]) {
                queueIdx[p] = i;
                Lq.set(p, 0, L[i]);
                miq.set(0, p, mi[i]);
                p++;
            }
        }
        Matrix Zm = new Matrix(1, 1);
        Zm.set(0, 0, Z);

        for (int k = 1; k <= n; k++) {
            if (pk[k] <= GlobalConstants.Zero) {
                continue;
            }
            Matrix N = new Matrix(1, 1);
            N.set(0, 0, k);
            // Pfqn_mva's `mi` is NOT a server count -- it enters only as the additive
            // term of C(i,s)=L(i,s)*(mi(i)+Qarv), so passing the real multiplicity
            // INFLATES the residence time instead of adding servers. Pfqn_mvams
            // forwards to Pfqn_mva when every station is a single server and to the
            // load-dependent recursion with mu(i,n)=min(n,S(i)) when one is not. Its
            // U is per STATION on that branch and reports 1-P(0) rather than the
            // [0,1] load, so utilization is recomputed here as U=X*L/S.
            // See _kb/03-api-layer.md (pfqn_mva: mi is not S).
            Matrix lambdaZero = new Matrix(1, 1);
            Ret.pfqnMVA out = Pfqn_mvams.pfqn_mvams(lambdaZero, Lq, N, Zm, null, miq);
            double Xk = out.X.get(0);
            for (int j = 0; j < nq; j++) {
                double srv = Math.max(1.0, miq.get(0, j));
                QN[queueIdx[j]] += pk[k] * out.Q.get(j, 0);
                UN[queueIdx[j]] += pk[k] * (Xk * Lq.get(j, 0) / srv);
                XN[queueIdx[j]] += pk[k] * Xk;
            }
            for (int i = 0; i < M; i++) {
                if (isDelay[i]) {
                    QN[i] += pk[k] * Xk * L[i];
                    UN[i] += pk[k] * Xk * L[i];
                    XN[i] += pk[k] * Xk;
                }
            }
        }

        for (int i = 0; i < M; i++) {
            RN[i] = (XN[i] > GlobalConstants.Zero) ? QN[i] / XN[i] : 0;
        }
        return new FesMapDeaggregateResult(QN, UN, XN, RN);
    }
}
