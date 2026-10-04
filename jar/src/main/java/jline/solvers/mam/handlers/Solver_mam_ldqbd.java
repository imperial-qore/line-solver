/*
 * Copyright (c) 2012-2026, QORE Lab, Imperial College London
 * All rights reserved.
 */
package jline.solvers.mam.handlers;

import jline.api.mam.Ldqbd;
import jline.api.mam.Qbd_setupdelayoff_closed;
import jline.api.mam.LdqbdOptions;
import jline.api.mam.LdqbdResult;
import jline.lang.NetworkStruct;
import jline.solvers.SolverOptions;
import jline.solvers.mam.MAMResult;
import jline.util.matrix.Matrix;

/**
 * Solver for single-class two-station queueing networks using a Level-Dependent QBD.
 *
 * <p>Two regimes share one block-tridiagonal generator, differing only in the
 * per-level arrival rate and the top level:
 * <ul>
 *   <li>CLOSED: one Delay (INF) + one Queue (FCFS), finite population N. Level n
 *       is the number of jobs at the queue, and the arrival rate out of level n
 *       is the finite-source rate (N-n)*lambda_eff, which vanishes at n = N and
 *       closes the chain by itself.</li>
 *   <li>OPEN: one Source (EXT) + one Queue, Poisson arrivals at a constant
 *       lambda_eff, truncated at options.cutoff or at a level where the tail
 *       probability is negligible.</li>
 * </ul>
 *
 * <p>Exactness: exact for exponential service at any number of servers, and for
 * PH service at any number of servers. The multiserver PH chain carries the
 * MULTISET of the phases the min(n,c) busy servers sit in; the collapsed
 * single-phase approximation this solver used until 2026-08-18 is gone.
 *
 * <p>The blocks come from {@link Solver_mam_ldqbd_statevec#solver_mam_ldqbd_ld},
 * the one construction shared with the SolverENV state-vector analyzer. This
 * class kept a second, closed-only copy until 2026-08-18, which silently ignored
 * sn.lldscaling (though the featset declared LoadDependence for this method) and
 * refused the open regime that MATLAB, Python and C++ all served.
 */
public final class Solver_mam_ldqbd {
    private Solver_mam_ldqbd() {}

    public static MAMResult solver_mam_ldqbd(NetworkStruct sn, SolverOptions options) {
        final int M = sn.nstations;
        final int K = sn.nclasses;

        // Shape validation, the service process, the per-level arrival rate and the
        // per-level service factor all live in the shared builder.
        Solver_mam_ldqbd_statevec.Ld ld = Solver_mam_ldqbd_statevec.solver_mam_ldqbd_ld(sn, options);

        final int Nlev = ld.Nlev;
        Matrix pi_ldqbd = null;
        double mean_queue;
        double x_setup = 0.0;
        if (ld.hasSetup) {
            // SETUP AND DELAY-OFF, the closed vacation queue. The level-dependent
            // chain this needs is the one the builder assembled with two extra
            // phase families -- the setup above level 0 and the delay-off at
            // level 0 -- and Qbd_setupdelayoff_closed builds and solves exactly
            // that, so it is called rather than duplicated. Without it the blocks
            // describe a server that is ALWAYS warm and the answer is
            // byte-identical across any setup mean (BUG-78).
            Qbd_setupdelayoff_closed.Result cr =
                    Qbd_setupdelayoff_closed.qbd_setupdelayoff_closed(
                            ld.N, 1.0 / ld.lambda_eff, ld.mu, ld.alpharate, ld.alphascv,
                            ld.betarate, ld.betascv);
            mean_queue = cr.QN;
            x_setup = cr.XN;
        } else {
            LdqbdOptions ldqbdOptions = new LdqbdOptions(options.tol, options.iter_max, false);
            LdqbdResult ldqbdResult = Ldqbd.ldqbd(ld.Q0, ld.Q1, ld.Q2, ldqbdOptions);
            pi_ldqbd = ldqbdResult.getPi();   // per-level, already summed over phases

            mean_queue = 0.0;
            for (int n = 0; n <= Nlev; n++) {
                mean_queue += n * pi_ldqbd.get(0, n);
            }
        }

        // Utilization is the fraction of the station's PEAK capacity in use,
        // sum_n pi(n)*sf(n)/utilPeak, which is the work-based convention CTMC,
        // MVA, NC and serial SSA all report; utilPeak = max(c, max(alpha)) is
        // CTMC's own normalizer. Already per-server: do not rescale by c again.
        //
        // Without load dependence sf(n) = min(n,c) and utilPeak = c, giving the
        // mean fraction of the c servers in use, and at c = 1 that is sf(n) = 1
        // for every n >= 1, so the sum collapses to 1 - pi(0). It used to report
        // 1 - pi(0) under load dependence, i.e. P(busy), which reads a station
        // running alpha(n) times faster as no busier than one at its nominal
        // rate: 0.9587 against CTMC's 0.6612 on a 4-job closed model with
        // alpha = [1 1.5 2 2.5].
        double util_queue = 0.0;
        if (ld.hasSetup) {
            // With a setup the server is DELIVERING work only in the busy phase,
            // so the level occupancy over-counts it: a level is occupied during
            // the setup too. The utilization law gives the same work-based number
            // without the per-phase vector, X*E[S]/peak, which is what the level
            // sum reduces to without a setup.
            util_queue = x_setup * ld.mean_service / ld.utilPeak;
        } else {
            for (int n = 1; n <= Nlev; n++) {
                util_queue += (ld.sf[n - 1] / ld.utilPeak) * pi_ldqbd.get(0, n);
            }
        }

        MAMResult result = new MAMResult();
        result.QN = new Matrix(M, K);
        result.UN = new Matrix(M, K);
        result.RN = new Matrix(M, K);
        result.TN = new Matrix(M, K);
        result.CN = new Matrix(1, K);
        result.XN = new Matrix(1, K);

        final int qi = ld.queueIdx;
        final int ri = ld.refIdx;
        if (ld.isOpen) {
            // Served throughput = arrival rate less the truncation blocking.
            double X = ld.lambda_eff * (1.0 - pi_ldqbd.get(0, Nlev));
            double R_queue = (X > 0) ? mean_queue / X : 0.0;
            // Source station: pass-through, no queueing.
            result.TN.set(ri, 0, X);
            result.QN.set(qi, 0, mean_queue);
            result.UN.set(qi, 0, util_queue);
            result.RN.set(qi, 0, R_queue);
            result.TN.set(qi, 0, X);
            result.XN.set(0, 0, X);
            result.CN.set(0, 0, R_queue);
        } else {
            double mean_delay = ld.N - mean_queue;
            double X = mean_delay * ld.lambda_eff;
            double R_queue = (X > 0) ? mean_queue / X : 0.0;
            double R_delay = 1.0 / ld.delayRate;
            // see _kb/06-solver-catalog.md for rationale: the delay completes at
            // mean_delay*lambda_d, of which only the rt(delay,queue) fraction
            // proceeds to the queue, so its throughput is NOT the queue flow X.
            result.QN.set(ri, 0, mean_delay);
            result.UN.set(ri, 0, mean_delay);   // infinite server: U = Q
            result.RN.set(ri, 0, R_delay);
            result.TN.set(ri, 0, mean_delay * ld.delayRate);
            result.QN.set(qi, 0, mean_queue);
            result.UN.set(qi, 0, util_queue);
            result.RN.set(qi, 0, R_queue);
            result.TN.set(qi, 0, X);
            result.XN.set(0, 0, X);
            result.CN.set(0, 0, R_delay + R_queue);
        }

        result.iter = 1;   // LDQBD is a direct method
        result.method = "ldqbd";
        return result;
    }
}
