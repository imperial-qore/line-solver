/**
 * Exact per-class busy-server fraction read off the enumerated CTMC state space.
 *
 * For a class that a G-network signal can annihilate, the departure-based
 * estimator T*E[S]/c is exact only under exponential service: a job destroyed
 * mid-service leaves behind busy time with no completion, so with phase-type
 * service T*E[S]/c under-counts (M/Er2/1 with lambda+ = 0.5, lambda- = 0.4
 * gives 0.34941 against a true 0.37696). The in-service occupancy computed here
 * is exact for any service process and coincides with T*E[S]/c when service is
 * exponential.
 *
 * PS-like disciplines share the servers among all resident jobs, so class k gets
 * the weighted share n_k w_k / sum_j n_j w_j of the busy servers; the remaining
 * disciplines expose the in-service indicator directly through ToMarginal.
 *
 * Mirrors matlab/src/solvers/CTMC/ctmc_signal_busy.m.
 *
 * @since LINE 3.0
 */
package jline.solvers.ctmc.handlers;

import jline.lang.NetworkStruct;
import jline.lang.constant.SchedStrategy;
import jline.lang.state.State;
import jline.lang.state.ToMarginal;
import jline.util.matrix.Matrix;

public final class CtmcSignalBusy {
    private CtmcSignalBusy() {}

    /**
     * @param sn            network structure
     * @param ind           node index of the station
     * @param ist           station index
     * @param sched         scheduling strategy at the station
     * @param nservers      number of servers at the station
     * @param stateSpace    working state space (rows indexed by wset entries)
     * @param istSpaceShift per-station column offset into the state space
     * @param wset          reachable state rows
     * @param probSysState  stationary probability of each state row
     * @param K             number of classes
     * @return per-class busy-server fraction (length K)
     */
    public static double[] busyFraction(NetworkStruct sn, int ind, int ist, SchedStrategy sched,
                                        double nservers, Matrix stateSpace, Matrix istSpaceShift,
                                        Matrix wset, Matrix probSysState, int K) {
        double[] UNb = new double[K];
        boolean isPS = sched == SchedStrategy.PS || sched == SchedStrategy.DPS
                || sched == SchedStrategy.GPS || sched == SchedStrategy.LPS;
        int colStart = (int) istSpaceShift.get(ist);
        int colEnd = colStart + sn.space.get(sn.stateful.get(ist)).getNumCols();
        for (int index = 0; index < wset.getNumCols(); index++) {
            int st = (int) wset.get(index);
            double pst = probSysState.get(st);
            if (pst == 0) {
                continue;
            }
            State.StateMarginalStatistics stats = ToMarginal.toMarginal(sn, ind,
                    Matrix.extract(stateSpace, st, st + 1, colStart, colEnd),
                    null, null, null, null, null);
            Matrix ni = stats.ni;
            Matrix nir = stats.nir;
            Matrix sir = stats.sir;
            double totJobs = 0.0;
            for (int j = 0; j < ni.length(); j++) {
                totJobs += ni.get(j);
            }
            if (totJobs <= 0) {
                continue;
            }
            if (isPS) {
                double wtot = 0.0;
                for (int k = 0; k < K; k++) {
                    wtot += nir.get(k) * sn.schedparam.get(ist, k);
                }
                if (wtot > 0) {
                    double busy = Math.min(totJobs, nservers) / nservers;
                    for (int k = 0; k < K; k++) {
                        UNb[k] += pst * (nir.get(k) * sn.schedparam.get(ist, k) / wtot) * busy;
                    }
                }
            } else {
                for (int k = 0; k < K; k++) {
                    UNb[k] += pst * sir.get(k) / nservers;
                }
            }
        }
        return UNb;
    }
}
