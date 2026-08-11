/**
 * @file CTMC reward analyzer using value iteration with uniformization
 *
 * @since LINE 3.0
 */
package jline.solvers.ctmc.analyzers;

import java.util.ArrayList;
import java.util.HashMap;
import java.util.LinkedHashMap;
import java.util.List;
import java.util.Map;

import jline.GlobalConstants;
import jline.api.mc.Ctmc_solve;
import jline.lang.NetworkStruct;
import jline.solvers.SolverOptions;
import jline.solvers.ctmc.handlers.Solver_ctmc;
import jline.solvers.ctmc.ResultCTMC;
import jline.util.matrix.Matrix;

public final class Solver_ctmc_reward {
    private Solver_ctmc_reward() {}

    /**
     * Compute rewards via value iteration on uniformized CTMC.
     */
    public static RewardResult solver_ctmc_reward(NetworkStruct sn, SolverOptions options) {
        long startTime = System.nanoTime();

        if (sn.reward == null || sn.reward.isEmpty()) {
            throw new IllegalStateException(
                    "No rewards defined. Use model.setReward() before calling reward analysis.");
        }

        int Tmax = (options.rewardIterations != null) ? options.rewardIterations : 1000;

        ResultCTMC ctmcResult = Solver_ctmc.solver_ctmc(sn, options);
        Matrix Q = ctmcResult.getQ();
        Matrix stateSpace = ctmcResult.getStateSpace();
        Matrix stateSpaceAggr = ctmcResult.getStateSpaceAggr();

        int nstates = Q.getNumRows();
        Map<String, ?> rewards = sn.reward;
        List<String> rewardNames = new ArrayList<String>(rewards.keySet());

        Map<String, double[]> R = new HashMap<String, double[]>();
        for (String name : rewardNames) {
            Object rewardFn = rewards.get(name);
            double[] rv = new double[nstates];
            for (int s = 0; s < nstates; s++) {
                Matrix stateRow = stateSpaceAggr.getRow(s);
                rv[s] = ((jline.lang.reward.RewardFunction) rewardFn).compute(stateRow, sn);
            }
            R.put(name, rv);
        }

        double q = 0.0;
        for (int i = 0; i < nstates; i++) {
            double diagVal = Math.abs(Q.get(i, i));
            if (diagVal > q) {
                q = diagVal;
            }
        }
        if (q == 0.0) {
            q = 1.0;
        }

        Matrix P = new Matrix(nstates, nstates);
        for (int i = 0; i < nstates; i++) {
            for (int j = 0; j < nstates; j++) {
                double pij = Q.get(i, j) / q;
                if (i == j) {
                    pij += 1.0;
                }
                P.set(i, j, pij);
            }
        }

        Map<String, Matrix> V = new LinkedHashMap<String, Matrix>();
        for (String name : rewardNames) {
            Matrix vMatrix = new Matrix(Tmax + 1, nstates);
            vMatrix.zero();
            double[] r = R.get(name);

            double[] vPrev = new double[nstates];

            for (int k = 0; k < Tmax; k++) {
                double[] vNew = new double[nstates];
                for (int s = 0; s < nstates; s++) {
                    double sum = r[s];
                    for (int sp = 0; sp < nstates; sp++) {
                        sum += P.get(s, sp) * vPrev[sp];
                    }
                    vNew[s] = sum;
                    vMatrix.set(k + 1, s, sum);
                }
                vPrev = vNew;
            }
            V.put(name, vMatrix);
        }

        double[] time = new double[Tmax + 1];
        for (int it = 0; it <= Tmax; it++) {
            time[it] = (double) it / q;
        }

        // Steady-state expected rewards from the stationary distribution:
        // E[r] = sum_s pi(s) * r(s), matching the MATLAB solver_ctmc_reward
        Matrix pi = Ctmc_solve.ctmc_solve(Q);
        double[] piClean = new double[nstates];
        double piSum = 0.0;
        for (int s = 0; s < nstates; s++) {
            double piv = (pi.getNumRows() == 1) ? pi.get(0, s) : pi.get(s, 0);
            if (piv < GlobalConstants.Zero) {
                piv = 0.0;
            }
            piClean[s] = piv;
            piSum += piv;
        }
        Map<String, Double> steadyState = new LinkedHashMap<String, Double>();
        for (String name : rewardNames) {
            double[] r = R.get(name);
            double avgReward = 0.0;
            for (int s = 0; s < nstates; s++) {
                avgReward += (piClean[s] / piSum) * r[s];
            }
            steadyState.put(name, avgReward);
        }

        double runtime = (System.nanoTime() - startTime) / 1e9;

        return new RewardResult(V, time, rewardNames, stateSpaceAggr, steadyState, runtime);
    }
}
