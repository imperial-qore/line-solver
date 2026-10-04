/*
 * Copyright (c) 2012-2026, Imperial College London
 * All rights reserved.
 */
package jline.inference.api;

import jline.lang.Network;
import jline.lang.nodes.Queue;
import jline.util.matrix.Matrix;

public final class Infer_minps {
    private Infer_minps() {}

    /**
     * MINPS demand estimation method.
     *
     * Runs both MLPS and RPS estimators and selects the one with the
     * smaller mean demand estimate.
     *
     * @param model LINE Network model with delay rates set
     * @param node PS queue node
     * @param rt response time samples (column vector)
     * @param classVec class of each sample (0-based)
     * @param ql queue lengths at arrival (n x R matrix)
     * @return estimated demands (1 x R array)
     */
    public static double[] infer_minps(Network model, Queue node, double[] rt, int[] classVec, Matrix ql) {
        int V = node.getNumberOfServers();

        double[] demandEstMLPS = Infer_mlps.infer_mlps(model, node, rt, classVec, ql);
        double[] demandEstRPS = Infer_rps.infer_rps(rt, classVec, ql, V);

        double avgMLPS = 0.0;
        for (double v : demandEstMLPS) {
            avgMLPS += v;
        }
        if (demandEstMLPS.length > 0) {
            avgMLPS /= demandEstMLPS.length;
        }
        double avgRPS = 0.0;
        for (double v : demandEstRPS) {
            avgRPS += v;
        }
        if (demandEstRPS.length > 0) {
            avgRPS /= demandEstRPS.length;
        }

        return avgMLPS < avgRPS ? demandEstMLPS : demandEstRPS;
    }
}
