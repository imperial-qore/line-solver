package jline.api.pfqn.ld;

import java.util.HashMap;
import java.util.List;

import jline.io.Ret;
import jline.lang.constant.SchedStrategy;
import jline.util.matrix.Matrix;

/**
 * Akyildiz-Bolch (A/B) linearizer method for closed product-form queueing networks.
 */
public final class Pfqn_ab {
    private Pfqn_ab() {}

    public static Ret.pfqnAB pfqn_ab(Matrix D, Matrix N, Matrix S, List<SchedStrategy> sched) {
        long startTime = System.nanoTime();
        int M = D.getNumRows();
        int R = D.getNumCols();

        Matrix QN = Matrix.zeros(M, R);
        Matrix UN = Matrix.zeros(M, R);
        Matrix RN = Matrix.zeros(M, R);
        Matrix TN = Matrix.zeros(M, R);
        Matrix CN = Matrix.zeros(M, R);
        Matrix XN = Matrix.zeros(1, R);

        boolean hasMultiServerPS = false;
        for (int ist = 0; ist < M && !hasMultiServerPS; ist++) {
            if (sched.get(ist) == SchedStrategy.PS) {
                for (int r = 0; r < R; r++) {
                    if (S.get(ist, r) > 1) { hasMultiServerPS = true; break; }
                }
            }
        }
        if (!hasMultiServerPS) return fallbackToStandardMVA(D, N, S, sched, startTime);

        int maxIter = 1000;
        double tolerance = 1e-6;
        int iter = 0;

        Matrix Q_prev = Matrix.zeros(M, R);
        for (int ist = 0; ist < M; ist++) {
            for (int r = 0; r < R; r++) Q_prev.set(ist, r, N.get(r) / (double) M);
        }

        while (iter < maxIter) {
            iter++;
            @SuppressWarnings("unchecked") // generic array creation is an inherent Java limitation
            HashMap<Integer, Double>[][] marginalProbs = new HashMap[M][R];
            for (int ist = 0; ist < M; ist++) {
                if (sched.get(ist) == SchedStrategy.PS) {
                    for (int r = 0; r < R; r++) {
                        if (S.get(ist, r) > 1) {
                            marginalProbs[ist][r] = findMarginalProbs(Q_prev.get(ist, r), (int) S.get(ist, r), N, r, "ab");
                        } else {
                            marginalProbs[ist][r] = findMarginalProbs(Q_prev.get(ist, r), 1, N, r, "scat");
                        }
                    }
                }
            }
            for (int r = 0; r < R; r++) {
                Matrix N_reduced = new Matrix(N);
                N_reduced.set(r, N_reduced.get(r) - 1);
                if (N_reduced.get(r) < 0) continue;
                for (int ist = 0; ist < M; ist++) {
                    SchedStrategy s = sched.get(ist);
                    if (s == SchedStrategy.PS) {
                        if (S.get(ist, r) > 1) {
                            CN.set(ist, r, computeMultiServerPSWaitingTime(D.get(ist, r), (int) S.get(ist, r), marginalProbs[ist][r]));
                        } else {
                            CN.set(ist, r, D.get(ist, r) * (1 + Q_prev.get(ist, r)));
                        }
                    } else if (s == SchedStrategy.FCFS) {
                        CN.set(ist, r, D.get(ist, r) * (1 + Q_prev.get(ist, r)));
                    } else if (s == SchedStrategy.INF) {
                        CN.set(ist, r, D.get(ist, r));
                    } else {
                        CN.set(ist, r, D.get(ist, r) * (1 + Q_prev.get(ist, r)));
                    }
                }
                double systemResponseTime = CN.sumCols().get(r);
                XN.set(r, N.get(r) / systemResponseTime);
            }
            Matrix Q_new = Matrix.zeros(M, R);
            for (int ist = 0; ist < M; ist++) {
                for (int r = 0; r < R; r++) Q_new.set(ist, r, XN.get(r) * CN.get(ist, r));
            }
            double maxChange = 0.0;
            for (int ist = 0; ist < M; ist++) {
                for (int r = 0; r < R; r++) {
                    double change = Math.abs(Q_new.get(ist, r) - Q_prev.get(ist, r));
                    if (change > maxChange) maxChange = change;
                }
            }
            Q_prev = Q_new;
            if (maxChange < tolerance) break;
        }

        QN = Q_prev;
        for (int ist = 0; ist < M; ist++) {
            for (int r = 0; r < R; r++) {
                if (sched.get(ist) == SchedStrategy.PS && S.get(ist, r) > 1) {
                    UN.set(ist, r, XN.get(r) * D.get(ist, r) / S.get(ist, r));
                } else {
                    UN.set(ist, r, XN.get(r) * D.get(ist, r));
                }
            }
        }
        for (int ist = 0; ist < M; ist++) {
            for (int r = 0; r < R; r++) RN.set(ist, r, XN.get(r) > 0 ? QN.get(ist, r) / XN.get(r) : 0.0);
        }
        for (int ist = 0; ist < M; ist++) {
            for (int r = 0; r < R; r++) TN.set(ist, r, XN.get(r));
        }
        double runtime = (System.nanoTime() - startTime) / 1e9;
        return new Ret.pfqnAB(QN, UN, RN, TN, CN, XN, "ab", iter, runtime);
    }

    private static Matrix weightFun(Matrix population, double alpha, double beta) {
        int maxClassPopulation = 0;
        for (int i = 0; i < population.length(); i++) {
            if (population.get(i) > maxClassPopulation) maxClassPopulation = (int) population.get(i);
        }
        double[] scalingFun = new double[maxClassPopulation + 1];
        if (maxClassPopulation >= 1) {
            scalingFun[1] = alpha;
            for (int n = 2; n <= maxClassPopulation; n++) scalingFun[n] = beta * scalingFun[n - 1];
        }
        double[][] weightFun = new double[maxClassPopulation + 1][maxClassPopulation + 1];
        weightFun[0][0] = 1.0;
        for (int l = 1; l <= maxClassPopulation; l++) {
            for (int j = 0; j < l; j++) weightFun[l][j] = weightFun[l - 1][j] - (weightFun[l - 1][j] * scalingFun[l]) / 100.0;
            double sum = 0.0;
            for (int j = 0; j < l; j++) sum += weightFun[l][j];
            weightFun[l][l] = 1.0 - sum;
        }
        return new Matrix(weightFun);
    }

    private static HashMap<Integer, Double> findMarginalProbs(double avgJobs, int numServers, Matrix population, int classIdx, String method) {
        HashMap<Integer, Double> marginalProbs = new HashMap<Integer, Double>();
        if ("scat".equalsIgnoreCase(method)) {
            int floorV = (int) Math.floor(avgJobs);
            int ceilV = (int) Math.ceil(avgJobs);
            marginalProbs.put(floorV, ceilV - avgJobs);
            marginalProbs.put(ceilV, avgJobs - floorV);
            return marginalProbs;
        }
        double ALPHA = 45.0;
        double BETA = 0.7;
        Matrix w = weightFun(population, ALPHA, BETA);
        int floorV = (int) Math.floor(avgJobs);
        int ceiling = floorV + 1;
        int maxVal = Math.min((2 * floorV) + 1, numServers - 2);
        for (int j = 0; j <= maxVal; j++) {
            double prob;
            if (j <= floorV) {
                int lDist = floorV - j;
                int lowerVal = floorV - lDist;
                int upperVal = ceiling + lDist;
                if (lDist > 25) prob = 0.0;
                else if (floorV < population.get(classIdx)) {
                    prob = w.get(floorV, lDist) * ((upperVal - avgJobs) / (upperVal - lowerVal));
                } else prob = 0.0;
            } else {
                int uDist = j - ceiling;
                int lowerVal = floorV - uDist;
                int upperVal = ceiling + uDist;
                if (uDist > 25) prob = 0.0;
                else if (ceiling < population.get(classIdx)) {
                    prob = w.get(ceiling, uDist) * ((avgJobs - lowerVal) / (upperVal - lowerVal));
                } else prob = 0.0;
            }
            marginalProbs.put(j, prob);
        }
        return marginalProbs;
    }

    private static double computeMultiServerPSWaitingTime(double demand, int numServers, HashMap<Integer, Double> marginalProbs) {
        double expectedWaitingTime = 0.0;
        for (java.util.Map.Entry<Integer, Double> e : marginalProbs.entrySet()) {
            int jobs = e.getKey();
            double prob = e.getValue();
            double serviceRate = Math.min((double) jobs, (double) numServers);
            if (serviceRate > 0) expectedWaitingTime += prob * demand * (1 + jobs / serviceRate);
            else expectedWaitingTime += prob * demand;
        }
        return expectedWaitingTime;
    }

    private static Ret.pfqnAB fallbackToStandardMVA(Matrix D, Matrix N, Matrix S, List<SchedStrategy> sched, long startTime) {
        int M = D.getNumRows();
        int R = D.getNumCols();
        Matrix QN = Matrix.zeros(M, R);
        Matrix UN = Matrix.zeros(M, R);
        Matrix RN = Matrix.zeros(M, R);
        Matrix TN = Matrix.zeros(M, R);
        Matrix CN = Matrix.zeros(M, R);
        Matrix XN = Matrix.zeros(1, R);
        for (int r = 0; r < R; r++) {
            double totalDemand = D.sumCols().get(r);
            XN.set(r, N.get(r) / totalDemand);
            for (int ist = 0; ist < M; ist++) {
                CN.set(ist, r, D.get(ist, r));
                QN.set(ist, r, XN.get(r) * CN.get(ist, r));
                UN.set(ist, r, XN.get(r) * D.get(ist, r));
                RN.set(ist, r, CN.get(ist, r));
                TN.set(ist, r, XN.get(r));
            }
        }
        double runtime = (System.nanoTime() - startTime) / 1e9;
        return new Ret.pfqnAB(QN, UN, RN, TN, CN, XN, "ab-fallback", 1, runtime);
    }

    /** PFQN ab algorithms. */
    public static final class PfqnAbAlgo {}
}
