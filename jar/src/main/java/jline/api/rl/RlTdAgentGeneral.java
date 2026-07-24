/**
 * General Temporal Difference Learning Agent for Queueing Control.
 */
package jline.api.rl;

import jline.util.matrix.Matrix;
import jline.io.Ret.SampleResult;

import java.util.ArrayList;
import java.util.HashMap;
import java.util.Iterator;
import java.util.Map;
import java.util.Random;

/**
 * General TD learning agent for queueing network control.
 */
public class RlTdAgentGeneral {
    public final double lr;
    public double epsilon;
    public final double epsDecay;

    public double[] v = new double[]{0.0};
    public int[] vSize = new int[]{0};

    private final Random random = new Random();

    public RlTdAgentGeneral() {
        this(0.1, 1.0, 0.9999);
    }

    public RlTdAgentGeneral(double lr) {
        this(lr, 1.0, 0.9999);
    }

    public RlTdAgentGeneral(double lr, double epsilon) {
        this(lr, epsilon, 0.9999);
    }

    public RlTdAgentGeneral(double lr, double epsilon, double epsDecay) {
        this.lr = lr;
        this.epsilon = epsilon;
        this.epsDecay = epsDecay;
    }

    public void reset(RlEnvGeneral env) {
        v = new double[]{0.0};
        vSize = new int[]{0};
        env.reset();
    }

    public double[] getValueFunction() {
        return v.clone();
    }

    public double[] solveForFixedPolicy(RlEnvGeneral env, int numEpisodes) {
        reset(env);

        int nqueues = env.nqueues;
        int dimSize = env.stateSize + 1;
        vSize = new int[nqueues];
        for (int i = 0; i < nqueues; i++) vSize[i] = dimSize;
        int totalSize = 1;
        for (int dim : vSize) totalSize *= dim;
        v = new double[totalSize];

        double t = 0.0;
        double c = 0.0;
        double bigT = 0.0;
        double bigC = 0.0;
        int[] x = new int[nqueues];
        int[] n = new int[nqueues];

        int j = 0;
        while (j < numEpisodes) {
            if (j % 1000 == 0) {
                System.out.printf("running episode #%d%n", j);
            }

            GeneralSampleEvent sampleEvent = env.sample();
            double dt = sampleEvent.dt;
            int depNode = sampleEvent.depNode;
            int arvNode = sampleEvent.arvNode;
            SampleResult sampleResult = sampleEvent.sampleResult;
            t += dt;

            int sumX = 0;
            for (int xi : x) sumX += xi;
            c += sumX * dt;

            if (contains(env.idxOfQueueInNodes, depNode)) {
                int depServer = indexOf(env.idxOfQueueInNodes, depNode);
                if (depServer >= 0) {
                    x[depServer] = x[depServer] - 1;
                }
            }

            if (contains(env.idxOfQueueInNodes, arvNode)) {
                int arvServer = indexOf(env.idxOfQueueInNodes, arvNode);
                if (arvServer >= 0) {
                    x[arvServer] = x[arvServer] + 1;
                }
            }

            env.update(sampleResult);

            if (env.isInStateSpace(x)) {
                j++;
                bigT = env.gamma * bigT + t;
                bigC = env.gamma * bigC + c;
                double meanCostRate = (bigT > 0) ? bigC / bigT : 0.0;

                int[] prevLoc = new int[nqueues];
                int[] curLoc = new int[nqueues];
                for (int i = 0; i < nqueues; i++) {
                    prevLoc[i] = n[i] + 1;
                    curLoc[i] = x[i] + 1;
                }

                int prevState = RlTdAgent.getStateFromLoc(vSize, prevLoc);
                int curState = RlTdAgent.getStateFromLoc(vSize, curLoc);

                if (prevState >= 0 && prevState < v.length && curState >= 0 && curState < v.length) {
                    v[prevState] = (1 - lr) * v[prevState] + lr * (c - t * meanCostRate + v[curState]);

                    double v0 = v[0];
                    for (int i = 0; i < v.length; i++) {
                        v[i] -= v0;
                    }
                }

                t = 0.0;
                c = 0.0;
                n = x.clone();
            }
        }

        return v.clone();
    }

    public double[] solve(RlEnvGeneral env, int numEpisodes) {
        reset(env);

        int nqueues = env.nqueues;
        int dimSize = env.stateSize + 1;
        vSize = new int[nqueues];
        for (int i = 0; i < nqueues; i++) vSize[i] = dimSize;
        int totalSize = 1;
        for (int dim : vSize) totalSize *= dim;
        v = new double[totalSize];

        double t = 0.0;
        double c = 0.0;
        double bigT = 0.0;
        double bigC = 0.0;
        int[] x = new int[nqueues];
        int[] n = new int[nqueues];

        double eps = epsilon;

        int j = 0;
        while (j < numEpisodes) {
            if (j % 1000 == 0) {
                System.out.printf("running episode #%d.%n", j);
            }

            eps *= epsDecay;

            GeneralSampleEvent sampleEvent = env.sample();
            double dt = sampleEvent.dt;
            int depNode = sampleEvent.depNode;
            int arvNode = sampleEvent.arvNode;
            SampleResult sampleResult = sampleEvent.sampleResult;
            t += dt;

            int sumX = 0;
            for (int xi : x) sumX += xi;
            c += sumX * dt;

            if (contains(env.idxOfQueueInNodes, depNode)) {
                int depServer = indexOf(env.idxOfQueueInNodes, depNode);
                if (depServer >= 0) {
                    x[depServer] = Math.max(0, x[depServer] - 1);
                }
            }

            if (contains(env.idxOfActionNodes, depNode) && env.isInActionSpace(x)) {
                int[] actions = env.actionSpace.get(depNode);
                if (actions != null && actions.length != 0) {
                    double[] nextValues = genNextValues(env, x, actions);
                    double[] policy = RlTdAgent.createGreedyPolicy(nextValues, eps, actions.length);

                    double r = random.nextDouble();
                    double cumSum = 0.0;
                    int selectedIdx = actions.length - 1;
                    for (int a = 0; a < policy.length; a++) {
                        cumSum += policy[a];
                        if (r < cumSum) {
                            selectedIdx = a;
                            break;
                        }
                    }
                    arvNode = actions[selectedIdx];

                    Matrix eventMatrix = ((SampleResult) sampleResult).event;
                    if (eventMatrix != null && eventMatrix.getNumRows() > 1) {
                        for (int row = 0; row < eventMatrix.getNumRows(); row++) {
                            if (row == 1 || (row > 0 && eventMatrix.getNumCols() > 1)) {
                                eventMatrix.set(row, 1, (double) arvNode);
                                break;
                            }
                        }
                    }
                }
            }

            if (contains(env.idxOfQueueInNodes, arvNode)) {
                int arvServer = indexOf(env.idxOfQueueInNodes, arvNode);
                if (arvServer >= 0) {
                    x[arvServer] = x[arvServer] + 1;
                }
            }
            env.update(sampleResult);

            if (env.isInStateSpace(x)) {
                j++;
                bigT = env.gamma * bigT + t;
                bigC = env.gamma * bigC + c;
                double meanCostRate = (bigT > 0) ? bigC / bigT : 0.0;

                int[] prevLoc = new int[nqueues];
                int[] curLoc = new int[nqueues];
                for (int i = 0; i < nqueues; i++) {
                    prevLoc[i] = n[i] + 1;
                    curLoc[i] = x[i] + 1;
                }

                int prevState = RlTdAgent.getStateFromLoc(vSize, prevLoc);
                int curState = RlTdAgent.getStateFromLoc(vSize, curLoc);

                if (prevState >= 0 && prevState < v.length && curState >= 0 && curState < v.length) {
                    v[prevState] = (1 - lr) * v[prevState] + lr * (c - t * meanCostRate + v[curState]);

                    double v0 = v[0];
                    for (int i = 0; i < v.length; i++) {
                        v[i] -= v0;
                    }
                }

                t = 0.0;
                c = 0.0;
                n = x.clone();
            }
        }

        return v.clone();
    }

    public static class HashmapResult {
        public final Matrix X;
        public final Matrix Y;
        public HashmapResult(Matrix X, Matrix Y) {
            this.X = X;
            this.Y = Y;
        }
    }

    public HashmapResult solveByHashmap(RlEnvGeneral env, int numEpisodes) {
        reset(env);

        int nqueues = env.nqueues;

        HashMap<String, Double> pointValues = new HashMap<String, Double>();
        pointValues.put(intArrayToString(new int[nqueues]), 0.0);
        pointValues.put("external", 0.0);

        double t = 0.0;
        double c = 0.0;
        double bigT = 0.0;
        double bigC = 0.0;
        int[] x = new int[nqueues];
        int[] n = new int[nqueues];

        double eps = epsilon;

        int j = 0;
        while (j < numEpisodes) {
            if (j % 1000 == 0) {
                System.out.printf("running episode #%d.%n", j);
            }

            eps *= epsDecay;

            GeneralSampleEvent sampleEvent = env.sample();
            double dt = sampleEvent.dt;
            int depNode = sampleEvent.depNode;
            int arvNode = sampleEvent.arvNode;
            SampleResult sampleResult = sampleEvent.sampleResult;
            t += dt;

            int sumX = 0;
            for (int xi : x) sumX += xi;
            c += sumX * dt;

            if (contains(env.idxOfQueueInNodes, depNode)) {
                int depServer = indexOf(env.idxOfQueueInNodes, depNode);
                if (depServer >= 0) {
                    x[depServer] = Math.max(0, x[depServer] - 1);
                }
            }

            if (contains(env.idxOfActionNodes, depNode) && env.isInActionSpace(x)) {
                int[] actions = env.actionSpace.get(depNode);
                if (actions != null && actions.length != 0) {
                    double[] nextPointValues = new double[actions.length];
                    for (int actI = 0; actI < actions.length; actI++) {
                        int qIdx = indexOf(env.idxOfQueueInNodes, actions[actI]);
                        int[] tmpNextState = x.clone();
                        if (qIdx >= 0) {
                            tmpNextState[qIdx] = tmpNextState[qIdx] + 1;
                        }
                        String key = intArrayToString(tmpNextState);
                        nextPointValues[actI] = pointValues.containsKey(key)
                                ? pointValues.get(key) : pointValues.get("external");
                    }
                    double[] policy = RlTdAgent.createGreedyPolicy(nextPointValues, eps, actions.length);

                    double r = random.nextDouble();
                    double cumSum = 0.0;
                    int selectedIdx = actions.length - 1;
                    for (int a = 0; a < policy.length; a++) {
                        cumSum += policy[a];
                        if (r < cumSum) {
                            selectedIdx = a;
                            break;
                        }
                    }
                    arvNode = actions[selectedIdx];

                    Matrix eventMatrix = ((SampleResult) sampleResult).event;
                    if (eventMatrix != null && eventMatrix.getNumRows() > 1) {
                        for (int row = 0; row < eventMatrix.getNumRows(); row++) {
                            if (row == 1 || (row > 0 && eventMatrix.getNumCols() > 1)) {
                                eventMatrix.set(row, 1, (double) arvNode);
                                break;
                            }
                        }
                    }
                }
            }

            if (contains(env.idxOfQueueInNodes, arvNode)) {
                int arvServer = indexOf(env.idxOfQueueInNodes, arvNode);
                if (arvServer >= 0) {
                    x[arvServer] = x[arvServer] + 1;
                }
            }
            env.update(sampleResult);

            if (env.isInStateSpace(x)) {
                j++;
                bigT = env.gamma * bigT + t;
                bigC = env.gamma * bigC + c;
                double meanCostRate = (bigT > 0) ? bigC / bigT : 0.0;

                String nKey = intArrayToString(n);
                String xKey = intArrayToString(x);

                if (!pointValues.containsKey(nKey)) {
                    pointValues.put(nKey, pointValues.get("external"));
                }

                double curVal = pointValues.containsKey(xKey)
                        ? pointValues.get(xKey) : pointValues.get("external");

                pointValues.put(nKey, (1 - lr) * pointValues.get(nKey) + lr * (c - t * meanCostRate + curVal));

                boolean allZero = true;
                for (int ni : n) {
                    if (ni != 0) {
                        allZero = false;
                        break;
                    }
                }
                if (allZero) {
                    double subtractor = pointValues.get(nKey);
                    for (String key : new ArrayList<String>(pointValues.keySet())) {
                        pointValues.put(key, pointValues.get(key) - subtractor);
                    }
                }

                t = 0.0;
                c = 0.0;
                n = x.clone();
            }
        }

        pointValues.remove("external");

        int count = pointValues.size();
        Matrix resultX = new Matrix(count, 1 + nqueues);
        Matrix resultY = new Matrix(count, 1);
        int iterator = 0;
        for (Map.Entry<String, Double> entry : pointValues.entrySet()) {
            resultX.set(iterator, 0, 1.0);
            int[] stateVals = stringToIntArray(entry.getKey());
            for (int k = 0; k < stateVals.length; k++) {
                resultX.set(iterator, 1 + k, (double) stateVals[k]);
            }
            resultY.set(iterator, 0, entry.getValue());
            iterator++;
        }

        return new HashmapResult(resultX, resultY);
    }

    public static class ApproximationResult {
        public final Matrix X;
        public final Matrix Y;
        public final Matrix coefficients;
        public ApproximationResult(Matrix X, Matrix Y, Matrix coefficients) {
            this.X = X;
            this.Y = Y;
            this.coefficients = coefficients;
        }
    }

    public ApproximationResult solveByLinear(RlEnvGeneral env, int numEpisodes) {
        HashmapResult hashmapResult = solveByHashmap(env, numEpisodes);
        Matrix resultX = hashmapResult.X;
        Matrix resultY = hashmapResult.Y;

        Matrix xt = resultX.transpose();
        Matrix xtx = xt.mult(resultX);
        Matrix xty = xt.mult(resultY);

        Matrix coeff = Matrix.createLike(xty);
        Matrix.solve(xtx, xty, coeff);

        return new ApproximationResult(resultX, resultY, coeff);
    }

    public ApproximationResult solveByQuad(RlEnvGeneral env, int numEpisodes) {
        HashmapResult hashmapResult = solveByHashmap(env, numEpisodes);
        Matrix baseX = hashmapResult.X;
        Matrix resultY = hashmapResult.Y;

        int numRows = baseX.getNumRows();
        int baseCols = baseX.getNumCols();

        int quadCols = 0;
        for (int i = 1; i < baseCols; i++) {
            for (int j2 = i; j2 < baseCols; j2++) {
                quadCols++;
            }
        }

        int totalCols = baseCols + quadCols;
        Matrix augX = new Matrix(numRows, totalCols);

        for (int row = 0; row < numRows; row++) {
            for (int col = 0; col < baseCols; col++) {
                augX.set(row, col, baseX.get(row, col));
            }
        }

        int colIdx = baseCols;
        for (int i = 1; i < baseCols; i++) {
            for (int j2 = i; j2 < baseCols; j2++) {
                for (int row = 0; row < numRows; row++) {
                    augX.set(row, colIdx, baseX.get(row, i) * baseX.get(row, j2));
                }
                colIdx++;
            }
        }

        Matrix xt = augX.transpose();
        Matrix xtx = xt.mult(augX);
        Matrix xty = xt.mult(resultY);

        Matrix coeff = Matrix.createLike(xty);
        Matrix.solve(xtx, xty, coeff);

        return new ApproximationResult(augX, resultY, coeff);
    }

    private double[] genNextValues(RlEnvGeneral env, int[] curState, int[] actions) {
        double[] values = new double[actions.length];
        for (int actI = 0; actI < actions.length; actI++) {
            int qIdx = indexOf(env.idxOfQueueInNodes, actions[actI]);
            int[] tmpLoc = new int[curState.length];
            for (int i = 0; i < curState.length; i++) tmpLoc[i] = curState[i] + 1;
            if (qIdx >= 0) {
                tmpLoc[qIdx] = tmpLoc[qIdx] + 1;
            }
            int stateIdx = RlTdAgent.getStateFromLoc(vSize, tmpLoc);
            values[actI] = (stateIdx >= 0 && stateIdx < v.length) ? v[stateIdx] : 0.0;
        }
        return values;
    }

    private static boolean contains(int[] arr, int value) {
        for (int v : arr) {
            if (v == value) return true;
        }
        return false;
    }

    private static int indexOf(int[] arr, int value) {
        for (int i = 0; i < arr.length; i++) {
            if (arr[i] == value) return i;
        }
        return -1;
    }

    private static String intArrayToString(int[] arr) {
        StringBuilder sb = new StringBuilder();
        for (int i = 0; i < arr.length; i++) {
            if (i > 0) sb.append(' ');
            sb.append(arr[i]);
        }
        return sb.toString();
    }

    private static int[] stringToIntArray(String s) {
        String[] parts = s.trim().split("\\s+");
        int[] result = new int[parts.length];
        for (int i = 0; i < parts.length; i++) {
            result[i] = Integer.parseInt(parts[i]);
        }
        return result;
    }
}
