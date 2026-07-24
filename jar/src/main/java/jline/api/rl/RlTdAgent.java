/**
 * @file Temporal Difference Learning Agent for Queueing Routing.
 *
 * @since LINE 3.0
 */
package jline.api.rl;

import jline.GlobalConstants;

import java.util.ArrayList;
import java.util.List;
import java.util.Random;

/**
 * TD learning agent for queueing network routing decisions.
 */
public class RlTdAgent {

    public final double lr;
    public double epsilon;
    public final double epsDecay;

    public double[] v = new double[]{0.0};
    public double[] q = new double[]{0.0};
    public int[] vSize = new int[]{0};
    public int[] qSize = new int[]{0};

    private final Random random = new Random();

    public RlTdAgent() {
        this(0.05, 1.0, 0.99);
    }

    public RlTdAgent(double lr) {
        this(lr, 1.0, 0.99);
    }

    public RlTdAgent(double lr, double epsilon) {
        this(lr, epsilon, 0.99);
    }

    public RlTdAgent(double lr, double epsilon, double epsDecay) {
        this.lr = lr;
        this.epsilon = epsilon;
        this.epsDecay = epsDecay;
    }

    public void reset(RlEnv env) {
        v = new double[]{0.0};
        q = new double[]{0.0};
        vSize = new int[]{0};
        qSize = new int[]{0};
        env.reset();
    }

    public double[] getValueFunction() {
        return v.clone();
    }

    public double[] getQFunction() {
        return q.clone();
    }

    public void solve(RlEnv env) {
        reset(env);

        int actionSize = env.actionSize;
        int stateSizePlusPad = env.stateSize + 5;

        vSize = new int[actionSize];
        for (int i = 0; i < actionSize; i++) vSize[i] = stateSizePlusPad;
        int totalVSize = 1;
        for (int dim : vSize) totalVSize *= dim;
        v = new double[totalVSize];

        qSize = new int[actionSize + 1];
        for (int i = 0; i < actionSize; i++) qSize[i] = stateSizePlusPad;
        qSize[actionSize] = actionSize;
        int totalQSize = 1;
        for (int dim : qSize) totalQSize *= dim;
        q = new double[totalQSize];
        for (int i = 0; i < totalQSize; i++) q[i] = random.nextDouble();

        int[] x = new int[actionSize];
        int[] n = new int[actionSize];
        double t = 0.0;
        double c = 0.0;
        double bigT = 0.0;
        double bigC = 0.0;

        int numEpisodes = 10000;
        double eps = epsilon;
        int j = 0;

        while (j < numEpisodes) {
            if (j % 1000 == 0) {
                System.out.printf("[rl_td_agent] running episode #%d.%n", j);
            }
            eps *= epsDecay;
            SampleEvent sampleEvent = env.sample();
            double dt = sampleEvent.t;
            int depNode = sampleEvent.depNode;
            t += dt;

            int sumX = 0;
            for (int xi : x) sumX += xi;
            c += sumX * dt;

            if (contains(env.idxOfSourceInNodes, depNode)) {
                if (env.isInActionSpace(env.model.getNodes())) {
                    double[] nextValues = new double[actionSize];
                    for (int a = 0; a < actionSize; a++) {
                        int[] nextLoc = new int[actionSize];
                        for (int kk = 0; kk < actionSize; kk++) nextLoc[kk] = x[kk] + 1;
                        nextLoc[a] = nextLoc[a] + 1;
                        int idx = getStateFromLoc(vSize, nextLoc);
                        nextValues[a] = (idx >= 0 && idx < v.length) ? v[idx] : 0.0;
                    }
                    double[] policy = createGreedyPolicy(nextValues, eps, actionSize);
                    double r = random.nextDouble();
                    double cumSum = 0.0;
                    int action = actionSize - 1;
                    for (int a = 0; a < actionSize; a++) {
                        cumSum += policy[a];
                        if (r < cumSum) { action = a; break; }
                    }
                    x[action] = x[action] + 1;
                    env.update(x);
                } else {
                    int minVal = x[0];
                    List<Integer> minIndices = new ArrayList<Integer>();
                    minIndices.add(0);
                    for (int kk = 1; kk < actionSize; kk++) {
                        if (x[kk] < minVal) {
                            minVal = x[kk];
                            minIndices.clear();
                            minIndices.add(kk);
                        } else if (x[kk] == minVal) {
                            minIndices.add(kk);
                        }
                    }
                    int action = (minIndices.size() > 1)
                            ? minIndices.get(random.nextInt(minIndices.size()))
                            : minIndices.get(0);
                    x[action] = x[action] + 1;
                    env.update(x);
                }
            } else if (contains(env.idxOfQueueInNodes, depNode)) {
                int queueIdx = indexOf(env.idxOfQueueInNodes, depNode);
                if (queueIdx >= 0) x[queueIdx] = Math.max(0, x[queueIdx] - 1);
                env.update(x);
            }

            if (env.isInStateSpace(env.model.getNodes())) {
                j++;
                bigT = env.gamma * bigT + t;
                bigC = env.gamma * bigC + c;
                double meanCostRate = (bigT > 0) ? bigC / bigT : 0.0;
                int[] prevLoc = new int[actionSize];
                int[] curLoc = new int[actionSize];
                for (int i = 0; i < actionSize; i++) {
                    prevLoc[i] = n[i] + 1;
                    curLoc[i] = x[i] + 1;
                }
                int prevState = getStateFromLoc(vSize, prevLoc);
                int curState = getStateFromLoc(vSize, curLoc);
                if (prevState >= 0 && prevState < v.length && curState >= 0 && curState < v.length) {
                    v[prevState] = (1 - lr) * v[prevState] + lr * (c - t * meanCostRate + v[curState]);
                    double v0 = v[0];
                    for (int i = 0; i < v.length; i++) v[i] -= v0;
                }
                t = 0.0;
                c = 0.0;
                n = x.clone();
            }
        }
    }

    public static double[] createGreedyPolicy(double[] stateQ, double epsilon, int nA) {
        double[] policy = new double[nA];
        for (int i = 0; i < nA; i++) policy[i] = epsilon / nA;
        double minVal = Double.MAX_VALUE;
        for (double qVal : stateQ) {
            if (qVal < minVal) minVal = qVal;
        }
        List<Integer> argmin = new ArrayList<Integer>();
        for (int i = 0; i < stateQ.length; i++) {
            if (stateQ[i] - minVal < GlobalConstants.FineTol) argmin.add(i);
        }
        double exploitProb = (1 - epsilon) / argmin.size();
        for (Integer idx : argmin) policy[idx] += exploitProb;
        return policy;
    }

    public static int getStateFromLoc(int[] objSize, int[] loc) {
        if (objSize.length != loc.length) return 0;
        int s = 0;
        int stride = 1;
        for (int i = 0; i < objSize.length; i++) {
            int locVal = loc[i] - 1;
            if (locVal < 0 || locVal >= objSize[i]) {
                int clamped = Math.max(0, Math.min(locVal, objSize[i] - 1));
                s += clamped * stride;
            } else {
                s += locVal * stride;
            }
            stride *= objSize[i];
        }
        return s;
    }

    public static int[] getStateFromLocs(int[] objSize, int[][] locs) {
        int[] result = new int[locs.length];
        for (int i = 0; i < locs.length; i++) {
            result[i] = getStateFromLoc(objSize, locs[i]);
        }
        return result;
    }

    private static boolean contains(int[] arr, int value) {
        for (int v : arr) if (v == value) return true;
        return false;
    }

    private static int indexOf(int[] arr, int value) {
        for (int i = 0; i < arr.length; i++) if (arr[i] == value) return i;
        return -1;
    }
}
