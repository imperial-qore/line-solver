package jline.lib.perm;

import java.util.ArrayList;
import java.util.HashMap;
import java.util.List;
import java.util.Map;

import jline.util.Pair;
import jline.util.matrix.Matrix;

/**
 * Queueing network model with think time.
 */
public class NetworkThink {

    private final int numberOfQueues;
    private final int numberOfClasses;
    private final int[] numberPerClass;
    private final double[] thinkTime;
    private final double[][] meanServiceDemand;
    private final boolean progress;

    public NetworkThink(int numberOfQueues, int numberOfClasses, int[] numberPerClass,
                        double[] thinkTime, double[][] meanServiceDemand, boolean progress) {
        this.numberOfQueues = numberOfQueues;
        this.numberOfClasses = numberOfClasses;
        this.numberPerClass = numberPerClass;
        this.thinkTime = thinkTime;
        this.meanServiceDemand = meanServiceDemand;
        this.progress = progress;
    }

    public NetworkThink(int numberOfQueues, int numberOfClasses, int[] numberPerClass,
                        double[] thinkTime, double[][] meanServiceDemand) {
        this(numberOfQueues, numberOfClasses, numberPerClass, thinkTime, meanServiceDemand, true);
    }

    public double joint(int[][] state) {
        double termMSD = 1.0;
        double invFactorial = 1.0;
        double factorialTerm = 1.0;
        double thinkTimeTerm = 1.0;

        for (int i = 0; i < numberOfQueues; i++) {
            for (int j = 0; j < numberOfClasses; j++) {
                termMSD *= Math.pow(meanServiceDemand[i][j], state[i][j]);
            }
        }

        for (int i = 0; i < numberOfQueues; i++) {
            for (int j = 0; j < numberOfClasses; j++) {
                invFactorial *= factorial(state[i][j]);
            }
            int queueTotal = 0;
            for (int v : state[i]) queueTotal += v;
            factorialTerm *= factorial(queueTotal);
        }

        for (int j = 0; j < numberOfClasses; j++) {
            int totalInClass = 0;
            for (int[] row : state) totalInClass += row[j];
            int thinkingJobs = numberPerClass[j] - totalInClass;
            thinkTimeTerm *= Math.pow(thinkTime[j], thinkingJobs);
        }

        return termMSD * factorialTerm * thinkTimeTerm / invFactorial;
    }

    public NetworkNoThink.MarginalResult marginal(PermSolver solver, int[] state, boolean preprocessing) {
        int totalJobs = 0;
        for (int v : state) totalJobs += v;

        double[][] matrixData = new double[totalJobs][totalJobs];

        int rowIndex = 0;
        for (int i = 0; i < numberOfQueues; i++) {
            for (int rep = 0; rep < state[i]; rep++) {
                int colIndex = 0;
                for (int j = 0; j < numberOfQueues; j++) {
                    for (int rep2 = 0; rep2 < state[j]; rep2++) {
                        matrixData[rowIndex][colIndex] = meanServiceDemand[i][j];
                        colIndex++;
                    }
                }
                int totalState = 0;
                for (int v : state) totalState += v;
                for (int j = 0; j < numberOfClasses; j++) {
                    int thinkingJobs = numberPerClass[j] - totalState;
                    for (int rep2 = 0; rep2 < thinkingJobs; rep2++) {
                        matrixData[rowIndex][colIndex] = thinkTime[j];
                        colIndex++;
                    }
                }
                rowIndex++;
            }
        }

        Matrix matrix = new Matrix(matrixData);

        Matrix processedMatrix;
        double rescalingFactor;
        if (preprocessing) {
            Pair<Matrix, Double> p = QueueingNetwork.preprocessingDS(matrix);
            processedMatrix = p.getLeft();
            rescalingFactor = p.getRight();
        } else {
            processedMatrix = matrix;
            rescalingFactor = 1.0;
        }

        PermSolver solverInstance;
        if (solver instanceof BethePermanent) {
            solverInstance = new BethePermanent(processedMatrix, 1e-3, 1000, true);
        } else if (solver instanceof NaivePermanent) {
            solverInstance = new NaivePermanent(processedMatrix, true);
        } else if (solver instanceof RyzerPermanent) {
            solverInstance = new RyzerPermanent(processedMatrix, "default", true);
        } else if (solver instanceof AdaPartSampler) {
            solverInstance = new AdaPartSampler(processedMatrix);
            solverInstance.solve();
        } else if (solver instanceof HuberLawSampler) {
            solverInstance = new HuberLawSampler(processedMatrix);
            solverInstance.solve();
        } else {
            throw new IllegalArgumentException("Unsupported solver type");
        }

        double factorialNormalization = 1.0;
        for (int jobs : state) {
            factorialNormalization *= factorial(jobs);
        }

        double probability = solverInstance.value * rescalingFactor / factorialNormalization;
        return new NetworkNoThink.MarginalResult(probability, solverInstance.time, solverInstance.memory);
    }

    public NetworkNoThink.MarginalResult marginal(PermSolver solver, int[] state) {
        return marginal(solver, state, false);
    }

    public Map<int[], NetworkNoThink.MarginalResult> generateMarginal(PermSolver solver, boolean preprocessing) {
        Map<int[], NetworkNoThink.MarginalResult> results = new HashMap<int[], NetworkNoThink.MarginalResult>();
        List<int[]> states = generateAllMarginalStates();
        for (int[] state : states) {
            NetworkNoThink.MarginalResult result = marginal(solver, state, preprocessing);
            results.put(state, result);
            if (progress) {
                System.out.println("Processed state: " + java.util.Arrays.toString(state));
            }
        }
        return results;
    }

    public Map<int[], NetworkNoThink.MarginalResult> generateMarginal(PermSolver solver) {
        return generateMarginal(solver, false);
    }

    private List<int[]> generateAllMarginalStates() {
        int totalJobs = 0;
        for (int v : numberPerClass) totalJobs += v;
        List<int[]> states = new ArrayList<int[]>();
        generateStatesRecursive(new int[numberOfQueues], 0, totalJobs, states);
        return states;
    }

    private void generateStatesRecursive(int[] currentState, int queueIndex, int remainingJobs, List<int[]> results) {
        if (queueIndex == numberOfQueues - 1) {
            currentState[queueIndex] = remainingJobs;
            results.add(currentState.clone());
            return;
        }
        for (int jobs = 0; jobs <= remainingJobs; jobs++) {
            currentState[queueIndex] = jobs;
            generateStatesRecursive(currentState, queueIndex + 1, remainingJobs - jobs, results);
        }
    }

    private static double factorial(int n) {
        if (n <= 1) return 1.0;
        double result = 1.0;
        for (int i = 2; i <= n; i++) result *= i;
        return result;
    }
}
