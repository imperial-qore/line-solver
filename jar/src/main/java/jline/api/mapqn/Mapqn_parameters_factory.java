package jline.api.mapqn;

import java.util.Map;

import jline.GlobalConstants;
import jline.lang.JobClass;
import jline.lang.NetworkStruct;
import jline.lang.nodes.Station;
import jline.util.matrix.Matrix;

/**
 * Factory class for creating Mapqn_parameters from NetworkStruct.
 */
public final class Mapqn_parameters_factory {
    private Mapqn_parameters_factory() {}

    public static Mapqn_parameters createFromNetworkStruct(NetworkStruct networkStruct) {
        networkStruct.validateStructuralConsistency();
        int M = networkStruct.nstations;
        int N = networkStruct.nclosedjobs;
        if (M <= 0) throw new IllegalArgumentException("Network must have at least one station");
        if (N <= 0) throw new IllegalArgumentException("Network must have closed jobs for MAPQN analysis");

        if (hasFiniteCapacity(networkStruct)) return createFiniteCapacityParameters(networkStruct);
        if (hasMultiplePhases(networkStruct)) return createLinearReductionParameters(networkStruct);
        return createBasicParameters(networkStruct);
    }

    public static LinearReductionParameters createLinearReductionParameters(NetworkStruct networkStruct) {
        int M = networkStruct.nstations;
        int N = networkStruct.nclosedjobs;
        int[] K = new int[M];
        for (int i = 0; i < M; i++) {
            K[i] = (networkStruct.phases != null) ? (int) networkStruct.phases.get(i, 0) : 1;
            if (K[i] < 1) K[i] = 1;
        }
        Matrix[] mu = new Matrix[M];
        for (int i = 0; i < M; i++) {
            Matrix muMat = null;
            if (networkStruct.mu != null && networkStruct.stations != null) {
                Station station = networkStruct.stations.get(i);
                Map<JobClass, Matrix> muMap = networkStruct.mu.get(station);
                if (muMap != null && !muMap.isEmpty()) {
                    JobClass firstClass = (networkStruct.jobclasses != null && !networkStruct.jobclasses.isEmpty())
                            ? networkStruct.jobclasses.get(0) : null;
                    Matrix muMatrix = muMap.get(firstClass);
                    if (muMatrix != null) muMat = convertToPhaseTransitionMatrix(muMatrix, K[i]);
                }
            }
            mu[i] = (muMat != null) ? muMat : Matrix.eye(K[i]);
        }
        Matrix r = extractRoutingMatrix(networkStruct);
        Matrix[] v = new Matrix[M];
        for (int i = 0; i < M; i++) v[i] = Matrix.zeros(K[i], K[i]);
        return new LinearReductionParameters(M, N, K, mu, r, v);
    }

    private static Mapqn_parameters createFiniteCapacityParameters(NetworkStruct networkStruct) {
        int M = networkStruct.nstations;
        int[] F = new int[M];
        for (int i = 0; i < M; i++) {
            F[i] = (networkStruct.cap != null) ? (int) networkStruct.cap.get(i, 0) : Integer.MAX_VALUE;
        }
        return createMapqn_qr_bounds_rsrd_parameters(networkStruct, F);
    }

    private static Mapqn_qr_bounds_rsrd_parameters createMapqn_qr_bounds_rsrd_parameters(NetworkStruct networkStruct, int[] F) {
        int M = networkStruct.nstations;
        int N = networkStruct.nclosedjobs;
        int[] K = new int[M];
        for (int i = 0; i < M; i++) {
            K[i] = (networkStruct.phases != null) ? (int) networkStruct.phases.get(i, 0) : 1;
            if (K[i] < 1) K[i] = 1;
        }
        Matrix[] mu = new Matrix[M];
        for (int i = 0; i < M; i++) mu[i] = extractServiceRateMatrix(networkStruct, i, K[i]);
        Matrix[] v = new Matrix[M];
        for (int i = 0; i < M; i++) v[i] = Matrix.zeros(K[i], K[i]);
        double[][] alpha = new double[M][N];
        for (int i = 0; i < M; i++) {
            for (int j = 0; j < N; j++) alpha[i][j] = 1.0;
        }
        if (networkStruct.cdscaling != null && networkStruct.stations != null) {
            for (int i = 0; i < M; i++) {
                Station station = networkStruct.stations.get(i);
                Object scalingFunc = networkStruct.cdscaling.get(station);
                if (scalingFunc != null) {
                    for (int n = 1; n <= N; n++) {
                        Matrix stateMatrix = new Matrix(1, 1);
                        stateMatrix.set(0, 0, (double) n);
                        try {
                            java.lang.reflect.Method m = scalingFunc.getClass().getMethod("apply", Matrix.class);
                            Object result = m.invoke(scalingFunc, stateMatrix);
                            if (result instanceof Number) alpha[i][n - 1] = ((Number) result).doubleValue();
                        } catch (Exception e) {
                            // ignore
                        }
                    }
                }
            }
        }
        Matrix r = extractRoutingMatrix(networkStruct);
        return new Mapqn_qr_bounds_rsrd_parameters(M, N, F, K, mu, v, alpha, r);
    }

    private static LinearReductionParameters createBasicParameters(NetworkStruct networkStruct) {
        int M = networkStruct.nstations;
        int N = networkStruct.nclosedjobs;
        int[] K = new int[M];
        for (int i = 0; i < M; i++) K[i] = 1;
        Matrix[] mu = new Matrix[M];
        for (int i = 0; i < M; i++) {
            double rate = (networkStruct.rates != null) ? networkStruct.rates.get(i, 0) : 1.0;
            mu[i] = new Matrix(new double[][]{{rate}});
        }
        Matrix r = extractRoutingMatrix(networkStruct);
        Matrix[] v = new Matrix[M];
        for (int i = 0; i < M; i++) v[i] = new Matrix(new double[][]{{0.0}});
        return new LinearReductionParameters(M, N, K, mu, r, v);
    }

    private static boolean hasFiniteCapacity(NetworkStruct networkStruct) {
        if (networkStruct.cap == null) return false;
        for (int i = 0; i < networkStruct.nstations; i++) {
            double capacity = networkStruct.cap.get(i, 0);
            if (capacity > 0 && capacity < GlobalConstants.Inf) return true;
        }
        return false;
    }

    private static boolean hasMultiplePhases(NetworkStruct networkStruct) {
        if (networkStruct.phases == null) return false;
        for (int i = 0; i < networkStruct.nstations; i++) {
            if (networkStruct.phases.get(i, 0) > 1) return true;
        }
        return false;
    }

    private static Matrix extractRoutingMatrix(NetworkStruct networkStruct) {
        int M = networkStruct.nstations;
        if (networkStruct.rt != null) {
            int nclasses = networkStruct.nclasses;
            Matrix r = new Matrix(M, M);
            for (int i = 0; i < M; i++) {
                for (int j = 0; j < M; j++) {
                    double sum = 0.0;
                    for (int c = 0; c < nclasses; c++) {
                        int fromIdx = i * nclasses + c;
                        int toIdx = j * nclasses + c;
                        sum += networkStruct.rt.get(fromIdx, toIdx);
                    }
                    r.set(i, j, sum / nclasses);
                }
            }
            for (int i = 0; i < M; i++) {
                double rowSum = 0.0;
                for (int j = 0; j < M; j++) rowSum += r.get(i, j);
                if (rowSum > 0) {
                    for (int j = 0; j < M; j++) r.set(i, j, r.get(i, j) / rowSum);
                }
            }
            return r;
        }
        return createDefaultRoutingMatrix(networkStruct);
    }

    private static Matrix createDefaultRoutingMatrix(NetworkStruct networkStruct) {
        int M = networkStruct.nstations;
        Matrix r = new Matrix(M, M);
        for (int i = 0; i < M; i++) {
            for (int j = 0; j < M; j++) r.set(i, j, j == (i + 1) % M ? 1.0 : 0.0);
        }
        return r;
    }

    private static Matrix extractServiceRateMatrix(NetworkStruct networkStruct, int stationIdx, int numPhases) {
        if (networkStruct.mu != null && networkStruct.stations != null) {
            Station station = networkStruct.stations.get(stationIdx);
            Map<JobClass, Matrix> muMap = networkStruct.mu.get(station);
            if (muMap != null && !muMap.isEmpty()) {
                JobClass firstClass = (networkStruct.jobclasses != null && !networkStruct.jobclasses.isEmpty())
                        ? networkStruct.jobclasses.get(0) : null;
                Matrix muMatrix = muMap.get(firstClass);
                if (muMatrix != null) return convertToPhaseTransitionMatrix(muMatrix, numPhases);
            }
        }
        if (networkStruct.rates != null) {
            double rate = networkStruct.rates.get(stationIdx, 0);
            if (rate > 0) {
                Matrix mu = Matrix.zeros(numPhases, numPhases);
                for (int k = 0; k < numPhases; k++) mu.set(k, k, rate);
                return mu;
            }
        }
        return Matrix.eye(numPhases);
    }

    private static Matrix convertToPhaseTransitionMatrix(Matrix serviceRates, int numPhases) {
        if (serviceRates.getNumRows() == numPhases && serviceRates.getNumCols() == numPhases) return serviceRates.copy();
        if (serviceRates.getNumRows() == 1 && serviceRates.getNumCols() == 1) {
            double rate = serviceRates.get(0, 0);
            Matrix mu = Matrix.zeros(numPhases, numPhases);
            for (int k = 0; k < numPhases; k++) mu.set(k, k, rate);
            return mu;
        }
        if ((serviceRates.getNumRows() == numPhases && serviceRates.getNumCols() == 1)
                || (serviceRates.getNumRows() == 1 && serviceRates.getNumCols() == numPhases)) {
            Matrix mu = Matrix.zeros(numPhases, numPhases);
            for (int k = 0; k < numPhases; k++) {
                double rate = (serviceRates.getNumRows() == 1) ? serviceRates.get(0, k) : serviceRates.get(k, 0);
                mu.set(k, k, rate);
            }
            return mu;
        }
        return Matrix.eye(numPhases);
    }
}
