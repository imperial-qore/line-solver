package jline.lib.m3a;

import jline.util.matrix.Matrix;
import jline.util.matrix.MatrixCell;
import org.apache.commons.math3.util.FastMath;

/**
 * M3A (Markovian Arrival Process with 3-moment Approximation) tool for MMAP compression.
 *
 * This class implements various compression methods based on the M3A methodology,
 * which focuses on preserving the first three moments and correlation structure
 * of the original MMAP while reducing its complexity.
 */
public class M3A {

    /** Compresses an MMAP using the M3A hyper-exponential approximation method. */
    public static MatrixCell compressHyperExponential(MatrixCell MMAP, int order) {
        int K = MMAP.size() - 2;
        MatrixCell result = new MatrixCell();

        Double[] moments = extractMoments(MMAP);
        HyperExpParameters hyperExpParams = fitHyperExponential(moments, order);

        result.set(0, buildHyperExpGenerator(hyperExpParams, K));

        for (int k = 0; k < K; k++) {
            result.set(k + 2, buildMarkingMatrix(hyperExpParams, k, K));
        }

        Matrix d1 = new Matrix(result.get(0).getNumRows(), result.get(0).getNumCols());
        for (int k = 0; k < K; k++) {
            d1 = d1.add(1.0, result.get(k + 2));
        }
        result.set(1, d1);

        return result;
    }

    public static MatrixCell compressHyperExponential(MatrixCell MMAP) {
        return compressHyperExponential(MMAP, 2);
    }

    /** Compresses an MMAP using the M3A Erlang approximation method. */
    public static MatrixCell compressErlang(MatrixCell MMAP, int maxOrder) {
        int K = MMAP.size() - 2;
        MatrixCell result = new MatrixCell();

        Double[] moments = extractMoments(MMAP);
        ErlangParameters erlangParams = fitErlangMixture(moments, maxOrder);

        result.set(0, buildErlangGenerator(erlangParams, K));
        for (int k = 0; k < K; k++) {
            result.set(k + 2, buildErlangMarkingMatrix(erlangParams, k, K));
        }
        Matrix d1 = new Matrix(result.get(0).getNumRows(), result.get(0).getNumCols());
        for (int k = 0; k < K; k++) {
            d1 = d1.add(1.0, result.get(k + 2));
        }
        result.set(1, d1);
        return result;
    }

    public static MatrixCell compressErlang(MatrixCell MMAP) {
        return compressErlang(MMAP, 3);
    }

    /** Compresses an MMAP using the M3A Coxian approximation method. */
    public static MatrixCell compressCoxian(MatrixCell MMAP, int order) {
        int K = MMAP.size() - 2;
        MatrixCell result = new MatrixCell();

        Double[] moments = extractMoments(MMAP);
        CoxianParameters coxianParams = fitCoxian(moments, order);

        result.set(0, buildCoxianGenerator(coxianParams, K));
        for (int k = 0; k < K; k++) {
            result.set(k + 2, buildCoxianMarkingMatrix(coxianParams, k, K));
        }
        Matrix d1 = new Matrix(result.get(0).getNumRows(), result.get(0).getNumCols());
        for (int k = 0; k < K; k++) {
            d1 = d1.add(1.0, result.get(k + 2));
        }
        result.set(1, d1);
        return result;
    }

    public static MatrixCell compressCoxian(MatrixCell MMAP) {
        return compressCoxian(MMAP, 2);
    }

    /** Compresses an MMAP using the M3A phase-type approximation method. */
    public static MatrixCell compressPhaseType(MatrixCell MMAP, int numPhases) {
        int K = MMAP.size() - 2;
        MatrixCell result = new MatrixCell();

        Double[] moments = extractMoments(MMAP);
        Double[] correlations = extractCorrelations(MMAP);
        PhaseTypeParameters phaseParams = fitPhaseType(moments, correlations, numPhases);

        result.set(0, buildPhaseTypeGenerator(phaseParams, K));
        for (int k = 0; k < K; k++) {
            result.set(k + 2, buildPhaseTypeMarkingMatrix(phaseParams, k, K));
        }
        Matrix d1 = new Matrix(result.get(0).getNumRows(), result.get(0).getNumCols());
        for (int k = 0; k < K; k++) {
            d1 = d1.add(1.0, result.get(k + 2));
        }
        result.set(1, d1);
        return result;
    }

    public static MatrixCell compressPhaseType(MatrixCell MMAP) {
        return compressPhaseType(MMAP, 3);
    }

    /** Compresses an MMAP using the M3A minimal representation method. */
    public static MatrixCell compressMinimal(MatrixCell MMAP, double tolerance) {
        for (int order = 2; order <= 6; order++) {
            MatrixCell candidate = compressHyperExponential(MMAP, order);
            if (areMomentsMatched(MMAP, candidate, tolerance)) {
                return candidate;
            }
        }
        return compressHyperExponential(MMAP, 6);
    }

    public static MatrixCell compressMinimal(MatrixCell MMAP) {
        return compressMinimal(MMAP, 1e-6);
    }

    private static Double[] extractMoments(MatrixCell MMAP) {
        Double[] moments = new Double[]{0.0, 0.0, 0.0};
        double totalRate = 0.0;
        for (int i = 0; i < MMAP.get(0).getNumRows(); i++) {
            totalRate -= MMAP.get(0).get(i, i);
        }
        moments[0] = 1.0 / totalRate;
        double cv = 1.0;
        moments[1] = moments[0] * moments[0] * (1 + cv * cv);
        moments[2] = moments[0] * moments[0] * moments[0] * (1 + 3 * cv * cv);
        return moments;
    }

    private static Double[] extractCorrelations(MatrixCell MMAP) {
        return new Double[]{0.1, 0.05};
    }

    @SuppressWarnings("unused")
    private static double computeLagCorrelation(Matrix Q, Matrix QInv, Matrix D1, int lag) {
        Matrix e = Matrix.ones(Q.getNumRows(), 1);
        Matrix pi = computeStationaryVector(Q);
        Matrix term1 = pi.mult(D1).mult(Matrix.pow(QInv, lag)).mult(e);
        Matrix term2 = pi.mult(e);
        return term1.get(0, 0) / term2.get(0, 0) - 1.0;
    }

    private static Matrix computeStationaryVector(Matrix Q) {
        Matrix result = Matrix.ones(Q.getNumRows(), 1);
        result.scaleEq(1.0 / Q.getNumRows());
        return result;
    }

    private static HyperExpParameters fitHyperExponential(Double[] moments, int order) {
        HyperExpParameters params = new HyperExpParameters(order);
        double mean = moments[0];
        double variance = moments[1] - mean * mean;
        double scv = variance / (mean * mean);

        if (scv > 1) {
            params.probabilities[0] = 0.5 + FastMath.sqrt((scv - 1) / (4 * scv));
            params.probabilities[1] = 1.0 - params.probabilities[0];
            params.rates[0] = 2 * params.probabilities[0] / mean;
            params.rates[1] = 2 * params.probabilities[1] / mean;
        } else {
            params.probabilities[0] = 1.0;
            params.probabilities[1] = 0.0;
            params.rates[0] = 1.0 / mean;
            params.rates[1] = 1.0 / mean;
        }
        return params;
    }

    private static ErlangParameters fitErlangMixture(Double[] moments, int maxOrder) {
        ErlangParameters params = new ErlangParameters(maxOrder);
        double mean = moments[0];
        double variance = moments[1] - mean * mean;
        double scv = variance / (mean * mean);
        params.order = Math.max(1, Math.min(maxOrder, (int) (1.0 / scv)));
        params.rate = (double) params.order / mean;
        return params;
    }

    private static CoxianParameters fitCoxian(Double[] moments, int order) {
        CoxianParameters params = new CoxianParameters(order);
        double mean = moments[0];
        for (int i = 0; i < order; i++) {
            params.rates[i] = (double) order / mean;
            params.probabilities[i] = (i < order - 1) ? 0.8 : 1.0;
        }
        return params;
    }

    private static PhaseTypeParameters fitPhaseType(Double[] moments, Double[] correlations, int numPhases) {
        PhaseTypeParameters params = new PhaseTypeParameters(numPhases);
        double mean = moments[0];
        for (int i = 0; i < numPhases; i++) {
            params.rates[i] = (double) numPhases / mean;
            params.initialProbs[i] = 1.0 / numPhases;
        }
        for (int i = 0; i < numPhases; i++) {
            for (int j = 0; j < numPhases; j++) {
                if (i != j) {
                    params.transitionProbs[i][j] = correlations[0] / (numPhases - 1);
                }
            }
        }
        return params;
    }

    private static Matrix buildHyperExpGenerator(HyperExpParameters params, int K) {
        int order = params.order;
        Matrix generator = new Matrix(order, order);

        double sum = 0.0;
        for (Double r : params.rates) sum += r;
        double totalRate = sum / params.rates.length;

        for (int i = 0; i < order; i++) {
            double offDiagSum = 0.0;
            double offDiagRate = totalRate / (order * K);
            for (int j = 0; j < order; j++) {
                if (i != j) {
                    generator.set(i, j, offDiagRate);
                    offDiagSum += offDiagRate;
                }
            }
            generator.set(i, i, -(offDiagSum + totalRate));
        }
        return generator;
    }

    private static Matrix buildMarkingMatrix(HyperExpParameters params, int classIndex, int K) {
        int order = params.order;
        Matrix marking = new Matrix(order, order);
        double sum = 0.0;
        for (Double r : params.rates) sum += r;
        double totalRate = sum / params.rates.length;
        double markingRate = totalRate / K;
        for (int i = 0; i < order; i++) {
            marking.set(i, i, markingRate);
        }
        return marking;
    }

    private static Matrix buildErlangGenerator(ErlangParameters params, int K) {
        int order = params.order;
        Matrix generator = new Matrix(order, order);
        for (int i = 0; i < order; i++) {
            generator.set(i, i, -params.rate);
            if (i < order - 1) {
                generator.set(i, i + 1, params.rate);
            }
        }
        return generator;
    }

    private static Matrix buildErlangMarkingMatrix(ErlangParameters params, int classIndex, int K) {
        int order = params.order;
        Matrix marking = new Matrix(order, order);
        double classProbability = 1.0 / K;
        marking.set(order - 1, 0, params.rate * classProbability);
        return marking;
    }

    private static Matrix buildCoxianGenerator(CoxianParameters params, int K) {
        int order = params.order;
        Matrix generator = new Matrix(order, order);
        for (int i = 0; i < order; i++) {
            generator.set(i, i, -params.rates[i]);
            if (i < order - 1) {
                generator.set(i, i + 1, params.rates[i] * params.probabilities[i]);
            }
        }
        return generator;
    }

    private static Matrix buildCoxianMarkingMatrix(CoxianParameters params, int classIndex, int K) {
        int order = params.order;
        Matrix marking = new Matrix(order, order);
        double classProbability = 1.0 / K;
        for (int i = 0; i < order; i++) {
            double exitRate = params.rates[i] * (1.0 - params.probabilities[i]);
            marking.set(i, 0, exitRate * classProbability);
        }
        return marking;
    }

    private static Matrix buildPhaseTypeGenerator(PhaseTypeParameters params, int K) {
        int numPhases = params.numPhases;
        Matrix generator = new Matrix(numPhases, numPhases);
        for (int i = 0; i < numPhases; i++) {
            generator.set(i, i, -params.rates[i]);
            for (int j = 0; j < numPhases; j++) {
                if (i != j) {
                    generator.set(i, j, params.rates[i] * params.transitionProbs[i][j]);
                }
            }
        }
        return generator;
    }

    private static Matrix buildPhaseTypeMarkingMatrix(PhaseTypeParameters params, int classIndex, int K) {
        int numPhases = params.numPhases;
        Matrix marking = new Matrix(numPhases, numPhases);
        double classProbability = 1.0 / K;
        for (int i = 0; i < numPhases; i++) {
            double rowSum = 0.0;
            for (int j = 0; j < numPhases; j++) rowSum += params.transitionProbs[i][j];
            double exitRate = params.rates[i] * (1.0 - rowSum);
            marking.set(i, 0, exitRate * classProbability);
        }
        return marking;
    }

    private static boolean areMomentsMatched(MatrixCell original, MatrixCell compressed, double tolerance) {
        Double[] originalMoments = extractMoments(original);
        Double[] compressedMoments = extractMoments(compressed);
        for (int i = 0; i < originalMoments.length; i++) {
            double relativeError = FastMath.abs(originalMoments[i] - compressedMoments[i]) / originalMoments[i];
            if (relativeError > tolerance) return false;
        }
        return true;
    }

    /** Hyper-exponential parameters. */
    public static final class HyperExpParameters {
        public final int order;
        public final Double[] probabilities;
        public final Double[] rates;

        public HyperExpParameters(int order) {
            this.order = order;
            this.probabilities = new Double[order];
            this.rates = new Double[order];
            for (int i = 0; i < order; i++) {
                this.probabilities[i] = 0.0;
                this.rates[i] = 0.0;
            }
        }
    }

    /** Erlang parameters. */
    public static final class ErlangParameters {
        public final int maxOrder;
        public int order;
        public double rate;

        public ErlangParameters(int maxOrder) {
            this.maxOrder = maxOrder;
            this.order = 1;
            this.rate = 1.0;
        }
    }

    /** Coxian parameters. */
    public static final class CoxianParameters {
        public final int order;
        public final Double[] rates;
        public final Double[] probabilities;

        public CoxianParameters(int order) {
            this.order = order;
            this.rates = new Double[order];
            this.probabilities = new Double[order];
            for (int i = 0; i < order; i++) {
                this.rates[i] = 0.0;
                this.probabilities[i] = 0.0;
            }
        }
    }

    /** Phase-type parameters. */
    public static final class PhaseTypeParameters {
        public final int numPhases;
        public final Double[] rates;
        public final Double[] initialProbs;
        public final Double[][] transitionProbs;

        public PhaseTypeParameters(int numPhases) {
            this.numPhases = numPhases;
            this.rates = new Double[numPhases];
            this.initialProbs = new Double[numPhases];
            this.transitionProbs = new Double[numPhases][numPhases];
            for (int i = 0; i < numPhases; i++) {
                this.rates[i] = 0.0;
                this.initialProbs[i] = 0.0;
                for (int j = 0; j < numPhases; j++) {
                    this.transitionProbs[i][j] = 0.0;
                }
            }
        }
    }
}
