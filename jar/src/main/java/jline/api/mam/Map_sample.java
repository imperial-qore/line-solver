/**
 * @file Markovian Arrival Process sample generation
 *
 * Generates random samples from MAP distributions for simulation and empirical analysis.
 * Essential for stochastic simulation, Monte Carlo methods, and model validation studies.
 *
 * @since LINE 3.0
 */
package jline.api.mam;

import jline.util.matrix.Matrix;
import jline.util.matrix.MatrixCell;

import java.util.Random;

public final class Map_sample {
    private Map_sample() {}

    /**
     * Generates samples of inter-arrival times from a MAP using a specified number of samples and a random generator.
     *
     * @param MAP    the MatrixCell representing the MAP, containing the D0 and D1 matrices
     * @param n      the number of samples to generate
     * @param random the random number generator to use
     * @return an array of doubles containing the generated samples
     */
    public static double[] map_sample(MatrixCell MAP, long n, Random random) {
        return map_sample(MAP.get(0), MAP.get(1), n, random);
    }

    /**
     * Generates samples of inter-arrival times from a MAP using a specified number of samples and a random generator.
     *
     * @param D0     the hidden transition matrix of the MAP
     * @param D1     the visible transition matrix of the MAP
     * @param n      the number of samples to generate
     * @param random the random number generator to use
     * @return an array of doubles containing the generated samples
     */
    public static double[] map_sample(Matrix D0, Matrix D1, long n, Random random) {
        double[] samples = new double[(int) n];
        if (D0.getNumElements() == 1) { // exponential distribution
            double lambda = D1.value();
            for (long i = 0; i < n; i++) {
                samples[(int) i] = -Math.log(random.nextDouble()) / lambda;
            }
        } else {
            int nphases = D0.getNumCols();
            Matrix pie = Map_pie.map_pie(D0, D1);
            int currentState = pie.getNumRows();
            double sum = 0.0;
            double r = random.nextDouble();
            for (int i = 0; i < nphases; i++) {
                sum += pie.get(0, i);
                if (r < sum) {
                    currentState = i;
                    break;
                }
            }

            double[] row = new double[2 * nphases];
            for (long i = 0; i < n; i++) {
                boolean continueSample = true;
                samples[(int) i] = 0.0;
                while (continueSample) {
                    double rate = -D0.get(currentState, currentState);
                    double time = -Math.log(random.nextDouble()) / rate;
                    samples[(int) i] += time;
                    for (int k = 0; k < nphases; k++) {
                        row[k] = D0.get(currentState, k);
                        row[nphases + k] = D1.get(currentState, k);
                    }
                    row[currentState] = 0.0;
                    int nextState = 2 * nphases - 1;
                    sum = 0.0;
                    r = random.nextDouble();
                    for (int j = 0; j < 2 * nphases; j++) {
                        sum += row[j] / rate;
                        if (r < sum) {
                            if (j >= nphases) {
                                nextState = j - nphases;
                                continueSample = false;
                                break;
                            } else {
                                nextState = j;
                                continueSample = true;
                                break;
                            }
                        }
                    }
                    currentState = nextState;
                }
            }
        }
        return samples;
    }

    /**
     * Stateful single-event MAP sampler that carries the modulating phase across calls.
     *
     * <p>Successive calls to {@link #next(Random)} continue the MAP's underlying CTMC
     * from the phase reached at the previous event, reproducing the process
     * autocorrelation. This differs from repeatedly invoking
     * {@link #map_sample(Matrix, Matrix, long, Random)} with {@code n=1}, which
     * re-initialises the phase from the stationary vector {@code map_pie} on every
     * call and therefore yields an i.i.d. (renewal) sequence with the correct
     * marginal but zero autocorrelation. Use this sampler when a MAP/MMPP2/RAP is
     * an arrival or service process whose burstiness must be preserved across
     * events (e.g. an M/MAP/1 queue); the renewal path is correct only for
     * phase-type renewal distributions whose {@code D1} restarts from a fixed entry
     * vector.</p>
     *
     * <p>The phase persists across idle periods (it is only advanced when
     * {@link #next(Random)} is called), matching the MAP/MAP/1 convention in which
     * the service process is frozen while the server is idle.</p>
     */
    public static final class MapSampler {
        private final Matrix D0;
        private final Matrix D1;
        private final int nphases;
        private final boolean exponential;
        private final double lambda;
        private int currentState = -1; // -1 => draw initial phase from map_pie on first call

        public MapSampler(Matrix D0, Matrix D1) {
            this.D0 = D0;
            this.D1 = D1;
            this.exponential = D0.getNumElements() == 1;
            this.nphases = D0.getNumCols();
            this.lambda = this.exponential ? D1.value() : 0.0;
        }

        /** Draws one inter-event time, advancing (and retaining) the modulating phase. */
        public double next(Random random) {
            if (exponential) {
                return -Math.log(random.nextDouble()) / lambda;
            }
            if (currentState < 0) {
                Matrix pie = Map_pie.map_pie(D0, D1);
                double sum = 0.0;
                double r = random.nextDouble();
                currentState = nphases - 1;
                for (int i = 0; i < nphases; i++) {
                    sum += pie.get(0, i);
                    if (r < sum) {
                        currentState = i;
                        break;
                    }
                }
            }
            double sample = 0.0;
            double[] row = new double[2 * nphases];
            boolean continueSample = true;
            while (continueSample) {
                double rate = -D0.get(currentState, currentState);
                sample += -Math.log(random.nextDouble()) / rate;
                for (int k = 0; k < nphases; k++) {
                    row[k] = D0.get(currentState, k);
                    row[nphases + k] = D1.get(currentState, k);
                }
                row[currentState] = 0.0;
                int nextState = 2 * nphases - 1;
                double sum = 0.0;
                double r = random.nextDouble();
                for (int j = 0; j < 2 * nphases; j++) {
                    sum += row[j] / rate;
                    if (r < sum) {
                        if (j >= nphases) {
                            nextState = j - nphases;
                            continueSample = false;
                        } else {
                            nextState = j;
                            continueSample = true;
                        }
                        break;
                    }
                }
                currentState = nextState;
            }
            return sample;
        }
    }

    /**
     * Stateful sampler for a BMAP (Batch Markovian Arrival Process) that retains
     * the modulating phase between draws.
     *
     * <p>Successive calls to {@link #next(Random)} continue the BMAP's underlying
     * CTMC from the phase reached at the previous batch arrival, reproducing both
     * the inter-arrival and batch-size autocorrelation. This differs from
     * repeatedly invoking {@link #bmap_sample(MatrixCell, long, Random)} with
     * {@code n=1}, which re-initialises the phase from the stationary vector
     * {@code map_pie} on every call and therefore yields an i.i.d. (renewal)
     * batch-arrival sequence with the correct marginals but zero autocorrelation.
     * Use this sampler when a BMAP is an arrival process whose burstiness must be
     * preserved across events.</p>
     *
     * <p>The phase persists across idle periods (it is only advanced when
     * {@link #next(Random)} is called), matching the MAP/MAP/1 convention in which
     * the arrival process is not frozen while downstream stations are busy.</p>
     */
    public static final class BmapSampler {
        private final Matrix D0;
        private final Matrix D1_total;
        private final Matrix[] Dk;
        private final int nphases;
        private final int maxBatchSize;
        private final boolean exponential;
        private final double lambda;
        private final double[] batchRates;
        private final double totalRate;
        private int currentState = -1; // -1 => draw initial phase from map_pie on first call

        public BmapSampler(MatrixCell bmap) {
            this.D0 = bmap.get(0);
            this.D1_total = bmap.get(1);
            this.maxBatchSize = bmap.size() - 2;
            this.Dk = new Matrix[maxBatchSize];
            for (int k = 0; k < maxBatchSize; k++) {
                this.Dk[k] = bmap.get(k + 2);
            }
            this.nphases = D0.getNumCols();
            this.exponential = D0.getNumElements() == 1;
            this.lambda = this.exponential ? D1_total.value() : 0.0;
            this.batchRates = new double[maxBatchSize];
            double sumRate = 0.0;
            if (this.exponential) {
                for (int k = 0; k < maxBatchSize; k++) {
                    this.batchRates[k] = Dk[k].value();
                    sumRate += this.batchRates[k];
                }
            }
            this.totalRate = sumRate;
        }

        /**
         * Draws one BMAP event (inter-arrival time and batch size), advancing (and
         * retaining) the modulating phase.
         */
        public BmapSample next(Random random) {
            if (exponential) {
                double interarrivalTime = -Math.log(random.nextDouble()) / lambda;
                int batchSize = 1;
                if (totalRate > 0) {
                    double r = random.nextDouble() * totalRate;
                    double cumSum = 0.0;
                    for (int k = 0; k < maxBatchSize; k++) {
                        cumSum += batchRates[k];
                        if (r < cumSum) {
                            batchSize = k + 1;
                            break;
                        }
                    }
                }
                return new BmapSample(interarrivalTime, batchSize);
            }

            if (currentState < 0) {
                Matrix pie = Map_pie.map_pie(D0, D1_total);
                double sum = 0.0;
                double r = random.nextDouble();
                currentState = nphases - 1;
                for (int i = 0; i < nphases; i++) {
                    sum += pie.get(0, i);
                    if (r < sum) {
                        currentState = i;
                        break;
                    }
                }
            }

            int totalCols = nphases + maxBatchSize * nphases;
            double[] row = new double[totalCols];
            boolean continuesSample = true;
            double interarrivalTime = 0.0;
            int batchSize = 1;

            while (continuesSample) {
                double rate = -D0.get(currentState, currentState);
                interarrivalTime += -Math.log(random.nextDouble()) / rate;

                for (int k = 0; k < nphases; k++) {
                    row[k] = D0.get(currentState, k);
                }
                row[currentState] = 0.0; // No self-transition
                for (int b = 0; b < maxBatchSize; b++) {
                    for (int k = 0; k < nphases; k++) {
                        row[nphases + b * nphases + k] = Dk[b].get(currentState, k);
                    }
                }

                int nextState = currentState;
                double sum = 0.0;
                double r = random.nextDouble();
                for (int j = 0; j < totalCols; j++) {
                    sum += row[j] / rate;
                    if (r < sum) {
                        if (j < nphases) {
                            nextState = j;
                            continuesSample = true;
                        } else {
                            int idx = j - nphases;
                            batchSize = (idx / nphases) + 1;
                            nextState = idx % nphases;
                            continuesSample = false;
                        }
                        break;
                    }
                }
                currentState = nextState;
            }

            return new BmapSample(interarrivalTime, batchSize);
        }
    }

    /**
     * Generates samples from a BMAP (Batch Markovian Arrival Process).
     */
    public static BmapSample[] bmap_sample(MatrixCell bmap, long n, Random random) {
        Matrix D0 = bmap.get(0);
        Matrix D1_total = bmap.get(1);
        int maxBatchSize = bmap.size() - 2;

        // Collect individual Dk matrices
        Matrix[] Dk = new Matrix[maxBatchSize];
        for (int k = 0; k < maxBatchSize; k++) {
            Dk[k] = bmap.get(k + 2);
        }

        return bmap_sample(D0, D1_total, Dk, n, random);
    }

    /**
     * Generates samples from a BMAP using D0, D1_total, and individual Dk matrices.
     */
    public static BmapSample[] bmap_sample(Matrix D0, Matrix D1_total, Matrix[] Dk, long n, Random random) {
        int nphases = D0.getNumCols();
        int maxBatchSize = Dk.length;

        // Handle degenerate case: single phase (exponential)
        if (D0.getNumElements() == 1) {
            double lambda = D1_total.value();
            // For single phase, determine batch size distribution
            double[] batchRates = new double[maxBatchSize];
            double totalRate = 0.0;
            for (int k = 0; k < maxBatchSize; k++) {
                batchRates[k] = Dk[k].value();
                totalRate += batchRates[k];
            }

            BmapSample[] result = new BmapSample[(int) n];
            for (int i = 0; i < n; i++) {
                double interarrivalTime = -Math.log(random.nextDouble()) / lambda;

                // Select batch size proportionally to rates
                int batchSize = 1;
                if (totalRate > 0) {
                    double r = random.nextDouble() * totalRate;
                    double cumSum = 0.0;
                    for (int k = 0; k < maxBatchSize; k++) {
                        cumSum += batchRates[k];
                        if (r < cumSum) {
                            batchSize = k + 1;
                            break;
                        }
                    }
                }
                result[i] = new BmapSample(interarrivalTime, batchSize);
            }
            return result;
        }

        // Multi-phase BMAP sampling
        BmapSample[] samples = new BmapSample[(int) n];
        for (int i = 0; i < n; i++) samples[i] = new BmapSample(0.0, 1);

        // Initialize state from stationary distribution
        Matrix pie = Map_pie.map_pie(D0, D1_total);
        int currentState = nphases - 1;
        double sum = 0.0;
        double r = random.nextDouble();
        for (int i = 0; i < nphases; i++) {
            sum += pie.get(0, i);
            if (r < sum) {
                currentState = i;
                break;
            }
        }

        // Build combined transition row: [D0 transitions | D1 | D2 | ... | Dk]
        // Total columns = nphases + maxBatchSize * nphases
        int totalCols = nphases + maxBatchSize * nphases;
        double[] row = new double[totalCols];

        for (int i = 0; i < (int) n; i++) {
            boolean continuesSample = true;
            double interarrivalTime = 0.0;
            int batchSize = 1;

            while (continuesSample) {
                double rate = -D0.get(currentState, currentState);
                double time = -Math.log(random.nextDouble()) / rate;
                interarrivalTime += time;

                // Build probability row
                for (int k = 0; k < nphases; k++) {
                    row[k] = D0.get(currentState, k);
                }
                row[currentState] = 0.0;  // No self-transition

                // Next entries: Dk transitions for each batch size
                for (int b = 0; b < maxBatchSize; b++) {
                    for (int k = 0; k < nphases; k++) {
                        row[nphases + b * nphases + k] = Dk[b].get(currentState, k);
                    }
                }

                // Select next state and batch size
                int nextState = currentState;
                sum = 0.0;
                r = random.nextDouble();

                for (int j = 0; j < totalCols; j++) {
                    sum += row[j] / rate;
                    if (r < sum) {
                        if (j < nphases) {
                            // D0 transition - no arrival, continue sampling
                            nextState = j;
                            continuesSample = true;
                        } else {
                            // Dk transition - arrival with batch size
                            int idx = j - nphases;
                            batchSize = (idx / nphases) + 1;
                            nextState = idx % nphases;
                            continuesSample = false;
                        }
                        break;
                    }
                }
                currentState = nextState;
            }

            samples[i] = new BmapSample(interarrivalTime, batchSize);
        }

        return samples;
    }
}
