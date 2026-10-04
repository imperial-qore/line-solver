/**
 * @file Marked Markovian Arrival Process sample generation
 *
 * @since LINE 3.0
 */
package jline.api.mam;

import java.util.Random;

import jline.io.Ret;
import jline.lang.processes.DiscreteSampler;
import jline.util.matrix.Matrix;
import jline.util.matrix.MatrixCell;

public final class Mmap_sample {
    private Mmap_sample() {}

    /**
     * Generates samples of inter-arrival times and event types from a MMAP.
     */
    public static Ret.mamMMAPSample mmap_sample(MatrixCell MMAP, long n, Random random) {
        Matrix D0 = MMAP.get(0);
        Matrix D1 = MMAP.get(1);
        int C = MMAP.size() - 2;
        int order = D0.getNumRows();
        DiscreteSampler[][] typeSampler = new DiscreteSampler[order][order];
        Matrix x = new Matrix(1, C);
        for (int c = 0; c < C; c++) {
            x.set(c, (double) c);
        }
        for (int i = 0; i < order; i++) {
            for (int j = 0; j < order; j++) {
                if (D1.get(i, j) > 0) {
                    Matrix pij = new Matrix(1, C);
                    for (int c = 0; c < C; c++) {
                        pij.set(c, MMAP.get(2 + c).get(i, j));
                    }
                    pij.scaleEq(1 / pij.elementSum());
                    if (pij.elementSum() > 0) {
                        typeSampler[i][j] = new DiscreteSampler(pij, x);
                    }
                }
            }
        }
        double[] samples = new double[(int) n];
        int[] types = new int[(int) n];
        int[] states = new int[(int) n];
        if (D0.getNumElements() == 1) {
            double lambda = D1.value();
            for (long i = 0; i < n; i++) {
                samples[(int) i] = -Math.log(random.nextDouble()) / lambda;
                types[(int) i] = 1;
                states[(int) i] = 0;
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
                boolean continue_sample = true;
                samples[(int) i] = 0.0;
                states[(int) i] = currentState;
                while (continue_sample) {
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
                                continue_sample = false;
                                types[(int) i] = (int) typeSampler[currentState][nextState].sample();
                                break;
                            } else {
                                nextState = j;
                                continue_sample = true;
                                break;
                            }
                        }
                    }
                    currentState = nextState;
                }
            }
        }
        return new Ret.mamMMAPSample(samples, C, types, states);
    }

    public static Ret.mamMMAPSample mmap_sample(MatrixCell MMAP, long n) {
        return mmap_sample(MMAP, n, new Random());
    }

    /** One marked draw: inter-arrival time plus the 1-based mark of the arrival. */
    public static final class MarkedSample {
        public final double interarrivalTime;
        public final int mark; // 1-based mark index (1..K)

        MarkedSample(double interarrivalTime, int mark) {
            this.interarrivalTime = interarrivalTime;
            this.mark = mark;
        }
    }

    /**
     * Stateful sampler for an MMAP {D0, D1_agg, D11..D1K} (M3A layout) that
     * retains the modulating phase between draws, so successive calls
     * reproduce the inter-arrival autocorrelation (mirrors
     * Map_sample.MapSampler / BmapSampler). Each draw also returns the
     * 1-based mark of the arrival, sampled with probability
     * D1k(i,j)/D1_agg(i,j) on the firing transition.
     */
    public static final class MmapSampler {
        private final Matrix D0;
        private final Matrix D1;
        private final MatrixCell mmap;
        private final int C;
        private final int nphases;
        private final boolean exponential;
        private final double lambda;
        private int currentState = -1; // -1 => draw initial phase from map_pie on first call

        public MmapSampler(MatrixCell mmap) {
            this.mmap = mmap;
            this.D0 = mmap.get(0);
            this.D1 = mmap.get(1);
            this.C = Math.max(0, mmap.size() - 2);
            this.nphases = D0.getNumCols();
            this.exponential = D0.getNumElements() == 1;
            this.lambda = this.exponential ? D1.value() : 0.0;
        }

        /** Samples the 1-based mark of an arrival firing on transition (i,j). */
        private int sampleMark(int i, int j, Random random) {
            if (C <= 0) {
                return 1;
            }
            double total = 0.0;
            for (int c = 0; c < C; c++) {
                total += Math.max(0.0, mmap.get(2 + c).get(i, j));
            }
            if (total <= 0.0) {
                return 1;
            }
            double r = random.nextDouble() * total;
            double sum = 0.0;
            for (int c = 0; c < C; c++) {
                sum += Math.max(0.0, mmap.get(2 + c).get(i, j));
                if (r < sum) {
                    return c + 1;
                }
            }
            return C;
        }

        /** Draws one marked inter-event time, advancing (and retaining) the phase. */
        public MarkedSample next(Random random) {
            if (exponential) {
                double time = -Math.log(random.nextDouble()) / lambda;
                return new MarkedSample(time, sampleMark(0, 0, random));
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
            while (true) {
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
                        nextState = j;
                        break;
                    }
                }
                if (nextState >= nphases) {
                    int destPhase = nextState - nphases;
                    int mark = sampleMark(currentState, destPhase, random);
                    currentState = destPhase;
                    return new MarkedSample(sample, mark);
                }
                currentState = nextState;
            }
        }
    }
}
