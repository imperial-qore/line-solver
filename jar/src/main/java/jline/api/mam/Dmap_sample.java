/**
 * Sample generation for discrete-time MAPs.
 *
 * @since LINE 3.0
 */
package jline.api.mam;

import jline.lib.butools.dmap.SamplesFromDMAP;
import jline.util.matrix.Matrix;
import jline.util.matrix.MatrixCell;

import java.util.Random;

public final class Dmap_sample {
    private Dmap_sample() {}

    /**
     * Generates samples of inter-arrival times from a discrete-time MAP.
     * Inter-arrival times are integer-valued, returned as doubles.
     */
    public static double[] dmap_sample(Matrix D0, Matrix D1, int n, Random random) {
        Random rng = (random != null) ? random : new Random();
        int[] intSamples = SamplesFromDMAP.samplesFromDMAP(D0, D1, n, null, 1e-14, rng);
        double[] out = new double[intSamples.length];
        for (int i = 0; i < intSamples.length; i++) {
            out[i] = (double) intSamples[i];
        }
        return out;
    }

    public static double[] dmap_sample(MatrixCell DMAP, long n, Random random) {
        return dmap_sample(DMAP.get(0), DMAP.get(1), (int) n, random);
    }
}
