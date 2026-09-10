package jline.opt.solver.de;

/**
 * Bit-exact port of the subset of {@code numpy.random.RandomState} consumed by
 * scipy's {@code differential_evolution}. Backed by {@link MT19937}.
 *
 * <p>Implements {@code random_sample}, {@code uniform}, the default-dtype
 * {@code randint} (32-bit masked rejection), and {@code shuffle}/
 * {@code permutation} (Fisher-Yates driven by {@code random_interval}), each
 * reproducing numpy's exact draw sequence and word consumption. This lets the
 * differential-evolution optimizer follow the identical trajectory to the
 * native-Python line-opt for a given integer seed.</p>
 *
 * <p>Java 8 compatible.</p>
 */
public class NumpyRandomState {

    private final MT19937 mt;

    public NumpyRandomState(long seed) {
        this.mt = new MT19937(seed);
    }

    public NumpyRandomState(MT19937 mt) {
        this.mt = mt;
    }

    public MT19937 getGenerator() {
        return mt;
    }

    /** numpy random_sample(): a 53-bit double in [0, 1) from two 32-bit words. */
    public double randomSample() {
        long a = mt.nextUint32() >>> 5;   // 27 bits
        long b = mt.nextUint32() >>> 6;   // 26 bits
        return (a * 67108864.0 + b) / 9007199254740992.0;
    }

    /** numpy random_sample(size): C-order fill of length n. */
    public double[] randomSample(int n) {
        double[] out = new double[n];
        for (int i = 0; i < n; i++) {
            out[i] = randomSample();
        }
        return out;
    }

    /** numpy uniform(low, high): scalar draw. */
    public double uniform(double low, double high) {
        return low + (high - low) * randomSample();
    }

    /** numpy uniform(low, high, size=n): C-order fill of length n. */
    public double[] uniform(double low, double high, int n) {
        double[] out = new double[n];
        for (int i = 0; i < n; i++) {
            out[i] = low + (high - low) * randomSample();
        }
        return out;
    }

    /** numpy uniform(size=n) == uniform(0, 1, n). */
    public double[] uniform(int n) {
        return uniform(0.0, 1.0, n);
    }

    /** Smallest 2^m - 1 that is &gt;= v (32-bit bit-fill). */
    private static long fillMask(long v) {
        v |= v >>> 1;
        v |= v >>> 2;
        v |= v >>> 4;
        v |= v >>> 8;
        v |= v >>> 16;
        return v;
    }

    /**
     * numpy randint(low, high) with the default integer dtype: a uniform
     * integer in [low, high) via 32-bit masked rejection (one word per
     * rejection iteration), matching numpy's word consumption exactly.
     */
    public long randint(long low, long high) {
        long rng = high - 1 - low;      // inclusive upper range
        if (rng == 0) {
            return low;
        }
        long mask = fillMask(rng);
        while (true) {
            long val = mt.nextUint32() & mask;
            if (Long.compareUnsigned(val, rng) <= 0) {
                return low + val;
            }
        }
    }

    /** Convenience: randint(0, high). */
    public long randint(long high) {
        return randint(0L, high);
    }

    /**
     * numpy random_interval(max): a uniform integer in [0, max] inclusive via
     * 32-bit masked rejection. Used by shuffle/permutation.
     */
    public long randomInterval(long max) {
        if (max == 0) {
            return 0;
        }
        long mask = fillMask(max);
        while (true) {
            long value = mt.nextUint32() & mask;
            if (Long.compareUnsigned(value, max) <= 0) {
                return value;
            }
        }
    }

    /** numpy shuffle: in-place Fisher-Yates using random_interval. */
    public void shuffle(int[] arr) {
        for (int i = arr.length - 1; i > 0; i--) {
            int j = (int) randomInterval(i);
            int tmp = arr[i];
            arr[i] = arr[j];
            arr[j] = tmp;
        }
    }

    /** numpy permutation(range(n)): shuffled copy of [0, 1, ..., n-1]. */
    public int[] permutation(int n) {
        int[] arr = new int[n];
        for (int i = 0; i < n; i++) {
            arr[i] = i;
        }
        shuffle(arr);
        return arr;
    }
}
