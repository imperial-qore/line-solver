package jline.opt.solver.de;

/**
 * Bit-exact reimplementation of numpy's legacy MT19937 core, matching the
 * generator underlying {@code numpy.random.RandomState}.
 *
 * <p>This reproduces the exact 32-bit output stream numpy produces so that
 * {@link NumpyRandomState} and the differential-evolution optimizer generate
 * identical draws to the native-Python line-opt implementation for a given
 * seed. Seeding follows numpy's {@code _legacy_seeding}: a scalar seed that
 * fits in 32 bits uses {@code init_genrand}; otherwise the seed words are fed
 * to {@code init_by_array}.</p>
 *
 * <p>Java 8 compatible. Unsigned 32-bit words are held as {@code int}; callers
 * that need the unsigned value use {@code &amp; 0xFFFFFFFFL}.</p>
 */
public class MT19937 {

    private static final int N = 624;
    private static final int M = 397;
    private static final int MATRIX_A = 0x9908b0df;
    private static final int UPPER_MASK = 0x80000000;
    private static final int LOWER_MASK = 0x7fffffff;

    private final int[] mt = new int[N];
    private int mti;

    /** Seed with a single non-negative integer, matching RandomState(int). */
    public MT19937(long seed) {
        seed(seed);
    }

    /** Seed matching numpy legacy behavior for a scalar integer seed. */
    public void seed(long seed) {
        if (seed < 0) {
            throw new IllegalArgumentException("seed must be non-negative");
        }
        if (seed <= 0xFFFFFFFFL) {
            initGenrand(seed & 0xFFFFFFFFL);
        } else {
            // Split into little-endian 32-bit words and init_by_array.
            java.util.List<Long> words = new java.util.ArrayList<Long>();
            long s = seed;
            while (s > 0) {
                words.add(s & 0xFFFFFFFFL);
                s >>>= 32;
            }
            long[] key = new long[words.size()];
            for (int i = 0; i < key.length; i++) {
                key[i] = words.get(i);
            }
            initByArray(key);
        }
    }

    /** numpy/Matsumoto init_genrand. */
    public void initGenrand(long s) {
        mt[0] = (int) (s & 0xFFFFFFFFL);
        for (int i = 1; i < N; i++) {
            long prev = mt[i - 1] & 0xFFFFFFFFL;
            long val = (1812433253L * (prev ^ (prev >>> 30)) + i) & 0xFFFFFFFFL;
            mt[i] = (int) val;
        }
        mti = N;
    }

    /** numpy/Matsumoto init_by_array (used for seeds &gt; 32 bits). */
    public void initByArray(long[] initKey) {
        initGenrand(19650218L);
        int i = 1;
        int j = 0;
        int k = Math.max(N, initKey.length);
        for (; k > 0; k--) {
            long prev = mt[i - 1] & 0xFFFFFFFFL;
            long val = ((mt[i] & 0xFFFFFFFFL)
                    ^ ((prev ^ (prev >>> 30)) * 1664525L)) + initKey[j] + j;
            mt[i] = (int) (val & 0xFFFFFFFFL);
            i++;
            j++;
            if (i >= N) {
                mt[0] = mt[N - 1];
                i = 1;
            }
            if (j >= initKey.length) {
                j = 0;
            }
        }
        for (k = N - 1; k > 0; k--) {
            long prev = mt[i - 1] & 0xFFFFFFFFL;
            long val = ((mt[i] & 0xFFFFFFFFL)
                    ^ ((prev ^ (prev >>> 30)) * 1566083941L)) - i;
            mt[i] = (int) (val & 0xFFFFFFFFL);
            i++;
            if (i >= N) {
                mt[0] = mt[N - 1];
                i = 1;
            }
        }
        mt[0] = 0x80000000;
        mti = N;
    }

    /** Draw the next tempered 32-bit word as an unsigned value in [0, 2^32). */
    public long nextUint32() {
        int y;
        if (mti >= N) {
            int kk;
            for (kk = 0; kk < N - M; kk++) {
                y = (mt[kk] & UPPER_MASK) | (mt[kk + 1] & LOWER_MASK);
                mt[kk] = mt[kk + M] ^ (y >>> 1) ^ ((y & 1) != 0 ? MATRIX_A : 0);
            }
            for (; kk < N - 1; kk++) {
                y = (mt[kk] & UPPER_MASK) | (mt[kk + 1] & LOWER_MASK);
                mt[kk] = mt[kk + (M - N)] ^ (y >>> 1) ^ ((y & 1) != 0 ? MATRIX_A : 0);
            }
            y = (mt[N - 1] & UPPER_MASK) | (mt[0] & LOWER_MASK);
            mt[N - 1] = mt[M - 1] ^ (y >>> 1) ^ ((y & 1) != 0 ? MATRIX_A : 0);
            mti = 0;
        }

        y = mt[mti++];
        y ^= (y >>> 11);
        y ^= (y << 7) & 0x9d2c5680;
        y ^= (y << 15) & 0xefc60000;
        y ^= (y >>> 18);
        return y & 0xFFFFFFFFL;
    }

    /** Copy of the internal 624-word key (for state validation/tests). */
    public int[] getStateKey() {
        return mt.clone();
    }

    /** Current position into the 624-word key (numpy get_state pos). */
    public int getPos() {
        return mti;
    }
}
