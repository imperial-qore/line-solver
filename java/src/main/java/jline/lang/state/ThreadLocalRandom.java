package jline.lang.state;

import org.apache.commons.math3.random.MersenneTwister;

import java.util.Random;
import java.util.concurrent.atomic.AtomicReference;

public class ThreadLocalRandom {
    private static ThreadLocal<MersenneTwister> threadLocalRandom = ThreadLocal.withInitial(() -> new MersenneTwister(0));

    public static void setSeed(final int seed) {
        threadLocalRandom.get().setSeed(seed);
    }

    public static double random() {
        return threadLocalRandom.get().nextDouble();
    }
}
