package jline.lang.state;

import org.apache.commons.math3.random.MersenneTwister;

import java.util.Random;

public class ThreadLocalRandom {
    private static ThreadLocal<MersenneTwister> threadLocalRandom = ThreadLocal.withInitial(() -> new MersenneTwister(1));

    public static void setSeed(long seed) {
        threadLocalRandom.set(new MersenneTwister(1));
    }

    public static double random() {
        return threadLocalRandom.get().nextDouble();
    }
}