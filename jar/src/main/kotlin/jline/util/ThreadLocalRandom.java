/*
 * Copyright (c) 2012-2026, QORE Lab, Imperial College London
 * All rights reserved.
 */

package jline.util;

/**
 * A class for random number generation in separate threads.
 * This class now delegates to RandomManager for centralized, reproducible random number generation.
 * 
 * @deprecated Use RandomManager directly for new code. This class is maintained for backward compatibility.
 */

public class ThreadLocalRandom {

    /**
     * Generates a random double value in [0.0, 1.0) using the thread-specific MersenneTwister.
     * 
     * @return random double value
     */
    public static double random() {
        return RandomManager.nextDouble();
    }

    /**
     * Sets the seed for the current thread's random number generator.
     * For backward compatibility, this only affects the current thread.
     * 
     * @param seed the seed value
     */
    public static void setSeed(final int seed) {
        // For backward compatibility, set seed only for current thread
        RandomManager.getThreadRandom().setSeed(seed);
    }
}
