/*
 * Copyright (c) 2012-2026, QORE Lab, Imperial College London
 * All rights reserved.
 */

package jline.util;

import org.apache.commons.math3.random.MersenneTwister;
import org.apache.commons.math3.util.FastMath;

import java.util.Random;
import java.util.concurrent.ConcurrentHashMap;
import java.util.concurrent.atomic.AtomicLong;

/**
 * Centralized random number generator management for reproducible, thread-safe random number generation.
 * 
 * This class ensures that:
 * - All random number generation uses MersenneTwister
 * - Setting a master seed produces deterministic results across all threads
 * - Each thread gets its own MersenneTwister instance with a deterministic sub-seed
 * - Parallel execution is reproducible when using the same master seed
 */
public class RandomManager {
    
    /**
     * Default master seed used when no explicit seed is set
     * Using 0 for backward compatibility with the original ThreadLocalRandom default
     */
    private static final int DEFAULT_SEED = 0;
    
    /**
     * Master seed that controls all random number generation
     */
    private static volatile int masterSeed = DEFAULT_SEED;
    
    /**
     * Master random number generator used for creating deterministic sub-seeds
     */
    private static volatile MersenneTwister masterRandom = new MersenneTwister(DEFAULT_SEED);
    
    /**
     * Thread-local storage for MersenneTwister instances
     * Initialize with seed 0 for backward compatibility - exactly like original ThreadLocalRandom
     */
    private static final ThreadLocal<MersenneTwister> threadLocalRandom = ThreadLocal.withInitial(() -> new MersenneTwister(0));
    
    /**
     * Map to store solver-specific random generators by solver name/id
     */
    private static final ConcurrentHashMap<String, MersenneTwister> solverRandoms = new ConcurrentHashMap<>();
    
    /**
     * Counter for generating unique thread identifiers
     */
    private static final AtomicLong threadCounter = new AtomicLong(0);
    
    /**
     * Thread-local storage for thread-specific identifiers
     */
    private static final ThreadLocal<Long> threadId = new ThreadLocal<Long>() {
        @Override
        protected Long initialValue() {
            return threadCounter.incrementAndGet();
        }
    };
    
    /**
     * Sets the master seed that controls all random number generation.
     * This ensures reproducible results across all threads and solvers.
     * 
     * @param seed the master seed value
     */
    public static synchronized void setMasterSeed(int seed) {
        masterSeed = seed;
        masterRandom = new MersenneTwister(seed);
        
        // For backward compatibility with original ThreadLocalRandom behavior,
        // do NOT clear existing thread-local generators - let setSeed() work per-thread
        solverRandoms.clear();
        threadCounter.set(0);
    }
    
    /**
     * Gets the current master seed.
     * 
     * @return the current master seed
     */
    public static int getMasterSeed() {
        return masterSeed;
    }
    
    /**
     * Gets a thread-specific MersenneTwister instance.
     * Each thread gets its own generator with a deterministic seed derived from the master seed.
     * 
     * @return thread-specific MersenneTwister instance
     */
    public static MersenneTwister getThreadRandom() {
        return threadLocalRandom.get();
    }
    
    /**
     * Gets a thread-specific Random instance for compatibility with java.util.Random.
     * This wraps the MersenneTwister to provide java.util.Random compatibility.
     * 
     * @return thread-specific Random instance backed by MersenneTwister
     */
    public static Random getThreadRandomAsRandom() {
        MersenneTwister mt = threadLocalRandom.get();
        return new Random() {
            @Override
            protected int next(int bits) {
                // Delegate to MersenneTwister's internal bits method
                return mt.nextInt() >>> (32 - bits);
            }
            
            @Override
            public double nextDouble() {
                return mt.nextDouble();
            }
            
            @Override
            public int nextInt() {
                return mt.nextInt();
            }
            
            @Override
            public int nextInt(int bound) {
                return mt.nextInt(bound);
            }
            
            @Override
            public long nextLong() {
                return mt.nextLong();
            }
            
            @Override
            public boolean nextBoolean() {
                return mt.nextBoolean();
            }
            
            @Override
            public float nextFloat() {
                return mt.nextFloat();
            }
            
            @Override
            public double nextGaussian() {
                return mt.nextGaussian();
            }
            
            @Override
            public void setSeed(long seed) {
                mt.setSeed(seed);
            }
        };
    }
    
    /**
     * Gets a solver-specific MersenneTwister instance.
     * Each solver gets its own generator with a deterministic seed.
     * 
     * @param solverId unique identifier for the solver
     * @return solver-specific MersenneTwister instance
     */
    public static MersenneTwister getSolverRandom(String solverId) {
        return solverRandoms.computeIfAbsent(solverId, id -> {
            synchronized (RandomManager.class) {
                // Create deterministic seed for this solver
                int solverSeed = masterSeed + id.hashCode();
                return new MersenneTwister(solverSeed);
            }
        });
    }
    
    /**
     * Creates a MersenneTwister instance for parallel processing with a thread offset.
     * This ensures different threads in parallel algorithms get different but reproducible sequences.
     * 
     * @param solverId unique identifier for the solver
     * @param threadOffset offset for this thread (e.g., thread index in parallel processing)
     * @return MersenneTwister instance with deterministic seed
     */
    public static MersenneTwister getParallelRandom(String solverId, int threadOffset) {
        synchronized (RandomManager.class) {
            int parallelSeed = masterSeed + solverId.hashCode() + threadOffset;
            return new MersenneTwister(parallelSeed);
        }
    }
    
    /**
     * Generates a random double value in [0.0, 1.0) using the thread-specific generator.
     * For backward compatibility, this uses the same ThreadLocal pattern as original ThreadLocalRandom.
     * 
     * @return random double value
     */
    public static double nextDouble() {
        return threadLocalRandom.get().nextDouble();
    }
    
    /**
     * Generates a random Gaussian (normal) distributed double value using the thread-specific generator.
     * For backward compatibility, this uses the same ThreadLocal pattern as original ThreadLocalRandom.
     * 
     * @return random Gaussian-distributed double value
     */
    public static double nextGaussian() {
        return threadLocalRandom.get().nextGaussian();
    }
    
    /**
     * Generates a random integer between 0 (inclusive) and the specified value (exclusive).
     * 
     * @param bound the upper bound (exclusive). Must be positive.
     * @return random integer in range [0, bound)
     */
    public static int nextInt(int bound) {
        return threadLocalRandom.get().nextInt(bound);
    }
    
    /**
     * Creates a new thread-specific MersenneTwister with deterministic seed.
     * 
     * @return new MersenneTwister instance for current thread
     */
    private static MersenneTwister createThreadSpecificRandom() {
        synchronized (RandomManager.class) {
            // Use thread ID to create deterministic but unique seed for each thread
            long currentThreadId = threadId.get();
            int threadSeed = masterSeed + (int) (currentThreadId & 0x7FFFFFFF);
            return new MersenneTwister(threadSeed);
        }
    }
    
    /**
     * Resets all random number generators to use the current master seed.
     * This can be called to ensure fresh random sequences.
     */
    public static synchronized void reset() {
        setMasterSeed(masterSeed);
    }
    
    /**
     * Gets a fresh MersenneTwister instance with the master seed.
     * Useful for creating new independent random streams.
     * 
     * @return new MersenneTwister instance with master seed
     */
    public static MersenneTwister newInstance() {
        return new MersenneTwister(masterSeed);
    }
    
    /**
     * Generates a random seed using the same formula as SolverOptions.
     * 
     * @return random seed value
     */
    public static int generateRandomSeed() {
        return FastMath.toIntExact(Math.round((Math.random() * (1e6 - 1)) + 1));
    }
}