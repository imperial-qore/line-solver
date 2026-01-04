package jline.api.fes

import jline.util.matrix.Matrix

/**
 * Options for Flow-Equivalent Server (FES) aggregation
 *
 * Copyright (c) 2012-2026, Imperial College London
 * All rights reserved.
 */
data class FESOptions(
    /** Solver to use for throughput computation ('mva' default) */
    val solver: String = "mva",

    /** Per-class population cutoffs (default: null, uses total jobs per class) */
    val cutoffs: Matrix? = null,

    /** Enable verbose output */
    val verbose: Boolean = false
) {
    companion object {
        /**
         * Create default options
         */
        @JvmStatic
        fun defaults(): FESOptions {
            return FESOptions()
        }

        /**
         * Create options with specified solver
         */
        @JvmStatic
        fun withSolver(solver: String): FESOptions {
            return FESOptions(solver = solver)
        }

        /**
         * Create options with specified cutoffs
         */
        @JvmStatic
        fun withCutoffs(cutoffs: Matrix): FESOptions {
            return FESOptions(cutoffs = cutoffs)
        }
    }
}
