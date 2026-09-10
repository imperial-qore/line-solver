/*
 * Copyright (c) 2012-2026, QORE Lab, Imperial College London
 * All rights reserved.
 */

package jline.solvers.tr;

import jline.solvers.SolverOptions;

/**
 * Method name table for {@link TransformSolve}, mirroring MATLAB
 * {@code matlab/src/solvers/TR/transform_method.m} and python
 * {@code transform_driver.transform_method}.
 *
 * <p>A method name names a model transformation: the strategy rewrites the model into
 * subproblems, those are solved by a real solver, and the metrics are mapped
 * back. An unrecognised token is an ERROR rather than a silent {@code none},
 * because a mistyped transform that quietly solved the untransformed model
 * would report a plausible number for the wrong problem.
 *
 * <p>NOT every transformation in the tree is reachable here. The fork-join tag
 * augmentation ({@link FJTagTransform}) is selected STRUCTURALLY, from the
 * presence of a Fork or Join node, and never by a token; the MMT/HT fixed point
 * keeps its own {@code options.config.fork_join}.
 */
public final class TransformMethod {

    /** No transformation; the caller solves its own model directly. */
    public static final String NONE = "none";

    /** Collapse every chain onto a single class, then deaggregate through alpha. */
    public static final String CHAINS = "chains";

    /** Load concealment (Birman-Kogan Algorithm 2), the first ITERATED strategy. */
    public static final String LC = "lc";

    private TransformMethod() {
    }

    /**
     * Normalises a transform method name onto its canonical name.
     *
     * @param token the raw {@code options.config.transform} value, possibly null
     * @return {@link #NONE} or {@link #CHAINS}
     */
    public static String canonical(Object token) {
        if (token == null) {
            return NONE;
        }
        String t = String.valueOf(token).toLowerCase();
        if (t.isEmpty() || NONE.equals(t) || "off".equals(t)) {
            return NONE;
        }
        if (CHAINS.equals(t) || "chain".equals(t) || "chainaggr".equals(t)
                || "chain_aggregation".equals(t)) {
            return CHAINS;
        }
        if (LC.equals(t) || "loadconceal".equals(t) || "load_concealment".equals(t)
                || "thinning".equals(t)) {
            return LC;
        }
        throw new RuntimeException("options.config.transform='" + token + "' is not a known model "
                + "transformation. Use 'none', 'chains' or 'lc'.");
    }

    /** The strategy implementing a canonical method name, or null for {@link #NONE}. */
    public static TransformSolve.Strategy strategy(String canonical) {
        if (NONE.equals(canonical)) {
            return null;
        }
        if (CHAINS.equals(canonical)) {
            return new ChainsStrategy();
        }
        if (LC.equals(canonical)) {
            return new LcStrategy();
        }
        throw new RuntimeException("no strategy for transform '" + canonical + "'.");
    }

    /** The nesting depth carried in {@code options.config.transform_depth}. */
    public static int depth(SolverOptions options) {
        Object v = options.config.get("transform_depth");
        if (v instanceof Number) {
            return ((Number) v).intValue();
        }
        return 0;
    }
}
