package jline.api.mdd;

/**
 * Next-state function of a structured model, over local-index tuples.
 *
 * <p>Given a global state as a K-tuple of 0-based local values, returns the
 * rows of its successor states. Used by {@link Mdd_reachset} to generate the
 * reachable set.</p>
 */
public interface MddNextState {

    /**
     * @param state a K-tuple of 0-based local values
     * @return an (m x K) array whose rows are the successors of state
     */
    int[][] next(int[] state);
}
