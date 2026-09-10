package jline.api.mdd;

/** Result of the Miner-Ciardo-Donatelli level aggregation. */
public class MddMcdResult {

    /** Mean occupancy per station (or place), in station order. */
    public double[] QLen;
    /** Per-station throughput; null when the descriptor carries no queueing parameters. */
    public double[] X;
    /** Per-station utilization; null when the descriptor carries no queueing parameters. */
    public double[] U;
    /** pik[k] is the level-k stationary vector over M_k, in paper orientation. */
    public double[][] pik;
    /** Mrows[k][r] = {node id, local value} of row r of M_k. */
    public int[][][] Mrows;
    /** |M_k| per paper level. */
    public int[] levelSizes;
    /** Fixed-point iterations performed. */
    public int iters;
    /**
     * max |A(p)| per paper level: the largest number of distinct root-to-node
     * paths at that level. 1 means no node there is shared, so conditioning on
     * the node equals conditioning on the whole path above it.
     */
    public double[] pathsPerLevel;
    /**
     * True certifies the result is EXACT with no reference solve needed; false
     * means "not certified by this test", never "approximate" -- a product-form
     * model is exact however much its diagram shares.
     */
    public boolean noAggregation;
}
