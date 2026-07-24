package jline.lang.nodeparam;

import jline.lang.NodeParam;
import jline.lang.constant.ReplacementStrategy;
import jline.util.matrix.Matrix;

import java.util.List;
import java.util.Map;
import java.util.Set;

/**
 * Parameter container for cache nodes in queueing networks.
 * 
 * <p>This class encapsulates all the parameters needed to configure a cache node,
 * including cache capacity, item costs, hit/miss class mappings, replacement strategies,
 * and access probabilities. Cache nodes model caching systems where items can be
 * stored temporarily to reduce access latencies.</p>
 * 
 * <p>Key cache characteristics managed by this class:
 * <ul>
 *   <li>Item capacity constraints and access costs</li>
 *   <li>Hit and miss class routing for cache hits/misses</li>
 *   <li>Replacement strategies (LRU, FIFO, etc.) when cache is full</li>
 *   <li>Read access probabilities for different items</li>
 *   <li>Actual hit/miss probabilities for performance analysis</li>
 * </ul>
 * </p>
 * 
 * @see jline.lang.nodes.Cache
 * @see ReplacementStrategy
 * @since 1.0
 */
public class CacheNodeParam extends NodeParam {
    /** Access cost matrix for cache items by class [items x classes x servers] */
    public Matrix[][] accost;
    
    /** Job class routing matrix for cache hits [items x classes] */
    public Matrix hitclass;
    
    /** Capacity matrix specifying maximum number of each item type [items x 1] */
    public Matrix itemcap;
    
    /** Job class routing matrix for cache misses [items x classes] */
    public Matrix missclass;

    /** Matrix containing the retrieval class for each item [items x classes] */
    public Matrix retrievalClasses;

    /** Set of indices for retrieval classes */
    public Set<Integer> retrievalClassIndices;

    /** Total number of distinct item types in the cache */
    public int nitems;

    /** Total cache capacity */
    public int totalCacheCapacity;

    /** Maximum number of items that can be retrieved simultaneously, 0 if no retrieval system, nItems - cache capacity
     * if a retrieval system exists
     */
    public int retrievalSystemCapacity;

    /** Retrieval System queue indices for each job class [jobClass -> node indices] */
    public Map<Integer, List<Integer>> retrievalSystemQueueIndices;

    /** Read access probabilities for each item by server [server -> list of probabilities by item] */
    public Map<Integer, List<Double>> pread;
    
    /** Replacement strategy used when cache is full (LRU, FIFO, etc.) */
    public ReplacementStrategy replacestrat;

    /** q-LRU admission probability on a miss (1.0 = always admit). Only used for QLRU. */
    public double qlru = 1.0;
    
    /** Actual hit probabilities computed during analysis [items x classes] */
    public Matrix actualhitprob;
    
    /** Actual miss probabilities computed during analysis [items x classes] */
    public Matrix actualmissprob;

    /** Actual delayed-hit fractions computed during analysis (retrieval system) [1 x classes] */
    public Matrix actualdelayedhitprob;

    /** Actual per-list (per-level) hit fractions computed during analysis [classes x lists] */
    public Matrix actualhitproblist;

    /** Actual expected latency computed during analysis [items x classes] */
    public Matrix actualresidt;

    /**
     * Checks if this cache parameter container is empty (no parameters are set).
     *
     * @return true if all cache parameters are null or default values, false otherwise
     */
    @Override
    public boolean isEmpty() {
        return accost == null &&
                hitclass == null &&
                itemcap == null &&
                missclass == null &&
                retrievalClasses == null &&
                (retrievalClassIndices == null || retrievalClassIndices.isEmpty()) &&
                nitems == 0 &&
                totalCacheCapacity == 0 &&
                (retrievalSystemQueueIndices == null || retrievalSystemQueueIndices.isEmpty()) &&
                (pread == null || pread.isEmpty()) &&
                replacestrat == null &&
                actualhitprob == null &&
                actualmissprob == null &&
                actualdelayedhitprob == null &&
                actualhitproblist == null &&
                actualresidt == null;
    }
}
