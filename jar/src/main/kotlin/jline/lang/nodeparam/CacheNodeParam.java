package jline.lang.nodeparam;

import jline.lang.NodeParam;
import jline.lang.constant.ReplacementStrategy;
import jline.util.matrix.Matrix;

import java.util.List;
import java.util.Map;

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
    
    /** Total number of distinct item types in the cache */
    public int nitems;
    
    /** Read access probabilities for each item by server [server -> list of probabilities by item] */
    public Map<Integer, List<Double>> pread;
    
    /** Replacement strategy used when cache is full (LRU, FIFO, etc.) */
    public ReplacementStrategy replacestrat;
    
    /** Actual hit probabilities computed during analysis [items x classes] */
    public Matrix actualhitprob;
    
    /** Actual miss probabilities computed during analysis [items x classes] */
    public Matrix actualmissprob;

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
                nitems == 0 &&
                (pread == null || pread.isEmpty()) &&
                replacestrat == null &&
                actualhitprob == null &&
                actualmissprob == null;
    }
}
