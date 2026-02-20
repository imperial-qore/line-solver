/**
 * @file NPFQN Traffic Merging
 * 
 * Implements traffic merging algorithms for non-product-form queueing networks.
 * Provides methods for combining multiple arrival streams in NPFQN analysis
 * while preserving important distributional properties.
 * 
 * @since LINE 3.0
 */
package jline.api.npfqn

import jline.api.mam.*
import jline.io.line_warning
import jline.util.matrix.MatrixCell

/**
 * NPFQN traffic merging algorithms.
 * 
 * Provides methods for merging MMAP traffic flows in non-product-form queueing networks
 * using various aggregation strategies. Supports multiple merge configurations to handle
 * different network topologies and performance requirements.
 *
 * @since LINE 3.0
 */

/**
 * Merges MMAP traffic flows with specified configurations.
 *
 * @param MMAPa            Map of MMAP traffic flows to be merged
 * @param config_merge_    Configuration for merging ("default", "super", "mixture", "interpos")
 * @param config_compress_ Configuration for compressing ("none", "default")
 * @return Merged and normalized MMAP traffic flow
 * @throws RuntimeException If unsupported configuration for merge is provided
 */
fun npfqn_traffic_merge(MMAPa: Map<Int?, MatrixCell?>, config_merge_: String, config_compress_: String?): MatrixCell? {
    // Filter out empty MMAPs
    val MMAP = MMAPa.values.filterNotNull().filter { !it.isEmpty() }
    val n = MMAP.size
    
    if (n == 0) {
        return mmap_normalize(MatrixCell())
    }
    
    if (n == 1) {
        return mmap_normalize(MMAP[0])
    }
    
    // Set merge configuration
    var merge = "default"
    if (config_merge_ == "super" || config_merge_ == "mixture" || config_merge_ == "interpos") {
        merge = config_merge_
    } else if (config_merge_ != "default") {
        throw RuntimeException("Unsupported configuration for merge: $config_merge_")
    }
    
    // Set compression configuration
    var compress = "default"
    if (config_compress_ == "none") {
        compress = "none"
    }
    
    // Perform merging based on configuration
    var SMMAP: MatrixCell? = when (merge) {
        "default", "super" -> {
            var result = MMAP[0]
            for (j in 1 until n) {
                result = mmap_super(result, MMAP[j], "match")!!
            }
            result
        }
        
        "mixture" -> {
            var result = MMAP[0]
            for (j in 1 until n) {
                result = mmap_super(result, MMAP[j], "match")!!
                if (mmap_isfeasible(result)) {
                    result = mmap_mixture_fit_mmap(result).MMAP
                }
            }
            result
        }
        
        "interpos" -> {
            // Note: This requires m3pp2m functions that are not yet ported
            // For now, fall back to super merge with a warning
            line_warning("npfqn_traffic_merge", "'interpos' merge mode requires m3pp2m functions that are not yet ported. Using 'super' merge instead.")
            var result = MMAP[0]
            for (j in 1 until n) {
                result = mmap_super(result, MMAP[j], "match")!!
            }
            result
        }
        
        else -> throw RuntimeException("Unsupported configuration for merge: $merge")
    }
    
    // Apply compression if requested
    if (compress == "default") {
        // SMMAP = mmap_compress(SMMAP!!)
    }
    
    return mmap_normalize(SMMAP!!)
}