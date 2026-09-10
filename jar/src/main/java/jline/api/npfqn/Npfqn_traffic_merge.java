/**
 * @file NPFQN Traffic Merging
 *
 * Implements traffic merging algorithms for non-product-form queueing networks.
 * Provides methods for combining multiple arrival streams in NPFQN analysis
 * while preserving important distributional properties.
 *
 * @since LINE 3.0
 */
package jline.api.npfqn;

import java.util.ArrayList;
import java.util.List;
import java.util.Map;

import jline.api.mam.Mmap_isfeasible;
import jline.api.mam.Mmap_mixture_fit_mmap;
import jline.api.mam.Mmap_normalize;
import jline.api.mam.Mmap_super;
import jline.io.InputOutput;
import jline.util.matrix.MatrixCell;

public final class Npfqn_traffic_merge {
    private Npfqn_traffic_merge() {}

    /**
     * Merges MMAP traffic flows with specified configurations.
     *
     * @param MMAPa            Map of MMAP traffic flows to be merged
     * @param config_merge_    Configuration for merging ("default", "super", "mixture", "interpos")
     * @param config_compress_ Configuration for compressing ("none", "default")
     * @return Merged and normalized MMAP traffic flow
     * @throws RuntimeException If unsupported configuration for merge is provided
     */
    public static MatrixCell npfqn_traffic_merge(Map<Integer, MatrixCell> MMAPa, String config_merge_, String config_compress_) {
        // Filter out empty MMAPs
        List<MatrixCell> MMAP = new ArrayList<MatrixCell>();
        for (MatrixCell c : MMAPa.values()) {
            if (c != null && !c.isEmpty()) {
                MMAP.add(c);
            }
        }
        int n = MMAP.size();

        if (n == 0) {
            return Mmap_normalize.mmap_normalize(new MatrixCell());
        }

        if (n == 1) {
            return Mmap_normalize.mmap_normalize(MMAP.get(0));
        }

        // Set merge configuration
        String merge = "default";
        if ("super".equals(config_merge_) || "mixture".equals(config_merge_) || "interpos".equals(config_merge_)) {
            merge = config_merge_;
        } else if (!"default".equals(config_merge_)) {
            throw new RuntimeException("Unsupported configuration for merge: " + config_merge_);
        }

        // Set compression configuration
        String compress = "default";
        if ("none".equals(config_compress_)) {
            compress = "none";
        }

        // Perform merging based on configuration
        MatrixCell SMMAP;
        if ("default".equals(merge) || "super".equals(merge)) {
            MatrixCell result = MMAP.get(0);
            for (int j = 1; j < n; j++) {
                result = Mmap_super.mmap_super(result, MMAP.get(j), "match");
            }
            SMMAP = result;
        } else if ("mixture".equals(merge)) {
            MatrixCell result = MMAP.get(0);
            for (int j = 1; j < n; j++) {
                result = Mmap_super.mmap_super(result, MMAP.get(j), "match");
                if (Mmap_isfeasible.mmap_isfeasible(result)) {
                    result = Mmap_mixture_fit_mmap.mmap_mixture_fit_mmap(result).MMAP;
                }
            }
            SMMAP = result;
        } else if ("interpos".equals(merge)) {
            // Note: This requires m3pp2m functions that are not yet ported
            // For now, fall back to super merge with a warning
            InputOutput.line_warning("npfqn_traffic_merge", "'interpos' merge mode requires m3pp2m functions that are not yet ported. Using 'super' merge instead.");
            MatrixCell result = MMAP.get(0);
            for (int j = 1; j < n; j++) {
                result = Mmap_super.mmap_super(result, MMAP.get(j), "match");
            }
            SMMAP = result;
        } else {
            throw new RuntimeException("Unsupported configuration for merge: " + merge);
        }

        // Apply compression if requested
        if ("default".equals(compress)) {
            // SMMAP = mmap_compress(SMMAP);
        }

        return Mmap_normalize.mmap_normalize(SMMAP);
    }
}
