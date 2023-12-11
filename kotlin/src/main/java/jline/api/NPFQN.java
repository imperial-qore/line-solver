package jline.api;

import jline.util.Matrix;

import java.util.HashMap;
import java.util.Map;

import static jline.lib.M3A.*;

/**
 * APIs for approximating non-product-form queueing networks
 */
public class NPFQN {
    /**
     * Merges MMAP traffic flows
     */
    // TODO: different style used here to specify input configuration
    public static Map<Integer, Matrix> npfqn_traffic_merge(Map<Integer, Map<Integer, Matrix>> MMAPa, String config_merge_, String config_compress_) {
        String merge = "default";
        String compress = "default";
        if (config_merge_.equals("super") || config_merge_.equals("mixture") || config_merge_.equals("interpos")) {
            merge = config_merge_;
        } else {
            throw new RuntimeException("Unsupported configuration for merge.");
        }

        int n = MMAPa.size();
        Map<Integer, Matrix> SMMAP = new HashMap<>();
        //TODO: incomplete implementation
        return mmap_normalize(SMMAP);
    }

    /**
     * Merges MMAP traffic flows with class switching with default parameters
     */
    public static Map<Integer, Matrix> npfqn_traffic_merge_cs(Map<Integer, Map<Integer, Matrix>> MMAPs, Matrix prob) {
        return npfqn_traffic_merge_cs(MMAPs, prob, "default");
    }

    /**
     * Merges MMAP traffic flows with class switching
     */
    public static Map<Integer, Matrix> npfqn_traffic_merge_cs(Map<Integer, Map<Integer, Matrix>> MMAPs, Matrix prob, String config) {
        int n = MMAPs.size();
        int R = prob.numCols;
        Map<Integer, Map<Integer, Matrix>> MMAP_copy = new HashMap<>();
        for (Map.Entry<Integer, Map<Integer, Matrix>> entry : MMAPs.entrySet()) {
            Integer key = entry.getKey();
            Map<Integer, Matrix> innerMap = entry.getValue();

            Map<Integer, Matrix> deepCopyInnerMap = new HashMap<>();

            for (Map.Entry<Integer, Matrix> innerEntry : innerMap.entrySet()) {
                Integer innerKey = innerEntry.getKey();
                Matrix innerValue = innerEntry.getValue();

                // Perform a deep copy of innerValue if needed

                deepCopyInnerMap.put(innerKey, innerValue);
            }

            MMAP_copy.put(key, deepCopyInnerMap);
        }
        for (int i = 0; i < n; i++) {
            Matrix P = new Matrix(R, R, R * R);
            for (int r = 0; r < R; r++) {
                for (int s = 0; s < R; s++) {
                    P.set(r, s, prob.get((i - 1) * R + r, s));
                }
            }
            MMAP_copy.put(i, mmap_mark(MMAP_copy.get(i), P));
        }
        Map<Integer, Matrix> SMMAP = new HashMap<>();
        if (n == 1) {
            SMMAP = MMAPs.get(0);
        } else {
            if (config.equals("default") || config.equals("super")) {
                SMMAP = MMAPs.get(0);
                for (int j = 1; j < n; j++) {
                    SMMAP = mmap_super(SMMAP, MMAP_copy.get(j), "match");
                }
            }
        }
        return SMMAP;
    }

    /**
     * Splits MMAP traffic flows with class switching
     * @param MMAP - a MMAP to be split
     * @param P - a class switching matrix
     * @return a split MMAP
     */
    public static Map<Integer, Map<Integer, Matrix>> npfqn_traffic_split_cs(Map<Integer, Matrix> MMAP, Matrix P) {
        int n = MMAP.size();
        int R = P.numRows;
        int J = P.numCols;
        int M = Math.round(J / (float) R);
        Map<Integer, Map<Integer, Matrix>> SMMAP = new HashMap<>();
        for (int j = 0; j < M; j++) {
            SMMAP.put(j, new HashMap<>());
            SMMAP.get(j).put(0, MMAP.get(0).add(1, MMAP.get(1)));
            SMMAP.get(j).put(1, new Matrix(MMAP.get(1).numRows, MMAP.get(1).numCols));
            for (int s = 0; s < R; s++) {
                SMMAP.get(j).put(2 + s, new Matrix(SMMAP.get(j).get(0).numRows, SMMAP.get(j).get(0).numCols));
                for (int r = 0; r < R; r++) {
                    Matrix a = MMAP.get(2 + r).clone();
                    a.scale(P.get(r, (j - 1) * R + s));
                    SMMAP.get(j).put(2 + s, SMMAP.get(j).get(2 + s).add(1, a));
                    SMMAP.get(j).put(1, SMMAP.get(j).get(1).add(1, a));
                    SMMAP.get(j).put(0, SMMAP.get(j).get(0).add(1, a));
                }
            }
            SMMAP.put(j, mmap_normalize(SMMAP.get(j)));
        }
        return SMMAP;
    }

}
