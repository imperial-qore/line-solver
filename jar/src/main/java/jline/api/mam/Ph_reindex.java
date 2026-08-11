/**
 * @file Phase-type distribution reindexing utilities
 *
 * Reindexes phase-type distribution maps for network models using integer station and class indices.
 * Essential utility for integrating MAM algorithms with LINE network modeling infrastructure.
 *
 * @since LINE 3.0
 */
package jline.api.mam;

import java.util.HashMap;
import java.util.Map;

import jline.lang.NetworkStruct;
import jline.util.matrix.MatrixCell;

public final class Ph_reindex {
    private Ph_reindex() {}

    /**
     * Reindexes phase-type (PH) distributions for a network model based on station and
     * job class indices.
     *
     * @param sn the NetworkStruct object containing the network structure, indexed PH
     *           distributions (sn.proc), and internal indexing information
     * @return a reindexed map where the key is an integer station index and the value is
     *         another map with integer job class indices and MatrixCell values
     */
    public static Map<Integer, Map<Integer, MatrixCell>> ph_reindex(NetworkStruct sn) {
        Map<Integer, Map<Integer, MatrixCell>> result = new HashMap<Integer, Map<Integer, MatrixCell>>();

        for (int i = 0; i < sn.nstations; i++) {
            for (int j = 0; j < sn.nclasses; j++) {
                if (j == 0) {
                    result.put(i, new HashMap<Integer, MatrixCell>());
                }
                result.get(i).put(j, sn.proc.get(sn.stations.get(i)).get(sn.jobclasses.get(j)));
            }
        }
        return result;
    }
}
