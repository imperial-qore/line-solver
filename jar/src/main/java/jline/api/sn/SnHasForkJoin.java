package jline.api.sn;

import jline.lang.NetworkStruct;

public final class SnHasForkJoin {
    private SnHasForkJoin() {}

    /**
     * Checks if the network uses fork and/or join nodes
     *
     * @param sn - NetworkStruct object for the queueing network model
     * @return boolean
     */
    public static boolean snHasForkJoin(NetworkStruct sn) {
        for (int i = 0; i < sn.fj.getNumRows(); i++) {
            for (int j = 0; j < sn.fj.getNumCols(); j++) {
                if (sn.fj.get(i, j) > 0) return true;
            }
        }
        return false;
    }
}
