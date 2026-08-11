package jline.api.sn;

import jline.lang.NetworkStruct;

public final class SnHasPriorities {
    private SnHasPriorities() {}

    /**
     * Checks if the network uses class priorities
     *
     * @param sn - NetworkStruct object for the queueing network model
     * @return boolean
     */
    public static boolean snHasPriorities(NetworkStruct sn) {
        for (int i = 0; i < sn.classprio.getNumRows(); i++) {
            for (int j = 0; j < sn.classprio.getNumCols(); j++) {
                if (sn.classprio.get(i, j) > 0) return true;
            }
        }
        return false;
    }
}
