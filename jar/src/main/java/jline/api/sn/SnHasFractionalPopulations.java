package jline.api.sn;

import jline.lang.NetworkStruct;

public final class SnHasFractionalPopulations {
    private SnHasFractionalPopulations() {}

    /**
     * Checks if the network has closed classes with non-integer populations
     *
     * @param sn NetworkStruct object for the queueing network model
     * @return boolean
     */
    public static boolean snHasFractionalPopulations(NetworkStruct sn) {
        return !sn.njobs.isInteger();
    }
}
