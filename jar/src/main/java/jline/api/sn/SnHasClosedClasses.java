package jline.api.sn;

import jline.lang.NetworkStruct;

public final class SnHasClosedClasses {
    private SnHasClosedClasses() {}

    /**
     * Checks if the network has one or more closed classes
     *
     * @param sn NetworkStruct object for the queueing network model
     * @return boolean
     */
    public static boolean snHasClosedClasses(NetworkStruct sn) {
        return sn.njobs.hasFinite();
    }
}
