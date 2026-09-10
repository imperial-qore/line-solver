package jline.api.sn;

import jline.lang.NetworkStruct;

public final class SnHasOpenClasses {
    private SnHasOpenClasses() {}

    /**
     * Checks if the network has one or more open classes
     *
     * @param sn NetworkStruct object for the queueing network model
     * @return boolean
     */
    public static boolean snHasOpenClasses(NetworkStruct sn) {
        return sn.njobs.hasInfinite();
    }
}
