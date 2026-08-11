package jline.api.sn;

import jline.lang.NetworkStruct;

public final class SnHasMixedClasses {
    private SnHasMixedClasses() {}

    /**
     * Checks if the network has both open and closed classes
     *
     * @param sn NetworkStruct object for the queueing network model
     * @return boolean
     */
    public static boolean snHasMixedClasses(NetworkStruct sn) {
        return sn.njobs.hasFinite() && sn.njobs.hasInfinite();
    }
}
