package jline.api.sn;

import jline.lang.NetworkStruct;

public final class SnHasLoadDependence {
    private SnHasLoadDependence() {}

    /**
     * Checks if the network has a station with load-dependent service process
     *
     * @param sn NetworkStruct object for the queueing network model
     * @return boolean
     */
    public static boolean snHasLoadDependence(NetworkStruct sn) {
        return sn.lldscaling.getNumCols() > 0;
    }
}
