package jline.api.sn;

import jline.lang.NetworkStruct;

public final class SnHasClassSwitching {
    private SnHasClassSwitching() {}

    /**
     * Checks if the network uses class-switching
     *
     * @param sn NetworkStruct object for the queueing network model
     * @return boolean
     */
    public static boolean snHasClassSwitching(NetworkStruct sn) {
        return sn.nclasses != sn.nchains;
    }
}
