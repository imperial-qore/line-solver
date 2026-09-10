package jline.api.sn;

import static jline.api.sn.SnHasMixedClasses.snHasMixedClasses;

import jline.lang.NetworkStruct;

public final class SnIsMixedModel {
    private SnIsMixedModel() {}

    /**
     * Checks if the network is a mixed model.
     *
     * @param sn the NetworkStruct object for the queueing network model
     * @return true if the network has mixed classes, false otherwise
     */
    public static boolean snIsMixedModel(NetworkStruct sn) {
        return snHasMixedClasses(sn);
    }
}
