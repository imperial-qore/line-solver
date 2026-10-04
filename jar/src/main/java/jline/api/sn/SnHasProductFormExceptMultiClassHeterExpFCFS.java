package jline.api.sn;

import jline.lang.NetworkStruct;

/**
 * Legacy name of {@link SnHasProductFormNotHetFCFS}, kept as a delegating
 * alias exactly as the python twin does: the pre-rework body here still
 * admitted LCFSPR and lacked the FCFS equal-service-means check, so it could
 * disagree with the predicate the AMVA dispatch actually uses.
 */
public final class SnHasProductFormExceptMultiClassHeterExpFCFS {
    private SnHasProductFormExceptMultiClassHeterExpFCFS() {}

    /**
     * Checks if the network satisfies product-form assumptions except multiclass heterogeneous FCFS
     *
     * @param sn - NetworkStruct object for the queueing network model
     * @return boolean
     */
    public static boolean snHasProductFormExceptMultiClassHeterExpFCFS(NetworkStruct sn) {
        return SnHasProductFormNotHetFCFS.snHasProductFormNotHetFCFS(sn);
    }
}
