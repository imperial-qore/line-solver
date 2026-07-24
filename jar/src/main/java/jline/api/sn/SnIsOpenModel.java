package jline.api.sn;

import jline.lang.NetworkStruct;

/**
 * Stochastic network model type classifier for open models.
 *
 * <p>Identifies open queueing network models with external arrivals and infinite
 * job populations. Open models are essential for analyzing systems with external
 * traffic sources and unlimited capacity for job creation.
 *
 * @since LINE 3.0
 */
public final class SnIsOpenModel {
    private SnIsOpenModel() {}

    /**
     * Checks if the network is an open model.
     *
     * @param sn the NetworkStruct object for the queueing network model
     * @return true if the network has infinite jobs and no finite jobs, false otherwise
     */
    public static boolean snIsOpenModel(NetworkStruct sn) {
        return !sn.njobs.hasFinite() && sn.njobs.hasInfinite();
    }
}
