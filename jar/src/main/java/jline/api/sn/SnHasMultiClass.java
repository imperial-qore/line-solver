/**
 * @file Stochastic network multi-class detection utility
 *
 * Identifies queueing networks with multiple job classes, which require specialized
 * analysis algorithms that account for class-dependent service parameters, routing
 * probabilities, and scheduling policies in multi-class queueing systems.
 *
 * @since LINE 3.0
 */
package jline.api.sn;

import jline.lang.NetworkStruct;

public final class SnHasMultiClass {
    private SnHasMultiClass() {}

    public static boolean snHasMultiClass(NetworkStruct sn) {
        return sn.nclasses > 1;
    }
}
