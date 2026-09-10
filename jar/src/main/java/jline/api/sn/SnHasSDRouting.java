/**
 * @file Stochastic network state-dependent routing checker
 *
 * Determines if a queueing network has state-dependent routing strategies that
 * violate product-form assumptions. State-dependent routing includes Round-Robin,
 * Weighted Round-Robin, Join Shortest Queue, Power of K Choices, and Reinforcement Learning.
 *
 * @since LINE 3.0
 */
package jline.api.sn;

import jline.lang.NetworkStruct;
import jline.lang.JobClass;
import jline.lang.nodes.Node;
import jline.lang.constant.RoutingStrategy;

import java.util.Map;

public final class SnHasSDRouting {
    private SnHasSDRouting() {}

    /**
     * Checks if the network has state-dependent routing strategies.
     *
     * <p>Product-form requires state-independent (Markovian) routing.
     * PROB and RAND are product-form compatible.
     * RROBIN, WRROBIN, JSQ, SQ are state-dependent and violate product-form.
     *
     * @param sn NetworkStruct object for the queueing network model
     * @return true if the network has state-dependent routing, false otherwise
     */
    public static boolean snHasSDRouting(NetworkStruct sn) {
        if (sn.routing == null || sn.routing.isEmpty()) {
            return false;
        }

        for (Map.Entry<Node, Map<JobClass, RoutingStrategy>> nodeEntry : sn.routing.entrySet()) {
            Map<JobClass, RoutingStrategy> classMap = nodeEntry.getValue();
            if (classMap != null) {
                for (Map.Entry<JobClass, RoutingStrategy> classEntry : classMap.entrySet()) {
                    RoutingStrategy strategy = classEntry.getValue();
                    if (strategy == RoutingStrategy.RROBIN ||
                            strategy == RoutingStrategy.WRROBIN ||
                            strategy == RoutingStrategy.JSQ ||
                            strategy == RoutingStrategy.SQ) {
                        return true;
                    }
                }
            }
        }
        return false;
    }
}
