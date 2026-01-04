package jline.api.sn

import jline.lang.NetworkStruct
import jline.lang.constant.NodeType
import jline.lang.constant.RoutingStrategy

/**
 * Prints the routing matrix of the network, optionally for a specific job class.
 *
 * @param sn            the NetworkStruct object for the queueing network model
 * @param onlyClassName the name of the specific job class to filter the routing matrix output, or null for all classes
 *
 * Prints the routing matrix of the network for all classes if onlyClassName is null.
 *
 * @param sn the NetworkStruct object for the queueing network model
 */
fun snPrintRoutingMatrix(sn: NetworkStruct, onlyClassName: String? = null) {
    // Node and class details
    val nodeNames = sn.nodenames
    val classNames = sn.classnames
    val rtNodes = sn.rtnodes
    val nNodes = sn.nnodes
    val nClasses = sn.nclasses

    // Iterate through all nodes and classes
    for (i in 0..<nNodes) {
        for (r in 0..<nClasses) {
            for (j in 0..<nNodes) {
                for (s in 0..<nClasses) {
                    if (rtNodes[i * nClasses + r, j * nClasses + s] > 0) {
                        val pr = if (sn.nodetype[i] == NodeType.Cache) {
                            "state-dependent"
                        } else if (sn.nodetype[i] == NodeType.Sink) {
                            continue
                        } else {
                            if (sn.routing[sn.nodes[i]]!![sn.jobclasses[r]] == RoutingStrategy.DISABLED) {
                                continue
                            } else {
                                String.format("%f", rtNodes[i * nClasses + r, j * nClasses + s])
                            }
                        }

                        if (onlyClassName == null) {
                            System.out.printf("\n%s [%s] => %s [%s] : Pr=%s",
                                nodeNames[i],
                                classNames[r],
                                nodeNames[j],
                                classNames[s],
                                pr)
                        } else {
                            if (classNames[r].equals(onlyClassName, ignoreCase = true) || classNames[s].equals(
                                    onlyClassName,
                                    ignoreCase = true)) {
                                System.out.printf("\n%s [%s] => %s [%s] : Pr=%s",
                                    nodeNames[i],
                                    classNames[r],
                                    nodeNames[j],
                                    classNames[s],
                                    pr)
                            }
                        }
                    }
                }
            }
        }
    }
    print("\n")
}
/**
 * Stochastic network PrintRoutingMatrix algorithms
 */
@Suppress("unused")
class SnprintroutingmatrixAlgo {
    companion object {
        // Class documentation marker for Dokka
    }
}