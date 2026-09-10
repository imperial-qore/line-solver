package jline.api.sn;

import jline.lang.NetworkStruct;
import jline.lang.constant.NodeType;
import jline.lang.constant.RoutingStrategy;

public final class SnPrintRoutingMatrix {
    private SnPrintRoutingMatrix() {}

    /**
     * Prints the routing matrix of the network for all classes.
     *
     * @param sn the NetworkStruct object for the queueing network model
     */
    public static void snPrintRoutingMatrix(NetworkStruct sn) {
        snPrintRoutingMatrix(sn, null);
    }

    /**
     * Prints the routing matrix of the network, optionally for a specific job class.
     *
     * @param sn            the NetworkStruct object for the queueing network model
     * @param onlyClassName the name of the specific job class to filter the routing matrix output, or null for all classes
     */
    public static void snPrintRoutingMatrix(NetworkStruct sn, String onlyClassName) {
        // Node and class details
        java.util.List<String> nodeNames = sn.nodenames;
        java.util.List<String> classNames = sn.classnames;
        jline.util.matrix.Matrix rtNodes = sn.rtnodes;
        int nNodes = sn.nnodes;
        int nClasses = sn.nclasses;

        for (int i = 0; i < nNodes; i++) {
            for (int r = 0; r < nClasses; r++) {
                for (int j = 0; j < nNodes; j++) {
                    for (int s = 0; s < nClasses; s++) {
                        if (rtNodes.get(i * nClasses + r, j * nClasses + s) > 0) {
                            String pr;
                            if (sn.nodetype.get(i) == NodeType.Cache) {
                                pr = "state-dependent";
                            } else if (sn.nodetype.get(i) == NodeType.Sink) {
                                continue;
                            } else {
                                if (sn.routing.get(sn.nodes.get(i)).get(sn.jobclasses.get(r)) == RoutingStrategy.DISABLED) {
                                    continue;
                                } else {
                                    pr = String.format("%f", rtNodes.get(i * nClasses + r, j * nClasses + s));
                                }
                            }

                            if (onlyClassName == null) {
                                System.out.printf("\n%s [%s] => %s [%s] : Pr=%s",
                                        nodeNames.get(i),
                                        classNames.get(r),
                                        nodeNames.get(j),
                                        classNames.get(s),
                                        pr);
                            } else {
                                if (classNames.get(r).equalsIgnoreCase(onlyClassName)
                                        || classNames.get(s).equalsIgnoreCase(onlyClassName)) {
                                    System.out.printf("\n%s [%s] => %s [%s] : Pr=%s",
                                            nodeNames.get(i),
                                            classNames.get(r),
                                            nodeNames.get(j),
                                            classNames.get(s),
                                            pr);
                                }
                            }
                        }
                    }
                }
            }
        }
        System.out.print("\n");
    }
}
