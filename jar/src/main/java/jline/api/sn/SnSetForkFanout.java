/**
 * @file Fork Fanout Modification for NetworkStruct
 *
 * @since LINE 3.0
 */
package jline.api.sn;

import java.util.ArrayList;
import java.util.HashMap;
import java.util.List;
import java.util.Map;

import jline.lang.NetworkStruct;
import jline.lang.nodes.Node;
import jline.lang.constant.NodeType;
import jline.lang.nodeparam.ForkNodeParam;
import jline.lang.NodeParam;

public final class SnSetForkFanout {
    private SnSetForkFanout() {}

    /**
     * Sets the fork fanout (number of tasks per output link) for a Fork node.
     */
    public static NetworkStruct snSetForkFanout(NetworkStruct sn, int forkNodeIdx, int fanOut,
                                                ModifyMode mode, ValidationLevel validation) {
        if (validation != ValidationLevel.NONE) {
            List<String> errors = new ArrayList<String>();

            String err = SnValidate.snValidateNodeIndex(sn, forkNodeIdx);
            if (err != null) {
                errors.add(err);
            }

            if (validation == ValidationLevel.FULL) {
                String err2 = SnValidate.snValidateNodeType(sn, forkNodeIdx, NodeType.Fork);
                if (err2 != null) {
                    errors.add(err2);
                }

                if (fanOut < 1) {
                    errors.add("fanOut=" + fanOut + " must be >= 1");
                }
            } else if (validation == ValidationLevel.MINIMAL) {
                if (forkNodeIdx >= 0 && forkNodeIdx < sn.nodetype.size()) {
                    NodeType nodeType = sn.nodetype.get(forkNodeIdx);
                    if (nodeType != NodeType.Fork) {
                        errors.add("Node " + forkNodeIdx + " is " + nodeType + ", expected Fork");
                    }
                }
            }

            if (!errors.isEmpty()) {
                throw new SnValidationException("snSetForkFanout validation failed", errors);
            }
        }

        NetworkStruct snWork = (mode == ModifyMode.COPY) ? sn.<NetworkStruct>copy() : sn;

        if (snWork.nodes == null || forkNodeIdx < 0 || forkNodeIdx >= snWork.nodes.size()) {
            return snWork;
        }
        Node node = snWork.nodes.get(forkNodeIdx);
        if (node == null) {
            return snWork;
        }

        NodeParam nodeParam = (snWork.nodeparam != null) ? snWork.nodeparam.get(node) : null;
        if (!(nodeParam instanceof ForkNodeParam)) {
            nodeParam = new ForkNodeParam();
            if (snWork.nodeparam == null) {
                snWork.nodeparam = new HashMap<Node, NodeParam>();
            }
            snWork.nodeparam.put(node, nodeParam);
        }

        ((ForkNodeParam) nodeParam).fanOut = (double) fanOut;

        return snWork;
    }

    public static NetworkStruct snSetForkFanout(NetworkStruct sn, int forkNodeIdx, int fanOut, ModifyMode mode) {
        return snSetForkFanout(sn, forkNodeIdx, fanOut, mode, ValidationLevel.MINIMAL);
    }

    public static NetworkStruct snSetForkFanout(NetworkStruct sn, int forkNodeIdx, int fanOut) {
        return snSetForkFanout(sn, forkNodeIdx, fanOut, ModifyMode.IN_PLACE, ValidationLevel.MINIMAL);
    }

    /**
     * Sets fork fanout for multiple Fork nodes in a single operation.
     */
    public static NetworkStruct snSetForkFanoutBatch(NetworkStruct sn, Map<Integer, Integer> fanOuts,
                                                     ModifyMode mode, ValidationLevel validation) {
        if (validation != ValidationLevel.NONE) {
            List<String> errors = new ArrayList<String>();

            for (Map.Entry<Integer, Integer> entry : fanOuts.entrySet()) {
                int nodeIdx = entry.getKey();
                int fanOut = entry.getValue();

                String err = SnValidate.snValidateNodeIndex(sn, nodeIdx);
                if (err != null) {
                    errors.add(err);
                }

                if (nodeIdx >= 0 && nodeIdx < sn.nodetype.size()) {
                    NodeType nodeType = sn.nodetype.get(nodeIdx);
                    if (nodeType != NodeType.Fork) {
                        errors.add("Node " + nodeIdx + " is " + nodeType + ", expected Fork");
                    }
                }

                if (validation == ValidationLevel.FULL && fanOut < 1) {
                    errors.add("fanOut for node " + nodeIdx + " = " + fanOut + " must be >= 1");
                }
            }

            if (!errors.isEmpty()) {
                throw new SnValidationException("snSetForkFanoutBatch validation failed", errors);
            }
        }

        NetworkStruct snWork = (mode == ModifyMode.COPY) ? sn.<NetworkStruct>copy() : sn;

        if (snWork.nodeparam == null) {
            snWork.nodeparam = new HashMap<Node, NodeParam>();
        }

        for (Map.Entry<Integer, Integer> entry : fanOuts.entrySet()) {
            int nodeIdx = entry.getKey();
            int fanOut = entry.getValue();

            if (snWork.nodes == null || nodeIdx < 0 || nodeIdx >= snWork.nodes.size()) {
                continue;
            }
            Node node = snWork.nodes.get(nodeIdx);
            if (node == null) {
                continue;
            }

            NodeParam nodeParam = snWork.nodeparam.get(node);
            if (!(nodeParam instanceof ForkNodeParam)) {
                nodeParam = new ForkNodeParam();
                snWork.nodeparam.put(node, nodeParam);
            }

            ((ForkNodeParam) nodeParam).fanOut = (double) fanOut;
        }

        return snWork;
    }

    public static NetworkStruct snSetForkFanoutBatch(NetworkStruct sn, Map<Integer, Integer> fanOuts, ModifyMode mode) {
        return snSetForkFanoutBatch(sn, fanOuts, mode, ValidationLevel.MINIMAL);
    }

    public static NetworkStruct snSetForkFanoutBatch(NetworkStruct sn, Map<Integer, Integer> fanOuts) {
        return snSetForkFanoutBatch(sn, fanOuts, ModifyMode.IN_PLACE, ValidationLevel.MINIMAL);
    }

    /**
     * Gets the current fork fanout for a Fork node.
     */
    public static double snGetForkFanout(NetworkStruct sn, int forkNodeIdx) {
        if (forkNodeIdx < 0 || forkNodeIdx >= sn.nnodes) {
            return Double.NaN;
        }

        if (sn.nodes == null || forkNodeIdx >= sn.nodes.size()) {
            return Double.NaN;
        }
        Node node = sn.nodes.get(forkNodeIdx);
        if (node == null) {
            return Double.NaN;
        }
        NodeParam nodeParam = (sn.nodeparam != null) ? sn.nodeparam.get(node) : null;
        if (nodeParam == null) {
            return Double.NaN;
        }

        if (nodeParam instanceof ForkNodeParam) {
            return ((ForkNodeParam) nodeParam).fanOut;
        } else {
            return Double.NaN;
        }
    }
}
