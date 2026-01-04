/**
 * @file Fork Fanout Modification for NetworkStruct
 *
 * Provides functions to directly modify fork fanout (tasksPerLink)
 * in a NetworkStruct without rebuilding the full Network object.
 *
 * Mirrors MATLAB implementation patterns for fork parameter modification.
 *
 * @since LINE 3.0
 */
package jline.api.sn

import jline.lang.NetworkStruct
import jline.lang.constant.NodeType
import jline.lang.nodeparam.ForkNodeParam

/**
 * Sets the fork fanout (number of tasks per output link) for a Fork node.
 *
 * Updates the `fanOut` field in the ForkNodeParam stored in `nodeparam`.
 *
 * @param sn NetworkStruct to modify
 * @param forkNodeIdx Node index of the Fork node (0-based)
 * @param fanOut Number of tasks to create per output link (must be >= 1)
 * @param mode IN_PLACE modifies sn directly, COPY returns a modified copy
 * @param validation Validation level to apply
 * @return Modified NetworkStruct
 * @throws SnValidationException if validation fails
 */
fun snSetForkFanout(
    sn: NetworkStruct,
    forkNodeIdx: Int,
    fanOut: Int,
    mode: ModifyMode = ModifyMode.IN_PLACE,
    validation: ValidationLevel = ValidationLevel.MINIMAL
): NetworkStruct {
    // Validation
    if (validation != ValidationLevel.NONE) {
        val errors = mutableListOf<String>()

        snValidateNodeIndex(sn, forkNodeIdx)?.let { errors.add(it) }

        if (validation == ValidationLevel.FULL) {
            snValidateNodeType(sn, forkNodeIdx, NodeType.Fork)?.let { errors.add(it) }

            if (fanOut < 1) {
                errors.add("fanOut=$fanOut must be >= 1")
            }
        } else if (validation == ValidationLevel.MINIMAL) {
            // At least check it's a Fork node
            if (forkNodeIdx >= 0 && forkNodeIdx < sn.nodetype.size) {
                val nodeType = sn.nodetype[forkNodeIdx]
                if (nodeType != NodeType.Fork) {
                    errors.add("Node $forkNodeIdx is $nodeType, expected Fork")
                }
            }
        }

        if (errors.isNotEmpty()) {
            throw SnValidationException("snSetForkFanout validation failed", errors)
        }
    }

    // Get working copy if needed
    val snWork = if (mode == ModifyMode.COPY) sn.copy<NetworkStruct>() else sn

    // Get the node object
    val node = snWork.nodes.getOrNull(forkNodeIdx) ?: return snWork

    // Update or create ForkNodeParam
    var nodeParam = snWork.nodeparam?.get(node)
    if (nodeParam == null || nodeParam !is ForkNodeParam) {
        nodeParam = ForkNodeParam()
        if (snWork.nodeparam == null) {
            snWork.nodeparam = HashMap()
        }
        snWork.nodeparam[node] = nodeParam
    }

    // Set fanout
    (nodeParam as ForkNodeParam).fanOut = fanOut.toDouble()

    return snWork
}

/**
 * Sets fork fanout for multiple Fork nodes in a single operation.
 *
 * More efficient than calling snSetForkFanout multiple times.
 *
 * @param sn NetworkStruct to modify
 * @param fanOuts Map from fork node index to fanout value
 * @param mode IN_PLACE modifies sn directly, COPY returns a modified copy
 * @param validation Validation level to apply
 * @return Modified NetworkStruct
 * @throws SnValidationException if validation fails
 */
fun snSetForkFanoutBatch(
    sn: NetworkStruct,
    fanOuts: Map<Int, Int>,
    mode: ModifyMode = ModifyMode.IN_PLACE,
    validation: ValidationLevel = ValidationLevel.MINIMAL
): NetworkStruct {
    // Validation
    if (validation != ValidationLevel.NONE) {
        val errors = mutableListOf<String>()

        for ((nodeIdx, fanOut) in fanOuts) {
            snValidateNodeIndex(sn, nodeIdx)?.let { errors.add(it) }

            if (nodeIdx >= 0 && nodeIdx < sn.nodetype.size) {
                val nodeType = sn.nodetype[nodeIdx]
                if (nodeType != NodeType.Fork) {
                    errors.add("Node $nodeIdx is $nodeType, expected Fork")
                }
            }

            if (validation == ValidationLevel.FULL && fanOut < 1) {
                errors.add("fanOut for node $nodeIdx = $fanOut must be >= 1")
            }
        }

        if (errors.isNotEmpty()) {
            throw SnValidationException("snSetForkFanoutBatch validation failed", errors)
        }
    }

    // Get working copy if needed
    val snWork = if (mode == ModifyMode.COPY) sn.copy<NetworkStruct>() else sn

    // Initialize nodeparam map if needed
    if (snWork.nodeparam == null) {
        snWork.nodeparam = HashMap()
    }

    // Update each fork node
    for ((nodeIdx, fanOut) in fanOuts) {
        val node = snWork.nodes.getOrNull(nodeIdx) ?: continue

        var nodeParam = snWork.nodeparam?.get(node)
        if (nodeParam == null || nodeParam !is ForkNodeParam) {
            nodeParam = ForkNodeParam()
            snWork.nodeparam[node] = nodeParam
        }

        (nodeParam as ForkNodeParam).fanOut = fanOut.toDouble()
    }

    return snWork
}

/**
 * Gets the current fork fanout for a Fork node.
 *
 * @param sn NetworkStruct to query
 * @param forkNodeIdx Node index of the Fork node (0-based)
 * @return Current fanout value, or NaN if not set
 */
fun snGetForkFanout(sn: NetworkStruct, forkNodeIdx: Int): Double {
    if (forkNodeIdx < 0 || forkNodeIdx >= sn.nnodes) {
        return Double.NaN
    }

    val node = sn.nodes.getOrNull(forkNodeIdx) ?: return Double.NaN
    val nodeParam = sn.nodeparam?.get(node) ?: return Double.NaN

    return if (nodeParam is ForkNodeParam) {
        nodeParam.fanOut
    } else {
        Double.NaN
    }
}

/**
 * Stochastic network SetForkFanout algorithms
 */
@Suppress("unused")
class SnsetforkfanoutAlgo {
    companion object {
        // Class documentation marker for Dokka
    }
}
