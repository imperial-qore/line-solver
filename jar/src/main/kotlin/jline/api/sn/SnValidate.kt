/**
 * @file NetworkStruct Validation Utilities
 *
 * Provides validation functions for NetworkStruct fields to ensure consistency
 * when modifying parameters directly via SN API methods.
 *
 * Mirrors MATLAB implementation patterns for validation.
 *
 * @since LINE 3.0
 */
package jline.api.sn

import jline.lang.NetworkStruct
import jline.lang.constant.NodeType
import jline.util.matrix.Matrix

/**
 * Validates NetworkStruct consistency at the specified validation level.
 *
 * @param sn Network structure to validate
 * @param level Validation level to apply
 * @return List of validation error messages (empty if valid)
 */
fun snValidate(sn: NetworkStruct, level: ValidationLevel = ValidationLevel.FULL): List<String> {
    if (level == ValidationLevel.NONE) {
        return emptyList()
    }

    val errors = mutableListOf<String>()

    // Minimal validation - always run
    errors.addAll(snValidateDimensions(sn))

    if (level == ValidationLevel.FULL) {
        errors.addAll(snValidateRates(sn))
        errors.addAll(snValidatePopulation(sn))
        errors.addAll(snValidateRouting(sn))
        errors.addAll(snValidateServers(sn))
    }

    return errors
}

/**
 * Validates that all matrix dimensions are consistent with nstations and nclasses.
 *
 * @param sn Network structure to validate
 * @return List of dimension-related validation errors
 */
fun snValidateDimensions(sn: NetworkStruct): List<String> {
    val errors = mutableListOf<String>()

    val M = sn.nstations
    val K = sn.nclasses

    // Check rates matrix
    if (sn.rates != null) {
        if (sn.rates.numRows != M || sn.rates.numCols != K) {
            errors.add("rates matrix dimensions (${sn.rates.numRows}x${sn.rates.numCols}) do not match (nstations=$M x nclasses=$K)")
        }
    }

    // Check scv matrix
    if (sn.scv != null) {
        if (sn.scv.numRows != M || sn.scv.numCols != K) {
            errors.add("scv matrix dimensions (${sn.scv.numRows}x${sn.scv.numCols}) do not match (nstations=$M x nclasses=$K)")
        }
    }

    // Check nservers matrix
    if (sn.nservers != null) {
        if (sn.nservers.numRows != M || sn.nservers.numCols != 1) {
            errors.add("nservers matrix dimensions (${sn.nservers.numRows}x${sn.nservers.numCols}) do not match (nstations=$M x 1)")
        }
    }

    // Check njobs matrix
    if (sn.njobs != null) {
        val totalElements = sn.njobs.numRows * sn.njobs.numCols
        if (totalElements != K) {
            errors.add("njobs matrix total elements ($totalElements) does not match nclasses=$K")
        }
    }

    // Check classprio matrix
    if (sn.classprio != null) {
        val totalElements = sn.classprio.numRows * sn.classprio.numCols
        if (totalElements != K) {
            errors.add("classprio matrix total elements ($totalElements) does not match nclasses=$K")
        }
    }

    // Check phases matrix
    if (sn.phases != null) {
        if (sn.phases.numRows != M || sn.phases.numCols != K) {
            errors.add("phases matrix dimensions (${sn.phases.numRows}x${sn.phases.numCols}) do not match (nstations=$M x nclasses=$K)")
        }
    }

    return errors
}

/**
 * Validates service and arrival rates for NaN/Inf and positive values.
 *
 * @param sn Network structure to validate
 * @return List of rate-related validation errors
 */
fun snValidateRates(sn: NetworkStruct): List<String> {
    val errors = mutableListOf<String>()

    if (sn.rates == null) {
        return errors
    }

    for (i in 0 until sn.rates.numRows) {
        for (j in 0 until sn.rates.numCols) {
            val rate = sn.rates.get(i, j)
            if (rate.isNaN()) {
                // NaN is allowed for undefined service (e.g., class not served at this station)
                continue
            }
            if (rate < 0) {
                errors.add("rates[$i,$j] = $rate is negative")
            }
        }
    }

    if (sn.scv != null) {
        for (i in 0 until sn.scv.numRows) {
            for (j in 0 until sn.scv.numCols) {
                val scvVal = sn.scv.get(i, j)
                if (scvVal.isNaN()) {
                    continue
                }
                if (scvVal < 0) {
                    errors.add("scv[$i,$j] = $scvVal is negative (SCV must be >= 0)")
                }
            }
        }
    }

    return errors
}

/**
 * Validates job population for closed classes.
 *
 * @param sn Network structure to validate
 * @return List of population-related validation errors
 */
fun snValidatePopulation(sn: NetworkStruct): List<String> {
    val errors = mutableListOf<String>()

    if (sn.njobs == null) {
        return errors
    }

    for (k in 0 until sn.nclasses) {
        val njob = getPopulation(sn.njobs, k)
        if (njob.isNaN()) {
            errors.add("njobs for class $k is NaN")
        } else if (njob.isFinite() && njob < 0) {
            errors.add("njobs for class $k = $njob is negative")
        }
        // Inf is valid for open classes
    }

    return errors
}

/**
 * Validates routing matrix is stochastic (row sums = 1).
 *
 * @param sn Network structure to validate
 * @return List of routing-related validation errors
 */
fun snValidateRouting(sn: NetworkStruct): List<String> {
    val errors = mutableListOf<String>()

    if (sn.rt == null) {
        return errors
    }

    val tolerance = 1e-6

    for (i in 0 until sn.rt.numRows) {
        var rowSum = 0.0
        var hasNonZero = false
        for (j in 0 until sn.rt.numCols) {
            val val_ = sn.rt.get(i, j)
            if (val_.isNaN()) {
                errors.add("rt[$i,$j] is NaN")
                continue
            }
            if (val_ < 0) {
                errors.add("rt[$i,$j] = $val_ is negative")
            }
            if (val_ > 0) {
                hasNonZero = true
            }
            rowSum += val_
        }
        // Only check row sum if row has non-zero entries
        if (hasNonZero && kotlin.math.abs(rowSum - 1.0) > tolerance) {
            errors.add("rt row $i sum = $rowSum (expected 1.0)")
        }
    }

    return errors
}

/**
 * Validates server counts are positive.
 *
 * @param sn Network structure to validate
 * @return List of server-related validation errors
 */
fun snValidateServers(sn: NetworkStruct): List<String> {
    val errors = mutableListOf<String>()

    if (sn.nservers == null) {
        return errors
    }

    for (i in 0 until sn.nservers.numRows) {
        val nserv = sn.nservers.get(i, 0)
        if (nserv.isNaN()) {
            errors.add("nservers[$i] is NaN")
        } else if (nserv <= 0 && !nserv.isInfinite()) {
            // Check node type - Source and Sink may have 0 servers
            val nodeType = if (i < sn.nodetype.size) sn.nodetype[i.toInt()] else null
            if (nodeType != NodeType.Source && nodeType != NodeType.Sink) {
                errors.add("nservers[$i] = $nserv must be positive")
            }
        }
    }

    return errors
}

/**
 * Validates a station index is within bounds.
 *
 * @param sn Network structure
 * @param stationIdx Station index to validate
 * @param paramName Parameter name for error message
 * @return Error message or null if valid
 */
fun snValidateStationIndex(sn: NetworkStruct, stationIdx: Int, paramName: String = "stationIdx"): String? {
    if (stationIdx < 0 || stationIdx >= sn.nstations) {
        return "$paramName=$stationIdx is out of bounds [0, ${sn.nstations - 1}]"
    }
    return null
}

/**
 * Validates a class index is within bounds.
 *
 * @param sn Network structure
 * @param classIdx Class index to validate
 * @param paramName Parameter name for error message
 * @return Error message or null if valid
 */
fun snValidateClassIndex(sn: NetworkStruct, classIdx: Int, paramName: String = "classIdx"): String? {
    if (classIdx < 0 || classIdx >= sn.nclasses) {
        return "$paramName=$classIdx is out of bounds [0, ${sn.nclasses - 1}]"
    }
    return null
}

/**
 * Validates a node index is within bounds.
 *
 * @param sn Network structure
 * @param nodeIdx Node index to validate
 * @param paramName Parameter name for error message
 * @return Error message or null if valid
 */
fun snValidateNodeIndex(sn: NetworkStruct, nodeIdx: Int, paramName: String = "nodeIdx"): String? {
    if (nodeIdx < 0 || nodeIdx >= sn.nnodes) {
        return "$paramName=$nodeIdx is out of bounds [0, ${sn.nnodes - 1}]"
    }
    return null
}

/**
 * Validates a node is of the expected type.
 *
 * @param sn Network structure
 * @param nodeIdx Node index to check
 * @param expectedType Expected node type
 * @return Error message or null if valid
 */
fun snValidateNodeType(sn: NetworkStruct, nodeIdx: Int, expectedType: NodeType): String? {
    if (nodeIdx < 0 || nodeIdx >= sn.nodetype.size) {
        return "nodeIdx=$nodeIdx is out of bounds"
    }
    val actualType = sn.nodetype[nodeIdx]
    if (actualType != expectedType) {
        return "Node $nodeIdx is $actualType, expected $expectedType"
    }
    return null
}

/**
 * Helper function to get population from njobs matrix.
 * Handles both (1 x K) and (K x 1) layouts.
 */
private fun getPopulation(njobs: Matrix, classIdx: Int): Double {
    return if (njobs.numRows == 1) {
        njobs.get(0, classIdx)
    } else {
        njobs.get(classIdx, 0)
    }
}

/**
 * Stochastic network Validate algorithms
 */
@Suppress("unused")
class SnvalidateAlgo {
    companion object {
        // Class documentation marker for Dokka
    }
}
