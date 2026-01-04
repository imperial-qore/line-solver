package jline.api.sn

import jline.lang.NetworkStruct
import jline.lang.nodes.Node
import jline.lang.nodes.StatefulNode
import jline.lang.nodes.Station
import jline.lang.JobClass
import jline.util.matrix.Matrix
import jline.util.matrix.MatrixCell
import jline.GlobalConstants
import jline.lang.Event
import jline.lang.Sync

/**
 * Prints comprehensive information about a NetworkStruct.
 * This function displays all fields, matrices, lists, and maps in a formatted manner
 * useful for debugging and inspection.
 *
 * @param sn the NetworkStruct object to print
 */
fun snPrint(sn: NetworkStruct) {
    // Basic integer fields
    println("nstations: ${sn.nstations}")
    println("nstateful: ${sn.nstateful}")
    println("nnodes: ${sn.nnodes}")
    println("nclasses: ${sn.nclasses}")
    println("nclosedjobs: ${sn.nclosedjobs}")
    println("nchains: ${sn.nchains}")
    
    // All matrix fields
    printMatrix("refstat", sn.refstat)
    printMatrix("njobs", sn.njobs)
    printMatrix("nservers", sn.nservers)
    printMatrix("connmatrix", sn.connmatrix)
    printMatrix("scv", sn.scv)
    printMatrix("isstation", sn.isstation)
    printMatrix("isstateful", sn.isstateful)
    printMatrix("isstatedep", sn.isstatedep)
    printMatrix("nodeToStateful", sn.nodeToStateful)
    printMatrix("nodeToStation", sn.nodeToStation)
    printMatrix("stationToNode", sn.stationToNode)
    printMatrix("stationToStateful", sn.stationToStateful)
    printMatrix("statefulToStation", sn.statefulToStation)
    printMatrix("statefulToNode", sn.statefulToNode)
    printMatrix("rates", sn.rates)
    printMatrix("classprio", sn.classprio)
    printMatrix("phases", sn.phases)
    printMatrix("phasessz", sn.phasessz)
    printMatrix("phaseshift", sn.phaseshift)
    printMatrix("schedparam", sn.schedparam)
    printMatrix("chains", sn.chains)
    printMatrix("rt", sn.rt)
    printMatrix("nvars", sn.nvars)
    printMatrix("rtnodes", sn.rtnodes)
    printMatrix("csmask", sn.csmask)
    printMatrix("isslc", sn.isslc)
    printMatrix("cap", sn.cap)
    printMatrix("classcap", sn.classcap)
    printMatrix("refclass", sn.refclass)
    printMatrix("lldscaling", sn.lldscaling)
    printMatrix("fj", sn.fj)
    
    // List fields with content
    printList("nodetype", sn.nodetype)
    printList("classnames", sn.classnames)
    printList("nodenames", sn.nodenames)
    
    // Map fields with detailed contents
    printMapContents("rtorig", sn.rtorig)
    printMapContents("lst", sn.lst)
    printMapContents("state", sn.state)
    printMapContents("stateprior", sn.stateprior)
    printMapContents("space", sn.space)
    printMapContents("routing", sn.routing)
    printMapContents("procid", sn.procid)
    printMapContents("mu", sn.mu)
    printMapContents("phi", sn.phi)
    printMapContents("proc", sn.proc)
    printMapContents("pie", sn.pie)
    printMapContents("sched", sn.sched)
    printMapContents("inchain", sn.inchain)
    printMapContents("visits", sn.visits)
    printMapContents("nodevisits", sn.nodevisits)
    printMapContents("droprule", sn.droprule)
    printMapContents("nodeparam", sn.nodeparam)
    printMapContents("sync", sn.sync)
    printMapContents("gsync", sn.gsync)
    printMapContents("cdscaling", sn.cdscaling)
    
    // Object list contents
    sn.stations?.let { stations ->
        val stationNames = stations.map { it.name }
        printList("stations", stationNames)
    } ?: printList("stations", null)
    
    sn.stateful?.let { stateful ->
        val statefulNames = stateful.map { it.name }
        printList("stateful", statefulNames)
    } ?: printList("stateful", null)
    
    sn.jobclasses?.let { jobclasses ->
        val jobClassNames = jobclasses.map { it.name }
        printList("jobclasses", jobClassNames)
    } ?: printList("jobclasses", null)
    
    sn.nodes?.let { nodes ->
        val nodeNames = nodes.map { it.name }
        printList("nodes", nodeNames)
    } ?: printList("nodes", null)
}

/**
 * Helper function to print matrix in compact format for maps
 */
private fun printMatrixCompact(matrix: Matrix?): String {
    return when {
        matrix == null -> "null"
        matrix.isEmpty -> "[]"
        else -> {
            val sb = StringBuilder("[")
            val numRows = matrix.numRows
            val numCols = matrix.numCols
            for (i in 0 until numRows) {
                if (i > 0) sb.append("; ")
                for (j in 0 until numCols) {
                    if (j > 0) sb.append(" ")
                    val value = matrix[i, j]
                    when {
                        value.isNaN() -> sb.append("NaN")
                        value == GlobalConstants.MaxInt.toDouble() -> sb.append("Inf")
                        value.isInfinite() -> sb.append("Inf")
                        matrix.isInteger -> sb.append(Math.round(value).toInt())
                        value == Math.floor(value) && !value.isNaN() -> sb.append(value.toInt())
                        else -> sb.append(value)
                    }
                }
            }
            sb.append("]")
            sb.toString()
        }
    }
}

/**
 * Helper function to print lists with proper formatting
 */
private fun printList(name: String, list: List<*>?) {
    print("$name: ")
    if (list.isNullOrEmpty()) {
        println("[]")
    } else {
        print("[")
        list.forEachIndexed { index, item ->
            if (index > 0) print(", ")
            when (item) {
                is String -> print("\"$item\"")
                else -> print(item)
            }
        }
        println("]")
    }
}

/**
 * Helper function to print matrix with all values
 */
private fun printMatrix(name: String, matrix: Matrix?) {
    when {
        matrix == null -> println("$name: null")
        matrix.isEmpty -> println("$name: []")
        else -> {
            print("$name: ")
            print(printMatrixCompact(matrix))
            println()
        }
    }
}

/**
 * Helper function to print map contents with detailed structure
 */
private fun printMapContents(name: String, map: Map<*, *>?) {
    when {
        map == null -> println("$name: null")
        map.isEmpty() -> println("$name: {}")
        else -> {
            print("$name: {")
            var first = true
            // Sort entries by key name for consistent output
            val sortedEntries = map.entries.sortedBy { entry ->
                when (val key = entry.key) {
                    is String -> key
                    is Int -> key.toString()
                    null -> "null"
                    else -> {
                        try {
                            val getName = key.javaClass.getMethod("getName")
                            getName.invoke(key).toString()
                        } catch (e: Exception) {
                            key.javaClass.simpleName
                        }
                    }
                }
            }
            sortedEntries.forEach { (key, value) ->
                if (!first) print(", ")
                first = false
                
                // Print key
                when {
                    key is String -> print("\"$key\"")
                    key is Int -> print(key)
                    key != null -> {
                        try {
                            val getName = key.javaClass.getMethod("getName")
                            print("\"${getName.invoke(key)}\"")
                        } catch (e: Exception) {
                            print("\"${key.javaClass.simpleName}\"")
                        }
                    }
                    else -> print("null")
                }
                
                print(": ")
                
                // Print value
                printValue(value)
            }
            println("}")
        }
    }
}

/**
 * Helper function to print individual values with appropriate formatting
 */
private fun printValue(value: Any?) {
    when (value) {
        null -> print("null")
        is Map<*, *> -> {
            if (value.isEmpty()) {
                print("{}")
            } else {
                print("{")
                var innerFirst = true
                // Sort inner map entries for consistent output
                val sortedInnerEntries = value.entries.sortedBy { entry ->
                    when (val key = entry.key) {
                        is String -> key
                        is Number -> key.toString()
                        null -> "null"
                        else -> {
                            try {
                                val getName = key.javaClass.getMethod("getName")
                                getName.invoke(key).toString()
                            } catch (e: Exception) {
                                key.javaClass.simpleName
                            }
                        }
                    }
                }
                sortedInnerEntries.forEach { (innerKey, innerValue) ->
                    if (!innerFirst) print(", ")
                    innerFirst = false
                    
                    // Print inner key
                    when {
                        innerKey is String -> print("\"$innerKey\"")
                        innerKey is Int -> print(innerKey)
                        innerKey is Long -> print(innerKey)
                        innerKey != null -> {
                            try {
                                val getName = innerKey.javaClass.getMethod("getName")
                                print("\"${getName.invoke(innerKey)}\"")
                            } catch (e: Exception) {
                                // For numeric types, just print the value
                                if (innerKey is Number) {
                                    print(innerKey)
                                } else {
                                    print("${innerKey.javaClass.simpleName}")
                                }
                            }
                        }
                        else -> print("null")
                    }
                    
                    print(": ")
                    
                    // Print inner value
                    printValue(innerValue)
                }
                print("}")
            }
        }
        is Matrix -> print(printMatrixCompact(value))
        is String -> print("\"$value\"")
        is Boolean -> print(value)
        is Int -> if (value == GlobalConstants.MaxInt) print("Inf") else print(value)
        is Double -> when {
            value.isNaN() -> print("NaN")
            value == GlobalConstants.MaxInt.toDouble() -> print("Inf")
            value.isInfinite() -> print("Inf")
            else -> print(value)
        }
        is Float -> when {
            value.isNaN() -> print("NaN")
            value == GlobalConstants.MaxInt.toFloat() -> print("Inf")
            value.isInfinite() -> print("Inf")
            else -> print(value)
        }
        is Long -> if (value == GlobalConstants.MaxInt.toLong()) print("Inf") else print(value)
        is Enum<*> -> print(value.toString())
        is List<*> -> {
            print("[")
            value.forEachIndexed { index, item ->
                if (index > 0) print(", ")
                when (item) {
                    is Matrix -> print(printMatrixCompact(item))
                    else -> printValue(item)
                }
            }
            print("]")
        }
        else -> {
            val str = value.toString()
            when {
                // Handle function handles and lambdas
                str.contains("Lambda") && str.contains("$") -> print("<Function>")
                // Handle MatrixCell objects (for proc field)
                value is MatrixCell -> {
                    print("{")
                    for (i in 0 until value.size()) {
                        if (i > 0) print(", ")
                        print("[${i}]: ")
                        print(printMatrixCompact(value.get(i)))
                    }
                    print("}")
                }
                // Handle Event objects
                value is Event -> {
                    print("\"(${value.event.name}: node: ${value.node}, class: ${value.jobClass}")
                    if (!value.prob.isNaN() && value.prob != 1.0) {
                        print(", prob: ${value.prob}")
                    }
                    if (!value.t.isNaN()) {
                        print(", t: ${value.t}")
                    }
                    if (!value.job.isNaN()) {
                        print(", job: ${value.job}")
                    }
                    print(")\"")
                }
                // Handle Sync objects
                value is Sync -> {
                    print("{")
                    print("\"active\": ")
                    printValue(value.active)
                    print(", \"passive\": ")
                    printValue(value.passive)
                    print("}")
                }
                // Handle object references with @ notation
                str.contains("@") && str.contains(".") -> {
                    // Try to get object's name first
                    try {
                        val getName = value.javaClass.getMethod("getName")
                        print(getName.invoke(value))
                    } catch (e: Exception) {
                        // Fallback to simple class name
                        print(value.javaClass.simpleName)
                    }
                }
                else -> print(str)
            }
        }
    }
}
/**
 * Stochastic network Print algorithms
 */
@Suppress("unused")
class SnprintAlgo {
    companion object {
        // Class documentation marker for Dokka
    }
}