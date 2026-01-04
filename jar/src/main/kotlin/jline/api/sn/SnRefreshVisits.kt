/**
 * Stochastic Network Visit Ratio Calculator
 * 
 * Computes visit ratios for each node and station in a queueing network by solving
 * the traffic equations using DTMC steady-state analysis. Visit ratios are fundamental
 * to queueing network analysis, determining relative utilization and throughput.
 *
 * @since LINE 3.0
 */
package jline.api.sn

import jline.api.mc.dtmc_solve
import jline.api.mc.dtmc_solve_reducible
import jline.lang.NetworkStruct
import jline.lang.constant.NodeType
import jline.util.matrix.Matrix

/**
 * Reculate the visits to each node
 *
 * @param sn      - NetworkStruct object for the queueing network model. The object will be updated by the method.
 * @param chains  - chain membership for each class
 * @param rt      - routing table among stations
 * @param rtnodes - routing table among nodes
 * @return updated sn structure
 */
fun snRefreshVisits(sn: NetworkStruct, chains: Matrix?, rt: Matrix, rtnodes: Matrix): NetworkStruct {
    val I = sn.nnodes
    val M = sn.nstateful
    val K = sn.nclasses
    val refstat = sn.refstat.copy() // Work on a copy to avoid modifying the original
    val nchains = sn.nchains

    /* Obtain chain characteristics */
    val inchain = sn.inchain
    for (c in 0..<nchains) {
        val inchain_c = inchain[c]
        val refstatValue = refstat[inchain_c!!.value().toInt(), 0]
        for (col in 1..<inchain_c.numCols) {
            val row = inchain_c[0, col].toInt()
            if (refstatValue != refstat[row, 0]) refstat[row, 0] = refstatValue
            //throw new RuntimeException("Classes within chain have different reference station");
        }
    }

    /* Transfer inchain to List<Integer> in order to reduce the time of type conversion (double -> int) which is time consuming) */
    val new_inchain: MutableMap<Int, List<Int>> = HashMap()
    for (c in 0..<nchains) {
        val inchain_c = inchain[c]
        val inchain_c_list: MutableList<Int> = ArrayList()
        for (i in 0..<inchain_c!!.numCols) inchain_c_list.add(inchain_c[i].toInt())
        new_inchain[c] = inchain_c_list
    }

    // Check once if network has Fork nodes (for routing matrix normalization)
    val hasFork = sn.nodetype.any { it == NodeType.Fork }

    /* Generate visits */
    val visits: MutableMap<Int, Matrix> = HashMap()
    for (c in 0..<nchains) {
        val inchain_c = new_inchain[c]!!
        val cols: MutableList<Int> =
            ArrayList() //If use JLineMatrix, there would be more data type transfer in Pchain creation
        for (ist in 0..<M) {
            for (ik in inchain_c.indices) {
                cols.add(ist * K + inchain_c[ik])
            }
        }

        //Pchain = rt(cols,cols);
        val Pchain = Matrix(cols.size, cols.size)
        for (row in cols.indices) {
            for (col in cols.indices) {
                Pchain[row, col] = rt[cols[row], cols[col]]
            }
        }

        // Normalize routing matrix for Fork-containing models
        // Fork nodes have row sums > 1 (sending to all branches with prob 1 each)
        // which causes dtmc_solve to fail. Normalize to make stochastic.
        if (hasFork) {
            for (row in 0..<Pchain.numRows) {
                val rs = Pchain.sumRows(row)
                if (rs > 1e-8) {
                    for (col in 0..<Pchain.numCols) {
                        Pchain[row, col] = Pchain[row, col] / rs
                    }
                }
            }
        }

        //visited = sum(Pchain,2) > 0;
        val visited = Matrix(Pchain.numRows, 1)
        var countTrue = 0
        for (row in 0..<Pchain.numRows) {
            if (Pchain.sumRows(row) > 0) {
                countTrue++
                visited[row, 0] = 1.0
            }
        }

        //alpha_visited = dtmc_solve(Pchain(visited,visited));
        val input = Matrix(countTrue, countTrue)
        var row_input = 0
        for (row in 0..<visited.numRows) {
            if (visited[row, 0] > 0) {
                var col_input = 0
                for (col in 0..<visited.numRows) {
                    if (visited[col, 0] > 0) {
                        input[row_input, col_input] = Pchain[row, col]
                        col_input++
                    }
                }
                row_input++
            }
        }

        // Use dtmc_solve as primary, fallback to dtmc_solve_reducible for chains with transient states
        var alpha_visited: Matrix
        try {
            alpha_visited = dtmc_solve(input)
            // Fallback to dtmc_solve_reducible if dtmc_solve fails (e.g., reducible chain)
            if (alpha_visited.elementMax() == 0.0 || alpha_visited.hasNaN()) {
                alpha_visited = dtmc_solve_reducible(input).first
            }
        } catch (e: Exception) {
            // Fallback to dtmc_solve_reducible if dtmc_solve throws an exception
            alpha_visited = dtmc_solve_reducible(input).first
        }

        //alpha = zeros(1,M*K); alpha(visited) = alpha_visited;
        val alpha = Matrix(1, M * K)
        var idx = 0
        for (row in 0..<visited.numRows) {
            if (visited[row, 0] > 0) alpha[0, row] = alpha_visited[0, idx++]
        }
        
        if (alpha.elementMax() >= 1.0 - 1e-12) {
            // Disabled because a self-looping customer is an absorbing chain
            // throw RuntimeException("One chain has an absorbing state.")
        }

        val visits_c = Matrix(M, K)
        for (ist in 0..<M) {
            for (k in inchain_c.indices) {
                visits_c[ist, inchain_c[k]] = alpha[0, ist * inchain_c.size + k]
            }
        }

        //visits{c} = visits{c} / sum(visits{c}(refstat(inchain{c}(1)),inchain{c}));
        var sum = 0.0
        val row = sn.stationToStateful[refstat[inchain_c[0]].toInt()].toInt()
        for (i in inchain_c.indices) {
            sum += visits_c[row, inchain_c[i]]
        }
        val visits_c_divide = Matrix(0, 0)
        visits_c.divide(sum, visits_c_divide, true)

        visits_c_divide.absEq()
        visits[c] = visits_c_divide
    }

    /* Generate node visits */
    val nodeVisits: MutableMap<Int, Matrix> = HashMap()
    for (c in 0..<nchains) {
        val inchain_c = new_inchain[c]!!
        val nodes_cols: MutableList<Int> =
            ArrayList() //If use JLineMatrix, there would be more data type transfer in Pchain creation
        for (ind in 0..<I) {
            for (ik in inchain_c.indices) {
                nodes_cols.add(ind * K + inchain_c[ik])
            }
        }

        val nodes_Pchain = Matrix(nodes_cols.size, nodes_cols.size)
        for (row in nodes_cols.indices) {
            for (col in nodes_cols.indices) {
                nodes_Pchain[row, col] = rtnodes[nodes_cols[row], nodes_cols[col]]
            }
        }

        // Normalize routing matrix for Fork-containing models
        if (hasFork) {
            for (row in 0..<nodes_Pchain.numRows) {
                val rs = nodes_Pchain.sumRows(row)
                if (rs > 1e-8) {
                    for (col in 0..<nodes_Pchain.numCols) {
                        nodes_Pchain[row, col] = nodes_Pchain[row, col] / rs
                    }
                }
            }
        }

        val nodes_visited = Matrix(nodes_Pchain.numRows, 1)
        var countTrue = 0
        for (row in 0..<nodes_Pchain.numRows) {
            if (nodes_Pchain.sumRows(row) > 0) {
                countTrue++
                nodes_visited[row, 0] = 1.0
            }
        }

        val input = Matrix(countTrue, countTrue)
        var row_input = 0
        for (row in 0..<nodes_visited.numRows) {
            if (nodes_visited[row, 0] > 0) {
                var col_input = 0
                for (col in 0..<nodes_visited.numRows) {
                    if (nodes_visited[col, 0] > 0) {
                        input[row_input, col_input] = nodes_Pchain[row, col]
                        col_input++
                    }
                }
                row_input++
            }
        }
        // Use dtmc_solve as primary, fallback to dtmc_solve_reducible for chains with transient states
        var nodes_alpha_visited: Matrix
        try {
            nodes_alpha_visited = dtmc_solve(input)
            // Fallback to dtmc_solve_reducible if dtmc_solve fails (e.g., reducible chain)
            if (nodes_alpha_visited.elementMax() == 0.0 || nodes_alpha_visited.hasNaN()) {
                nodes_alpha_visited = dtmc_solve_reducible(input).first
            }
        } catch (e: Exception) {
            // Fallback to dtmc_solve_reducible if dtmc_solve throws an exception
            nodes_alpha_visited = dtmc_solve_reducible(input).first
        }

        val nodes_alpha = Matrix(1, I * K)
        var idx = 0
        for (row in 0..<nodes_visited.numRows) {
            if (nodes_visited[row, 0] > 0) nodes_alpha[0, row] = nodes_alpha_visited[0, idx++]
        }

        val node_visits_c = Matrix(I, K)
        for (ind in 0..<I) {
            for (k in inchain_c.indices) {
                node_visits_c[ind, inchain_c[k]] = nodes_alpha[0, ind * inchain_c.size + k]
            }
        }

        val ref = refstat[inchain_c[0]].toInt()
        var sum = 0.0
        for (k in inchain_c.indices) {
            sum += node_visits_c[sn.statefulToNode[ref].toInt(), inchain_c[k]]
        }
        val node_visits_c_divide = Matrix(0, 0)
        node_visits_c.divide(sum, node_visits_c_divide, true)

        // Remove small numerical perturbations (negative values)
        for (i in 0..<node_visits_c_divide.numRows) {
            for (j in 0..<node_visits_c_divide.numCols) {
                if (node_visits_c_divide[i, j] < 0) {
                    node_visits_c_divide[i, j] = 0.0
                }
            }
        }
        
        // Remove NaN values  
        for (i in 0..<node_visits_c_divide.numRows) {
            for (j in 0..<node_visits_c_divide.numCols) {
                if (node_visits_c_divide[i, j].isNaN()) {
                    node_visits_c_divide[i, j] = 0.0
                }
            }
        }
        nodeVisits[c] = node_visits_c_divide
    }

    sn.visits = visits
    sn.nodevisits = nodeVisits
    sn.isslc = Matrix(sn.nclasses, 1)
    sn.refstat = refstat

    return sn
}
/**
 * Stochastic network RefreshVisits algorithms
 */
@Suppress("unused")
class SnrefreshvisitsAlgo {
    companion object {
        // Class documentation marker for Dokka
    }
}
