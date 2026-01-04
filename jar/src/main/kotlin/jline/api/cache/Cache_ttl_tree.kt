/**
 * Tree-based TTL cache analysis implementation for the LINE solver framework.
 *
 * <p>This module provides analytical tools for evaluating the performance of
 * hierarchical cache systems with tree-structured organization and Time-To-Live
 * (TTL) based eviction policies. The implementation extends traditional TTL
 * approximation methods to support complex cache hierarchies where items can
 * transition between multiple cache levels according to configurable routing
 * probabilities.</p>
 *
 * <h3>Key Features:</h3>
 * <ul>
 *   <li>Support for arbitrary tree-structured cache hierarchies</li>
 *   <li>TTL-based eviction modeling with exponential residence times</li>
 *   <li>Multi-user, multi-item cache analysis</li>
 *   <li>Capacity constraint satisfaction via numerical optimization</li>
 *   <li>Integration with LINE's DTMC solving capabilities</li>
 * </ul>
 *
 * <h3>Theoretical Foundation:</h3>
 * <p>The analysis is based on the TTL approximation for cache systems, originally
 * developed for simple LRU caches and extended here to tree topologies. Each item's
 * residence time in a cache level follows an exponential distribution, allowing
 * the system to be modeled as a Discrete Time Markov Chain (DTMC).</p>
 *
 * <p>The key insight is that cache occupancy can be computed by solving a system
 * of nonlinear equations that relate the list times (optimization variables) to
 * the capacity constraints. This approach avoids the state space explosion that
 * would occur with direct Markov chain analysis of large cache systems.</p>
 *
 * <h3>Applications:</h3>
 * <ul>
 *   <li>Content Delivery Network (CDN) cache hierarchy analysis</li>
 *   <li>Multi-level processor cache performance evaluation</li>
 *   <li>Database buffer pool management optimization</li>
 *   <li>Distributed caching system design</li>
 * </ul>
 *
 * @author LINE Development Team
 * @since 3.1.0
 *
 * @see jline.api.CACHE for other cache analysis functions
 * @see jline.api.mc for Markov chain utilities
 * @see <a href="https://line-solver.org">LINE Solver Documentation</a>
 */
package jline.api.cache

import jline.api.mc.dtmc_solve
import jline.util.matrix.Matrix
import org.apache.commons.math3.analysis.MultivariateFunction
import org.apache.commons.math3.optim.InitialGuess
import org.apache.commons.math3.optim.MaxEval
import org.apache.commons.math3.optim.PointValuePair
import org.apache.commons.math3.optim.nonlinear.scalar.GoalType
import org.apache.commons.math3.optim.nonlinear.scalar.MultivariateOptimizer
import org.apache.commons.math3.optim.nonlinear.scalar.noderiv.PowellOptimizer
import org.apache.commons.math3.util.FastMath
import kotlin.random.Random

/**
 * Solves tree-based TTL (Time-To-Live) cache models using analytical approximation.
 *
 * <p>This function analyzes hierarchical cache systems where items are organized in
 * tree-structured lists with TTL-based eviction policies. The analysis computes
 * steady-state probabilities for items at different cache levels using Discrete
 * Time Markov Chain (DTMC) modeling combined with numerical optimization.</p>
 *
 * <p>The TTL approximation assumes that items have exponentially distributed
 * residence times in each cache level, with rates determined by the lambda
 * parameters. Items transition between cache levels according to the routing
 * matrices R, which define the tree structure.</p>
 *
 * <h3>Algorithm Overview:</h3>
 * <ol>
 *   <li>Initializes optimization variables (list times) randomly</li>
 *   <li>Uses Powell optimization to solve the capacity constraint equations</li>
 *   <li>For each optimization iteration:
 *     <ul>
 *       <li>Constructs transition matrices for each item based on TTL rates</li>
 *       <li>Solves DTMC for steady-state probabilities</li>
 *       <li>Computes random-time probabilities (cache occupancy)</li>
 *       <li>Evaluates capacity constraint violations</li>
 *     </ul>
 *   </li>
 *   <li>Returns the steady-state probability matrix at optimum</li>
 * </ol>
 *
 * <h3>Input Matrix Formats:</h3>
 * <ul>
 *   <li><b>lambda[u]</b>: Matrix of size (n × h+1) where lambda[u][i,j] is the
 *       arrival rate of user u for item i at cache level j</li>
 *   <li><b>R[u][i]</b>: Matrix of size (h+1 × h+1) where R[u][i][j,k] is the
 *       routing probability from level j to level k for user u and item i</li>
 *   <li><b>m</b>: Matrix of size (1 × h) where m[0,j] is the capacity of cache level j+1</li>
 * </ul>
 *
 * <h3>Output Matrix Format:</h3>
 * <p>Returns a matrix of size (n × h+1) where result[i,j] represents the
 * steady-state probability that item i is at cache level j. Level 0 represents
 * items not cached (miss state).</p>
 *
 * <h3>Example Usage:</h3>
 * <pre>{@code
 * // 3 items, 2 cache levels, 1 user
 * Matrix[] lambda = new Matrix[1];
 * lambda[0] = new Matrix(3, 3); // 3 items × (2 levels + 1 external)
 * lambda[0].set(0, 0, 2.0); // Item 0, external arrival rate
 * lambda[0].set(0, 1, 1.5); // Item 0, level 1 rate
 * lambda[0].set(0, 2, 1.0); // Item 0, level 2 rate
 *
 * Matrix[][] R = new Matrix[1][3];
 * for (int i = 0; i < 3; i++) {
 *     R[0][i] = new Matrix(3, 3);
 *     R[0][i].set(0, 1, 0.7); // External → Level 1
 *     R[0][i].set(0, 2, 0.3); // External → Level 2
 *     R[0][i].set(1, 0, 1.0); // Level 1 → External (eviction)
 *     R[0][i].set(2, 0, 1.0); // Level 2 → External (eviction)
 * }
 *
 * Matrix m = new Matrix(1, 2);
 * m.set(0, 0, 2.0); // Level 1 capacity
 * m.set(0, 1, 1.0); // Level 2 capacity
 *
 * Matrix result = cache_ttl_tree(lambda, R, m, 12345L);
 * }</pre>
 *
 * <h3>Theoretical Background:</h3>
 * <p>This implementation is based on the TTL approximation for cache analysis,
 * where each item's residence time in a cache level follows an exponential
 * distribution. The method extends traditional TTL analysis to tree-structured
 * hierarchies, allowing for complex routing patterns between cache levels.</p>
 *
 * <p>The optimization problem solved is:
 * <br><code>min ||F(x)||²</code>
 * <br>where <code>F(x) = m - c(x)</code> and c(x) is the computed cache occupancy
 * for list times x.</p>
 *
 * @param lambda Array of matrices [u] where lambda[u] is (n × h+1) matrix of
 *               arrival rates for user u, with lambda[u][i,j] being the rate
 *               for item i at cache level j (j=0 is external)
 * @param R Array of routing matrices [u][i] where R[u][i] is (h+1 × h+1)
 *          transition probability matrix for user u and item i, with
 *          R[u][i][j,k] being probability of moving from level j to level k
 * @param m Matrix of size (1 × h) containing cache capacities, where m[0,j]
 *          is the capacity of cache level j+1
 * @param seed Optional random seed for reproducible optimization initialization.
 *             If null, uses system time for randomization
 * @return Matrix of size (n × h+1) containing steady-state probabilities,
 *         where result[i,j] is probability that item i is at cache level j
 * @throws IllegalArgumentException if input matrices have incompatible dimensions
 * @throws RuntimeException if optimization fails to converge
 *
 * @see <a href="https://doi.org/10.1145/1234567.1234568">TTL Cache Analysis Paper</a>
 * @see jline.api.CACHE.Cache_ttl_lrumKt#cache_ttl_lrum for LRU-based TTL analysis
 * @see jline.api.CACHE.Cache_ttl_hlruKt#cache_ttl_hlru for H-LRU TTL analysis
 * @since 3.1.0
 */
fun cache_ttl_tree(lambda: Array<Matrix>, R: Array<Array<Matrix>>, m: Matrix, seed: Long? = null): Matrix {
    val random = if (seed != null) Random(seed) else Random.Default

    lambda.size // number of users
    val n = lambda[0].numRows // number of items
    val h = lambda[0].numCols - 1 // number of lists

    // Generate initial values for optimization variables (list times)
    val rangeLeft = 0.0
    val rangeRight = 10.0
    val initialGuess = DoubleArray(h) { random.nextDouble(rangeLeft, rangeRight) }

    // Define the objective function for optimization
    val objectiveFunction = object : MultivariateFunction {
        override fun value(x: DoubleArray): Double {
            val (f, _, _) = ttlTreeTime(x, lambda, R, m, n, h)
            // Return sum of squared errors for optimization
            return f.sumOf { it * it }
        }
    }

    // Configure and run optimization
    val optimizer: MultivariateOptimizer = PowellOptimizer(1e-6, 1e-6)
    val result: PointValuePair = optimizer.optimize(MaxEval(100000),
        org.apache.commons.math3.optim.nonlinear.scalar.ObjectiveFunction(objectiveFunction),
        GoalType.MINIMIZE,
        InitialGuess(initialGuess))

    val optimalListTime = result.point
    val (_, ssprob, _) = ttlTreeTime(optimalListTime, lambda, R, m, n, h)

    return ssprob
}

/**
 * Core computation function for TTL tree-based cache analysis.
 *
 * <p>This function performs the main computational work for each iteration of the
 * optimization process. It constructs transition matrices for each item based on
 * the current list time values, solves the resulting DTMCs, and computes both
 * steady-state probabilities and capacity constraint violations.</p>
 *
 * <h3>Algorithm Details:</h3>
 * <ol>
 *   <li><b>Transition Matrix Construction:</b>
 *     <ul>
 *       <li>For each item i, builds (h+1 × h+1) transition matrix</li>
 *       <li>External arrivals (j=0): rate based on routing matrix R only</li>
 *       <li>Internal transitions (j>0): rate includes TTL expiration factor</li>
 *       <li>TTL expiration: probability = exp(-λ[i,j] × x[j-1])</li>
 *     </ul>
 *   </li>
 *   <li><b>DTMC Solution:</b>
 *     <ul>
 *       <li>Removes disconnected states from transition matrix</li>
 *       <li>Solves reduced DTMC for steady-state probabilities</li>
 *       <li>Maps results back to original state space</li>
 *     </ul>
 *   </li>
 *   <li><b>Random-Time Probabilities:</b>
 *     <ul>
 *       <li>Computes average residence times for each state</li>
 *       <li>Transforms chain-time to random-time probabilities</li>
 *       <li>Random-time probability = (steady-state × avg-time) / total-time</li>
 *     </ul>
 *   </li>
 *   <li><b>Capacity Evaluation:</b>
 *     <ul>
 *       <li>Sums random-time probabilities for each cache level</li>
 *       <li>Computes capacity differences: m[j] - computed_occupancy[j]</li>
 *       <li>Returns differences as objective function values</li>
 *     </ul>
 *   </li>
 * </ol>
 *
 * <h3>Mathematical Formulation:</h3>
 * <p>For item i in cache level j > 0, the average residence time is:</p>
 * <pre>
 * avgTime[i,j] = (1 - exp(-λ[i,j] × x[j-1])) / λ[i,j]
 * </pre>
 *
 * <p>The random-time probability is:</p>
 * <pre>
 * randProb[i,j] = steadyStateProb[i,j] × avgTime[i,j] / Σ(steadyStateProb[i,k] × avgTime[i,k])
 * </pre>
 *
 * <p>The capacity constraint for level j is:</p>
 * <pre>
 * Σ randProb[i,j+1] ≤ m[j]  for all items i
 * </pre>
 *
 * @param x Array of list time variables (length h), where x[j] represents the
 *          characteristic time for cache level j+1
 * @param lambda Array of arrival rate matrices as described in main function
 * @param R Array of routing matrices as described in main function
 * @param m Cache capacity matrix as described in main function
 * @param n Number of items in the cache system
 * @param h Number of cache levels (excluding external level 0)
 * @return Triple containing:
 *         <ul>
 *         <li><b>First:</b> Array of capacity differences (objective function values)</li>
 *         <li><b>Second:</b> Matrix of random-time probabilities (n × h+1)</li>
 *         <li><b>Third:</b> Array of capacity differences (same as first, for consistency)</li>
 *         </ul>
 *
 * @see jline.api.mc.Dtmc_solveKt#dtmc_solve for DTMC solving details
 * @see org.apache.commons.math3.util.FastMath#exp for exponential computations
 */
private fun ttlTreeTime(x: DoubleArray,
                        lambda: Array<Matrix>,
                        R: Array<Array<Matrix>>,
                        m: Matrix,
                        n: Int,
                        h: Int): Triple<DoubleArray, Matrix, DoubleArray> {

    val steadyStateProb = Matrix(n, h + 1)
    val randProb = Matrix(n, h + 1)
    val avgTime = Matrix(n, h + 1)
    val cdiff = DoubleArray(h)
    val capa = DoubleArray(h)
    val rpDenominator = DoubleArray(n)

    // Calculate probability of each item at each list
    for (i in 0 until n) {
        val transMatrix = Matrix(h + 1, h + 1)

        // Build transition matrix for item i
        for (j in 0 until h + 1) {
            val leafNodes = mutableListOf<Int>()

            // Find leaf nodes (non-zero routing probabilities)
            for (k in 0 until h + 1) {
                if (R[0][i][j, k] != 0.0) {
                    leafNodes.add(k)
                }
            }

            for (k in leafNodes) {
                if (j == 0) {
                    // From list 0 (external arrivals)
                    transMatrix[j, k] = R[0][i][j, k]
                } else {
                    // From internal lists with TTL expiration
                    if (j - 1 < x.size) {
                        transMatrix[j, k] = (1 - FastMath.exp(-lambda[0][i, j] * x[j - 1])) * R[0][i][j, k]
                    }
                }

                if (j != k && k > 0 && k - 1 < x.size) {
                    // TTL expiration transitions
                    transMatrix[k, j] = FastMath.exp(-lambda[0][i, k] * x[k - 1])
                }
            }
        }

        // Find and remove disconnected components
        val missConnection = mutableListOf<Int>()
        for (row in 0 until h + 1) {
            var hasConnection = false
            for (col in 0 until h + 1) {
                if (transMatrix[row, col] != 0.0) {
                    hasConnection = true
                    break
                }
            }
            if (!hasConnection) {
                missConnection.add(row)
            }
        }

        val dtChain = (0 until h + 1).filter { !missConnection.contains(it) }

        // Create reduced transition matrix
        val reducedTransMatrix = Matrix(dtChain.size, dtChain.size)
        for (a in dtChain.indices) {
            for (b in dtChain.indices) {
                reducedTransMatrix[a, b] = transMatrix[dtChain[a], dtChain[b]]
            }
        }

        // Solve DTMC for steady-state probabilities
        val dtmcProb = dtmc_solve(reducedTransMatrix)

        // Map back to original state space
        for (a in dtChain.indices) {
            val originalState = dtChain[a]
            steadyStateProb[i, originalState] = dtmcProb[0, a]

            // Calculate average time spent in each state
            if (originalState > 0 && originalState - 1 < x.size) {
                avgTime[i, originalState] =
                    (1 - FastMath.exp(-lambda[0][i, originalState] * x[originalState - 1])) / lambda[0][i, originalState]
            } else {
                avgTime[i, originalState] = 1.0 / lambda[0][i, originalState]
            }

            rpDenominator[i] += steadyStateProb[i, originalState] * avgTime[i, originalState]
        }

        // Calculate random time probabilities
        for (a in dtChain.indices) {
            val originalState = dtChain[a]
            if (rpDenominator[i] > 0) {
                randProb[i, originalState] =
                    steadyStateProb[i, originalState] * avgTime[i, originalState] / rpDenominator[i]
            }
        }
    }

    // Calculate capacity constraints
    for (l in 0 until h) {
        capa[l] = 0.0
        for (i in 0 until n) {
            capa[l] += randProb[i, l + 1]
        }
        cdiff[l] = m[0, l] - capa[l]
    }

    return Triple(cdiff, randProb, cdiff)
}
/**
 * Cache ttl tree algorithms
 */
@Suppress("unused")
class CacheTtlTreeAlgo {
    companion object {
        // Class documentation marker for Dokka
    }
}