/**
 * @file Auxiliary EC terms computation for load-dependent mixed MVA
 * 
 * Provides auxiliary functionality for computing EC terms used in load-dependent Mean Value
 * Analysis for mixed queueing networks. Handles the complex mathematical calculations required
 * for analyzing systems with both open and closed classes and state-dependent service rates.
 * 
 * @since LINE 3.0
 */
package jline.api.pfqn.ld

import jline.io.Ret
import jline.util.matrix.Matrix
import org.apache.commons.math3.util.FastMath

/**
 * Auxiliary function used by pfqn_mvaldmx to compute the EC terms
 */

fun pfqn_mvaldmx_ec(lambda: Matrix, D: Matrix?, mu: Matrix): Ret.pfqnMVALDMXEC {
    val M = mu.numRows
    mu.numCols
    val Lo = Matrix(M, 1)
    for (i in 0..<M) {
        Lo[i] = lambda.mult(Matrix.extractRows(D, i, i + 1, null).transpose())[0]
    }
    val b = Matrix(M, 1) // Limited load dependence level
    for (i in 0..<M) {
        var idx = 0
        while (idx < mu.numCols && mu[i, idx] != mu[i, mu.numCols - 1]) {
            idx++
        }
        b[i] = idx.toDouble()
    }
    val Nt = mu.numCols // Compute extra elements if present
    val oldEnd = mu.numCols - 1
    mu.expandMatrix(mu.numRows,
        mu.numCols + 2 + b.elementMax().toInt(),
        mu.nonZeroLength + (2 + b.elementMax().toInt()) * Matrix.extractColumn(mu, oldEnd, null).nonZeroLength)
    for (i in 0..<mu.numRows) {
        for (j in oldEnd + 1..<mu.numCols) {
            mu[i, j] = mu[i, oldEnd]
        }
    }
    val C = Matrix(mu.numRows, mu.numCols)
    mu.divide(1.0, C, false)
    val EC = Matrix(M, Nt)
    val E = Matrix(M, 1 + Nt)
    val Eprime = Matrix(M, 1 + Nt)
    for (i in 0..<M) {
        val E1 = Matrix(1 + Nt, 1 + Nt)
        val E2 = Matrix(1 + Nt, 1 + Nt)
        val E3 = Matrix(1 + Nt, 1 + Nt)
        val F2 = Matrix(1 + Nt, 2 + b[i].toInt() - 2)
        val F3 = Matrix(1 + Nt, 2 + b[i].toInt() - 2)

        val E2prime = Matrix(1 + Nt, 1 + Nt)
        val F2prime = Matrix(1 + Nt, 2 + b[i].toInt() - 2)
        for (n in 0..Nt) {
            if (n >= b[i] + 1) {
                E[i, n] = 1.0 / FastMath.pow(1 - Lo[i] * C[i, b[i].toInt()], n + 1)
                Eprime[i, n] = C[i, b[i].toInt()] * E[i, n]
            } else {
                // Compute E1
                if (n == 0) {
                    E1[n] = 1 / (1 - Lo[i] * C[i, b[i].toInt()])
                    var j = 0
                    while (j < b[i] - 1 + 1) {
                        E1[n] = E1[n] * C[i, j] / C[i, b[i].toInt()]
                        j++
                    }
                } else {
                    E1[n] = 1 / (1 - Lo[i] * C[i, b[i].toInt()]) * C[i, b[i].toInt()] / C[i, n - 1] * E1[n - 1]
                }

                // Compute F2
                run {
                    var n0 = 0
                    while (n0 <= b[i] - 2 + 1) {
                        if (n0 == 0) {
                            F2[n, n0] = 1
                        } else {
                            F2[n, n0] = ((n.toDouble() + n0) / n0 * Lo[i] * C[i, n + n0 - 1] * F2[n, n0 - 1])
                        }
                        n0++
                    }
                }

                // Compute E2
                var sumf2 = 0.0
                run {
                    var k = -1
                    while (k < b[i] - 2 + 1) {
                        sumf2 += F2[n, 1 + k]
                        k++
                    }
                }
                E2[n] = sumf2

                // Compute F3
                run {
                    var n0 = 0
                    while (n0 <= b[i] - 2 + 1) {
                        if (n == 0 && n0 == 0) {
                            F3[n, n0] = 1
                            var j = 0
                            while (j < b[i] - 1 + 1) {
                                F3[n, n0] = F3[n, n0] * C[i, j] / C[i, b[i].toInt()]
                                j++
                            }
                        } else if (n > 0 && n0 == 0) {
                            F3[n, n0] = C[i, b[i].toInt()] / C[i, n - 1] * F3[n - 1, 0]
                        } else {
                            F3[n, n0] = (n.toDouble() + n0) / n0 * Lo[i] * C[i, b[i].toInt()] * F3[n, n0 - 1]
                        }
                        n0++
                    }
                }

                // Compute E3
                var sumf3 = 0.0
                run {
                    var k = -1
                    while (k < b[i] - 2 + 1) {
                        sumf3 += F3[n, 1 + k]
                        k++
                    }
                }
                E3[n] = sumf3

                // Compute F2prime
                var n0 = 0
                while (n0 <= b[i] - 2 + 1) {
                    if (n0 == 0) {
                        F2prime[n, n0] = C[i, n]
                    } else {
                        F2prime[n, n0] = (n.toDouble() + n0) / n0 * Lo[i] * C[i, n + n0] * F2prime[n, n0 - 1]
                    }
                    n0++
                }

                // Compute E2prime
                var sumf2p = 0.0
                var k = -1
                while (k < b[i] - 2 + 1) {
                    sumf2p += F2prime[n, 1 + k]
                    k++
                }
                E2prime[n] = sumf2p

                // Compute E, Eprime and EC
                E[i, n] = E1[n] + E2[n] - E3[n]
                if (n < b[i] - 1 + 1) {
                    Eprime[i, n] = C[i, b[i].toInt()] * E1[n] + E2prime[n] - C[i, b[i].toInt()] * E3[n]
                } else {
                    Eprime[i, n] = C[i, b[i].toInt()] * E[i, n]
                }
            }
        }
        for (n in 0..<Nt) {
            EC[i, n] = C[i, n] * E[i, n + 1] / E[i, n]
        }
    }
    return Ret.pfqnMVALDMXEC(EC, E, Eprime, Lo)
}
/**
 * PFQN mvaldmx ec algorithms
 */
@Suppress("unused")
class PfqnMvaldmxEcAlgo {
    companion object {
        // Class documentation marker for Dokka
    }
}