/**
 * @file Monte Carlo Integration (MCI) methods for normalizing constant computation
 * 
 * Implements Monte Carlo integration approaches including Importance Monte Carlo Integration
 * (IMCI) and standard MCI for computing normalizing constants in closed product-form queueing
 * networks. Uses variance reduction techniques with adaptive importance sampling.
 * 
 * @since LINE 3.0
 */
package jline.api.pfqn.nc

import jline.api.pfqn.mva.pfqn_bs
import jline.io.Ret
import jline.util.Maths
import jline.util.matrix.Matrix
import org.apache.commons.math3.util.FastMath
import java.util.*
import kotlin.math.ln

fun pfqn_mci(D: Matrix, N: Matrix, Z: Matrix, I: Int, variant: String): Ret.pfqnNc {
    val M = D.numRows
    val R = D.numCols

    val lGn: Double
    if (D.isEmpty || D.elementSum() < 1e-4) {
        lGn = -N.factln().elementSum() + N.scale(ln(Z.elementSum())).elementSum()
        val G = FastMath.exp(lGn)
        return Ret.pfqnNc(G, lGn)
    }

    var tput = Matrix(1, R)
    val util: Matrix
    val gamma = Matrix(1, M)

    when (variant.lowercase(Locale.getDefault())) {
        "imci" -> {
            tput = pfqn_bs(D, N, Z).X
            util = D.mult(tput.transpose())
            var i = 0
            while (i < M) {
                gamma[i] = FastMath.max(0.01, 1 - util[i])
                i++
            }
        }

        "mci" -> {
            tput = pfqn_bs(D, N, Z).X
            util = D.mult(tput.transpose())
            var i = 0
            while (i < M) {
                if (util[i] > 0.9) {
                    gamma[i] = 1 / FastMath.sqrt(N.elementMax())
                } else {
                    gamma[i] = 1 - util[i]
                }
                i++
            }
        }

        "rm" -> {
            var r = 0
            while (r < R) {
                val Dr = D.getColumn(r)
                tput[r] = N[r] / (Dr.elementSum() + Z[r] + Dr.elementMax() * (N.elementSum() - 1))
                r++
            }
            util = D.mult(tput.transpose())
            var i = 0
            while (i < M) {
                if (util[i] > 0.9) {
                    gamma[i] = 1 / FastMath.sqrt(N.elementMax())
                } else {
                    gamma[i] = 1 - util[i]
                }
                i++
            }
        }
    }

    val logfact = Matrix(1, R)
    for (r in 0..<R) {
        var n = 0
        while (n < N[r]) {
            logfact[r] = FastMath.log((1 + n).toDouble())
            n++
        }
    }

    val rand = Random()
    val VL = Matrix(I, M)
    for (i in 0..<I) {
        for (j in 0..<M) {
            VL[i, j] = FastMath.log(rand.nextDouble())
        }
    }

    val V = Matrix(I, M)
    for (i in 0..<I) {
        for (j in 0..<M) {
            V[i, j] = -1.0 / gamma[j] * VL[i, j]
        }
    }

    val ZI = Z.repmat(I, 1)
    val ones = Matrix.ones(1, M)
    val lZ = ones.sub(gamma).scale(-1.0).mult(V.transpose())
    lZ.subEq(gamma.log().elementSum())
    lZ.subEq(logfact.elementSum())
    lZ.addEq(N.mult(V.mult(D).add(ZI).log().transpose()))
    val lG = Maths.logmeanexp(lZ)

    return Ret.pfqnNc(FastMath.exp(lG), lG)
}
/**
 * PFQN mci algorithms
 */
@Suppress("unused")
class PfqnMciAlgo {
    companion object {
        // Class documentation marker for Dokka
    }
}