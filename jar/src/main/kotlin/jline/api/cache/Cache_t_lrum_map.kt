/**
 * @file LRUM-MAP Cache Response Time Analysis
 * 
 * Analyzes response times for LRUM cache systems with Markovian Arrival
 * Process (MAP) input streams. Combines MAP-based arrival modeling with
 * LRUM cache replacement policies for accurate performance evaluation.
 * 
 * @since LINE 3.0
 */
package jline.api.cache

import de.xypron.jcobyla.Calcfc
import de.xypron.jcobyla.Cobyla
import jline.api.mam.map_piq
import jline.util.Maths
import jline.util.matrix.Matrix
import jline.util.matrix.MatrixCell
import kotlin.math.floor
import kotlin.math.log10
import kotlin.math.pow

/**
 * APIs for stochastic models of caches
 */


fun cache_t_lrum_map(D0: Array<MatrixCell>, D1: Array<MatrixCell>, m: Matrix): Matrix {
    val n = D0.size
    val h = D0[0].size()

    val x = DoubleArray(h)
    for (i in 0..<h) {
        x[i] = 1.0
    }


    val rhobeg = 0.5
    val rhoend = 1.0e-6
    val maxFunEvals =
        500 * m.length() * 10.0.pow(floor(log10(n.toDouble()) / 4.0)).toInt() // maximum number of function evaluations

    val objectiveFunction = Calcfc { comn, comm, comx, comcon ->
        val result = lrummapTime(comx, D0, D1, m, n, h)
        val objVal = result.norm()
        objVal
    }

    Cobyla.findMinimum(objectiveFunction,
        h,
        h,
        x,
        rhobeg,
        rhoend,
        0,
        maxFunEvals) //1: print result; 0: no print result

    // Extract the optimized characteristic times
    val t = Matrix(1, h)
    for (i in 0..<h) {
        t[0, i] = x[i]
    }

    return t
}

fun lrummapTime(x: DoubleArray,
                D0Matrix: Array<MatrixCell>,
                D1Matrix: Array<MatrixCell>,
                m: Matrix,
                n: Int,
                h: Int): Matrix {
    val F = Matrix(1, h)

    val d = D0Matrix[0][0].numRows
    val Nh = arrayOfNulls<MatrixCell>(n)
    val N0 = arrayOfNulls<MatrixCell>(n)


    val expD0 = arrayOfNulls<MatrixCell>(n)
    val trans = arrayOfNulls<MatrixCell>(n)
    val pih = arrayOfNulls<MatrixCell>(n)
    val pi0 = arrayOfNulls<MatrixCell>(n)
    val R = arrayOfNulls<MatrixCell>(n)
    val neg_D0 = arrayOfNulls<MatrixCell>(n)

    val Q = arrayOfNulls<MatrixCell>(n)


    for (i in 0..<n) {
        expD0[i] = MatrixCell(h)
        trans[i] = MatrixCell(h)
        pih[i] = MatrixCell(h)
        pi0[i] = MatrixCell(1)
        R[i] = MatrixCell(h)
        Nh[i] = MatrixCell(h)
        N0[i] = MatrixCell(1)
        neg_D0[i] = MatrixCell(h)
        Q[i] = MatrixCell(h)
        pi0[i]!![0] = map_piq(D0Matrix[i][0], D1Matrix[i][0])
    }



    for (k in 0..<n) {
        for (l in 0..<h) {
            val index = Matrix.scaleMult(D0Matrix[k][l], x[l])
            val expD0kl = Maths.matrixExp(index)
            expD0[k]!![l] = expD0kl
            val temp = (Matrix.negative(D0Matrix[k][l])).inv()
            Nh[k]!![l] = Matrix.oneMinusMatrix(expD0kl).mult(temp)
            Nh[k]!![l]
            trans[k]!![l] = Matrix.oneMinusMatrix(expD0kl).mult(temp).mult(D1Matrix[k][l])
        }


        for (l in h - 1 downTo 0) {
            if (l == h - 1) {
                R[k]!![l] = (trans[k]!![l - 1]).mult(Matrix.oneMinusMatrix(trans[k]!![l]).inv())
            } else if (l == 0) {
                val temp = (Matrix.negative(D0Matrix[k][l])).inv().mult(D1Matrix[k][l]) //(A0)
                R[k]!![l] = temp.mult(Matrix.oneMinusMatrix(R[k]!![l + 1].mult(expD0[k]!![l + 1])).inv())
            } else {
                R[k]!![l] = trans[k]!![l - 1].mult(Matrix.oneMinusMatrix(R[k]!![l + 1].mult(expD0[k]!![l + 1])).inv())
            }
        }

        N0[k]!![0] = (Matrix.negative(D0Matrix[k][0])).inv()
    }

    val prod_R = arrayOfNulls<MatrixCell>(n)
    val sum_prod_R = MatrixCell()
    for (k in 0..<n) {
        sum_prod_R[k] = R[k]!![0]
        prod_R[k] = MatrixCell(h + 1)
        prod_R[k]!![0] = R[k]!![0]
        for (l in 1..<h) {
            prod_R[k]!![l] = prod_R[k]!![l - 1].mult(R[k]!![l])
            sum_prod_R[k] = sum_prod_R[k].add(prod_R[k]!![l])
        }
    }


    /*
    for (int k = 0; k < n; k++) {
        for (int l = 0; l < h; l++) {
            pih[k].set(l,map_pie(D0Matrix[k].get(l), D1Matrix[k].get(l)));
        }
    }
    */
    for (k in 0..<n) {
        pih[k]!![0] = pi0[k]!![0].mult(R[k]!![0])
        for (l in 1..<h) {
            pih[k]!![l] = pih[k]!![l - 1].mult(R[k]!![l])
        }
    }


    val e = Matrix(d, 1)
    e.ones()
    val production = Array(n) { DoubleArray(h + 1) }

    //double[][] productionD = new double[n][h + 1];
    for (k in 0..<n) {
        production[k][0] = (pi0[k]!![0]).mult(N0[k]!![0]).mult(e)[0, 0]
        //productionD[k][0] = (pi0[k].get(0)).mult(N0[k].get(0)).mult(e).get(0, 0);
    }
    for (l in 1..<h + 1) {
        for (k in 0..<n) {
            production[k][l] = (pih[k]!![l - 1]).mult(Nh[k]!![l - 1]).mult(e)[0, 0]
            //productionD[k][l] = (pih[k].get(l - 1)).mult(Nh[k].get(l - 1)).mult(D1Matrix[k].get(l-1)).mult(e).get(0, 0);
        }
    }

    val m_matrix = DoubleArray(h + 1)
    // double[] m_matrixD = new double[h + 1];
    for (l in 1..<h + 1) {
        m_matrix[l] = 0.0
        for (k in 0..<n) {
            var down_sum = 0.0
            for (j in 0..<h + 1) {
                down_sum += production[k][j]
            }
            m_matrix[l] += production[k][l] / down_sum
            //m_matrixD[l] += productionD[k][l] / down_sum;
        }
    }

    /*
    Matrix e = new Matrix(d, 1);
    e.ones();
    double[][] production = new double[n][h + 1];
    double[][] productionD = new double[n][h + 1];
    for (int k = 0; k < n; k++) {
        production[k][0] = (pi0[k].get(0)).mult(N0[k].get(0)).mult(e).get(0, 0);
        productionD[k][0] = (pi0[k].get(0)).mult(N0[k].get(0)).mult(e).get(0, 0);
    }
    for (int l = 1; l < h + 1; l++) {
        for (int k = 0; k < n; k++) {
            production[k][l] = (pih[k].get(l - 1)).mult(Nh[k].get(l - 1)).mult(e).get(0, 0);
            productionD[k][l] = (pih[k].get(l - 1)).mult(Nh[k].get(l - 1)).mult(D1Matrix[k].get(l-1)).mult(e).get(0, 0);
        }
    }

    double[] m_matrix = new double[h + 1];
    double[] m_matrixD = new double[h + 1];
    for (int l = 1; l < h + 1; l++) {
        m_matrix[l] = 0;
        for (int k = 0; k < n; k++) {
            double down_sum = 0;
            for (int j = 0; j < h + 1; j++) {
                down_sum += production[k][j];
            }
            m_matrix[l] += production[k][l] / down_sum;
            m_matrixD[l] += productionD[k][l] / down_sum;
        }
    }
    */
    for (l in 0..<h) {
        F[0, l] = m[0, l] - m_matrix[l + 1]
    }


    return F
}
/**
 * Cache t lrum map algorithms
 */
@Suppress("unused")
class CacheTLrumMapAlgo {
    companion object {
        // Class documentation marker for Dokka
    }
}