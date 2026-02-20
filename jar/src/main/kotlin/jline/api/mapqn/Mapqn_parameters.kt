/**
 * @file MAPQN Parameters Base Class
 * 
 * Defines the base parameter structure for MAP queueing network analysis.
 * Provides common interface and data structures for various MAPQN analysis
 * algorithms including bounds computation and performance evaluation.
 * 
 * @since LINE 3.0
 */
package jline.api.mapqn

import jline.lang.NetworkStruct
import jline.util.matrix.Matrix

/**
 * Base class for MAPQN model parameters
 */
abstract class Mapqn_parameters {
    abstract val M: Int  // number of queues
    abstract val N: Int  // population
    
    /**
     * Validate parameters
     */
    open fun validate() {
        require(M > 0) { "M must be positive" }
        require(N > 0) { "N must be positive" }
    }
    
    companion object {
        /**
         * Creates appropriate Mapqn_parameters from a NetworkStruct.
         * This factory method analyzes the network structure and returns
         * the appropriate Mapqn_parameters subclass.
         * 
         * @param networkStruct The network structure to convert
         * @return Appropriate Mapqn_parameters instance
         * @throws IllegalArgumentException if the network structure is not supported
         */
        @JvmStatic
        fun fromNetworkStruct(networkStruct: NetworkStruct): Mapqn_parameters {
            return Mapqn_parameters_factory.createFromNetworkStruct(networkStruct)
        }
    }
}

/**
 * Parameters for linear reduction models with phases
 */
class LinearReductionParameters(
    override val M: Int,
    override val N: Int,
    val K: IntArray,  // number of phases for each queue
    val mu: Array<Matrix>,  // completion transition rates [i][k][h]
    val r: Matrix,  // routing probabilities [i][j]
    val v: Array<Matrix>   // background transition rates [i][k][h]
) : Mapqn_parameters() {
    
    override fun validate() {
        super.validate()
        require(K.size == M) { "K array size must equal M" }
        require(mu.size == M) { "mu array size must equal M" }
        require(r.numRows == M && r.numCols == M) { "r must be MxM matrix" }
        require(v.size == M) { "v array size must equal M" }
        
        for (i in 0 until M) {
            require(K[i] > 0) { "K[$i] must be positive" }
            require(mu[i].numRows == K[i] && mu[i].numCols == K[i]) { 
                "mu[$i] must be ${K[i]}x${K[i]} matrix" 
            }
            require(v[i].numRows == K[i] && v[i].numCols == K[i]) { 
                "v[$i] must be ${K[i]}x${K[i]} matrix" 
            }
            // Check non-negativity
            for (row in 0 until r.numRows) {
                for (col in 0 until r.numCols) {
                    require(r.get(row, col) >= 0) { "Routing probabilities must be non-negative" }
                }
            }
            for (row in 0 until mu[i].numRows) {
                for (col in 0 until mu[i].numCols) {
                    require(mu[i].get(row, col) >= 0) { "Service rates must be non-negative" }
                }
            }
            for (row in 0 until v[i].numRows) {
                for (col in 0 until v[i].numCols) {
                    require(v[i].get(row, col) >= 0) { "Background rates must be non-negative" }
                }
            }
        }
    }
    
    /**
     * Compute q parameter as defined in AMPL model
     */
    fun q(i: Int, j: Int, k: Int, h: Int): Double {
        return if (j != i) {
            r.get(i, j) * mu[i].get(k, h)
        } else {
            v[i].get(k, h) + r.get(i, i) * mu[i].get(k, h)
        }
    }
}

/**
 * Parameters for MVA version models
 */
class MVAVersionParameters(
    override val M: Int,
    override val N: Int,
    val K: Int,  // number of levels (scalar in MVA version)
    val muM: DoubleArray,  // service rate for queues 1 to M-1
    val muMAP: Matrix,  // service rate for MAP queue [k][h]
    val r: Matrix,  // routing probabilities [i][j]
    val v: Matrix   // level change rates [k][h]
) : Mapqn_parameters() {
    
    override fun validate() {
        super.validate()
        require(K > 0) { "K must be positive" }
        require(muM.size == M - 1) { "muM array size must equal M-1" }
        require(muMAP.numRows == K && muMAP.numCols == K) { "muMAP must be KxK matrix" }
        require(r.numRows == M && r.numCols == M) { "r must be MxM matrix" }
        require(v.numRows == K && v.numCols == K) { "v must be KxK matrix" }
        
        require(muM.all { it >= 0 }) { "Service rates must be non-negative" }
        // Check non-negativity for matrices
        for (row in 0 until muMAP.numRows) {
            for (col in 0 until muMAP.numCols) {
                require(muMAP.get(row, col) >= 0) { "MAP service rates must be non-negative" }
            }
        }
        for (row in 0 until r.numRows) {
            for (col in 0 until r.numCols) {
                require(r.get(row, col) >= 0) { "Routing probabilities must be non-negative" }
            }
        }
        for (row in 0 until v.numRows) {
            for (col in 0 until v.numCols) {
                require(v.get(row, col) >= 0) { "Level change rates must be non-negative" }
            }
        }
    }
    
    /**
     * Compute q parameter as defined in AMPL model
     */
    fun q(i: Int, j: Int, k: Int, h: Int): Double {
        return if (i < M - 1) {
            if (k == h) r.get(i, j) * muM[i] else 0.0
        } else {  // i == M-1 (MAP queue)
            if (j < M - 1) {
                r.get(M - 1, j) * muMAP.get(k, h)
            } else {
                if (k != h) {
                    v.get(k, h) + r.get(M - 1, M - 1) * muMAP.get(k, h)
                } else {
                    0.0
                }
            }
        }
    }
}

/**
 * Parameters for QR Bounds BAS (Blocking After Service) model
 */
class Mapqn_qr_bounds_bas_parameters(
    override val M: Int,  // number of queues
    override val N: Int,  // population
    val MR: Int,  // number of independent blocking configurations
    val f: Int,  // finite capacity queue
    val K: IntArray,  // number of phases for each queue
    val F: IntArray,  // capacity for each queue
    val MM: Matrix,  // blocking order [m][i]
    val MM1: Matrix,  // blocking order [m][i]
    val ZZ: IntArray,  // nonzeros in independent blocking configurations
    val BB: Matrix,  // blocking state [m][i]
    val mu: Array<Matrix>,  // completion transition rates [i][k][h]
    val v: Array<Matrix>,  // background transition rates [i][k][h]
    val r: Matrix  // routing probabilities [i][j]
) : Mapqn_parameters() {
    
    val ZM: Int = ZZ.maxOrNull() ?: 0  // max of ZZ
    
    override fun validate() {
        super.validate()
        require(MR > 0) { "MR must be positive" }
        require(f in 1..M) { "f must be between 1 and M" }
        require(K.size == M) { "K array size must equal M" }
        require(F.size == M) { "F array size must equal M" }
        require(MM.numRows == MR && MM.numCols == 2) { "MM must be MRx2 matrix (blocking order pairs)" }
        require(MM1.numRows == MR && MM1.numCols == M) { "MM1 must be MRxM matrix" }
        require(ZZ.size == MR) { "ZZ array size must equal MR" }
        require(BB.numRows == MR && BB.numCols == M) { "BB must be MRxM matrix" }
        require(mu.size == M) { "mu array size must equal M" }
        require(v.size == M) { "v array size must equal M" }
        require(r.numRows == M && r.numCols == M) { "r must be MxM matrix" }
        
        for (i in 0 until M) {
            require(K[i] > 0) { "K[$i] must be positive" }
            require(F[i] > 0) { "F[$i] must be positive" }
            require(mu[i].numRows == K[i] && mu[i].numCols == K[i]) { 
                "mu[$i] must be ${K[i]}x${K[i]} matrix" 
            }
            require(v[i].numRows == K[i] && v[i].numCols == K[i]) { 
                "v[$i] must be ${K[i]}x${K[i]} matrix" 
            }
        }
        
        // Check non-negativity and constraints
        for (row in 0 until MM.numRows) {
            for (col in 0 until MM.numCols) {
                require(MM.get(row, col) >= 0) { "MM values must be non-negative" }
            }
        }
        require(ZZ.all { it >= 0 }) { "ZZ values must be non-negative" }
        for (row in 0 until BB.numRows) {
            for (col in 0 until BB.numCols) {
                val value = BB.get(row, col).toInt()
                require(value in 0..1) { "BB values must be 0 or 1" }
            }
        }
        for (row in 0 until r.numRows) {
            for (col in 0 until r.numCols) {
                require(r.get(row, col) >= 0) { "Routing probabilities must be non-negative" }
            }
        }
        for (i in 0 until M) {
            for (row in 0 until mu[i].numRows) {
                for (col in 0 until mu[i].numCols) {
                    require(mu[i].get(row, col) >= 0) { "Service rates must be non-negative" }
                }
            }
            for (row in 0 until v[i].numRows) {
                for (col in 0 until v[i].numCols) {
                    require(v[i].get(row, col) >= 0) { "Background rates must be non-negative" }
                }
            }
        }
    }
    
    /**
     * Compute q parameter as defined in AMPL model (scalar version)
     */
    fun q(i: Int, j: Int, k: Int, h: Int): Double {
        return if (j != i) {
            r.get(i, j) * mu[i].get(k, h)
        } else {
            v[i].get(k, h) + r.get(i, i) * mu[i].get(k, h)
        }
    }

    /**
     * Compute transition rates q(i,j) as 2D array of matrices
     * q[i][j] is a K[i] x K[i] matrix
     * Matches MATLAB qrf_bas.m computation
     */
    fun q(): Array<Array<Matrix>> {
        val result = Array(M) { i ->
            Array(M) { j ->
                val mat = Matrix(K[i], K[i])
                for (ki in 0 until K[i]) {
                    for (hi in 0 until K[i]) {
                        mat.set(ki, hi, q(i, j, ki, hi))
                    }
                }
                mat
            }
        }
        return result
    }
}

/**
 * Parameters for QR Bounds RSRD (Repetitive Service Random Destination) model
 */
class Mapqn_qr_bounds_rsrd_parameters(
    override val M: Int,  // number of queues
    override val N: Int,  // population
    val F: IntArray,  // capacity of each queue
    val K: IntArray,  // number of phases for each queue
    val mu: Array<Matrix>,  // completion transition rates [i][k][h]
    val v: Array<Matrix>,  // background transition rates [i][k][h]
    val alpha: Array<DoubleArray>,  // load dependent rates [i][n], default 1
    val r: Matrix  // routing probabilities [i][j]
) : Mapqn_parameters() {
    
    override fun validate() {
        super.validate()
        require(F.size == M) { "F array size must equal M" }
        require(K.size == M) { "K array size must equal M" }
        require(mu.size == M) { "mu array size must equal M" }
        require(v.size == M) { "v array size must equal M" }
        require(alpha.size == M) { "alpha array size must equal M" }
        require(r.numRows == M && r.numCols == M) { "r must be MxM matrix" }
        
        for (i in 0 until M) {
            require(F[i] > 0) { "F[$i] must be positive" }
            require(K[i] > 0) { "K[$i] must be positive" }
            require(alpha[i].size == N) { "alpha[$i] array size must equal N" }
            require(mu[i].numRows == K[i] && mu[i].numCols == K[i]) { 
                "mu[$i] must be ${K[i]}x${K[i]} matrix" 
            }
            require(v[i].numRows == K[i] && v[i].numCols == K[i]) { 
                "v[$i] must be ${K[i]}x${K[i]} matrix" 
            }
            require(alpha[i].all { it >= 0 }) { "alpha[$i] values must be non-negative" }
        }
        
        // Check non-negativity
        for (row in 0 until r.numRows) {
            for (col in 0 until r.numCols) {
                require(r.get(row, col) >= 0) { "Routing probabilities must be non-negative" }
            }
        }
        for (i in 0 until M) {
            for (row in 0 until mu[i].numRows) {
                for (col in 0 until mu[i].numCols) {
                    require(mu[i].get(row, col) >= 0) { "Service rates must be non-negative" }
                }
            }
            for (row in 0 until v[i].numRows) {
                for (col in 0 until v[i].numCols) {
                    require(v[i].get(row, col) >= 0) { "Background rates must be non-negative" }
                }
            }
        }
    }
    
    /**
     * Compute q parameter as defined in AMPL model
     */
    fun q(i: Int, j: Int, k: Int, h: Int, n: Int): Double {
        return if (n == 0) {
            0.0
        } else {
            if (j != i) {
                r.get(i, j) * mu[i].get(k, h) * alpha[i][n - 1]
            } else {
                v[i].get(k, h) * alpha[i][n - 1] + r.get(i, i) * mu[i].get(k, h) * alpha[i][n - 1]
            }
        }
    }
}