package jline.lib.m3a

import jline.util.matrix.Matrix
import jline.util.matrix.MatrixCell
import jline.util.Maths
import jline.api.mam.*
import org.apache.commons.math3.util.FastMath
import de.xypron.jcobyla.Cobyla
import de.xypron.jcobyla.Calcfc
import kotlin.math.*

/**
 * Utility functions for M3A (Markovian Arrival Process with 3-moment Approximation) compression.
 * 
 * This class provides supporting functions for moment computation, optimization,
 * and quality assessment for MMAP compression algorithms.
 */
class M3aUtils {
    
    companion object {
        
        /**
         * Computes the autocorrelation function of an MMAP up to the specified lag.
         * Uses the first two matrices (D0, D1) of the MMAP which form a valid MAP.
         * 
         * @param MMAP the MMAP to analyze
         * @param maxLag maximum lag to compute
         * @return array of autocorrelation values for lags 1 to maxLag
         */
        fun computeAutocorrelation(MMAP: MatrixCell, maxLag: Int): Array<Double> {
            // Extract D0 and D1 matrices which form a valid MAP
            val D0 = MMAP[0]
            val D1 = MMAP[1]
            
            // Create a 1 x maxLag matrix with lag values 1, 2, ..., maxLag
            val lags = Matrix(1, maxLag, maxLag)
            for (i in 0..<maxLag) {
                lags.set(0, i, (i + 1).toDouble())
            }
            
            val acfMatrix = map_acf(D0, D1, lags)
            return Array(maxLag) { acfMatrix.get(it) }
        }
        
        /**
         * Computes the index of dispersion for counts (IDC) of an MMAP.
         * Uses the first two matrices (D0, D1) of the MMAP which form a valid MAP.
         * 
         * @param MMAP the MMAP to analyze
         * @param timeWindow the time window for IDC computation (ignored, uses asymptotic IDC)
         * @return the index of dispersion for counts
         */
        fun computeIDC(MMAP: MatrixCell, timeWindow: Double): Double {
            // Extract D0 and D1 matrices which form a valid MAP
            val D0 = MMAP[0]
            val D1 = MMAP[1]
            
            return map_idc(D0, D1)
        }
        
        /**
         * Computes the coefficient of variation of an MMAP.
         * Uses the first two matrices (D0, D1) of the MMAP which form a valid MAP.
         * 
         * @param MMAP the MMAP to analyze
         * @return the coefficient of variation
         */
        fun computeCoeffVar(MMAP: MatrixCell): Double {
            // Extract D0 and D1 matrices which form a valid MAP
            val D0 = MMAP[0]
            val D1 = MMAP[1]
            
            return sqrt(map_scv(D0, D1))
        }
        
        /**
         * Computes the first n moments of an MMAP.
         * Uses the first two matrices (D0, D1) of the MMAP which form a valid MAP.
         * 
         * @param MMAP the MMAP to analyze
         * @param n number of moments to compute
         * @return array of the first n moments
         */
        fun computeMoments(MMAP: MatrixCell, n: Int): Array<Double> {
            // Extract D0 and D1 matrices which form a valid MAP
            val D0 = MMAP[0]
            val D1 = MMAP[1]
            
            val moments = Array(n) { 0.0 }
            for (k in 1..n) {
                moments[k - 1] = map_moment(D0, D1, k)
            }
            
            return moments
        }
        
        /**
         * Computes the spectral gap of an MMAP generator matrix.
         * 
         * @param MMAP the MMAP to analyze
         * @return the spectral gap (difference between largest and second-largest eigenvalue)
         */
        fun computeSpectralGap(MMAP: MatrixCell): Double {
            val Q = MMAP[0].add(1.0, MMAP[1])
            val eigenResult = Q.eigval()
            val eigenvalues = eigenResult.values
            
            // Convert eigenvalues to list and sort in descending order
            val eigenList = mutableListOf<Double>()
            for (i in 0 until eigenvalues.numRows) {
                eigenList.add(eigenvalues.get(i, 0))
            }
            val sortedEigenvalues = eigenList.sortedDescending()
            
            // Spectral gap is difference between largest and second-largest
            return if (sortedEigenvalues.size >= 2) {
                sortedEigenvalues[0] - sortedEigenvalues[1]
            } else {
                0.0
            }
        }
        
        /**
         * Optimizes MMAP parameters using COBYLA optimization.
         * 
         * @param initialParams initial parameter vector
         * @param objectiveFunction objective function to minimize
         * @param constraints constraint functions
         * @param tolerance optimization tolerance
         * @return optimized parameter vector
         */
        fun optimizeParameters(
            initialParams: Array<Double>,
            objectiveFunction: (Array<Double>) -> Double,
            constraints: Array<(Array<Double>) -> Double>,
            tolerance: Double = 1e-6
        ): Array<Double> {
            
            // Convert to DoubleArray for COBYLA
            val x = initialParams.toDoubleArray()
            
            val calcfc = Calcfc { n, m, xVars, con ->
                val xArray = Array(n) { xVars[it] }
                
                // Set constraint values (constraints are >= 0)
                for (i in 0 until constraints.size) {
                    con[i] = constraints[i](xArray)
                }
                
                // Return objective function value
                objectiveFunction(xArray)
            }
            
            Cobyla.findMinimum(
                calcfc,
                initialParams.size,
                constraints.size,
                x,
                0.5, // rhobeg
                tolerance, // rhoend
                0, // print level (0 = no output)
                1000 // max iterations
            )
            
            // Convert back to Array<Double>
            return Array(x.size) { x[it] }
        }
        
        /**
         * Computes the Kullback-Leibler divergence between two MMAPs.
         * 
         * @param MMAP1 first MMAP
         * @param MMAP2 second MMAP
         * @param numSamples number of samples for estimation
         * @param seed random seed for reproducible results
         * @return KL divergence estimate
         */
        fun computeKLDivergence(MMAP1: MatrixCell, MMAP2: MatrixCell, numSamples: Int = 10000, seed: Int = 23000): Double {
            val samples1 = sampleInterArrivalTimes(MMAP1, numSamples, seed)
            val samples2 = sampleInterArrivalTimes(MMAP2, numSamples, seed)
            
            // Compute empirical distributions
            val hist1 = computeHistogram(samples1, 50)
            val hist2 = computeHistogram(samples2, 50)
            
            // Compute KL divergence
            var kl = 0.0
            for (i in hist1.indices) {
                if (hist1[i] > 0 && hist2[i] > 0) {
                    kl += hist1[i] * log(hist1[i] / hist2[i], kotlin.math.E)
                }
            }
            
            return kl
        }
        
        /**
         * Validates that a matrix represents a valid MMAP.
         * 
         * @param MMAP the MMAP to validate
         * @return true if the MMAP is valid, false otherwise
         */
        fun validateMMAP(MMAP: MatrixCell): Boolean {
            if (MMAP.size() < 2) return false
            
            val D0 = MMAP[0]
            val D1 = MMAP[1]
            
            // Check dimensions
            if (D0.numRows != D0.numCols || D1.numRows != D1.numCols) return false
            if (D0.numRows != D1.numRows) return false
            
            // Check that D0 has non-positive diagonal elements
            for (i in 0..<D0.numRows) {
                if (D0.get(i, i) > 0) return false
            }
            
            // Check that off-diagonal elements of D0 are non-negative
            for (i in 0..<D0.numRows) {
                for (j in 0..<D0.numCols) {
                    if (i != j && D0.get(i, j) < 0) return false
                }
            }
            
            // Check that D1 has non-negative elements
            for (i in 0..<D1.numRows) {
                for (j in 0..<D1.numCols) {
                    if (D1.get(i, j) < 0) return false
                }
            }
            
            // Check that rows of D0 + D1 sum to zero
            val Q = D0.add(1.0, D1)
            val e = Matrix.ones(Q.numRows, 1)
            val rowSums = Q.mult(e)
            
            for (i in 0..<rowSums.numRows) {
                if (abs(rowSums.get(i, 0)) > 1e-10) return false
            }
            
            return true
        }
        
        // Helper methods
        private fun computeStationaryDistribution(Q: Matrix): Matrix {
            val n = Q.numRows
            val A = Q.transpose().add(1.0, Matrix.ones(n, n))
            val b = Matrix.ones(n, 1)
            val result = Matrix.zeros(n, 1)
            
            Matrix.solve(A, b, result)
            return result
        }
        
        private fun factorial(n: Int): Long {
            return if (n <= 1) 1 else n * factorial(n - 1)
        }
        
        private fun sampleInterArrivalTimes(MMAP: MatrixCell, numSamples: Int, seed: Int = 23000): Array<Double> {
            val samples = Array(numSamples) { 0.0 }
            
            // Simple simulation approach with fixed seed for reproducibility
            val random = kotlin.random.Random(seed)
            var state = 0
            val Q = MMAP[0].add(1.0, MMAP[1])
            val rates = Array(Q.numRows) { -Q.get(it, it) }
            
            for (i in 0..<numSamples) {
                // Generate exponential random variable
                val u = random.nextDouble()
                samples[i] = -log(u, kotlin.math.E) / rates[state]
                
                // Transition to next state (simplified)
                state = (state + 1) % Q.numRows
            }
            
            return samples
        }
        
        private fun computeHistogram(samples: Array<Double>, numBins: Int): Array<Double> {
            val hist = Array(numBins) { 0.0 }
            val minVal = samples.minOrNull() ?: 0.0
            val maxVal = samples.maxOrNull() ?: 1.0
            val binWidth = (maxVal - minVal) / numBins
            
            for (sample in samples) {
                val binIndex = min(numBins - 1, ((sample - minVal) / binWidth).toInt())
                hist[binIndex] += 1.0
            }
            
            // Normalize
            val total = hist.sum()
            for (i in hist.indices) {
                hist[i] /= total
            }
            
            return hist
        }
    }
}
