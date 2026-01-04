package jline.lib.butools

import jline.api.mc.ctmc_solve
import jline.util.Maths
import jline.util.matrix.Matrix
import jline.util.matrix.MatrixCell
import org.apache.commons.math3.util.CombinatoricsUtils

/**
 * Analyzes multi-class MMAP[K]/PH[K]/1 queue with preemptive-resume (PRPR) priority.
 *
 * This function computes the moments of queue length and sojourn time distributions
 * for a multi-class Markovian Arrival Process (MMAP) with class-dependent Phase-Type
 * (PH) service time distributions in a single-server queue with preemptive-resume
 * priority scheduling.
 *
 * In preemptive-resume scheduling, when a higher-priority job arrives, it preempts the
 * current job in service. The preempted job resumes its service later from where it
 * was interrupted.
 *
 * Reference: G. Horvath, "Efficient analysis of the MMAP[K]/PH[K]/1 priority queue",
 * European Journal of Operational Research, 246(1), 128-139, 2015.
 *
 * @param D The MMAP arrival process matrices (D0...DK). D1 corresponds to lowest priority,
 *          DK to highest priority.
 * @param sigma The initial phase probability vectors for service distributions (indexed 0..K-1)
 * @param S The transient generators of PH service distributions (indexed 0..K-1)
 * @param numOfQLMoms Number of queue length moments to compute (default: null)
 * @param numOfQLProbs Number of queue length probability values (default: null)
 * @param numOfSTMoms Number of sojourn time moments (default: null)
 * @param stCdfPoints Matrix of points for sojourn time CDF evaluation (default: null)
 * @param prec Numerical precision for iterations (default: 1e-14)
 * @param erlMaxOrder_ Maximum order for Erlang approximation (default: 200)
 * @param classes_ Indices of classes to analyze (default: all classes, 0-indexed)
 * @return Map of performance measures including "ncMoms", "ncDistr", "stMoms", "stDistr"
 */
fun MMAPPH1PRPR(D: MatrixCell,
                sigma: MatrixCell,
                S: MatrixCell,
                numOfQLMoms: Int?,
                numOfQLProbs: Int?,
                numOfSTMoms: Int?,
                stCdfPoints: Matrix?,
                prec: Double?,
                erlMaxOrder_: Int?,
                classes_: Matrix?): Map<String, MutableMap<Int, Matrix>> {

    val K = D.size() - 1

    // Parse options
    val erlMaxOrder = erlMaxOrder_ ?: 200
    val precision = prec ?: 1e-14

    var classes = Matrix(1, K, K)
    for (i in 0 until K) {
        classes[i] = i.toDouble()
    }
    if (classes_ != null) {
        classes = classes_
    }

    // Some preparation
    val D0 = D[0]
    val N = D0.numRows
    val I = Matrix.eye(N)

    // sD = sum of all D matrices
    var sD = Matrix(N, N, N * N)
    for (i in 0..K) {
        sD = sD.add(1.0, D[i])
    }

    // s{i} = -S{i} * ones (exit rate vectors)
    val s: MutableMap<Int, Matrix> = HashMap()
    // M(i) = number of phases in service distribution for class i
    val M = Matrix(1, K, K)
    for (i in 0 until K) {
        val neg_Si = S[i].copy()
        neg_Si.scaleEq(-1.0)
        s[i] = neg_Si.sumCols()
        M[i] = sigma[i].length().toDouble()
    }

    // Initialize result maps
    val result: MutableMap<String, MutableMap<Int, Matrix>> = HashMap()
    if (numOfSTMoms != null) {
        result["stMoms"] = HashMap()
    }
    if (stCdfPoints != null) {
        result["stDistr"] = HashMap()
    }
    if (numOfQLMoms != null) {
        result["ncMoms"] = HashMap()
    }
    if (numOfQLProbs != null) {
        result["ncDistr"] = HashMap()
    }

    // Process each class
    for (g in 0 until classes.length()) {
        val k = classes[g].toInt()

        // ======================================================
        // Step 1. Solution of the workload process of the system
        // ======================================================

        // sM = sum of phases for classes k to K-1
        var sM = 0
        for (i in k until K) {
            sM += M[i].toInt()
        }

        // Qwmm = D0 + D1 + ... + Dk (arrivals of classes 0..k-1 don't affect workload)
        var Qwmm = D0.copy()
        for (i in 0 until k) {
            Qwmm = Qwmm.add(1.0, D[i + 1])
        }

        // Build Qwpm, Qwmp, Qwpp for classes k to K-1
        val Qwpm = Matrix(N * sM, N, N * sM * N)
        val Qwmp = Matrix(N, N * sM, N * N * sM)
        val Qwpp = Matrix(N * sM, N * sM, N * sM * N * sM)

        var kix = 0
        for (i in k until K) {
            val bs = N * M[i].toInt()
            // Qwmp(:, kix:kix+bs-1) = kron(D{i+1}, sigma{i})
            Qwmp.insertSubMatrix(0, kix, N, kix + bs, D[i + 1].kron(sigma[i]))
            // Qwpm(kix:kix+bs-1, :) = kron(I, s{i})
            Qwpm.insertSubMatrix(kix, 0, kix + bs, N, I.kron(s[i]!!))
            // Qwpp(kix:kix+bs-1, kix:kix+bs-1) = kron(I, S{i})
            Qwpp.insertSubMatrix(kix, kix, kix + bs, kix + bs, I.kron(S[i]))
            kix += bs
        }

        // Calculate fundamental matrices
        val fluidResult = FluidFundamentalMatrices(Qwpp, Qwpm, Qwmp, Qwmm, precision, null, null)
        val Psiw = fluidResult["P"]!!
        val Kw = fluidResult["K"]!!
        val Uw = fluidResult["U"]!!

        // Calculate boundary vector pm
        val neg_Kw = Kw.copy()
        neg_Kw.scaleEq(-1.0)
        val Ua = Matrix.ones(N, 1).add(2.0, Qwmp.mult(neg_Kw.pinv()).sumRows())

        // Solve overdetermined system using least-squares (MATLAB linsolve equivalent)
        val A_ls = Matrix.concatColumns(Uw, Ua, null).transpose()
        val b_ls = Matrix.concatColumns(Matrix(1, N, 0), Matrix.singleton(1.0), null).transpose()
        var pm = A_ls.leftMatrixDivide(b_ls).transpose()

        // Bw matrix - exit from class k service
        val Bw = Matrix(N * sM, N, N * sM * N)
        Bw.insertSubMatrix(0, 0, N * M[k].toInt(), N, I.kron(s[k]!!))

        // kappa = initial distribution for sojourn time analysis
        val denom = pm.mult(Qwmp).mult(neg_Kw.pinv()).mult(Bw).elementSum()
        val kappa = pm.mult(Qwmp)
        kappa.scaleEq(1.0 / denom)

        if (k < K - 1) {
            // ====================================================================
            // Step 2. Construct fluid model for the remaining sojourn time process
            // ====================================================================

            // Qsmm = D0 + D1 + ... + D{k+1}
            var Qsmm = D0.copy()
            for (i in 0..k) {
                Qsmm = Qsmm.add(1.0, D[i + 1])
            }

            val Np = Kw.numRows
            // Size for higher priority classes (k+1 to K-1)
            var sizeHigh = 0
            for (i in k + 1 until K) {
                sizeHigh += N * M[i].toInt()
            }

            val Qspm = Matrix(Np + sizeHigh, N, (Np + sizeHigh) * N)
            val Qsmp = Matrix(N, Np + sizeHigh, N * (Np + sizeHigh))
            val Qspp = Matrix(Np + sizeHigh, Np + sizeHigh, (Np + sizeHigh) * (Np + sizeHigh))

            // Insert Kw and Bw for the workload process
            Qspp.insertSubMatrix(0, 0, Np, Np, Kw)
            Qspm.insertSubMatrix(0, 0, Np, N, Bw)

            // Add higher priority classes (k+1 to K-1)
            kix = Np
            for (i in k + 1 until K) {
                val bs = N * M[i].toInt()
                Qsmp.insertSubMatrix(0, kix, N, kix + bs, D[i + 1].kron(sigma[i]))
                Qspm.insertSubMatrix(kix, 0, kix + bs, N, I.kron(s[i]!!))
                Qspp.insertSubMatrix(kix, kix, kix + bs, kix + bs, I.kron(S[i]))
                kix += bs
            }

            // Initial distribution for sojourn time fluid model
            val inis = Matrix.concatColumns(kappa, Matrix(1, sizeHigh, 0), null)

            // Solve for Psis
            val Psis = FluidFundamentalMatrices(Qspp, Qspm, Qsmp, Qsmm, precision, null, null)["P"]!!

            // ==========================================
            // Step 3. Calculate the performance measures
            // ==========================================

            // MOMENTS OF THE SOJOURN TIME
            if (numOfSTMoms != null) {
                val Pn: MutableMap<Int, Matrix> = HashMap()
                Pn[0] = Psis
                val rtMoms = Matrix(1, numOfSTMoms, numOfSTMoms)

                for (n in 1..numOfSTMoms) {
                    val A = Qspp.add(1.0, Psis.mult(Qsmp))
                    val B = Qsmm.add(1.0, Qsmp.mult(Psis))
                    var C = Pn[n - 1]!!.copy()
                    C.scaleEq((-2.0 * n))

                    var bino = 1
                    for (i in 1 until n) {
                        bino = bino * (n - i + 1) / i
                        C = C.add(bino.toDouble(), Pn[i]!!.mult(Qsmp).mult(Pn[n - i]!!))
                    }

                    val P = Matrix.lyap(A, B, C, null)
                    Pn[n] = P
                    rtMoms[n - 1] = inis.mult(P).elementSum() * Math.pow(-1.0, n.toDouble()) / Math.pow(2.0, n.toDouble())
                }
                result["stMoms"]!![k] = rtMoms
            }

            // DISTRIBUTION OF THE SOJOURN TIME
            if (stCdfPoints != null) {
                var res = Matrix(1, stCdfPoints.length(), stCdfPoints.length())
                for (o in 0 until stCdfPoints.length()) {
                    val t = stCdfPoints[o]
                    val L = erlMaxOrder
                    val lambda = L / t / 2.0

                    val Psie = FluidFundamentalMatrices(
                        Qspp.add(-lambda, Matrix.eye(Qspp.numRows)),
                        Qspm,
                        Qsmp,
                        Qsmm.add(-lambda, Matrix.eye(Qsmm.numRows)),
                        precision, null, null
                    )["P"]!!

                    val Pn: MutableMap<Int, Matrix> = HashMap()
                    Pn[0] = Psie
                    var pr = inis.mult(Psie).elementSum()

                    for (n in 1 until L) {
                        val A = Qspp.add(1.0, Psie.mult(Qsmp)).add(-lambda, Matrix.eye(Qspp.numRows))
                        val B = Qsmm.add(1.0, Qsmp.mult(Psie)).add(-lambda, Matrix.eye(Qsmm.numRows))
                        var C = Pn[n - 1]!!.copy()
                        C.scaleEq(2.0 * lambda)

                        for (i in 1 until n) {
                            C = C.add(1.0, Pn[i]!!.mult(Qsmp).mult(Pn[n - i]!!))
                        }

                        val P = Matrix.lyap(A, B, C, null)
                        Pn[n] = P
                        pr += inis.mult(P).elementSum()
                    }
                    res[o] = pr
                }
                result["stDistr"]!![k] = res
            }

            // MOMENTS OF THE NUMBER OF JOBS
            if (numOfQLMoms != null) {
                // First calculate at departure instants
                val QLDPn: MutableMap<Int, Matrix> = HashMap()
                QLDPn[0] = Psis
                val dqlMoms = Matrix(1, numOfQLMoms, numOfQLMoms)

                for (n in 1..numOfQLMoms) {
                    val A = Qspp.add(1.0, Psis.mult(Qsmp))
                    val B = Qsmm.add(1.0, Qsmp.mult(Psis))
                    var C = QLDPn[n - 1]!!.mult(D[k + 1])
                    C.scaleEq(n.toDouble())

                    var bino = 1
                    for (i in 1 until n) {
                        bino = bino * (n - i + 1) / i
                        C = C.add(bino.toDouble(), QLDPn[i]!!.mult(Qsmp).mult(QLDPn[n - i]!!))
                    }

                    val P = Matrix.lyap(A, B, C, null)
                    QLDPn[n] = P
                    dqlMoms[n - 1] = inis.mult(P).elementSum()
                }
                val dqlMomsConverted = MomsFromFactorialMoms(dqlMoms)

                // Now calculate at random time instance
                val pi = ctmc_solve(sD)
                val lambdak = pi.mult(D[k + 1]).elementSum()
                val QLPn: MutableMap<Int, Matrix> = HashMap()
                QLPn[0] = pi
                val qlMoms = Matrix(1, numOfQLMoms, numOfQLMoms)

                val iTerm = Matrix.ones(N, 1).mult(pi).add(-1.0, sD).pinv()

                for (n in 1..numOfQLMoms) {
                    val sumP = inis.mult(QLDPn[n]!!).elementSum() +
                            n * (inis.mult(QLDPn[n - 1]!!).add(-1.0 / lambdak, QLPn[n - 1]!!.mult(D[k + 1])))
                                .mult(iTerm).mult(D[k + 1].sumRows())[0]
                    val P = Matrix.scaleMult(pi, sumP).add(
                        n.toDouble(),
                        QLPn[n - 1]!!.mult(D[k + 1]).add(-lambdak, inis.mult(QLDPn[n - 1]!!)).mult(iTerm)
                    )
                    QLPn[n] = P
                    qlMoms[n - 1] = P.elementSum()
                }
                result["ncMoms"]!![k] = MomsFromFactorialMoms(qlMoms)
            }

            // DISTRIBUTION OF THE NUMBER OF JOBS
            if (numOfQLProbs != null) {
                var sDk = D0.copy()
                for (i in 0 until k) {
                    sDk = sDk.add(1.0, D[i + 1])
                }

                // First calculate at departure instants
                val Psid = FluidFundamentalMatrices(Qspp, Qspm, Qsmp, sDk, precision, null, null)["P"]!!
                val Pn: MutableMap<Int, Matrix> = HashMap()
                Pn[0] = Psid
                var dqlProbs = inis.mult(Psid)

                for (n in 1 until numOfQLProbs) {
                    val A = Qspp.add(1.0, Psid.mult(Qsmp))
                    val B = sDk.add(1.0, Qsmp.mult(Psid))
                    var C = Pn[n - 1]!!.mult(D[k + 1])

                    for (i in 1 until n) {
                        C = C.add(1.0, Pn[i]!!.mult(Qsmp).mult(Pn[n - i]!!))
                    }

                    val P = Matrix.lyap(A, B, C, null)
                    Pn[n] = P
                    dqlProbs = Matrix.concatRows(dqlProbs, inis.mult(P), null)
                }

                // Now calculate at random time instance
                val pi = ctmc_solve(sD)
                val lambdak = pi.mult(D[k + 1]).elementSum()
                val iTerm = Matrix.negative(sD.add(-1.0, D[k + 1])).pinv()

                var qlProbs = Matrix.scaleMult(Matrix.extractRows(dqlProbs, 0, 1, null).mult(iTerm), lambdak)
                for (n in 1 until numOfQLProbs) {
                    val P = Matrix.extractRows(qlProbs, n - 1, n, null).mult(D[k + 1])
                        .add(lambdak, Matrix.extractRows(dqlProbs, n, n + 1, null)
                            .add(-1.0, Matrix.extractRows(dqlProbs, n - 1, n, null)))
                        .mult(iTerm)
                    qlProbs = Matrix.concatRows(qlProbs, P, null)
                }
                result["ncDistr"]!![k] = qlProbs.sumRows().transpose()
            }
        } else if (k == K - 1) {
            // ==========================================
            // Step 3. Performance measures for highest priority class (k == K-1)
            // ==========================================

            // MOMENTS OF THE SOJOURN TIME
            if (numOfSTMoms != null) {
                val neg_Kw_inv = neg_Kw.pinv()
                val rtMoms = Matrix(1, numOfSTMoms, numOfSTMoms)
                for (i in 1..numOfSTMoms) {
                    rtMoms[i - 1] = CombinatoricsUtils.factorial(i) *
                            kappa.mult(Matrix.pow(neg_Kw_inv, i + 1)).mult(Bw.sumRows())[0]
                }
                result["stMoms"]!![k] = rtMoms
            }

            // DISTRIBUTION OF THE SOJOURN TIME
            if (stCdfPoints != null) {
                val neg_Kw_inv = neg_Kw.pinv()
                val rtDistr = Matrix(1, stCdfPoints.length(), stCdfPoints.length())
                for (o in 0 until stCdfPoints.length()) {
                    val t = stCdfPoints[o]
                    val Kwt = Kw.copy()
                    Kwt.scaleEq(t)
                    rtDistr[o] = kappa.mult(neg_Kw_inv)
                        .mult(Matrix.eye(Kw.numRows).add(-1.0, Maths.matrixExp(Kwt)))
                        .mult(Bw.sumRows())[0]
                }
                result["stDistr"]!![k] = rtDistr
            }

            // MOMENTS AND DISTRIBUTION OF THE NUMBER OF JOBS (using QBD approach)
            if (numOfQLMoms != null || numOfQLProbs != null) {
                // Build QBD matrices
                val Mk = M[k].toInt()
                val qbdSize = N * Mk

                // L = kron(sD - D{k+1}, eye(M(k))) + kron(eye(N), S{k})
                val L = sD.add(-1.0, D[k + 1]).kron(Matrix.eye(Mk)).add(1.0, I.kron(S[k]))
                // B = kron(eye(N), s{k} * sigma{k})
                val B = I.kron(s[k]!!.mult(sigma[k]))
                // F = kron(D{k+1}, eye(M(k)))
                val F = D[k + 1].kron(Matrix.eye(Mk))
                // L0 = kron(sD - D{k+1}, eye(M(k)))
                val L0 = sD.add(-1.0, D[k + 1]).kron(Matrix.eye(Mk))

                val R = QBDFundamentalMatrices(B, L, F, precision, null, null, null)["R"]!!

                // p0 = CTMCSolve(L0 + R*B)
                var p0 = ctmc_solve(L0.add(1.0, R.mult(B)))
                // Normalize: p0 = p0 / sum(p0 * inv(I - R))
                val IminusR_inv = Matrix.eye(R.numRows).add(-1.0, R).pinv()
                p0 = Matrix.scaleMult(p0, 1.0 / p0.mult(IminusR_inv).elementSum())

                if (numOfQLMoms != null) {
                    val qlMoms = Matrix(1, numOfQLMoms, numOfQLMoms)
                    for (i in 1..numOfQLMoms) {
                        qlMoms[i - 1] = CombinatoricsUtils.factorial(i) *
                                p0.mult(Matrix.pow(R, i)).mult(Matrix.pow(IminusR_inv, i + 1)).elementSum()
                    }
                    result["ncMoms"]!![k] = MomsFromFactorialMoms(qlMoms)
                }

                if (numOfQLProbs != null) {
                    var qlProbs = p0.copy()
                    for (i in 1 until numOfQLProbs) {
                        qlProbs = Matrix.concatRows(qlProbs, p0.mult(Matrix.pow(R, i)), null)
                    }
                    result["ncDistr"]!![k] = qlProbs.sumRows().transpose()
                }
            }
        }
    }

    return result
}
