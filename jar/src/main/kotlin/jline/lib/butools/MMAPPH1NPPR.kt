package jline.lib.butools

import jline.api.mc.ctmc_solve
import jline.util.Maths
import jline.util.matrix.Matrix
import jline.util.matrix.MatrixCell
import org.apache.commons.math3.util.CombinatoricsUtils
import org.apache.commons.math3.util.FastMath
import kotlin.math.pow

fun MMAPPH1NPPR(D: MatrixCell,
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
    var erlMaxOrder = 200
    if (erlMaxOrder_ != null) {
        erlMaxOrder = erlMaxOrder_
    }
    var precision = 1e-14
    if (prec != null) {
        precision = prec
    }
    var classes = Matrix(1, K, K)
    for (i in 0..<K) {
        classes[i] = i.toDouble()
    }
    if (classes_ != null) {
        classes = classes_
    }

    val D0 = D[0]
    val N = D0.numRows
    val I = Matrix.eye(N)
    var sD = Matrix(N, N, N * N)

    for (i in 0..<K + 1) {
        sD = sD.add(1.0, D[i])
    }

    val M = Matrix(1, K, K)
    val s: MutableMap<Int, Matrix> = HashMap()
    for (i in 0..<K) {
        val neg_si = S[i].copy()
        neg_si.scaleEq(-1.0)
        s[i] = neg_si.sumCols()
        M[i] = M[i] + sigma[i].length()
    }
    val QWMM = D0.copy()
    val QWPP = Matrix(N * M.elementSum().toInt(),
        N * M.elementSum().toInt(),
        N * M.elementSum().toInt() * N * M.elementSum().toInt())
    val QWMP = Matrix(N, N * M.elementSum().toInt(), N * N * M.elementSum().toInt())
    val QWPM = Matrix(N * M.elementSum().toInt(), N, N * N * M.elementSum().toInt())
    var kix = 0
    for (i in 0..<K) {
        val bs = N * M[i].toInt()
        QWPP.insertSubMatrix(kix, kix, kix + bs, kix + bs, Matrix.eye(N).kron(S[i]))
        QWMP.insertSubMatrix(0, kix, QWMP.numRows, kix + bs, D[i + 1].kron(sigma[i]))
        QWPM.insertSubMatrix(kix, 0, kix + bs, QWPM.numCols, Matrix.eye(N).kron(s[i]!!))
        kix = kix + bs
    }

    val FluidFundamentalMatrices = FluidFundamentalMatrices(QWPP, QWPM, QWMP, QWMM, precision, null, null)

    var Kw = FluidFundamentalMatrices["K"]
    val Uw = FluidFundamentalMatrices["U"]
    val neg_Kw = Kw!!.copy()
    neg_Kw.scaleEq(-1.0)

    val Ua = Matrix.ones(N, 1).add(2.0, QWMP.mult(neg_Kw.pinv()).sumRows())
    // Solve overdetermined system using least-squares (MATLAB linsolve equivalent)
    val A_ls = Matrix.concatColumns(Uw!!, Ua, null).transpose()
    val b_ls = Matrix.concatColumns(Matrix(1, N, 0), Matrix.singleton(1.0), null).transpose()
    var pm = A_ls.leftMatrixDivide(b_ls).transpose()

    val ro = ((1 - pm.elementSum()) / 2) / (pm.elementSum() + (1 - pm.elementSum()) / 2)
    val kappa = pm.copy()
    kappa.scaleEq(1 / pm.elementSum())

    val pi = ctmc_solve(sD)
    val lambda = Matrix(K, 1, 1)
    for (i in 0..<K) {
        lambda[i] = pi.mult(D[i + 1]).elementSum()
    }

    val Psiw: MutableMap<Int, Matrix?> = HashMap()
    val Qwmp: MutableMap<Int, Matrix> = HashMap()
    val Qwzp: MutableMap<Int, Matrix> = HashMap()
    val Qwpp: MutableMap<Int, Matrix> = HashMap()
    val Qwmz: MutableMap<Int, Matrix> = HashMap()
    val Qwpz: MutableMap<Int, Matrix> = HashMap()
    val Qwzz: MutableMap<Int, Matrix> = HashMap()
    val Qwmm: MutableMap<Int, Matrix> = HashMap()
    val Qwpm: MutableMap<Int, Matrix> = HashMap()
    val Qwzm: MutableMap<Int, Matrix> = HashMap()

    for (k in 0..<K) {
        val Mlo = if (k == 0) {
            0.0
        } else {
            M.sumSubMatrix(0, 1, 0, k)
        }

        val Mhi = M.elementSum() - Mlo
        val a = (N * Mlo * Mhi + N * Mhi).toInt()
        val Qkwpp = Matrix(a, a, a * a)
        val Qkwpz = Matrix(a, N * Mlo.toInt(), a * N * Mlo.toInt())
        val Qkwpm = Matrix(a, N, a * N)
        val Qkwmz = Matrix(N, N * Mlo.toInt())
        val Qkwmp = Matrix(N, a, a * N)
        var Dlo = D0.copy()
        for (i in 0..<k) {
            Dlo = Dlo.add(1.0, D[i + 1])
        }
        val Qkwmm = Dlo
        val Qkwzp = Matrix(N * Mlo.toInt(), a, N * Mlo.toInt() * a)
        val Qkwzm = Matrix(N * Mlo.toInt(), N, N * Mlo.toInt() * N)
        val Qkwzz = Matrix(N * Mlo.toInt(), N * Mlo.toInt(), N * Mlo.toInt() * N * Mlo.toInt())
        kix = 0
        for (i in k..<K) {
            var kix2 = 0
            for (j in 0..<k) {
                val bs = (N * M[j] * M[i]).toInt()
                val bs2 = N * M[j].toInt()
                Qkwpp.insertSubMatrix(kix,
                    kix,
                    kix + bs,
                    kix + bs,
                    Matrix.eye(N).kron(Matrix.eye(M[j].toInt()).kron(S[i])))
                Qkwpz.insertSubMatrix(kix,
                    kix2,
                    kix + bs,
                    kix2 + bs2,
                    Matrix.eye(N).kron(Matrix.eye(M[j].toInt()).kron(s[i]!!)))
                Qkwzp.insertSubMatrix(kix2,
                    kix,
                    kix2 + bs2,
                    kix + bs,
                    D[i + 1].kron(Matrix.eye(M[j].toInt()).kron(sigma[i])))
                kix = kix + bs
                kix2 = kix2 + bs2
            }
        }
        for (i in k..<K) {
            val bs = N * M[i].toInt()
            Qkwpp.insertSubMatrix(kix, kix, kix + bs, kix + bs, Matrix.eye(N).kron(S[i]))
            Qkwpm.insertSubMatrix(kix, 0, kix + bs, Qkwpm.numCols, Matrix.eye(N).kron(s[i]!!))
            Qkwmp.insertSubMatrix(0, kix, Qkwmp.numRows, kix + bs, D[i + 1].kron(sigma[i]))
            kix = kix + bs
        }
        kix = 0
        for (j in 0..<k) {
            val bs = N * M[j].toInt()
            Qkwzz.insertSubMatrix(kix,
                kix,
                kix + bs,
                kix + bs,
                Dlo.kron(Matrix.eye(M[j].toInt())).add(1.0, Matrix.eye(N).kron(S[j])))
            Qkwzm.insertSubMatrix(kix, 0, kix + bs, Qkwzm.numCols, Matrix.eye(N).kron(s[j]!!))
            kix = kix + bs
        }
        // When Mlo = 0, Qkwzz is empty (0x0), so skip the z-related terms
        val Psikw = if (Mlo > 0) {
            val neg_Qkwzz = Qkwzz.copy()
            neg_Qkwzz.scaleEq(-1.0)
            FluidFundamentalMatrices(Qkwpp.add(1.0, Qkwpz.mult(neg_Qkwzz.pinv()).mult(Qkwzp)),
                Qkwpm.add(1.0, Qkwpz.mult(neg_Qkwzz.pinv()).mult(Qkwzm)),
                Qkwmp,
                Qkwmm,
                precision,
                null,
                null)["P"]
        } else {
            // No lower priority classes, so Qkwzz-related terms are zero
            FluidFundamentalMatrices(Qkwpp, Qkwpm, Qkwmp, Qkwmm, precision, null, null)["P"]
        }
        Psiw[k] = Psikw
        Qwzp[k] = Qkwzp
        Qwmp[k] = Qkwmp
        Qwpp[k] = Qkwpp
        Qwmz[k] = Qkwmz
        Qwpz[k] = Qkwpz
        Qwzz[k] = Qkwzz
        Qwmm[k] = Qkwmm
        Qwpm[k] = Qkwpm
        Qwzm[k] = Qkwzm
    }

    val lambdaS = lambda.elementSum()
    val phi: MutableMap<Int, Matrix> = HashMap()
    val neg_D0 = D0.copy()
    neg_D0.scaleEq(-1.0)
    val phi0 = kappa.mult(neg_D0)
    phi0.scaleEq((1 - ro) / lambdaS)
    phi[0] = phi0

    val q0: MutableMap<Int, Matrix> = HashMap()
    val qL: MutableMap<Int, Matrix> = HashMap()
    q0[0] = Matrix(0, 0, 0)
    qL[0] = Matrix(0, 0, 0)

    for (k in 0..<K - 1) {
        var sDk = D[0]
        for (j in 0..k) {
            sDk = sDk.add(1.0, D[j + 1])
        }
        var pk = 0.0
        for (j in 0..k) {
            pk = pk + lambda[j]
        }
        pk = pk / lambdaS - (1 - ro) * kappa.mult(sDk.sumRows())[0] / lambdaS
        val Qwzpk = Qwzp[k + 1]
        var vix = 0
        val Ak: MutableMap<Int, Matrix> = HashMap()
        for (ii in 0..k) {
            val bs = (N * M[ii]).toInt()
            val V1 = Matrix.extractRows(Qwzpk!!, vix, vix + bs, null)
            val a = sDk.kron(Matrix.eye(M[ii].toInt()))
            a.scaleEq(-1.0)
            a.add(-1.0, I.kron(S[ii]))
            Ak[ii] = I.kron(sigma[ii]).mult(a.pinv()).mult(I.kron(s[ii]!!)).add(1.0, V1.mult(Psiw[k + 1]!!))
            vix = vix + bs
        }
        val Qwmpk = Qwmp[k + 1]
        val Bk = Qwmpk!!.mult(Psiw[k + 1]!!)
        // MATLAB: ztag = phi{1}*(inv(-D0)*D{k+1}*Ak{k} - Ak{1} + inv(-D0)*Bk)
        var ztag =
            phi[0]!!.mult(neg_D0.pinv().mult(D[k + 1]).mult(Ak[k]!!).add(-1.0, Ak[0]!!).add(1.0, neg_D0.pinv().mult(Bk)))
        // MATLAB: for i=1:k-1, ztag = ztag + phi{i+1}*(Ak{i}-Ak{i+1}) + phi{1}*inv(-D0)*D{i+1}*Ak{i}
        for (i in 0..<k) {
            ztag = ztag.add(1.0, phi[i + 1]!!.mult(Ak[i]!!.add(-1.0, Ak[i + 1]!!)))
                .add(1.0, phi[0]!!.mult(neg_D0.pinv()).mult(D[i + 1]).mult(Ak[i]!!))
        }
        val Mx = Matrix.eye(Ak[k]!!.numCols).add(-1.0, Ak[k]!!)
        Mx.insertSubMatrix(0, 0, Mx.numRows, 1, Matrix.ones(N, 1))
        phi[k + 1] = Matrix.concatColumns(Matrix.singleton(pk),
            Matrix.extractRows(ztag.transpose(), 1, ztag.length(), null).transpose(),
            null).mult(Mx.pinv())
        q0[k + 1] = phi[0]!!.mult(neg_D0.pinv())
        qL[k + 1] = Matrix(0, 0, 0)
        // MATLAB: for ii=1:k (processes classes 1 to k, which in 0-based is 0 to k)
        for (ii in 0..k) {
            val a = sDk.kron(Matrix.eye(M[ii].toInt()))
            a.scaleEq(-1.0)
            val qLii = phi[ii + 1]!!.add(-1.0, phi[ii]!!).add(1.0, phi[0]!!.mult(neg_D0.pinv()).mult(D[ii + 1]))
                .mult(I.kron(sigma[ii])).mult(a.add(-1.0, I.kron(S[ii])).pinv())
            qL[k + 1] = Matrix.concatColumns(qL[k + 1]!!, qLii, null)
        }
    }

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
    for (g in 0..<classes.length()) {
        val k = classes[g].toInt()
        var sD0k = D0
        for (i in 0..k - 1) {
            sD0k = sD0k.add(1.0, D[i + 1])
        }
        // MATLAB: if k<K (1-based), which is k=1,...,K-1, maps to 0-based k=0,...,K-2 -> k < K-1
        if (k < K - 1) {
            // Check if there are lower priority classes (Mlo > 0 means non-empty z-related matrices)
            val Mlo_k = if (k == 0) 0.0 else M.sumSubMatrix(0, 1, 0, k)
            Kw = if (Mlo_k > 0) {
                val neg_Qwzz_k = Qwzz[k]!!.copy()
                neg_Qwzz_k.scaleEq(-1.0)
                Qwpp[k]!!.add(1.0, Qwpz[k]!!.mult(neg_Qwzz_k.pinv()).mult(Qwzp[k]!!)).add(1.0, Psiw[k]!!.mult(Qwmp[k]!!))
            } else {
                // No lower priority classes, skip z-related terms
                Qwpp[k]!!.add(1.0, Psiw[k]!!.mult(Qwmp[k]!!))
            }
            // MATLAB: BM = []; CM = []; DM = []; for i=1:k-1, ... end
            var BM = Matrix(0, 0, 0)
            var CM = Matrix(0, N, 0)  // Empty matrix with N columns for vertical concat
            var DM = Matrix(0, 0, 0)
            for (i in 0..k - 1) {
                BM = BM.createBlockDiagonal(I.kron(S[i]))
                DM = DM.createBlockDiagonal(D[k + 1].kron(Matrix.eye(M[i].toInt())))
                CM = Matrix.concatRows(CM, I.kron(s[i]!!), null)
            }
            // Build Kwu, handling case where BM/DM are empty (k=0)
            val Kwu = if (Mlo_k > 0) {
                Matrix.concatRows(Matrix.concatColumns(Kw,
                    Qwpz[k]!!.add(1.0, Psiw[k]!!.mult(Qwmz[k]!!)).mult(Matrix.negative(Qwzz[k]!!).pinv()).mult(DM),
                    null), Matrix.concatColumns(Matrix(BM.numRows, Kw.numRows, 0), BM, null), null)
            } else {
                // No lower priority classes, Kwu is just Kw
                Kw
            }
            val Bwu = if (Mlo_k > 0) {
                Matrix.concatRows(Psiw[k]!!.mult(D[k + 1]), CM, null)
            } else {
                Psiw[k]!!.mult(D[k + 1])
            }
            val iniw: Matrix
            val pwu: Matrix
            if (k > 0) {
                iniw = Matrix.concatColumns(q0[k]!!.mult(Qwmp[k]!!).add(1.0, qL[k]!!.mult(Qwzp[k]!!)),
                    qL[k]!!.mult(DM),
                    null)
                pwu = q0[k]!!.mult(D[k + 1])
            } else {
                iniw = pm.mult(Qwmp[k]!!)
                pwu = pm.mult(D[k + 1])
            }
            val norm = pwu.elementSum() + iniw.mult(Matrix.negative(Kwu).pinv()).mult(Bwu).elementSum()
            pwu.scaleEq(1 / norm)
            iniw.scaleEq(1 / norm)
            val KN = Kwu.numRows
            // Sum M from column k+1 to K-1 (higher priority classes' phases)
            val sizeHigh = M.sumSubMatrix(0, 1, k + 1, K)
            val Qspp = Matrix((KN + N * sizeHigh).toInt(),
                (KN + N * sizeHigh).toInt(),
                (KN + N * sizeHigh).pow(2.0).toInt())
            val Qspm = Matrix((KN + N * sizeHigh).toInt(),
                N,
                N * (KN + N * sizeHigh).toInt())
            val Qsmp = Matrix(N,
                (KN + N * sizeHigh).toInt(),
                N * (KN + N * sizeHigh).toInt())
            val Qsmm = sD0k.add(1.0, D[k + 1])
            kix = 0
            for (i in k + 1..<K) {
                val bs = N * M[i].toInt()
                Qspp.insertSubMatrix(KN + kix, KN + kix, KN + kix + bs, KN + kix + bs, I.kron(S[i]))
                Qspm.insertSubMatrix(KN + kix, 0, KN + kix + bs, Qspm.numCols, I.kron(s[i]!!))
                Qsmp.insertSubMatrix(0, KN + kix, Qsmp.numRows, KN + kix + bs, D[i + 1].kron(sigma[i]))
                kix = kix + bs
            }
            Qspp.insertSubMatrix(0, 0, KN, KN, Kwu)
            Qspm.insertSubMatrix(0, 0, KN, Qspm.numCols, Bwu)
            val inis = Matrix.concatColumns(iniw, Matrix(1, (N * sizeHigh).toInt()), null)

            val Psis = FluidFundamentalMatrices(Qspp, Qspm, Qsmp, Qsmm, precision, null, null)["P"]

            if (numOfSTMoms != null) {
                // MATLAB: Pn = {Psis}; (1-based, so Pn{1}=Psis)
                val Pn: MutableMap<Int, Matrix?> = HashMap()
                Pn[0] = Psis  // Pn[0] corresponds to MATLAB Pn{1}
                val wtMoms = Matrix(1, numOfSTMoms, numOfSTMoms)
                // MATLAB: for n=1:numOfSTMoms
                for (n in 1..numOfSTMoms) {
                    val A = Qspp.add(1.0, Psis!!.mult(Qsmp))
                    val B = Qsmm.add(1.0, Qsmp.mult(Psis))
                    // MATLAB: C = -2*n*Pn{n}; MUST copy to avoid modifying Pn!
                    var C = Pn[n - 1]!!.copy()
                    C.scaleEq((-2 * n).toDouble())
                    var bino = 1
                    // MATLAB: for i=1:n-1
                    for (i in 1..n - 1) {
                        bino = bino * (n - i + 1) / i
                        // MATLAB: C = C + bino * Pn{i+1}*Qsmp*Pn{n-i+1}
                        // In 0-based: Pn[i]*Qsmp*Pn[n-i]
                        C = C.add(bino.toDouble(), Pn[i]!!.mult(Qsmp).mult(Pn[n - i]!!))
                    }
                    val P = Matrix.lyap(A, B, C, null)
                    Pn[n] = P
                    // MATLAB: wtMoms(n) = sum(inis*P*(-1)^n) / 2^n
                    wtMoms[n - 1] = inis.mult(P).elementSum() * (-1.0).pow(n.toDouble()) / 2.0.pow(n.toDouble())
                }
                // MATLAB: calculate RESPONSE time moments
                val Pnr: MutableMap<Int, Matrix?> = HashMap()
                // MATLAB: Pnr = {sum(inis*Pn{1})*sigma{k}}
                Pnr[0] = Matrix.scaleMult(sigma[k], inis.mult(Pn[0]!!).elementSum())
                val rtMoms = Matrix(1, numOfSTMoms, numOfSTMoms)
                val negSk = S[k].copy()
                negSk.scaleEq(-1.0)
                for (n in 1..numOfSTMoms) {
                    // MATLAB: P = n*Pnr{n}*inv(-S{k}) + (-1)^n*sum(inis*Pn{n+1})*sigma{k} / 2^n
                    val term1 = Pnr[n - 1]!!.mult(negSk.pinv())
                    term1.scaleEq(n.toDouble())
                    val coeff = (-1.0).pow(n.toDouble()) / 2.0.pow(n.toDouble()) * inis.mult(Pn[n]!!).elementSum()
                    val P = term1.add(coeff, sigma[k])
                    Pnr[n] = P
                    // MATLAB: rtMoms(n) = sum(P)+sum(pwu)*factorial(n)*sum(sigma{k}*inv(-S{k})^n)
                    rtMoms[n - 1] = P.elementSum() + pwu.elementSum() * CombinatoricsUtils.factorial(n) *
                            sigma[k].mult(Matrix.pow(negSk.pinv(), n)).elementSum()
                }
                result["stMoms"]!![k] = rtMoms
            }

            if (stCdfPoints != null) {
                var res = Matrix(1, 0, 0)
                for (o in 0..<stCdfPoints.length()) {
                    val t = stCdfPoints[o].toInt()
                    val L = erlMaxOrder
                    val lambdae = L / t.toDouble() / 2.0
                    val Psie = FluidFundamentalMatrices(Qspp.add(-lambdae, Matrix.eye(Qspp.numRows)),
                        Qspm,
                        Qsmp,
                        Qsmm.add(-lambdae, Matrix.eye(Qsmm.numRows)),
                        precision,
                        null,
                        null)["P"]
                    val Pn: MutableMap<Int, Matrix?> = HashMap()
                    Pn[0] = Psis
                    var pr = (pwu.elementSum() + inis.mult(Psie!!)
                        .elementSum()) * (1 - sigma[k].mult(Matrix.pow(Matrix.eye(S[k].numRows).pinv()
                        .add(-1 / 2.0 / lambdae, S[k]), L)).elementSum())
                    for (n in 0..L - 1) {
                        val A = Qspp.add(1.0, Psie.mult(Qsmp)).add(-lambdae, Matrix.eye(Qspp.length()))
                        val B = Qsmm.add(1.0, Qsmp.mult(Psie)).add(-lambdae, Matrix.eye(Qsmm.length()))
                        var C = Pn[n]
                        C!!.scaleEq(2 * lambdae)
                        for (i in 0..<n - 1) {
                            C = C!!.add(1.0, Pn[i + 1]!!.mult(Qsmp).mult(Pn[n - i + 1]!!))
                        }
                        val P = Matrix.lyap(A, B, C!!, null)
                        Pn[n + 1] = P
                        pr = pr + inis.mult(P).elementSum() * (1 - sigma[k].mult(Matrix.pow(Matrix.eye(S[k].length())
                            .pinv().add(-1 / 2.0 / lambdae, S[k]), L - n)).elementSum())
                    }
                    res = Matrix.concatColumns(res, Matrix.singleton(pr), null)
                }
                result["stDistr"]!![k] = res
            }
            if (numOfQLMoms != null || numOfQLProbs != null) {
                val W =
                    Matrix.negative(sD.add(-1.0, D[k + 1]).kron(Matrix.eye(M[k].toInt()))).add(-1.0, I.kron(S[k])).pinv()
                        .mult(D[k + 1].kron(Matrix.eye(M[k].toInt())))
                val iW = Matrix.eye(W.numCols).add(-1.0, W).pinv()
                val w = Matrix.eye(N).kron(sigma[k])
                val omega =
                    Matrix.negative(sD.add(-1.0, D[k + 1]).kron(Matrix.eye(M[k].toInt()))).add(-1.0, I.kron(S[k])).pinv()
                        .mult(I.kron(s[k]!!))
                if (numOfQLMoms != null) {
                    // MATLAB: Psii = {Psis}; QLDPn = {inis*Psii{1}*w*iW}; (1-based)
                    val Psii: MutableMap<Int, Matrix?> = HashMap()
                    Psii[0] = Psis  // Psii[0] = MATLAB Psii{1}
                    val QLDPn: MutableMap<Int, Matrix> = HashMap()
                    QLDPn[0] = inis.mult(Psii[0]!!).mult(w).mult(iW)  // QLDPn[0] = MATLAB QLDPn{1}
                    // MATLAB: for n=1:numOfQLMoms
                    for (n in 1..numOfQLMoms) {
                        val A = Qspp.add(1.0, Psis!!.mult(Qsmp))
                        val B = Qsmm.add(1.0, Qsmp.mult(Psis))
                        // MATLAB: C = n*Psii{n}*D{k+1}
                        var C = Psii[n - 1]!!.mult(D[k + 1])
                        C.scaleEq(n.toDouble())
                        var bino = 1
                        // MATLAB: for i=1:n-1
                        for (i in 1..n - 1) {
                            bino = bino * (n - i + 1) / i
                            // MATLAB: C = C + bino * Psii{i+1}*Qsmp*Psii{n-i+1}
                            // 0-based: Psii[i]*Qsmp*Psii[n-i]
                            C = C.add(bino.toDouble(), Psii[i]!!.mult(Qsmp).mult(Psii[n - i]!!))
                        }
                        val P = Matrix.lyap(A, B, C, null)
                        Psii[n] = P  // MATLAB: Psii{n+1} = P
                        // MATLAB: QLDPn{n+1} = n*QLDPn{n}*iW*W + inis*P*w*iW
                        val qlDTerm = QLDPn[n - 1]!!.mult(iW).mult(W)
                        qlDTerm.scaleEq(n.toDouble())
                        QLDPn[n] = qlDTerm.add(1.0, inis.mult(P).mult(w).mult(iW))
                    }
                    // MATLAB: for n=0:numOfQLMoms
                    for (n in 0..numOfQLMoms) {
                        // MATLAB: QLDPn{n+1} = (QLDPn{n+1} + pwu*w*iW^(n+1)*W^n)*omega
                        QLDPn[n] = QLDPn[n]!!.add(1.0,
                                pwu.mult(w).mult(Matrix.pow(iW, n + 1)).mult(Matrix.pow(W, n))).mult(omega)
                    }
                    // MATLAB: QLPn = {pi}
                    val QLPn: MutableMap<Int, Matrix> = HashMap()
                    QLPn[0] = pi  // QLPn[0] = MATLAB QLPn{1}
                    var qlMOms = Matrix(1, numOfQLMoms, numOfQLMoms)
                    val iTerm = Matrix.ones(N, 1).mult(pi).add(-1.0, sD).pinv()
                    // MATLAB: for n=1:numOfQLMoms
                    for (n in 1..numOfQLMoms) {
                        // MATLAB: sumP = sum(QLDPn{n+1}) + n*(QLDPn{n} - QLPn{n}*D{k+1}/lambda(k))*iTerm*sum(D{k+1},2)
                        val sumP =
                            QLDPn[n]!!.elementSum() + n * (QLDPn[n - 1]!!.add(-1.0 / lambda[k], QLPn[n - 1]!!.mult(D[k + 1]))
                                .mult(iTerm).mult(D[k + 1].sumRows()))[0]
                        // MATLAB: P = sumP*pi + n*(QLPn{n}*D{k+1} - QLDPn{n}*lambda(k))*iTerm
                        val P = Matrix.scaleMult(pi, sumP)
                            .add(n.toDouble(), QLPn[n - 1]!!.mult(D[k + 1]).add(-lambda[k], QLDPn[n - 1]!!).mult(iTerm))
                        QLPn[n] = P  // MATLAB: QLPn{n+1} = P
                        qlMOms[n - 1] = P.elementSum()  // MATLAB: qlMoms(n)
                    }
                    qlMOms = MomsFromFactorialMoms(qlMOms)
                    result["ncMoms"]!![k] = qlMOms
                }
                if (numOfQLProbs != null) {
                    // MATLAB: Psid = FluidFundamentalMatrices(...); Pn = {Psid}; (1-based)
                    val Psid =
                        FluidFundamentalMatrices(Qspp, Qspm, Qsmp, sD0k, precision, null, null)["P"]
                    val Pn: MutableMap<Int, Matrix?> = HashMap()
                    Pn[0] = Psid  // Pn[0] = MATLAB Pn{1}
                    var XDn = inis.mult(Psid!!).mult(w)
                    var dqlProbs = XDn.add(1.0, pwu.mult(w)).mult(omega)
                    // MATLAB: for n=1:numOfQLProbs-1
                    for (n in 1..<numOfQLProbs) {
                        val A = Qspp.add(1.0, Psid.mult(Qsmp))
                        val B = sD0k.add(1.0, Qsmp.mult(Psid))
                        // MATLAB: C = Pn{n}*D{k+1}
                        var C = Pn[n - 1]!!.mult(D[k + 1])
                        // MATLAB: for i=1:n-1
                        for (i in 1..n - 1) {
                            // MATLAB: C = C + Pn{i+1}*Qsmp*Pn{n-i+1}
                            // 0-based: Pn[i]*Qsmp*Pn[n-i]
                            C = C.add(1.0, Pn[i]!!.mult(Qsmp).mult(Pn[n - i]!!))
                        }
                        val P = Matrix.lyap(A, B, C, null)
                        Pn[n] = P  // MATLAB: Pn{n+1} = P
                        XDn = XDn.mult(W).add(1.0, inis.mult(P).mult(w))
                        // MATLAB: dqlProbs = [dqlProbs; (XDn+pwu*w*W^n)*omega]
                        dqlProbs = Matrix.concatRows(dqlProbs,
                            XDn.add(1.0, pwu.mult(w).mult(Matrix.pow(W, n))).mult(omega),
                            null)
                    }
                    // MATLAB: iTerm = inv(-(sD-D{k+1}))
                    val iTerm = Matrix.negative(sD.add(-1.0, D[k + 1])).pinv()
                    // MATLAB: qlProbs = lambda(k)*dqlProbs(1,:)*iTerm
                    var qlProbs = Matrix.scaleMult(Matrix.extractRows(dqlProbs, 0, 1, null).mult(iTerm), lambda[k])
                    // MATLAB: for n=1:numOfQLProbs-1
                    for (n in 1..<numOfQLProbs) {
                        // MATLAB: P = (qlProbs(n,:)*D{k+1}+lambda(k)*(dqlProbs(n+1,:)-dqlProbs(n,:)))*iTerm
                        val P = Matrix.extractRows(qlProbs, n - 1, n, null).mult(D[k + 1]).add(1.0,
                            Matrix.scaleMult(Matrix.extractRows(dqlProbs, n, n + 1, null)
                                .add(-1.0, Matrix.extractRows(dqlProbs, n - 1, n, null)), lambda[k])).mult(iTerm)
                        qlProbs = Matrix.concatRows(qlProbs, P, null)
                    }
                    result["ncDistr"]!![k] = qlProbs.sumRows().transpose()
                }
            }
        } else if (k == K - 1) {  // MATLAB: k==K (highest priority class, 0-based)
            if (numOfSTMoms != null || stCdfPoints != null) {
                Kw = Qwpp[k]!!.add(1.0, Qwpz[k]!!.mult(Matrix.negative(Qwzz[k]!!).pinv()).mult(Qwzp[k]!!))
                    .add(1.0, Psiw[k]!!.mult(Qwmp[k]!!))
                var AM = Matrix(0, 0, 0)
                var BM = Matrix(0, 0, 0)
                var CM = Matrix(0, s[0]!!.numCols, 0)
                var DM = Matrix(0, 0, 0)
                for (i in 0..<k) {
                    AM = AM.createBlockDiagonal(Matrix(N, 1).kron(Matrix.eye(M[i].toInt()).kron(s[k]!!)))
                    BM = BM.createBlockDiagonal(S[i])
                    CM = Matrix.concatRows(CM, s[i]!!, null)
                    DM = DM.createBlockDiagonal(D[k + 1].kron(Matrix.eye(M[i].toInt())))
                }
                val Z = Matrix.concatRows(Matrix.concatColumns(Kw,
                    Matrix.concatRows(AM, Matrix(N * M[k].toInt(), AM.numCols), null),
                    null), Matrix.concatColumns(Matrix(BM.numRows, Kw.numCols), BM, null), null)
                val z =
                    Matrix.concatRows(Matrix.concatRows(Matrix(AM.numRows, 1, 0), Matrix.ones(N, 1).kron(s[k]!!), null),
                        CM,
                        null)
                val iniw = Matrix.concatColumns(q0[k]!!.mult(Qwmp[k]!!).add(1.0, qL[k]!!.mult(Qwzp[k]!!)),
                    Matrix(1, BM.numRows),
                    null)
                val zeta = Matrix.scaleMult(iniw, 1 / iniw.mult(Matrix.negative(Z).pinv()).mult(z).elementSum())
                if (numOfSTMoms != null) {
                    // MATLAB: for i=1:numOfSTMoms, rtMomsH(i) = factorial(i)*zeta*inv(-Z)^(i+1)*z
                    val rtMomsH = Matrix(1, numOfSTMoms, numOfSTMoms)
                    for (i in 0..<numOfSTMoms) {
                        // MATLAB i=1,2,... maps to Kotlin i=0,1,...
                        // factorial(i) in MATLAB = factorial(i+1) in Kotlin
                        // inv(-Z)^(i+1) in MATLAB = pow(i+2) in Kotlin
                        rtMomsH[i] =
                            CombinatoricsUtils.factorial(i + 1) * zeta.mult(Matrix.pow(Matrix.negative(Z).pinv(), i + 2))
                                .mult(z).elementSum()
                    }
                    result["stMoms"]!![k] = rtMomsH
                }
                if (stCdfPoints != null) {
                    var rtDistr = zeta.mult(Matrix.negative(Z).pinv())
                        .mult(Matrix.eye(Z.numCols).add(-1.0, Z.mult(Maths.matrixExp(Matrix.scaleMult(Z, stCdfPoints[0])))))
                        .mult(z)
                    for (i in 1..<stCdfPoints.length()) {
                        rtDistr = Matrix.concatColumns(rtDistr,
                            zeta.mult(Matrix.negative(Z).pinv()).mult(Matrix.eye(Z.numCols)
                                .add(-1.0, Z.mult(Maths.matrixExp(Matrix.scaleMult(Z, stCdfPoints[i]))))).mult(z),
                            null)
                    }
                    result["stDistr"]!![k] = rtDistr
                }
            }

            // MATLAB: elseif strcmp(varargin{argIx},'ncMoms') || strcmp(varargin{argIx},'ncDistr')
            if (numOfQLMoms != null || numOfQLProbs != null) {
                val L = Matrix(N * M.elementSum().toInt(),
                    N * M.elementSum().toInt(),
                    FastMath.pow(N * M.elementSum(), 2).toInt())
                val B = Matrix(N * M.elementSum().toInt(),
                    N * M.elementSum().toInt(),
                    FastMath.pow(N * M.elementSum(), 2).toInt())
                val F = Matrix(N * M.elementSum().toInt(),
                    N * M.elementSum().toInt(),
                    FastMath.pow(N * M.elementSum(), 2).toInt())
                kix = 0
                // MATLAB: for i=1:K
                for (i in 0..<K) {
                    val bs = N * M[i].toInt()
                    F.insertSubMatrix(kix, kix, kix + bs, kix + bs, D[k + 1].kron(Matrix.eye(M[i].toInt())))
                    // MATLAB: L(...) = kron(sD0k,eye(M(i))) + kron(I,S{i})
                    L.insertSubMatrix(kix,
                        kix,
                        kix + bs,
                        kix + bs,
                        sD0k.kron(Matrix.eye(M[i].toInt())).add(1.0, I.kron(S[i])))
                    // MATLAB: if i<K (1-based), in 0-based: if i < K-1
                    if (i < K - 1) {
                        // MATLAB: L(...,N*sum(M(1:k-1))+1:end) - in 0-based, +1 becomes +0
                        L.insertSubMatrix(kix,
                            N * M.sumSubMatrix(0, 1, 0, k).toInt(),
                            kix + bs,
                            L.numCols,
                            I.kron(s[i]!!.mult(sigma[k])))
                    } else {
                        B.insertSubMatrix(kix,
                            N * M.sumSubMatrix(0, 1, 0, k).toInt(),
                            kix + bs,
                            B.numCols,
                            I.kron(s[i]!!.mult(sigma[k])))
                    }
                    kix = kix + bs
                }
                val R = QBDFundamentalMatrices(B, L, F, precision, null, null, null)["R"]
                val P0 = Matrix.concatColumns(qL[k]!!, q0[k]!!.mult(I.kron(sigma[k])), null)
                P0.scaleEq(1 / P0.mult(Matrix.eye(R!!.numRows).add(-1.0, R).pinv()).elementSum())

                if (numOfQLMoms != null) {
                    // MATLAB: for i=1:numOfQLMoms, qlMoms(i) = sum(factorial(i)*p0*R^i*inv(eye(size(R))-R)^(i+1))
                    val qlMoms = Matrix(1, numOfQLMoms, numOfQLMoms)
                    for (i in 0..<numOfQLMoms) {
                        // MATLAB i=1,2,... maps to Kotlin i=0,1,...
                        // factorial(i) in MATLAB = factorial(i+1) in Kotlin
                        // R^i in MATLAB = R^(i+1) in Kotlin
                        // inv(...)^(i+1) in MATLAB = inv(...)^(i+2) in Kotlin
                        qlMoms[i] = Matrix.scaleMult(P0.mult(Matrix.pow(R, i + 1))
                            .mult(Matrix.pow(Matrix.eye(R.numRows).add(-1.0, R), i + 2).pinv()),
                            CombinatoricsUtils.factorial(i + 1).toDouble()).elementSum()
                    }
                    result["ncMoms"]!![k] = MomsFromFactorialMoms(qlMoms)
                }

                if (numOfQLProbs != null) {
                    // MATLAB: qlProbs = p0; for i=1:numOfQLProbs-1, qlProbs = [qlProbs; p0*R^i]
                    var qlProbs = P0.copy()
                    for (i in 1..<numOfQLProbs) {
                        // MATLAB i=1,2,... uses R^i, Kotlin i=1,2,... also uses R^i
                        qlProbs = Matrix.concatRows(qlProbs, P0.mult(Matrix.pow(R, i)), null)
                    }
                    result["ncDistr"]!![k] = qlProbs.sumRows()
                }
            }
        }
    }
    return result
}