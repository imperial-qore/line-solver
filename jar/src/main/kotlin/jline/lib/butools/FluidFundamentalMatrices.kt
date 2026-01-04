package jline.lib.butools

import jline.util.matrix.Matrix
import org.apache.commons.math3.util.FastMath

fun FluidFundamentalMatrices(Fpp: Matrix,
                             Fpm: Matrix,
                             Fmp: Matrix,
                             Fmm: Matrix,
                             precision_: Double?,
                             maxNumIt_: Int?,
                             method_: String?): Map<String, Matrix> {
    var Psi = Matrix.singleton(0.0)
    val precision = precision_ ?: 1e-14
    val maxNumIt = maxNumIt_ ?: 150
    val method = method_ ?: "ADDA"
    var numit = 0
    if (Fpp.numRows == 0) {
        Psi = Matrix(0, Fmm.numRows)
    } else if (method == "CR") {
    } else if (method == "ADDA" || method == "SDA") {
        var A = Fpp.copy()
        A.scaleEq(-1.0)
        val B = Fpm.copy()
        val C = Fmp.copy()
        var D = Fmm.copy()
        D.scaleEq(-1.0)
        val diag_A = Matrix.singleton(0.0)
        Matrix.extractDiag(A, diag_A)
        var gamma1 = diag_A.elementMax()
        val diag_D = Matrix.singleton(0.0)
        Matrix.extractDiag(D, diag_D)
        var gamma2 = diag_D.elementMax()
        if (method == "SDA") {
            gamma1 = FastMath.max(gamma1, gamma2)
            gamma2 = gamma1
        }
        val sA = A.numRows
        val sD = D.numRows
        val IA = Matrix.eye(sA)
        val ID = Matrix.eye(sD)
        val gamma2IA = IA.copy()
        gamma2IA.scaleEq(gamma2)
        val gamma1ID = ID.copy()
        gamma1ID.scaleEq(gamma1)
        A = A.add(1.0, gamma2IA)
        D = D.add(1.0, gamma1ID)
        val Dginv = D.inv()
        var Vginv = D.add(-1.0, C.mult(A.inv()).mult(B)).inv()
        var Wginv = A.add(-1.0, B.mult(Dginv).mult(C)).inv()
        val gammaVginv = Vginv.copy()
        gammaVginv.scaleEq(gamma1 + gamma2)
        var Eg = ID.add(-1.0, gammaVginv)
        val gammaWginV = Wginv.copy()
        gammaWginV.scaleEq(gamma1 + gamma2)
        var Fg = IA.add(-1.0, gammaWginV)
        var Gg = Dginv.mult(C).mult(Wginv)
        Gg.scaleEq(gamma1 + gamma2)
        var Hg = Wginv.mult(B).mult(Dginv)
        Hg.scaleEq(gamma1 + gamma2)

        var diff = 1.0
        while (diff > precision && numit < maxNumIt) {
            Vginv = Eg.mult(ID.add(-1.0, Gg.mult(Hg)).inv())
            Wginv = Fg.mult(IA.add(-1.0, Hg.mult(Gg)).inv())
            Gg = Gg.add(1.0, Vginv.mult(Gg).mult(Fg))
            Hg = Hg.add(1.0, Wginv.mult(Hg).mult(Eg))
            Eg = Vginv.mult(Eg)
            Fg = Wginv.mult(Fg)
            val neg = Matrix.firstNorm(Eg)
            val nfg = Matrix.firstNorm(Fg)
            if (method == "ADDA") {
                val eta = FastMath.sqrt(nfg / neg)
                Eg.scaleEq(eta)
                Fg.scaleEq(1 / eta)
                diff = neg * nfg
            } else {
                diff = FastMath.min(neg, nfg)
            }
            numit++
        }
        Psi = Hg
    }


    val result: MutableMap<String, Matrix> = HashMap()
    result["P"] = Psi
    result["K"] = Fpp.add(1.0, Psi.mult(Fmp))
    result["U"] = Fmm.add(1.0, Fmp.mult(Psi))

    return result
}