package jline.lib.smc

import jline.util.matrix.Matrix
import kotlin.math.pow

fun QBD_NI(A0: Matrix,
           A1: Matrix,
           A2: Matrix,
           MaxNumIt_: Int?,
           Verbose_: Int?,
           Mode_: String?,
           RAPComp_: Int?): Map<String?, Matrix?> {
    var A1 = A1
    var Mode = "Sylvest"
    var MaxNumIt = 50
    var Verbose = false
    var RAPComp = false

    if (MaxNumIt_ != null) {
        MaxNumIt = MaxNumIt_
    }

    if (Mode_ != null) {
        if (Mode_ == "Sylvest" || Mode_ == "Estimat" || Mode_ == "DirectSum") {
            Mode = Mode_
        } else {
            throw RuntimeException("QBD_LR mode not recognized")
        }
    }

    if (Verbose_ != null) {
        if (Verbose_ == 1) {
            Verbose = true
        }
    }

    if (RAPComp_ != null) {
        if (RAPComp_ == 1) {
            RAPComp = true
        }
    }
    val A1_diag = Matrix(0, 0, 0)
    Matrix.extractDiag(A1, A1_diag)
    var continues = true
    val lamb = Matrix.negative(A1_diag).elementMax()
    val m = A1.numRows
    if (!RAPComp) {
        continues = false
        if (A1_diag.elementSum() < 0) {
            continues = true
            A0.scaleEq(1 / lamb)
            A1.scaleEq(1 / lamb)
            A1 = A1.add(1.0, Matrix.eye(m))
            A2.scaleEq(1 / lamb)
        }
        QBD_ParsePara(A0, A1, A2)
    } else {
        QBD_ParsePara(A0, A1, A2)
        A0.scaleEq(1 / lamb)
        A1.scaleEq(1 / lamb)
        A1 = A1.add(1.0, Matrix.eye(m))
        A2.scaleEq(1 / lamb)
    }
    var result = QBD_EG(A0, A1, A2, Verbose)
    val G = result["G"]
    if (G!!.length() == 0) {
        return result
    }
    var R = Matrix(m, m, m * m)
    var check = 1.0
    var numit = 0
    while (check > 10.0.pow(-12.0) && numit < MaxNumIt) {
        numit = numit + 1
        val YK: Matrix
        if (numit == 1) {
            YK = A2.mult(Matrix.eye(m).add(-1.0, A1).inv())
        } else {
            if (Mode == "Estimat") {
                val FRK = A2.add(1.0, R.mult(A1.add(-1.0, Matrix.eye(m)).add(1.0, R.mult(A0))))
                val ZK = FRK.mult(Matrix.eye(m).add(-1.0, A1)).inv()
                YK = FRK.add(1.0, FRK.add(1.0, ZK.mult(A1)).add(1.0, R.mult(ZK).add(1.0, ZK.mult(R)).mult(A0)))
            } else {
                val D = Matrix.scaleMult(A2.add(1.0, R.mult(A1.add(-1.0, Matrix.eye(m)).add(1.0, R.mult(A0)))), -1.0)
                val C = A1.add(1.0, R.mult(A0)).add(-1.0, Matrix.eye(m))
                if (Mode == "Sylvest") {
                    YK = QBD_NI_Sylvest(A0.transpose(), R.transpose(), C.transpose(), D.transpose())
                } else {
                    val D_reshape = D.copy()
                    D_reshape.reshape(m * m, 1)
                    YK = A0.transpose().krons(R).add(1.0, C.transpose().krons(Matrix.eye(m)).inv().mult(D_reshape))
                    YK.reshape(m, m)
                }
            }
        }
        R = R.add(1.0, YK)
        check = Matrix.infNorm(YK)
    }
    if (Verbose && numit == MaxNumIt && check > 10.0.pow(-12.0)) {
        println("Maximum Number of Iterations reached")
    }

    R = A2.mult(Matrix.eye(m).add(-1.0, A1.add(1.0, A2.mult(G))).inv())

    var U = A1.add(1.0, R.mult(A0))
    if (continues) {
        U = Matrix.scaleMult(U.add(-1.0, Matrix.eye(m)), lamb)
    }

    result = HashMap()
    result["G"] = G
    result["R"] = R
    result["U"] = U
    return result
}