package jline.lib.smc

import jline.util.matrix.Matrix
import org.apache.commons.math3.util.FastMath
import kotlin.math.pow

fun QBD_CR(A0: Matrix,
           A1: Matrix,
           A2: Matrix,
           MaxNumIt_: Int?,
           Verbose_: Int?,
           Mode_: String?,
           RAPComp_: Int?): Map<String?, Matrix?> {
    var A0 = A0
    var A1 = A1
    var A2 = A2
    var Mode = "Shift"
    var MaxNumIt = 50
    var Verbose = false
    var RAPComp = false

    if (MaxNumIt_ != null) {
        MaxNumIt = MaxNumIt_
    }

    if (Mode_ != null) {
        if (Mode_ == "Shift" || Mode_ == "Basic") {
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
    var G = result["G"]
    // QBD_EG only handles special cases (rank-1 A0 or A2)
    // For general cases, proceed with the iterative cyclic reduction method
    val theta = stat(A0.add(1.0, A1).add(1.0, A2))
    val drift = theta.mult(A0.sumRows())[0] - theta.mult(A2.sumRows())[0]
    val A2old = A2.copy()
    val uT = Matrix.scaleMult(Matrix.ones(1, m), 1.0 / m.toDouble())
    val A0old = A0.copy()
    if (Mode == "Shift") {
        if (drift < 0) {
            A2 = A2.add(-1.0, Matrix.ones(m, 1).mult(theta.mult(A2)))
            A1 = A1.add(1.0, Matrix.ones(m, 1).mult(theta.mult(A0)))
        } else {
            A0 = A0.add(-1.0, A0.sumRows().mult(uT))
            A1 = A1.add(1.0, A2.sumRows().mult(uT))
        }
    }
    var A = A1.copy()
    var B = A2.copy()
    var C = A0.copy()

    var Ahat = A.copy()
    var check = 1.0
    var numit = 0

    while (check > 10.0.pow(-14.0) && numit < MaxNumIt) {
        var Atemp = Matrix.eye(m).add(-1.0, A).inv()
        val BAtemp = B.mult(Atemp)
        Atemp = C.mult(Atemp)
        Ahat = A.add(1.0, BAtemp.mult(C))
        A = A.add(1.0, BAtemp.mult(A)).add(1.0, Atemp.mult(B))
        B = BAtemp.mult(B)
        C = Atemp.mult(C)
        numit = numit + 1
        check = FastMath.min(Matrix.infNorm(B), Matrix.infNorm(C))
    }
    if (Verbose && numit == MaxNumIt && check > 10.0.pow(-14.0)) {
        println("Maximum Number of Iterations reached")
    }

    G = Matrix.eye(m).add(-1.0, Ahat).inv().mult(A0)

    if (Mode == "Shift") {
        if (drift < 0) {
            A1 = A1.add(-1.0, Matrix.ones(m, 1).mult(theta).mult(A0))
            A2 = A2old.copy()
        } else {
            G = G.add(1.0, Matrix.ones(m, 1).mult(uT))
            A1 = A1.add(-1.0, A2.sumRows().mult(uT))
            A0 = A0old.copy()
        }
    }
    val R = A2.mult(Matrix.eye(m).add(-1.0, A1.add(1.0, A2.mult(G))).inv())

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