package jline.lib.smc

import jline.util.matrix.Matrix

fun QBD_EG(A0: Matrix, A1: Matrix, A2: Matrix, optVerbose: Boolean): MutableMap<String?, Matrix?> {
    var G = Matrix(0, 0, 0)
    var R = Matrix(0, 0, 0)
    var U = Matrix(0, 0, 0)
    val m = A1.numRows
    val theta = stat(A0.add(1.0, A1).add(1.0, A2))
    val drift = theta.mult(A0.sumRows())[0] - theta.mult(A2.sumRows())[0]
    if (drift > 0) {
        if (A0.rank() == 1) {
            // find index of first nonzero row
            var non_zero_row = 0
            val row_sum = A0.sumRows()
            for (i in 0..<A0.numRows) {
                if (row_sum[i] > 0) {
                    non_zero_row = i
                    break
                }
            }
            val beta = Matrix.scaleMult(Matrix.extractRows(A0, non_zero_row, non_zero_row + 1, null),
                1 / A0.sumRows(non_zero_row))
            G = Matrix.ones(m, 1).mult(beta)
            R = A2.mult(Matrix.eye(m).add(-1.0, A1.add(1.0, A2.mult(G))).inv())
        } else if (A2.rank() == 1) {
            val eta = QBD_CAUDAL(A0, A1, A2)
            R = A2.mult(Matrix.eye(m).add(-1.0, A1).add(-eta, A0).inv())
            G = Matrix.eye(m).add(-1.0, A1.add(1.0, R.mult(A0))).inv().mult(A0)
        }
    } else if (drift < 0) {
        if (A2.rank() == 1) {
            val alpha = A2.mult(Matrix.ones(m, 1))
            R = Matrix.scaleMult(alpha.mult(theta), theta.mult(alpha)[0])
            G = Matrix.eye(m).add(-1.0, A1.add(1.0, R.mult(A0))).inv().mult(A0)
        } else if (A0.rank() == 0) {
            val A0hat = Matrix.diag(*theta.elementPow(-1.0).toArray1D()).mult(A2.transpose())
                .mult(Matrix.diag(*theta.toArray1D()))
            val A1hat = Matrix.diag(*theta.elementPow(-1.0).toArray1D()).mult(A1.transpose())
                .mult(Matrix.diag(*theta.toArray1D()))
            val A2hat = Matrix.diag(*theta.elementPow(-1.0).toArray1D()).mult(A0.transpose())
                .mult(Matrix.diag(*theta.toArray1D()))
            val etahat = QBD_CAUDAL(A0hat, A1hat, A2hat)
            G = Matrix.diag(*theta.elementPow(-1.0).toArray1D())
                .mult(A2hat.mult(Matrix.eye(m).add(-1.0, A1hat).add(-etahat, A0hat).inv()).transpose())
                .mult(Matrix.diag(*theta.toArray1D()))
            R = A2.mult(Matrix.eye(m).add(-1.0, A1.add(1.0, A2.mult(G))).inv())
        }
    }

    if (!R.isEmpty) {
        U = A1.add(1.0, R.mult(A0))
    }

    val result: MutableMap<String?, Matrix?> = HashMap()
    result["G"] = G
    result["R"] = R
    result["U"] = U
    return result
}