package jline.lib.smc

import jline.util.matrix.Matrix
import kotlin.math.pow

fun QBD_ParsePara(A0: Matrix, A1: Matrix, A2: Matrix) {
    if (A0.numCols != A0.numRows) {
        throw RuntimeException("A0 is not a square matrix")
    }
    if (A1.numCols != A1.numRows) {
        throw RuntimeException("A1 is not a square matrix")
    }
    if (A2.numCols != A2.numRows) {
        throw RuntimeException("A2 is not a square matrix")
    }
    if (A0.numCols != A1.numCols) {
        throw RuntimeException("The matrices A0 and A1 do not have the same dimension")
    }
    if (A0.numCols != A2.numCols) {
        throw RuntimeException("The matrices A0 and A2 do not have the same dimension")
    }
    if (A0.elementMin() <= -(-10.0).pow(-14.0)) {
        throw RuntimeException("The matrix A0 contains negative data")
    }
    if (A1.elementMin() <= -(-10.0).pow(-14.0)) {
        throw RuntimeException("The matrix A1 contains negative data")
    }
    if (A2.elementMin() <= -(-10.0).pow(-14.0)) {
        throw RuntimeException("he matrix A2 contains negative data")
    }
    if (A0.add(1.0, A1).add(1.0, A2).sumRows().elementMax() > 1 + 10.0.pow(-14.0)) {
        throw RuntimeException("The matrix A0+A1+A2 has to be (sub)stochastic")
    }
}