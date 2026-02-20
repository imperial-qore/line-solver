/*
 * Copyright (c) 2012-2026, QORE Lab, Imperial College London
 * All rights reserved.
 */
package jline.lib.smc

import jline.util.matrix.Matrix
import jline.util.matrix.Matrix.scaleMult
import kotlin.math.pow

fun QBD_CAUDAL(A0: Matrix, A1: Matrix, A2: Matrix, Dual: Boolean = false): Double {
    var A0 = A0
    var A2 = A2
    if (Dual) {
        val A2old = A2
        A2 = A0
        A0 = A2old
    }

    var eta_min = 0.0
    var eta_max = 1.0
    var eta = 0.5
    while (eta_max - eta_min > 10.0.pow(-15.0)) {
        val A3 = A2.add(1.0, scaleMult(A1, eta)).add(1.0, scaleMult(A0, eta.pow(2.0)))
        val eigs = A3.eigval()
        val new_eta = eigs.values.elementMax()
        if (new_eta > eta) {
            eta_min = eta
        } else {
            eta_max = eta
        }
        eta = (eta_max + eta_min) / 2
    }
    return eta
}


