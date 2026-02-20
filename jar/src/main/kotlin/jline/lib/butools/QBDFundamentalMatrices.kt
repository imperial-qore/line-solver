package jline.lib.butools

import jline.lib.smc.QBD_FI
import jline.lib.smc.QBD_IS
import jline.lib.smc.QBD_LR
import jline.lib.smc.QBD_NI
import jline.util.matrix.Matrix

fun QBDFundamentalMatrices(B: Matrix?,
                           L: Matrix,
                           F: Matrix?,
                           precision_: Double?,
                           maxNumIt_: Int?,
                           method_: String?,
                           Verbose_: Int?): Map<String?, Matrix?> {
    var precision = 1e-14
    if (precision_ != null) {
        precision = precision_
    }

    var maxNumIt = 50
    if (maxNumIt_ != null) {
        maxNumIt = maxNumIt_
    }

    var method = "CR"
    if (method_ != null) {
        method = method_
    }

    var Verbose = 0
    if (Verbose_ != null) {
        Verbose = Verbose_
    }

    if (method == "LR") {
        return QBD_LR(B!!, L, F!!, maxNumIt, Verbose, null, null)
    }

    if (method == "NI") {
        return QBD_NI(B!!, L, F!!, maxNumIt, Verbose, null, null)
    }

    if (method == "IS") {
        return QBD_IS(B!!, L, F!!, maxNumIt, Verbose, null, null)
    }

    if (method == "FI") {
        return QBD_FI(B!!, L, F!!, maxNumIt, Verbose, null, null, null)
    }

    return jline.lib.smc.QBD_CR(B!!, L, F!!, maxNumIt, Verbose, null, null)
}