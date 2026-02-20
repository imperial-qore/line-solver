package jline.lib.butools

import jline.GlobalConstants.Inf
import jline.GlobalConstants.NegInf
import kotlin.math.sqrt

fun APH3rdMomentUpperBound(m1: Double, m2: Double, n: Int): Double {
    val n2 = m2 / m1 / m1
    return if (n2 < (n + 1.0) / n) {
        NegInf
    } else if (n2 <= n / (n - 1.0)) {
        m1 * m2 * (2.0 * (n - 2.0) * (n * n2 - n - 1.0) * sqrt(1.0 + (n * (n2 - 2.0)) / (n - 1.0)) + (n + 2.0) * (3.0 * n * n2 - 2.0 * n - 2.0)) / (n * n * n2)
    } else {
        Inf
    }
}