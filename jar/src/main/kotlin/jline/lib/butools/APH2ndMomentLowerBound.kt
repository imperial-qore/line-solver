package jline.lib.butools

fun APH2ndMomentLowerBound(m1: Double, n: Int): Double {
    return m1 * m1 * (n + 1) / n
}