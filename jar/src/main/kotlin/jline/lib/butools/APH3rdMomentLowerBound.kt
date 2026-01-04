package jline.lib.butools

import jline.GlobalConstants.Inf
import org.apache.commons.math3.complex.Complex
import kotlin.math.sqrt

fun APH3rdMomentLowerBound(m1: Double, m2: Double, n: Int): Double {
    val ni2 = m2 / m1 / m1
    if (ni2 < (n + 1.0) / n) {
        return Inf
    } else if (ni2 < (n + 4.0) / (n + 1.0)) {
        val n2 = Complex(ni2, 0.0)
        val three = Complex(3.0, 0.0)
        val negTwo = Complex(-2.0, 0.0)

        val numeratorP = (Complex(n + 1.0, 0.0)).multiply(n2.subtract(2.0))
        val denominatorP = three.multiply(n2).multiply(Complex(n - 1.0, 0.0))
        val sqrtPart =
            (negTwo.multiply(Complex(sqrt(n + 1.0), 0.0))).divide(Complex(sqrt(-3.0 * n * ni2 + 4.0 * n + 4.0), 0.0))
                .subtract(1.0)
        val p = numeratorP.divide(denominatorP).multiply(sqrtPart)

        val numeratorA = n2.subtract(2.0)
        val sqrtPartA = p.multiply(p)
            .add(p.multiply(Complex(n.toDouble(), 0.0)).multiply(n2.subtract(2.0)).divide(Complex(n - 1.0, 0.0))).sqrt()
        val denominatorA = p.multiply(Complex(1.0 - n2.real, 0.0)).add(sqrtPartA)
        val a = numeratorA.divide(denominatorA)

        val numeratorL1 = (Complex(3.0 + a.real, 0.0)).multiply(Complex(n - 1.0, 0.0)).add(a.multiply(2.0))
        val denominatorL1 = Complex((n - 1.0) * (1.0 + a.multiply(p).real), 0.0)
        val numeratorL2 = a.multiply(Complex(2.0 * (n + 1.0), 0.0))
        val denominatorL2 =
            Complex(2.0 * (n - 1.0) + a.multiply(p).multiply(Complex(n * a.real + 2.0 * n - 2.0, 0.0)).real, 0.0)
        val l = numeratorL1.divide(denominatorL1).subtract(numeratorL2.divide(denominatorL2))

        return l.real * m1 * m2
    } else {
        return (n + 1.0) / n * ni2 * m1 * m2
    }
}