package jline.lib.butools

fun NormMomsFromMoms(m: DoubleArray): DoubleArray {
    val length = m.size
    val nm = DoubleArray(length)
    nm[0] = m[0]

    for (i in 1..<length) {
        nm[i] = m[i] / (m[i - 1] * m[0])
    }

    return nm
}