package jline.api.mam

import jline.util.matrix.Matrix
import jline.util.matrix.MatrixCell
import java.util.*

/*
         * Safely combines multiple MMAPs into a single superposed MMAP while considering order constraints.
         *
         *
         * This method creates a superposed MMAP from a given set of MMAPs, ensuring that the combined MMAP
         * does not exceed a specified maximum order (`maxorder`). The MMAPs are combined in an order determined
         * by the squared coefficient of variation (SCV) of the unmarked processes. If combining two MMAPs
         * would exceed the maximum order, alternative methods are used, such as fitting to an exponential or
         * using a simplified MAP.
         *
         * @param MMAPS    a map of MMAPs to be combined
         * @param maxorder the maximum allowed order for the resulting superposed MMAP
         * @param method   the method for combining MMAPs; "default" or "match"
         * @return a MatrixCell representing the combined superposed MMAP
         */
fun mmap_super_safe(MMAPS: Map<Int?, MatrixCell>, maxorder: Int, method: String = "default"): MatrixCell? {
    // Handle empty hashmap case
    if (MMAPS.isEmpty()) {
        val lambda = Matrix(1, 1)
        lambda[0, 0] = 1.0
        return mmap_exponential(lambda)
    }
    
    var sup: MatrixCell? = MatrixCell()
    val scv_unmarked: MutableList<Double> = ArrayList()
    for (i in 0..<MMAPS.size) {
        scv_unmarked.add(map_scv(MMAPS[i]!![0], MMAPS[i]!![1]))
    }

    val indices = arrayOfNulls<Int>(scv_unmarked.size)
    for (i in indices.indices) {
        indices[i] = i
    }
    Arrays.sort(indices, Comparator.comparingDouble { index: Int? ->
        scv_unmarked[index!!]
    })
    val sortedIndices = Arrays.stream(indices).mapToInt { obj: Int? -> obj!!.toInt() }.toArray()

    for (i in sortedIndices.indices) {
        val smallest_value = sortedIndices[i]
        if (sup!!.isEmpty) {
            sup = MMAPS[smallest_value]
            if (maxorder == 1) {
                sup = mmap_exponential(mmap_lambda(MMAPS[smallest_value]!!))
            }
        } else {
            sup = if (sup[0].length() * MMAPS[smallest_value]!![0].length() > maxorder) {
                if (sup[0].length() * 2 < maxorder) {
                    mmap_super(sup, mamap2m_fit_gamma_fb_mmap(MMAPS[smallest_value]!!), method)
                } else {
                    mmap_super(sup, mmap_exponential(mmap_lambda(MMAPS[smallest_value]!!))!!, method)
                }
            } else {
                mmap_super(sup, MMAPS[smallest_value]!!, method)
            }
        }
    }


    return sup
}
/**
 * MMAP super safe algorithms
 */
@Suppress("unused")
class MmapSuperSafeAlgo {
    companion object {
        // Class documentation marker for Dokka
    }
}