/**
 * @file Multi-class trace mean computation
 * 
 * Computes per-class means for multi-class empirical trace data. Enables separate 
 * statistical characterization of different job classes or traffic types within 
 * a single measurement trace for multi-class queueing model parameterization.
 * 
 * @since LINE 3.0
 */
package jline.api.trace

import jline.util.matrix.Matrix

/**
 * Computes the mean of a trace, divided by types.
 *
 * @param trace  the array containing the trace data
 * @param ntypes the number of different types
 * @param type   an array indicating the type of each element in the trace
 * @return a matrix containing the mean values for each type
 */
fun mtrace_mean(trace: DoubleArray, ntypes: Int, type: IntArray): Matrix {
    val mean = DoubleArray(ntypes)
    for (c in 0..<ntypes) {
        var sum = 0.0
        var ctr = 0
        for (i in trace.indices) {
            if (type[i] == c) {
                sum += trace[i]
                ctr++
            }
        }
        mean[c] = if (ctr > 0) sum / ctr.toDouble() else Double.NaN
    }
    return Matrix(mean).transpose()
}

/**    public static void main(String[] args) {
Matrix D0 = new Matrix("[-2,1;0,-3]");
 * /        Matrix D1 = new Matrix("[1,0;2,1]");
 * /        Matrix D11 = new Matrix("[1,0;2,0]");
 * /        Matrix D12 = new Matrix("[0,0;0,1]");
 * /        MatrixCell MAP = new MatrixCell(D0, D1);
 * /        MAP.set(0, D0);
 * /        MAP.set(1, D1);
 * /        MatrixCell MMAP = new MatrixCell(D0, D1);
 * /        MMAP.set(0, D0);
 * /        MMAP.set(1, D1);
 * /        MMAP.set(2, D11);
 * /        MMAP.set(3, D12);
 * /        System.out .println(map_mean(MAP));
 * /        System.out .println(mmap_lambda(MMAP).reciprocal());
 * /        mmapSampleReturn samples = mmap_sample(MMAP, 100000);
 * /        System.out .println(mtrace_mean(samples.samples, samples.ntypes, samples.types)); */
/**
 * Mtrace Mean algorithms
 */
@Suppress("unused")
class MtraceMeanAlgo {
    companion object {
        // Class documentation marker for Dokka
    }
}