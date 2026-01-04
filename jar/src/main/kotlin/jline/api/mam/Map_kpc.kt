/**
 * @file Markovian Arrival Process Kronecker product composition
 * 
 * Computes Kronecker product composition of multiple MAPs for building complex arrival processes.
 * Used for modeling independent superposition and synchronization in multi-stream systems.
 * 
 * @since LINE 3.0
 */
package jline.api.mam

import jline.util.matrix.Matrix

/**
 * Kronecker product composition of MAPs.
 * 
 * @param MAPa First MAP as Array<Matrix> or List<Array<Matrix>> for multiple MAPs
 * @param MAPb Second MAP as Array<Matrix> (optional if MAPa contains multiple MAPs)
 * @return Composite MAP as Array<Matrix> where result[0] = D0 and result[1] = D1
 */
fun map_kpc(MAPa: Any, MAPb: Array<Matrix>? = null): Array<Matrix> {
    return when {
        // Case: MAPa is a list of MAPs, compose them all
        MAPa is List<*> && MAPa.isNotEmpty() && MAPa[0] is Array<*> -> {
            @Suppress("UNCHECKED_CAST")
            val mapList = MAPa as List<Array<Matrix>>
            
            if (mapList.size < 2) {
                throw IllegalArgumentException("Need at least 2 MAPs for composition")
            }
            
            // Start with first two MAPs
            var result = map_kpc(mapList[0], mapList[1])
            
            // Compose with remaining MAPs
            for (k in 2..<mapList.size) {
                result = map_kpc(result, mapList[k])
            }
            
            result
        }
        
        // Case: MAPa and MAPb are both single MAPs
        MAPa is Array<*> && MAPb != null -> {
            @Suppress("UNCHECKED_CAST")
            val mapA = MAPa as Array<Matrix>
            
            if (mapA.size != 2 || MAPb.size != 2) {
                throw IllegalArgumentException("Each MAP must have exactly 2 matrices [D0, D1]")
            }
            
            // MAP = {-kron(MAPa{1}, MAPb{1}), kron(MAPa{2}, MAPb{2})}
            val D0 = mapA[0].kron(MAPb[0])
            D0.scale(-1.0)  // Negate the result
            
            val D1 = mapA[1].kron(MAPb[1])
            
            arrayOf(D0, D1)
        }
        
        else -> {
            throw IllegalArgumentException("Invalid input: MAPa must be Array<Matrix> or List<Array<Matrix>>, and MAPb must be provided for single MAP case")
        }
    }
}

/**
 * Convenience function for composing exactly two MAPs.
 */
fun map_kpc(MAPa: Array<Matrix>, MAPb: Array<Matrix>): Array<Matrix> {
    return map_kpc(MAPa as Any, MAPb)
}

/**
 * Convenience function for composing a list of MAPs.
 */
fun map_kpc(maps: List<Array<Matrix>>): Array<Matrix> {
    return map_kpc(maps as Any, null)
}
/**
 * MAP kpc algorithms
 */
@Suppress("unused")
class MapKpcAlgo {
    companion object {
        // Class documentation marker for Dokka
    }
}