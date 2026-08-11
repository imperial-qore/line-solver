package jline.lib.mom.util;

import java.util.Arrays;
import java.util.Comparator;

/**
 * MOM-specific utility functions for combinatorial operations
 * Complements the existing jline.util.Maths class
 */
public final class MomUtils {
    private MomUtils() {}

    /**
     * Sort combinations by number of non-zeros and their positions.
     * This is used to optimize the order of processing in the solver.
     *
     * @param combinations Array of combinations to sort
     * @return Sorted array of combinations
     */
    public static int[][] sortByNnzPos(int[][] combinations) {
        int[][] result = new int[combinations.length][];
        for (int i = 0; i < combinations.length; i++) {
            result[i] = combinations[i];
        }
        Arrays.sort(result, new Comparator<int[]>() {
            @Override
            public int compare(int[] a, int[] b) {
                int nnzA = countNonZeros(a);
                int nnzB = countNonZeros(b);
                if (nnzA != nnzB) return Integer.compare(nnzA, nnzB);
                // Then by position of non-zeros (lexicographic)
                StringBuilder sa = new StringBuilder();
                for (int i = 0; i < a.length; i++) {
                    if (a[i] != 0) sa.append(i);
                }
                StringBuilder sb = new StringBuilder();
                for (int i = 0; i < b.length; i++) {
                    if (b[i] != 0) sb.append(i);
                }
                return sa.toString().compareTo(sb.toString());
            }
        });
        return result;
    }

    /**
     * Count non-zero elements in an array.
     *
     * @param array The array to count non-zeros in
     * @return Number of non-zero elements
     */
    public static int countNonZeros(int[] array) {
        int count = 0;
        for (int v : array) {
            if (v != 0) count++;
        }
        return count;
    }

    /**
     * Hash a population vector to a unique index.
     * This is used for efficient lookup of population states.
     *
     * @param population The population vector
     * @param maxPop Maximum population per class (for hash calculation)
     * @return Hash index
     */
    public static int hashPop(int[] population, int[] maxPop) {
        int hash = 0;
        int multiplier = 1;
        for (int i = 0; i < population.length; i++) {
            hash += population[i] * multiplier;
            multiplier *= (maxPop[i] + 1);
        }
        return hash;
    }
}
