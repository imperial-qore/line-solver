package jline.lib.perm;

import java.util.ArrayList;
import java.util.List;

import jline.util.matrix.Matrix;

/**
 * Implementation of the naive exact permanent computation.
 *
 * This solver computes the exact permanent of a matrix by iterating over all
 * possible permutations and summing the products of matrix elements according
 * to each permutation. While this gives the exact result, it has O(n!) complexity
 * and is only practical for small matrices.
 *
 * The permanent of an n×n matrix A is defined as:
 * perm(A) = Σ_π ∏_{i=1}^n a_{i,π(i)}
 * where the sum is over all permutations π of {1,2,...,n}.
 */
public class NaivePermanent extends PermSolver {

    public NaivePermanent(Matrix matrix) {
        this(matrix, false);
    }

    public NaivePermanent(Matrix matrix, boolean solve) {
        super(matrix);
        if (solve) {
            solve();
        }
    }

    @Override
    public void compute() {
        value = naive();
    }

    /**
     * Compute the exact permanent using the naive permutation-based method.
     */
    private double naive() {
        double permanent = 0.0;

        // Generate all permutations of column indices
        int[] indices = new int[n];
        for (int i = 0; i < n; i++) indices[i] = i;
        List<int[]> permutations = generatePermutations(indices);

        // Iterate over all permutations
        for (int[] permutation : permutations) {
            double product = 1.0;

            // Compute product for this permutation
            for (int i = 0; i < n; i++) {
                product *= matrix.get(i, permutation[i]);
            }

            permanent += product;
        }

        return permanent;
    }

    /**
     * Generate all permutations of the given array using Heap's algorithm.
     */
    private List<int[]> generatePermutations(int[] array) {
        List<int[]> result = new ArrayList<int[]>();
        heapPermutation(array.clone(), array.length, result);
        return result;
    }

    /**
     * Heap's algorithm for generating permutations.
     */
    private void heapPermutation(int[] array, int size, List<int[]> result) {
        if (size == 1) {
            result.add(array.clone());
            return;
        }

        for (int i = 0; i < size; i++) {
            heapPermutation(array, size - 1, result);

            // If size is odd, swap first and last element
            if (size % 2 == 1) {
                int temp = array[0];
                array[0] = array[size - 1];
                array[size - 1] = temp;
            } else {
                // If size is even, swap ith and last element
                int temp = array[i];
                array[i] = array[size - 1];
                array[size - 1] = temp;
            }
        }
    }
}
