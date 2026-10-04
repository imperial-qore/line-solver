/**
 * @file Cache Access Factor Computation
 *
 * @since LINE 3.0
 */
package jline.api.cache;

import java.util.ArrayDeque;
import java.util.ArrayList;
import java.util.Collections;
import java.util.Deque;
import java.util.List;

import jline.io.Ret;
import jline.util.matrix.Matrix;
import jline.util.matrix.MatrixCell;

public final class Cache_gamma {
    private Cache_gamma() {}

    /**
     * Computes access factors for the cache.
     */
    public static Ret.cacheGamma cache_gamma(MatrixCell lambda, MatrixCell R) {
        int u = lambda.size();
        int n = lambda.get(0).getNumRows();
        int h = lambda.get(0).getNumCols() - 1;
        Matrix gamma = new Matrix(n, h);

        for (int i = 0; i < n; i++) {
            for (int j = 0; j < h; j++) {
                MatrixCell RMatrixCell = (MatrixCell) (Object) R.get(0);
                Matrix graph = RMatrixCell.get(i);

                List<Integer> Pj = findShortestPath(graph, 0, j);

                if (Pj.isEmpty()) {
                    gamma.set(i, j, 0.0);
                } else {
                    double gammaValue = 0.0;
                    for (int v = 0; v < u; v++) {
                        gammaValue += lambda.get(v).get(i, 0);
                    }

                    for (int li = 1; li < Pj.size(); li++) {
                        double y = 0.0;
                        int l_1 = Pj.get(li - 1);
                        int l = Pj.get(li);
                        for (int v = 0; v < u; v++) {
                            MatrixCell RvMatrixCell = (MatrixCell) (Object) R.get(v);
                            y += lambda.get(v).get(i, l_1) * RvMatrixCell.get(i).get(l_1, l);
                        }
                        gammaValue *= y;
                    }
                    gamma.set(i, j, gammaValue);
                }
            }
        }

        return new Ret.cacheGamma(gamma, u, n, h);
    }

    private static List<Integer> findShortestPath(Matrix adjacencyMatrix, int source, int destination) {
        int n = adjacencyMatrix.getNumRows();
        if (source >= n || destination >= n || source < 0 || destination < 0) {
            return Collections.emptyList();
        }

        boolean[] visited = new boolean[n];
        int[] parent = new int[n];
        for (int i = 0; i < n; i++) parent[i] = -1;
        Deque<Integer> queue = new ArrayDeque<Integer>();

        queue.add(source);
        visited[source] = true;

        while (!queue.isEmpty()) {
            int current = queue.removeFirst();

            if (current == destination) {
                List<Integer> path = new ArrayList<Integer>();
                int node = destination;
                while (node != -1) {
                    path.add(0, node);
                    node = parent[node];
                }
                return path;
            }

            for (int next = 0; next < n; next++) {
                if (!visited[next] && adjacencyMatrix.get(current, next) > 0) {
                    visited[next] = true;
                    parent[next] = current;
                    queue.add(next);
                }
            }
        }
        return Collections.emptyList();
    }
}
