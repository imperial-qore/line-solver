package jline.util.graph;

import jline.util.matrix.Matrix;

import java.util.Comparator;
import java.util.Iterator;
import java.util.Set;
import java.util.TreeSet;

/**
 * An undirected graph data structure with weighted edges represented as an adjacency matrix.
 * 
 * <p>This class provides graph algorithms for undirected graphs including connected component
 * detection and graph analysis. The adjacency matrix is symmetric, with edges (i,j) and (j,i)
 * having the same weight.</p>
 * 
 * <p>Key features:
 * <ul>
 *   <li>Weakly connected component (WCC) detection</li>
 *   <li>Symmetric adjacency matrix representation</li>
 *   <li>Weighted edge support with normalization options</li>
 *   <li>Column filtering for selective analysis</li>
 * </ul>
 * </p>
 * 
 * @see DirectedGraph
 * @see Matrix
 */
public class UndirectedGraph {

    /** Number of vertices in the graph */
    private final int V;
    /** Symmetric adjacency matrix storing edge weights */
    private final Matrix adjacencyMatrix;  // Matrix to store the edges and weights
    /** Set of weakly connected components */
    private final Set<Set<Integer>> wcc;
    /** Set of column indices to ignore during analysis */
    private final Set<Integer> colsToIgnore;

    /**
     * Constructs an undirected graph with the given adjacency matrix and column filter.
     * Weights for bidirectional edges are summed.
     * 
     * @param param adjacency matrix where entry (i,j) represents edge weight
     * @param colsToIgnore set of column indices to ignore, can be null
     */
    public UndirectedGraph(Matrix param, Set<Integer> colsToIgnore) {
        this(param, colsToIgnore, false);
    }

    /**
     * Constructs an undirected graph with optional weight normalization.
     * 
     * @param param adjacency matrix where entry (i,j) represents edge weight
     * @param colsToIgnore set of column indices to ignore, can be null
     * @param normalize if true, weights are averaged for links (i,j) and (j,i);
     *                  if false, weights are summed
     */
    public UndirectedGraph(Matrix param, Set<Integer> colsToIgnore, boolean normalize) {
        this.V = param.getNumCols();
        this.adjacencyMatrix = param.copy().add(param.copy().transpose());
        if (normalize) {
            this.adjacencyMatrix.scaleEq(1.0 / 2);  // Initialize an empty VxV matrix
        }
        this.colsToIgnore = colsToIgnore;

        this.wcc = new TreeSet<>(new Comparator<Set<Integer>>() {
            @Override
            public int compare(Set<Integer> arg0, Set<Integer> arg1) {
                Iterator<Integer> it0 = arg0.iterator();
                Iterator<Integer> it1 = arg1.iterator();
                while (it0.hasNext()) {
                    if (it1.hasNext()) {
                        int val0 = it0.next();
                        int val1 = it1.next();
                        if (val0 != val1)
                            return val0 - val1;
                    } else {
                        return 1;
                    }
                }

                if (it1.hasNext())
                    return -1;
                else
                    return 0;
            }
        });
    }

    public UndirectedGraph(Matrix param) {
        this(param, null);
    }

    // Add an undirected edge between vertices `s` and `d` with the specified weight
    public void addEdge(int s, int d, double weight) {
        if (s >= V || d >= V) {
            throw new RuntimeException("The index of row or column is out of bounds");
        }
        adjacencyMatrix.set(s, d, weight);  // Set the edge weight in the matrix
        adjacencyMatrix.set(d, s, weight);  // Since undirected, set both ways
    }

    public void computeWeaklyConnectedComponents() {
        boolean[] visitedVertices = new boolean[V];

        for (int i = 0; i < V; i++) {
            if (!visitedVertices[i]) {
                Set<Integer> res = new TreeSet<>();
                findCC(i, visitedVertices, res);
                if (colsToIgnore != null) {
                    for (int num : colsToIgnore) {
                        res.remove(num);
                    }
                }
                wcc.add(res);
            }
        }
    }

    // Helper method to find connected components using DFS
    public void findCC(int i, boolean[] visited, Set<Integer> chain) {
        visited[i] = true;
        chain.add(i);

        for (int neighbor = 0; neighbor < V; neighbor++) {
            if (adjacencyMatrix.get(i, neighbor) > 0 && !visited[neighbor]) {  // If there's an edge and neighbor is not visited
                findCC(neighbor, visited, chain);
            }
        }
    }

    public Set<Set<Integer>> getWCC() {
        return this.wcc;
    }

    // Method to convert the graph to an adjacency matrix representation
    public Matrix toMatrix() {
        return this.adjacencyMatrix;  // Return the internal adjacency matrix
    }
}
