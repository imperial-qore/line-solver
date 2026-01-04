package jline.util.graph;

import jline.GlobalConstants;
import jline.util.matrix.Matrix;

import java.util.*;

/**
 * A directed graph data structure with weighted edges represented as an adjacency matrix.
 * 
 * <p>This class provides graph algorithms for directed graphs including cycle detection,
 * topological sorting, and connectivity analysis. It's particularly useful for routing
 * matrix analysis, dependency resolution, and network topology validation.</p>
 * 
 * <p>Key graph algorithms supported:
 * <ul>
 *   <li>Directed Acyclic Graph (DAG) detection</li>
 *   <li>Topological sorting using Kahn's algorithm</li>
 *   <li>Adjacency matrix-based representation</li>
 *   <li>Weighted edge support</li>
 *   <li>Column filtering for selective analysis</li>
 * </ul>
 * </p>
 * 
 * @see UndirectedGraph
 * @see Matrix
 * @since 1.0
 */
public class DirectedGraph {

    /** Number of vertices in the graph */
    private final int V;
    
    /** Adjacency matrix storing edge weights [V x V] */
    private final Matrix adjacencyMatrix;
    
    /** Set of column indices to ignore during analysis */
    private final Set<Integer> colsToIgnore;

    /**
     * Constructs a directed graph with the given adjacency matrix and column filter.
     * 
     * @param param adjacency matrix [V x V] where entry (i,j) represents edge weight from vertex i to j
     * @param colsToIgnore set of column indices to ignore during graph operations, can be null
     */
    public DirectedGraph(Matrix param, Set<Integer> colsToIgnore) {
        this.V = param.getNumCols();
        this.adjacencyMatrix = new Matrix(param);  // Initialize an empty VxV matrix
        this.colsToIgnore = colsToIgnore;
    }

    /**
     * Constructs a directed graph with the given adjacency matrix (no column filtering).
     * 
     * @param param adjacency matrix [V x V] where entry (i,j) represents edge weight from vertex i to j
     */
    public DirectedGraph(Matrix param) {
        this(param, null);
    }

    /**
     * Tests if the given adjacency matrix represents a Directed Acyclic Graph (DAG).
     * 
     * @param adj adjacency matrix [V x V] to test for cycles
     * @return true if the graph is acyclic, false if cycles are detected
     */
    public static boolean isDAG(Matrix adj) {
        return kahn(adj).getNumCols() == adj.getNumRows();
    }

    /**
     * Performs topological sorting using Kahn's algorithm.
     * 
     * <p>This method computes a topological ordering of vertices in a directed graph.
     * If the graph contains cycles, the returned ordering will be incomplete.</p>
     * 
     * @param adj adjacency matrix [V x V] representing the directed graph
     * @return matrix containing topologically sorted vertex indices, incomplete if graph has cycles
     */
    public static Matrix kahn(Matrix adj) {
        final double tol = GlobalConstants.FineTol;
        final int n = adj.getNumRows();

        int[] indegree = new int[n];
        for (int col = 0; col < n; col++) {
            for (int row = 0; row < n; row++) {
                if (Math.abs(adj.get(row, col)) > tol) indegree[col]++;
            }
        }

        List<Integer> order = new ArrayList<>(n);
        Deque<Integer> q = new ArrayDeque<>();
        for (int v = 0; v < n; v++) if (indegree[v] == 0) q.add(v);

        while (!q.isEmpty()) {
            int i = q.removeFirst();
            order.add(i);
            for (int j = 0; j < n; j++) {
                if (Math.abs(adj.get(i, j)) > tol && --indegree[j] == 0) q.addLast(j);
            }
        }

        Matrix res = new Matrix(1, order.size());
        for (int k = 0; k < order.size(); k++) res.set(0, k, order.get(k));
        return res;
    }

    // Add a directed edge from vertex `s` to vertex `d` with the specified weight
    public void addEdge(int s, int d, double weight) {
        if (s >= V || d >= V) {
            throw new RuntimeException("The index of row or column is out of bounds");
        }
        adjacencyMatrix.set(s, d, weight);  // Set the edge weight in the adjacency matrix
    }

    // Find all neighbors of a given node that have outgoing edges
    private List<Integer> find(int node) {
        List<Integer> neighbors = new ArrayList<>();
        for (int j = 0; j < V; j++) {
            if (adjacencyMatrix.get(node, j) > 0) {  // If there is an edge with weight > 0
                neighbors.add(j);
            }
        }
        return neighbors;
    }

    private boolean isInArray(int value, int[] array) {
        for (int i : array) {
            if (i == value) {
                return true;
            }
        }
        return false;
    }

    public SCCResult stronglyconncomp() {
        int idx = 0;
        List<int[]> SCC = new ArrayList<>();
        List<Integer> stk = new ArrayList<>();

        int[] v_idx = new int[V];
        int[] v_low = new int[V];
        boolean[] v_stk = new boolean[V];

        for (int i = 0; i < V; i++) {
            if (v_idx[i] == 0) {
                SCCAuxResult auxResult = stronglyconncomp_aux(i, v_idx, v_low, v_stk, SCC, stk, idx);
                v_idx = auxResult.v_idx;
                v_low = auxResult.v_low;
                v_stk = auxResult.v_stk;
                SCC = auxResult.SCC;
                stk = auxResult.stk;
                idx = auxResult.idx;
            }
        }

        // Sort SCCs by size
        SCC.sort((a, b) -> Integer.compare(b.length, a.length));

        int[] I = new int[V];
        for (int j = 0; j < SCC.size(); j++) {
            for (int node : SCC.get(j)) {
                I[node] = j + 1;
            }
        }

        boolean[] recurrent = new boolean[SCC.size()];
        for (int j = 0; j < SCC.size(); j++) {
            int[] scc = SCC.get(j);
            boolean is_recurrent = true;

            // Check if any node in the SCC has outgoing edges to nodes outside the SCC
            for (int node : scc) {
                List<Integer> out_edges = find(node);
                for (int out_node : out_edges) {
                    if (!isInArray(out_node, scc)) {
                        is_recurrent = false;
                        break;
                    }
                }
                if (!is_recurrent) {
                    break;
                }
            }

            recurrent[j] = is_recurrent;
        }

        return new SCCResult(I, recurrent);
    }

    private SCCAuxResult stronglyconncomp_aux(int i, int[] v_idx, int[] v_low, boolean[] v_stk, List<int[]> SCC, List<Integer> stk, int idx) {
        idx++;
        v_idx[i] = idx;
        v_low[i] = idx;
        stk.add(0, i);
        v_stk[i] = true;

        List<Integer> out_edges = find(i);

        for (int j : out_edges) {
            if (v_idx[j] == 0) {
                SCCAuxResult auxResult = stronglyconncomp_aux(j, v_idx, v_low, v_stk, SCC, stk, idx);
                v_idx = auxResult.v_idx;
                v_low = auxResult.v_low;
                v_stk = auxResult.v_stk;
                SCC = auxResult.SCC;
                stk = auxResult.stk;
                idx = auxResult.idx;
                v_low[i] = Math.min(v_low[i], v_low[j]);
            } else if (v_stk[j]) {
                v_low[i] = Math.min(v_low[i], v_idx[j]);
            }
        }

        if (v_low[i] == v_idx[i]) {
            int pos = stk.indexOf(i);
            int[] scc = new int[pos + 1];
            for (int k = 0; k <= pos; k++) {
                scc[k] = stk.get(k);
            }
            stk.subList(0, pos + 1).clear();
            for (int node : scc) {
                v_stk[node] = false;
            }
            SCC.add(scc);
        }

        return new SCCAuxResult(v_idx, v_low, v_stk, SCC, stk, idx);
    }

    // Updated toMatrix method to return the internal adjacency matrix
    public Matrix toMatrix() {
        return this.adjacencyMatrix;  // Return the internal adjacency matrix
    }

    public static class SCCAuxResult {
        public int[] v_idx;
        public int[] v_low;
        public boolean[] v_stk;
        public List<int[]> SCC;
        public List<Integer> stk;
        public int idx;

        public SCCAuxResult(int[] v_idx, int[] v_low, boolean[] v_stk, List<int[]> SCC, List<Integer> stk, int idx) {
            this.v_idx = v_idx;
            this.v_low = v_low;
            this.v_stk = v_stk;
            this.SCC = SCC;
            this.stk = stk;
            this.idx = idx;
        }
    }

    public static class SCCResult {
        public int[] I;
        public boolean[] recurrent;

        public SCCResult(int[] I, boolean[] recurrent) {
            this.I = I;
            this.recurrent = recurrent;
        }
    }
}
