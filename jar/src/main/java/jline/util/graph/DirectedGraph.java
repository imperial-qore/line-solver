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

    // Find all nodes with incoming edges to this node (matching MATLAB's find(e(:,i)) for SCC traversal)
    private List<Integer> findIncoming(int node) {
        List<Integer> neighbors = new ArrayList<>();
        for (int j = 0; j < V; j++) {
            if (adjacencyMatrix.get(j, node) > 0) {  // If there is an edge from j to node (column lookup)
                neighbors.add(j);
            }
        }
        return neighbors;
    }

    // Find all nodes with outgoing edges from this node (matching MATLAB's find(A(node,:)) for recurrence check)
    private List<Integer> findOutgoing(int node) {
        List<Integer> neighbors = new ArrayList<>();
        for (int j = 0; j < V; j++) {
            if (adjacencyMatrix.get(node, j) > 0) {  // If there is an edge from node to j (row lookup)
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

    /**
     * Strongly connected components, with the bottom ones (BSCCs) flagged.
     *
     * <p>{@code I[v]} is the 1-based component index of vertex v, components numbered by
     * decreasing size (ties broken by discovery order). {@code recurrent[c]} is true when no
     * edge leaves component c; for the transition graph of a Markov chain those components are
     * exactly its recurrent classes and the rest are transient.</p>
     *
     * <p>Tarjan's algorithm with an explicit depth-first stack. The recursion is unrolled
     * deliberately: this runs on CTMC/DTMC state spaces whose graphs routinely contain paths
     * far longer than the JVM stack can carry, and a recursive formulation overflows on them.</p>
     */
    public SCCResult stronglyconncomp() {
        int[] v_idx = new int[V];    // discovery index, 0 while unvisited
        int[] v_low = new int[V];    // lowlink
        boolean[] v_stk = new boolean[V]; // currently on the Tarjan stack
        int[] comp = new int[V];     // component index in discovery order
        Arrays.fill(comp, -1);
        int idx = 0;
        int ncomp = 0;

        int[] stk = new int[V];      // Tarjan stack
        int nstk = 0;

        int[] frameV = new int[V];   // DFS stack: vertex of each open frame
        int[] frameK = new int[V];   // out-neighbours of that vertex already consumed
        @SuppressWarnings("unchecked")
        List<Integer>[] frameN = (List<Integer>[]) new List[V]; // neighbour list, fetched once

        for (int root = 0; root < V; root++) {
            if (v_idx[root] != 0) {
                continue;
            }

            idx++;
            v_idx[root] = idx;
            v_low[root] = idx;
            stk[nstk++] = root;
            v_stk[root] = true;

            int nf = 1;
            frameV[0] = root;
            frameK[0] = 0;
            frameN[0] = findOutgoing(root);

            while (nf > 0) {
                int v = frameV[nf - 1];
                List<Integer> nbrs = frameN[nf - 1];
                if (frameK[nf - 1] < nbrs.size()) {
                    int w = nbrs.get(frameK[nf - 1]);
                    frameK[nf - 1]++;
                    if (v_idx[w] == 0) {
                        idx++;
                        v_idx[w] = idx;
                        v_low[w] = idx;
                        stk[nstk++] = w;
                        v_stk[w] = true;
                        frameV[nf] = w;
                        frameK[nf] = 0;
                        frameN[nf] = findOutgoing(w);
                        nf++;
                    } else if (v_stk[w]) {
                        v_low[v] = Math.min(v_low[v], v_idx[w]);
                    }
                } else {
                    // v is exhausted: close its component if it is a root, then hand its
                    // lowlink back to the parent frame
                    if (v_low[v] == v_idx[v]) {
                        while (true) {
                            int w = stk[--nstk];
                            v_stk[w] = false;
                            comp[w] = ncomp;
                            if (w == v) {
                                break;
                            }
                        }
                        ncomp++;
                    }
                    nf--;
                    if (nf > 0) {
                        int p = frameV[nf - 1];
                        v_low[p] = Math.min(v_low[p], v_low[v]);
                    }
                }
            }
        }

        // Renumber by decreasing component size; the sort is stable, so components of equal
        // size keep their discovery order
        int[] counts = new int[ncomp];
        for (int v = 0; v < V; v++) {
            counts[comp[v]]++;
        }
        List<Integer> order = new ArrayList<Integer>(ncomp);
        for (int c = 0; c < ncomp; c++) {
            order.add(c);
        }
        final int[] countsRef = counts;
        Collections.sort(order, new Comparator<Integer>() {
            @Override
            public int compare(Integer a, Integer b) {
                return Integer.compare(countsRef[b], countsRef[a]);
            }
        });
        int[] relabel = new int[ncomp];
        for (int rank = 0; rank < ncomp; rank++) {
            relabel[order.get(rank)] = rank + 1;
        }

        int[] I = new int[V];
        for (int v = 0; v < V; v++) {
            I[v] = relabel[comp[v]];
        }

        // A component is recurrent iff no edge leaves it
        boolean[] recurrent = new boolean[ncomp];
        Arrays.fill(recurrent, true);
        for (int v = 0; v < V; v++) {
            List<Integer> out = findOutgoing(v);
            for (int k = 0; k < out.size(); k++) {
                if (I[out.get(k)] != I[v]) {
                    recurrent[I[v] - 1] = false;
                    break;
                }
            }
        }

        return new SCCResult(I, recurrent);
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
