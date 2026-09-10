package jline.util;

/**
 * IndexedMinHeap
 *
 * A min-heap over a fixed set of integer keys 0..n-1, each associated
 * with a double priority. Supports O(log n) update of any key by index,
 * and O(1) peek at the minimum.
 *
 * Typical use in NRM: keys are reaction indices, priorities are absolute
 * firing times τ_k.
 */
public class IndexedMinHeap {

    private final int   n;
    private final int[] heap;   // heap[pos]  = reaction index at heap position pos
    private final int[] pos;    // pos[k]     = current position of reaction k in heap
    public  final double[] key; // key[k]     = current priority of reaction k

    /**
     * Construct a heap of size n. All keys initialised to +∞.
     * Bulk-set key[] then call buildHeap(), or call update() individually.
     */
    public IndexedMinHeap(int n) {
        this.n    = n;
        this.heap = new int[n];
        this.pos  = new int[n];
        this.key  = new double[n];
        for (int i = 0; i < n; i++) {
            heap[i] = i;
            pos[i]  = i;
            key[i]  = Double.POSITIVE_INFINITY;
        }
    }

    // ------------------------------------------------------------------
    // Public API
    // ------------------------------------------------------------------

    /** Return the index of the reaction with the smallest key. O(1). */
    public int peekMin() {
        return heap[0];
    }

    /** Return the smallest key value. O(1). */
    public double peekMinKey() {
        return key[heap[0]];
    }

    /** Update the priority of reaction k and restore heap order. O(log n). */
    public void update(int k, double newKey) {
        key[k] = newKey;
        siftUp(pos[k]);
        siftDown(pos[k]);
    }

    /** Read the current priority of reaction k. O(1). */
    public double getKey(int k) {
        return key[k];
    }

    /**
     * Establish heap order after bulk-writing into key[] directly.
     * Runs in O(n) — prefer this over n individual update() calls at init.
     */
    public void buildHeap() {
        for (int p = n / 2 - 1; p >= 0; p--) siftDown(p);
    }

    // ------------------------------------------------------------------
    // Internal helpers
    // ------------------------------------------------------------------

    private void swap(int i, int j) {
        int tmp  = heap[i]; heap[i]  = heap[j]; heap[j]  = tmp;
        pos[heap[i]] = i;   pos[heap[j]] = j;
    }

    private void siftUp(int p) {
        while (p > 0) {
            int parent = (p - 1) / 2;
            if (key[heap[parent]] > key[heap[p]]) {
                swap(parent, p);
                p = parent;
            } else break;
        }
    }

    private void siftDown(int p) {
        while (true) {
            int smallest = p;
            int l = 2 * p + 1;
            int r = 2 * p + 2;
            if (l < n && key[heap[l]] < key[heap[smallest]]) smallest = l;
            if (r < n && key[heap[r]] < key[heap[smallest]]) smallest = r;
            if (smallest != p) { swap(p, smallest); p = smallest; }
            else break;
        }
    }
}