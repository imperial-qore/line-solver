package jline.api.mdd;

/**
 * Plain-array export of an {@link MDD}, the input contract of {@link Mdd_mcd}.
 *
 * <p>Mirrors MDD.toStruct in MATLAB and MDD.to_struct in python.</p>
 */
public class MddStruct {

    /** Number of variable levels. */
    public final int K;
    /** domain[k] is the number of local states at level k. */
    public final int[] domain;
    /** Id of the top (level-0) node; MDD.TERM_FALSE for the empty set. */
    public final int root;
    /** nnodes[k] is the live node count at level k. */
    public final int[] nnodes;
    /**
     * node[k][p][v] is the child of arc v of level-k node id p+1: a level-(k+1)
     * node id when k &lt; K-1, or a terminal when k == K-1.
     */
    public final int[][][] node;

    public MddStruct(int K, int[] domain, int root, int[] nnodes, int[][][] node) {
        this.K = K;
        this.domain = domain;
        this.root = root;
        this.nnodes = nnodes;
        this.node = node;
    }
}
