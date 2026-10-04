package jline.api.mdd;

/** Storage description of the set held in an {@link MDD}. */
public class MddStats {

    /** Number of variable levels. */
    public final int levels;
    /** Reachable node count per level. */
    public final int[] nodesPerLevel;
    /** Reachable non-terminal nodes. */
    public final int numNodes;
    /** Nodes physically held in the tables, dead ones included. */
    public final int tableNodes;
    /** |S|. */
    public final long numStates;
    /** Integers in the reachable arc arrays, the diagram footprint. */
    public final long mddInts;
    /** Integers an explicit state list would need, |S| * K. */
    public final long explicitInts;

    public MddStats(int levels, int[] nodesPerLevel, int numNodes, int tableNodes,
                    long numStates, long mddInts, long explicitInts) {
        this.levels = levels;
        this.nodesPerLevel = nodesPerLevel;
        this.numNodes = numNodes;
        this.tableNodes = tableNodes;
        this.numStates = numStates;
        this.mddInts = mddInts;
        this.explicitInts = explicitInts;
    }

    /** Explicit footprint divided by the diagram footprint. */
    public double compression() {
        return (double) explicitInts / Math.max(mddInts, 1L);
    }
}
