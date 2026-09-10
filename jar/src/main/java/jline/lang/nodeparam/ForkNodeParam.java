package jline.lang.nodeparam;

import jline.lang.NodeParam;

/**
 * Parameter container for fork nodes in queueing networks.
 * 
 * <p>This class configures fork nodes that split incoming jobs into multiple
 * parallel tasks. Fork nodes are fundamental components in fork-join queueing
 * models where jobs require parallel processing across multiple servers or
 * subsystems before synchronization.</p>
 * 
 * <p>Fork node capabilities:
 * <ul>
 *   <li>Splitting jobs into parallel subtasks with configurable fan-out</li>
 *   <li>Supporting deterministic and probabilistic splitting</li>
 *   <li>Enabling parallel processing workflows</li>
 *   <li>Fork-join network modeling for multi-stage systems</li>
 * </ul>
 * </p>
 * 
 * @see jline.lang.nodes.Fork
 * @see JoinNodeParam
 * @since 1.0
 */
public class ForkNodeParam extends NodeParam {
    
    /**
     * Fan-out ratio specifying the number of parallel tasks created per job.
     * Default is NaN (not configured).
     */
    public double fanOut = Double.NaN;

    /**
     * Variable forking levels, all (nnodes x nclasses) and indexed by
     * DESTINATION NODE rather than link ordinal, so a relink cannot silently
     * permute them. {@code fanOut} above stays the scalar every existing
     * consumer reads.
     *
     * <p>{@code fanOutLink}: expected tasks sent to destination k for class r,
     * zero on a link the class does not take. {@code fanOutProb}: probability
     * the branch fires at all. {@code fanOutDist}: the jobs-per-link
     * distribution, indexed {@code [k][r]} with null meaning degenerate at
     * {@code fanOutLink}.</p>
     */
    public jline.util.matrix.Matrix fanOutLink;
    public jline.util.matrix.Matrix fanOutProb;
    public jline.lang.processes.DiscreteSampler[][] fanOutDist;

    /** FJ tag augmentation (ModelAdapter.fjtag): original classes forked here */
    public java.util.List<Integer> fjClasses;
    /** FJ tag augmentation: matched join node per forked class entry */
    public java.util.List<Integer> fjJoins;
    /** FJ tag augmentation: B x T auxiliary class indices per original class */
    public java.util.Map<Integer, jline.util.matrix.Matrix> fjAuxmatrix;
    /** FJ tag augmentation: branch head node indices per original class */
    public java.util.Map<Integer, int[]> fjBranchheads;

    /**
     * Checks if this fork parameter container is empty (no fan-out configured).
     * 
     * @return true if fan-out ratio is not set (NaN), false otherwise
     */
    @Override
    public boolean isEmpty() {
        return Double.isNaN(fanOut);
    }
}
