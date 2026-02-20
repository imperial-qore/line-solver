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
     * Checks if this fork parameter container is empty (no fan-out configured).
     * 
     * @return true if fan-out ratio is not set (NaN), false otherwise
     */
    @Override
    public boolean isEmpty() {
        return Double.isNaN(fanOut);
    }
}
