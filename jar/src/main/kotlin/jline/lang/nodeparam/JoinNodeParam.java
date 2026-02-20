package jline.lang.nodeparam;

import jline.lang.JobClass;
import jline.lang.NodeParam;
import jline.lang.constant.JoinStrategy;

import java.util.Map;

/**
 * Parameter container for join nodes in queueing networks.
 * 
 * <p>This class configures join nodes that synchronize parallel tasks back into
 * single jobs. Join nodes are the complement to fork nodes in fork-join queueing
 * models, handling task completion synchronization and job reconstruction.</p>
 * 
 * <p>Join node capabilities:
 * <ul>
 *   <li>Task synchronization with various join strategies</li>
 *   <li>Configurable fan-in ratios per job class</li>
 *   <li>Required task completion thresholds</li>
 *   <li>Class-specific joining policies</li>
 * </ul>
 * </p>
 * 
 * <p>Common join strategies include waiting for all tasks, partial synchronization,
 * or quorum-based completion policies.</p>
 * 
 * @see jline.lang.nodes.Join
 * @see ForkNodeParam
 * @see JoinStrategy
 * @since 1.0
 */
public class JoinNodeParam extends NodeParam {
    
    /** Join strategy configuration per job class (e.g., wait-for-all, partial) */
    public Map<JobClass, JoinStrategy> joinStrategy;
    
    /** Fan-in ratio specifying expected number of parallel tasks per job class */
    public Map<JobClass, Double> fanIn;
    
    /** Required number of tasks that must complete before job reconstruction per class */
    public Map<JobClass, Double> joinRequired;

    /**
     * Checks if this join parameter container is empty (no join configuration set).
     * 
     * @return true if no join parameters are configured, false otherwise
     */
    @Override
    public boolean isEmpty() {
        return (joinStrategy == null || joinStrategy.isEmpty()) &&
                (fanIn == null || fanIn.isEmpty()) &&
                (joinRequired == null || joinRequired.isEmpty());
    }
}
