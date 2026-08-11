package jline.lang;

import jline.lang.processes.Distribution;
import jline.util.matrix.Matrix;

import java.io.Serializable;
import java.util.HashMap;
import java.util.Map;

/**
 * Base class for node parameters in queueing network models.
 *
 * <p>NodeParam stores job class-dependent parameters that control the behavior
 * of network nodes (stations). This includes routing weights, outbound links,
 * memory-dependent parameters, and other node-specific configuration data.
 *
 * <p>The parameters are organized by job class, allowing different classes
 * to have different behaviors at the same node. This is essential for
 * multi-class queueing networks where job classes may have different
 * service requirements, routing patterns, or other characteristics.
 *
 * <p>Subclasses extend this base functionality for specific node types
 * such as caches, forks, joins, and service stations.
 *
 * @see jline.lang.nodeparam
 * @see JobClass
 * @see Matrix
 */
public class NodeParam implements Serializable {

    /**
     * Routing weights by job class - controls probabilistic routing decisions
     */
    public Map<JobClass, Matrix> weights;

    /**
     * Outbound link specifications by job class - defines outgoing connections
     */
    public Map<JobClass, Matrix> outlinks;

    /**
     * Weighted round-robin cycle by job class: the outlink node indices with
     * each destination repeated by its integer weight. The WRROBIN pointer
     * holds a POSITION (1..len) in this cycle rather than a node index, so
     * repeated destinations are visited in the correct proportion.
     */
    public Map<JobClass, Matrix> weightedOutlinks;

    /**
     * SQ(d): d, the number of sampled destinations, per job class
     * (mirrors the MATLAB nodeparam.d field).
     */
    public Map<JobClass, Integer> d;

    /**
     * Patience distributions by job class for customer abandonment modeling
     */
    public Map<JobClass, Distribution> patience;

    /**
     * RL value function by job class. For tabular RL (rlStateSize=0) this is a
     * flat representation of an N-dimensional table indexed by per-queue
     * occupancy; the corresponding shape is in {@link #rlValueFunctionShape}.
     * For linear approximation (rlStateSize&gt;0) this is the coefficient row
     * vector applied to a quadratic feature lift of the queue lengths.
     */
    public Map<JobClass, Matrix> rlValueFunction;

    /**
     * Per-class shape of the tabular RL value function table. Each int[] is the
     * size along each dimension; only used when rlStateSize=0.
     */
    public Map<JobClass, int[]> rlValueFunctionShape;

    /**
     * Per-class list of node indices that consult the RL value function. Nodes
     * not in this list fall back to JSQ-style routing.
     */
    public Map<JobClass, int[]> rlNodesNeedAction;

    /**
     * Per-class state-size hint. 0 selects tabular RL, &gt;0 selects linear
     * approximation, missing/negative selects JSQ fallback.
     */
    public Map<JobClass, Integer> rlStateSize;

    /**
     * Constructs an empty NodeParam with initialized parameter maps.
     *
     * <p>All parameter maps are initialized as empty HashMaps, ready to
     * store job class-specific configuration data.
     */
    public NodeParam() {
        weights = new HashMap<JobClass, Matrix>();
        outlinks = new HashMap<JobClass, Matrix>();
        weightedOutlinks = new HashMap<JobClass, Matrix>();
        d = new HashMap<JobClass, Integer>();
        patience = new HashMap<JobClass, Distribution>();
        rlValueFunction = new HashMap<JobClass, Matrix>();
        rlValueFunctionShape = new HashMap<JobClass, int[]>();
        rlNodesNeedAction = new HashMap<JobClass, int[]>();
        rlStateSize = new HashMap<JobClass, Integer>();
    }

    /**
     * Checks if this NodeParam contains no configuration data.
     *
     * @return true if all parameter maps are null or empty, false otherwise
     */
    public boolean isEmpty() {
        return (weights == null || weights.isEmpty()) &&
                (outlinks == null || outlinks.isEmpty()) &&
                (d == null || d.isEmpty()) &&
                (patience == null || patience.isEmpty());
    }
}
