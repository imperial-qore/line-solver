package jline.lang.nodeparam;

import jline.lang.NodeParam;
import jline.lang.constant.HeteroSchedPolicy;
import jline.lang.constant.ProcessType;
import jline.util.matrix.Matrix;
import jline.util.matrix.MatrixCell;

import java.util.List;
import java.util.Map;

/**
 * Parameter container for service nodes in queueing networks.
 *
 * <p>This class stores configuration parameters for service nodes that process
 * jobs according to externally defined service specifications. Service nodes
 * can reference external files or configurations that define their behavior.</p>
 *
 * <p>It is also the home of the per-station heterogeneous-server fields. These
 * are ragged, node-type-conditional structures (server-type lists, compatibility
 * matrices, per-type/per-class rates and processes) and therefore live in the
 * {@code nodeparam} container rather than as flat {@code NetworkStruct} root
 * fields. They are populated only for {@link jline.lang.nodes.Queue} stations
 * that declare server types; for all other stations {@link #nservertypes} stays
 * 0 (homogeneous).</p>
 *
 * @see jline.lang.nodes.ServiceNode
 * @since 1.0
 */
public class ServiceNodeParam extends NodeParam {

    /** List of file names or identifiers for external service definitions */
    public List<String> fileName;

    /** Number of server types at this station; 0 indicates a homogeneous queue. */
    public int nservertypes = 0;

    /** Server type names; servertypenames.get(t) is the name of server type t. */
    public List<String> servertypenames;

    /** Number of servers per server type: Matrix (nTypes x 1). */
    public Matrix serverspertype;

    /** Server-class compatibility matrix: Matrix (nTypes x K), 1.0 if type t serves class r. */
    public Matrix servercompat;

    /** Heterogeneous service rates: heterorates.get(serverTypeId).get(classId) -> rate. */
    public Map<Integer, Map<Integer, Double>> heterorates;

    /** Heterogeneous service processes: heteroproc.get(serverTypeId).get(classId) -> PH matrices. */
    public Map<Integer, Map<Integer, MatrixCell>> heteroproc;

    /** Heterogeneous process types: heteroprocid.get(serverTypeId).get(classId) -> ProcessType. */
    public Map<Integer, Map<Integer, ProcessType>> heteroprocid;

    /** Heterogeneous scheduling policy for this station. */
    public HeteroSchedPolicy heteroschedpolicy;

    /**
     * Servers seized at once by a job, per class: Matrix (1 x K), all ones unless
     * some class declares job parallelism. Null when the station declares none.
     */
    public Matrix serverparallelism;

    /**
     * Checks if this service parameter container is empty.
     *
     * @return true if no file names are specified, false otherwise
     */
    @Override
    public boolean isEmpty() {
        return (fileName == null);
    }
}
