package jline.lang.nodeparam;

import jline.lang.NodeParam;

import java.util.List;

/**
 * Parameter container for replayer nodes in queueing networks.
 * 
 * <p>This class configures replayer nodes that read and replay job traces from
 * external files. Replayer nodes are essential for trace-driven simulation,
 * workload reproduction, and validation against real-world data patterns.</p>
 * 
 * <p>Replay capabilities include:
 * <ul>
 *   <li>Reading job arrival traces from files</li>
 *   <li>Reproducing empirical interarrival time patterns</li>
 *   <li>Supporting multiple trace file formats</li>
 *   <li>Enabling workload characterization studies</li>
 *   <li>Validating analytical models against real data</li>
 * </ul>
 * </p>
 * 
 * @see jline.lang.processes.Replayer
 * @see jline.lang.nodes.Source
 * @since 1.0
 */
public class ReplayerNodeParam extends NodeParam {
    
    /** List of trace file names to replay */
    public List<String> fileName;
    
    /** Base file path for trace files */
    public String filePath;

    /**
     * Checks if this replayer parameter container is empty (no trace files configured).
     * 
     * @return true if no trace files are specified, false otherwise
     */
    @Override
    public boolean isEmpty() {
        return (fileName == null || fileName.isEmpty()) &&
                filePath == null;
    }
}
