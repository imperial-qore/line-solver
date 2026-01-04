package jline.lang.nodeparam;

import jline.lang.NodeParam;

import java.util.List;

/**
 * Parameter container for service nodes in queueing networks.
 * 
 * <p>This class stores configuration parameters for service nodes that process
 * jobs according to externally defined service specifications. Service nodes
 * can reference external files or configurations that define their behavior.</p>
 * 
 * @see jline.lang.nodes.ServiceNode
 * @since 1.0
 */
public class ServiceNodeParam extends NodeParam {
    
    /** List of file names or identifiers for external service definitions */
    public List<String> fileName;

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
