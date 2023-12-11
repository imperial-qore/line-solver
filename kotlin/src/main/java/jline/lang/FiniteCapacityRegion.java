package jline.lang;

import jline.lang.nodes.Node;

import java.util.List;
import java.util.Map;

/**
 * Collection of stations with constraints on the number of admitted jobs
 */
public class FiniteCapacityRegion {

    public static final int UNBOUNDED = -1;

    public List<Node> nodes;
    public List<JobClass> classes;
    private int globalMaxJobs;
    private int globalMaxMemory;
    public Map<JobClass, Integer> classMaxJobs;
    public Map<JobClass, Integer> classMaxMemory;
    public Map<JobClass, Boolean> dropRule;
    public Map<JobClass, Integer> classSize;

    public FiniteCapacityRegion(List<Node> nodes, List<JobClass> classes) {
        this.nodes = nodes;
        this.classes = classes;
        this.globalMaxJobs = UNBOUNDED;
        this.globalMaxMemory = UNBOUNDED;
    }

    public void setGlobalMaxJobs(int njobs) {
        this.globalMaxJobs = njobs;
    }

    public void setGlobalMaxMemory(int memlim) {
        this.globalMaxMemory = memlim;
    }

    public List<Node> getNodes() {
        return nodes;
    }

    public int getGlobalMaxJobs() {
        return globalMaxJobs;
    }

    public int getGlobalMaxMemory() {
        return globalMaxMemory;
    }
}

