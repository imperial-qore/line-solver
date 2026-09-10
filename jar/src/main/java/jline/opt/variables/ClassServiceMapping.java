package jline.opt.variables;

import jline.lang.JobClass;
import jline.lang.Network;
import jline.lang.RoutingMatrix;
import jline.lang.nodes.Node;
import jline.lang.nodes.Station;
import jline.util.matrix.Matrix;

import java.util.HashSet;
import java.util.List;
import java.util.Set;

/**
 * Optimize the class-to-station mapping by rerouting a job class through a
 * selected candidate station and bypassing the others, preserving default
 * routing of every other class. Per-node routing is reconstructed from the
 * connection matrix (uniform split over outgoing links); for the mapped class
 * the non-selected candidate columns are zeroed so no traffic enters them.
 * Mirrors native-Python {@code ClassServiceMapping}.
 */
public class ClassServiceMapping extends DecisionVariable {

    private final JobClass jobclass;
    private final List<Station> stations;

    public ClassServiceMapping(JobClass jobclass, List<Station> stations) {
        this(jobclass, stations, jobclass.getName() + "_mapping");
    }

    public ClassServiceMapping(JobClass jobclass, List<Station> stations, String name) {
        super(name);
        this.jobclass = jobclass;
        this.stations = stations;
    }

    public JobClass getJobClass() {
        return jobclass;
    }

    public List<Station> getStations() {
        return stations;
    }

    public double[][] getBounds() {
        return unitBounds(1);
    }

    public Object decode(double[] x) {
        int n = stations.size();
        int idx = (int) Math.floor(x[0] * n);
        return Math.min(idx, n - 1);
    }

    public void apply(Network model, Object value) {
        if (stations.isEmpty()) {
            return;
        }
        int v = ((Number) value).intValue();
        if (v < 0) {
            v = 0;
        }
        if (v > stations.size() - 1) {
            v = stations.size() - 1;
        }

        List<Node> nodes = model.getNodes();
        int n = nodes.size();
        Matrix conn = connectionMatrix(model);

        String selectedName = stations.get(v).getName();
        Set<Integer> blocked = new HashSet<Integer>();
        for (int k = 0; k < stations.size(); k++) {
            if (!stations.get(k).getName().equals(selectedName)) {
                int idx = indexOfNode(nodes, stations.get(k).getName());
                if (idx >= 0) {
                    blocked.add(idx);
                }
            }
        }

        JobClass mapped = resolveClass(model, jobclass);
        if (mapped == null) {
            return;
        }

        RoutingMatrix rt = model.initRoutingMatrix();
        for (JobClass c : model.getClasses()) {
            // per-class adjacency copy
            double[][] adj = new double[n][n];
            for (int i = 0; i < n; i++) {
                for (int j = 0; j < n; j++) {
                    adj[i][j] = conn.get(i, j);
                }
            }
            if (c.getName().equals(mapped.getName())) {
                for (int b : blocked) {
                    for (int i = 0; i < n; i++) {
                        adj[i][b] = 0.0;
                    }
                }
            }
            for (int i = 0; i < n; i++) {
                double outDegree = 0.0;
                for (int j = 0; j < n; j++) {
                    outDegree += adj[i][j];
                }
                if (outDegree <= 0) {
                    continue;
                }
                for (int j = 0; j < n; j++) {
                    if (adj[i][j] > 0) {
                        rt.set(c, c, nodes.get(i), nodes.get(j), adj[i][j] / outDegree);
                    }
                }
            }
        }
        model.link(rt);
    }

    public String getVariableType() {
        return "class_mapping";
    }
}
