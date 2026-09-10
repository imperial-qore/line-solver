package jline.opt.variables;

import jline.lang.Network;
import jline.lang.nodes.Node;
import jline.lang.nodes.Station;

/**
 * Optimize the number of identical copies (replicas) of a station. N replicas
 * of a load-balanced identical station are represented as a single multiserver
 * station with N times the base server count (an M/M/c-equivalent reduction),
 * keeping topology and node names fixed. Mirrors native-Python
 * {@code StationReplicas}.
 */
public class StationReplicas extends DecisionVariable {

    private final Station station;
    private final int minReplicas;
    private final int maxReplicas;

    public StationReplicas(Station station, int minReplicas, int maxReplicas) {
        this(station, minReplicas, maxReplicas, station.getName() + "_replicas");
    }

    public StationReplicas(Station station, int minReplicas, int maxReplicas, String name) {
        super(name);
        this.station = station;
        this.minReplicas = minReplicas;
        this.maxReplicas = maxReplicas;
    }

    public Station getStation() {
        return station;
    }

    public double[][] getBounds() {
        return unitBounds(1);
    }

    public Object decode(double[] x) {
        double continuous = minReplicas + x[0] * (maxReplicas - minReplicas);
        int v = (int) Math.round(continuous);
        if (v < minReplicas) {
            v = minReplicas;
        }
        if (v > maxReplicas) {
            v = maxReplicas;
        }
        return v;
    }

    public void apply(Network model, Object value) {
        int v = ((Number) value).intValue();
        for (Node node : model.getNodes()) {
            if (node.getName().equals(station.getName()) && node instanceof Station) {
                Station st = (Station) node;
                int base;
                try {
                    base = st.getNumberOfServers();
                } catch (RuntimeException e) {
                    base = 1;
                }
                if (base < 1 || base == Integer.MAX_VALUE) {
                    base = 1;
                }
                st.setNumberOfServers(base * v);
                return;
            }
        }
    }

    public String getVariableType() {
        return "station_replicas";
    }
}
