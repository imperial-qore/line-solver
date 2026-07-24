package jline.opt.variables;

import jline.lang.Network;
import jline.lang.nodes.Node;
import jline.lang.nodes.Station;

/**
 * Optimize the number of servers at a station. Encodes an integer server count
 * in [minServers, maxServers] via a rounded linear map. Mirrors native-Python
 * {@code ServerAllocation}.
 */
public class ServerAllocation extends DecisionVariable {

    private final Station station;
    private final int minServers;
    private final int maxServers;

    public ServerAllocation(Station station, int minServers, int maxServers) {
        this(station, minServers, maxServers, station.getName() + "_servers");
    }

    public ServerAllocation(Station station, int minServers, int maxServers, String name) {
        super(name);
        this.station = station;
        this.minServers = minServers;
        this.maxServers = maxServers;
    }

    public Station getStation() {
        return station;
    }

    public int getMinServers() {
        return minServers;
    }

    public int getMaxServers() {
        return maxServers;
    }

    public double[][] getBounds() {
        return unitBounds(1);
    }

    public Object decode(double[] x) {
        double continuous = minServers + x[0] * (maxServers - minServers);
        int v = (int) Math.round(continuous);
        if (v < minServers) {
            v = minServers;
        }
        if (v > maxServers) {
            v = maxServers;
        }
        return v;
    }

    public void apply(Network model, Object value) {
        int v = ((Number) value).intValue();
        for (Node node : model.getNodes()) {
            if (node.getName().equals(station.getName()) && node instanceof Station) {
                ((Station) node).setNumberOfServers(v);
                return;
            }
        }
    }

    public String getVariableType() {
        return "server_allocation";
    }
}
