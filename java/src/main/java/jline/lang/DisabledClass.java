package jline.lang;

import jline.lang.constant.JobClassType;
import jline.lang.constant.JoinStrategy;
import jline.lang.constant.RoutingStrategy;
import jline.lang.nodes.Join;
import jline.lang.nodes.Node;
import jline.lang.nodes.Station;

import java.io.Serializable;

/**
 * Class of jobs that perpetually loop at a given station
 */
public class DisabledClass extends JobClass implements Serializable {

    public DisabledClass(Network model, String name,  Station refstat) throws Exception {
        super("Disabled", name);
        this.type = JobClassType.DISABLED;
        this.priority = 0;
        model.addJobClass(this);
        this.setReferenceStation(refstat);
        for (int i=0; i<model.getNumberOfNodes(); i++){
            Node node_i = model.getNodeByIndex(i);
            node_i.setRouting(this, RoutingStrategy.DISABLED);
            if (node_i instanceof Join){
                ((Join) node_i).setStrategy(this, JoinStrategy.STD);
                ((Join) node_i).setRequired(this, -1);
            }
        }
    }

}
