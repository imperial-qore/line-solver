package jline.lang.sections;

import java.io.Serializable;
import java.util.ArrayList;
import java.util.List;
import jline.lang.*;
import jline.lang.constant.RoutingStrategy;
import jline.lang.distributions.*;
import jline.lang.nodes.*;
import jline.lang.sections.*;

/**
 * Output section that routes jobs
 */
public class Dispatcher extends OutputSection implements Serializable {
    public Dispatcher(List<JobClass> customerClasses) {
        super("Dispatcher");
    }

    public void initDispatcherJobClasses(List<JobClass> customerClasses){
        for(int r = 0; r < customerClasses.size(); r++){
            if(r < this.outputStrategies.size()){
                this.outputStrategies.set(r, new OutputStrategy(customerClasses.get(r), RoutingStrategy.RAND));
            } else {
                this.outputStrategies.add(new OutputStrategy(customerClasses.get(r), RoutingStrategy.RAND));
            }
        }
    }
}
