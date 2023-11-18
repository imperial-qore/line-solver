
package jline.lang.layerednetworks;

import jline.util.Matrix;

import java.util.HashMap;
import java.util.Map;

/**
 * Service exposed by a Task object
 */
public class Entry extends LayeredNetworkElement{
    protected Task parent;
    protected Map<Integer,String> replyActivity = new HashMap();//TODO:K/V type
    private final double openArrivalRate;
    protected Matrix scheduling = new Matrix(0,0,0);

    public Entry(LayeredNetwork model, String name) {
        super(name);
        this.openArrivalRate = 0.0;
        model.entries.put(model.entries.size(),this);
        this.model = model;
    }

    public  void on(Task newParent){//TODO
        newParent.addEntry(this);
        this.parent = newParent;
    }
}

