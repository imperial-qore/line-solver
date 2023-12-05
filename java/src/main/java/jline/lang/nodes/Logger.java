package jline.lang.nodes;

import java.io.Serializable;

import jline.lang.Network;
import jline.lang.sections.Forker;

/**
 * A node that logs passage of jobs
 */
public class Logger extends Node implements Serializable {

    protected Network model;
    protected String fileName;

    public Logger(Network model,  String name, String logfile) {
        // TODO: this class doesn't work.
        super("Logger");
        //this.setName(name);
        model.addNode(this);
        this.model = model;
        //this.fileName = logfile;
    }

    @Override
    public Network getModel() {
        return this.model;
    }

}
