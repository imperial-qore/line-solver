package jline.lang.nodes;

import java.io.Serializable;
import java.util.List;

import jline.lang.*;
import jline.lang.constant.JoinStrategy;
import jline.lang.sections.Dispatcher;
import jline.lang.sections.Joiner;
import jline.lang.sections.ServiceTunnel;

/**
 * A node that reassembles a set of sibling tasks into the original parent job
 */
public class Join extends Station implements Serializable {
	public Node joinOf;

    public Join(Network model){
        this(model, "Join");
    }

    public Join(Network model, String name){
        this(model, name, null);
    }

    public Join(Network model, String name, Node fork) {
        super(name);
        this.server = new ServiceTunnel();
        this.numberOfServers = Integer.MAX_VALUE;
        this.model = model;
        this.joinOf = fork;
        initJoinJobClasses();
        model.addNode(this);
    }

    public void initJoinJobClasses() {
        List<JobClass> classes = this.model.getClasses();
        this.input = new Joiner(classes);
        this.output = new Dispatcher(classes);
    }

    @Override
    public Network getModel() {
        return this.model;
    }

    public void setStrategy(JobClass jobClass, JoinStrategy joinStrategy){
        Joiner joiner = (Joiner) this.input;
        joiner.setStrategy(jobClass, joinStrategy);
        this.input = joiner;
    }

    public void setRequired(JobClass jobClass, double njobs){
        Joiner joiner = (Joiner) this.input;
        joiner.setRequired(jobClass, njobs);
        this.input = joiner;
    }
}
