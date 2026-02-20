/*
 * Copyright (c) 2012-2026, QORE Lab, Imperial College London
 * All rights reserved.
 */

package jline.lang.nodes;

import jline.lang.JobClass;
import jline.lang.Network;
import jline.lang.constant.JoinStrategy;
import jline.lang.sections.Dispatcher;
import jline.lang.sections.Joiner;
import jline.lang.sections.ServiceTunnel;

import java.io.Serializable;
import java.util.List;

/**
 * A node that reassembles a set of sibling tasks into the original parent job
 */
public class Join extends Station implements Serializable {
    public Node joinOf;

    public Join(Network model) {
        this(model, "Join");
    }

    public Join(Network model, String name) {
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

    @Override
    public Network getModel() {
        return this.model;
    }

    public void initJoinJobClasses() {
        List<JobClass> classes = this.model.getClasses();
        this.input = new Joiner(classes);
        this.output = new Dispatcher(classes);
    }

    public void setRequired(JobClass jobClass, double njobs) {
        Joiner joiner = (Joiner) this.input;
        joiner.setRequired(jobClass, njobs);
        this.input = joiner;
    }

    public void setStrategy(JobClass jobClass, JoinStrategy joinStrategy) {
        Joiner joiner = (Joiner) this.input;
        joiner.setStrategy(jobClass, joinStrategy);
        this.input = joiner;
    }
}
