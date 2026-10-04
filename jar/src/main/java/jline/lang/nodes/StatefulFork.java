/*
 * Copyright (c) 2012-2026, QORE Lab, Imperial College London
 * All rights reserved.
 */

package jline.lang.nodes;

import jline.lang.JobClass;
import jline.lang.Network;
import jline.lang.constant.SchedStrategy;
import jline.lang.sections.Buffer;
import jline.lang.sections.Forker;
import jline.lang.sections.ServiceTunnel;

import java.io.Serializable;
import java.util.List;

/**
 * A stateful Fork used only on FJ tag-augmented model copies (see
 * ModelAdapter.fjtag). It reports NodeType.Fork but extends StatefulNode
 * so that it can hold the parent job for one vanishing state before the
 * fork firing (sn.fjsync) emits the sibling tasks. Never used in
 * user-facing models.
 */
public class StatefulFork extends StatefulNode implements Serializable {
    private final int cap;
    private final SchedStrategy schedStrategy;

    public StatefulFork(Network model, String name) {
        super(name);
        List<JobClass> classes = model.getClasses();
        this.cap = Integer.MAX_VALUE;
        this.input = new Buffer(classes);
        this.schedStrategy = SchedStrategy.FORK;
        this.server = new ServiceTunnel();
        this.output = new Forker(classes);
        this.setModel(model);
        model.addNode(this);
    }

    @Override
    public Network getModel() {
        return this.model;
    }
}
