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

import static jline.io.InputOutput.line_error;
import static jline.io.InputOutput.line_warning;
import static jline.io.InputOutput.mfilename;

/**
 * A node that forks an incoming job into a set of sibling tasks.
 *
 * <p>Fork nodes split incoming jobs into parallel sibling tasks that are sent to multiple
 * output destinations. The number of tasks per output link can be configured using
 * {@link #setTasksPerLink(int)}.</p>
 *
 * <p><b>Solver compatibility for tasksPerLink &gt; 1:</b></p>
 * <ul>
 *   <li><b>SolverLDES</b>: Fully supported - correctly simulates multiple tasks per link</li>
 *   <li><b>SolverJMT</b>: Fully supported - simulation handles multiple tasks correctly</li>
 *   <li><b>SolverMVA (H-T method)</b>: Not supported - throws error when tasksPerLink &gt; 1</li>
 *   <li><b>SolverMVA (MMT method)</b>: Supported - analytical approximation</li>
 * </ul>
 *
 * @see Join
 * @see Forker
 */
public class Fork extends Node implements Serializable {
    private final int cap;
    private final SchedStrategy schedStrategy;

    /**
     * Creates a new fork node with default name "Fork".
     * 
     * @param model the network model to add this fork node to
     */
    public Fork(Network model) {
        this(model, "Fork");
    }

    /**
     * Creates a new fork node with the specified name.
     * Initializes the fork with buffer input, service tunnel server, and forker output sections.
     * 
     * @param model the network model to add this fork node to
     * @param name the name for this fork node
     */
    public Fork(Network model, String name) {
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

    /**
     * Sets the number of tasks sent out on each outgoing link.
     *
     * <p>By default, a Fork node sends exactly one task per outgoing link. This method allows
     * configuring the Fork to send multiple identical tasks on each link. The total number of
     * tasks created will be: (number of outgoing links) × tasksPerLink.</p>
     *
     * <p><b>Important:</b> Values greater than 1 are only fully supported by simulation-based
     * solvers (SolverLDES, SolverJMT). Analytical solvers may produce errors or inaccurate results.</p>
     *
     * @param nTasks the number of tasks to send on each outgoing link (default: 1)
     * @see Fork class documentation for solver compatibility details
     */
    public void setTasksPerLink(int nTasks) {
        if (nTasks != 1) {
            line_warning(mfilename(new Object() {}), "The setTasksPerLink feature is experimental and results may be inaccurate for analytical solvers.");
        }
        Forker f = (Forker) this.output;
        f.tasksPerLink = nTasks;
    }

    /**
     * Sets the tasks per link for one class only, leaving every other class on
     * the node-wide value.
     *
     * @param jobClass the class the count applies to
     * @param nTasks   the number of tasks emitted on each outgoing link
     */
    public void setTasksPerLink(JobClass jobClass, double nTasks) {
        if (nTasks != 1) {
            line_warning(mfilename(new Object() {}), "The setTasksPerLink feature is experimental and results may be inaccurate for analytical solvers.");
        }
        Forker f = (Forker) this.output;
        f.tasksPerLinkByDest.add(new Forker.ForkOverride("", jobClass.getIndex(), nTasks));
    }

    /**
     * Makes the number of tasks emitted on each outgoing link RANDOM.
     *
     * <p>The degree is redrawn independently for every link and every forked
     * job, which is the variable forking level of JMT's JobsPerLinkDis. Exact
     * under SolverJMT and SolverLDES, which draw it at the fork epoch; the
     * analytical solvers see E[dist].</p>
     *
     * @param jobClass the class the distribution applies to
     * @param dist     a DiscreteSampler over the tasks-per-link support
     */
    public void setTasksPerLinkDistribution(JobClass jobClass,
                                            jline.lang.processes.DiscreteSampler dist) {
        setTasksPerLinkDistribution(jobClass, dist, null);
    }

    /**
     * Makes the number of tasks emitted on ONE outgoing link random.
     *
     * @param jobClass the class the distribution applies to
     * @param dist     a DiscreteSampler over the tasks-per-link support
     * @param destNode the destination node of the link, or null for every link
     */
    public void setTasksPerLinkDistribution(JobClass jobClass,
                                            jline.lang.processes.DiscreteSampler dist,
                                            Node destNode) {
        if (dist == null) {
            line_error(mfilename(new Object() {}), "The tasks-per-link distribution must be a DiscreteSampler.");
        }
        Forker f = (Forker) this.output;
        String dest = (destNode == null) ? "" : destNode.getName();
        f.tasksPerLinkDist.add(new Forker.ForkOverride(dest, jobClass.getIndex(), dist.getMean(), dist));
    }

    /**
     * Activates an outgoing branch only with probability {@code prob}.
     *
     * <p>Branches are activated independently, so the number of siblings a job
     * produces is random even when the tasks per link are deterministic. The
     * matched Join must be told what to wait for: under a standard join a job
     * that skipped a branch would block forever, so any probability below one
     * requires JoinStrategy.PARTIAL or a quorum.</p>
     *
     * @param jobClass the class the probability applies to
     * @param destNode the destination node of the branch
     * @param prob     the activation probability, in [0,1]
     */
    public void setBranchProbability(JobClass jobClass, Node destNode, double prob) {
        if (prob < 0.0 || prob > 1.0) {
            line_error(mfilename(new Object() {}), "A branch activation probability must lie in [0,1].");
        }
        Forker f = (Forker) this.output;
        f.branchProb.add(new Forker.ForkOverride(destNode.getName(), jobClass.getIndex(), prob));
    }
}
