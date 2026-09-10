/*
 * Copyright (c) 2012-2026, QORE Lab, Imperial College London
 * All rights reserved.
 */
package jline.examples.java.advanced;

import java.util.ArrayList;
import java.util.Arrays;
import java.util.List;

import jline.api.pfqn.Pfqn_sdr;
import jline.lang.ClosedClass;
import jline.lang.Network;
import jline.lang.StateDepRouting;
import jline.lang.constant.SchedStrategy;
import jline.lang.nodes.Node;
import jline.lang.nodes.Queue;
import jline.lang.processes.Exp;
import jline.solvers.ctmc.SolverCTMC;
import jline.solvers.nc.SolverNC;

/**
 * Product-form state-dependent routing, Krzesinski (1987), Perform. Eval.
 * 7:125-143.
 *
 * <p>Section 2.5's central server with two peripheral centers and two levels of
 * subnetwork nesting. With C = (-1,-1), d_12 = 1, d_13 = 2 and d_23 = 3 the
 * routing probabilities out of the central server are Table 1 of the paper, and
 * the branch populations are capped at m_2 &lt;= 1 and m_3 &lt;= 3 by the
 * coefficients themselves. When both branches are full the customer is returned
 * to the central server and served again before retrying, which is the busy form
 * of waiting of Section 2.5.</p>
 *
 * <p>The routing violates the usual product form yet keeps one of its own, so
 * SolverNC solves it exactly and agrees with SolverCTMC to machine precision.</p>
 *
 * @see SdrMultiBranchExample for a branch holding several centers
 */
public class StateDepRoutingExample {

    public static void main(String[] args) {
        double[] mu = {1.0, 0.8, 0.5};
        int N = 3;

        Network model = new Network("sdrCentralServer");
        Queue cpu = new Queue(model, "CPU", SchedStrategy.FCFS);
        Queue disk1 = new Queue(model, "Disk1", SchedStrategy.FCFS);
        Queue disk2 = new Queue(model, "Disk2", SchedStrategy.FCFS);
        ClosedClass jobclass = new ClosedClass(model, "Class1", N, cpu, 0);
        cpu.setService(jobclass, new Exp(mu[0]));
        disk1.setService(jobclass, new Exp(mu[1]));
        disk2.setService(jobclass, new Exp(mu[2]));

        model.addLink(cpu, cpu);      // denied entry: the busy form of waiting
        model.addLink(cpu, disk1);
        model.addLink(cpu, disk2);
        model.addLink(disk1, cpu);
        model.addLink(disk2, cpu);
        disk1.setProbRouting(jobclass, cpu, 1.0);
        disk2.setProbRouting(jobclass, cpu, 1.0);

        // Branch index 1 is the complement M-V and stays empty
        List<List<Node>> branches = new ArrayList<List<Node>>();
        branches.add(new ArrayList<Node>());
        branches.add(new ArrayList<Node>(Arrays.asList((Node) disk1)));
        branches.add(new ArrayList<Node>(Arrays.asList((Node) disk2)));
        double[][] d = {{0, 1, 2}, {0, 0, 3}};
        cpu.setStateDepRouting(jobclass, cpu, branches, new int[]{0, 1, 2}, new double[]{-1, -1}, d);

        // Table 1 of the paper, read straight off eq. (10)
        StateDepRouting sdr = model.getStruct(true).sdr;
        Pfqn_sdr.Coeff c = Pfqn_sdr.pfqn_sdrcoeff(sdr);
        System.out.println("Krzesinski (1987) Table 1");
        int[][] table = {{0, 0}, {0, 1}, {1, 0}, {0, 2}, {1, 1}, {1, 2}};
        for (int t = 0; t < table.length; t++) {
            double[] n = {N - table[t][0] - table[t][1], table[t][0], table[t][1]};
            double[] P = Pfqn_sdr.pfqn_sdrprob(c, n);
            System.out.printf("  (m2,m3)=(%d,%d)  P12=%.4f  P13=%.4f  P11=%.4f%n",
                    table[t][0], table[t][1], P[1], P[2], Pfqn_sdr.pfqn_sdrped(P));
        }

        System.out.println();
        new SolverNC(model).getAvgTable().print();
        System.out.println("Section 4 MVA and convolution, same model:");
        jline.solvers.SolverOptions mvaOpt = new SolverNC(model).getOptions();
        mvaOpt.method = "sdr.mva";
        new SolverNC(model, mvaOpt).getAvgTable().print();
        new SolverCTMC(model).getAvgTable().print();
    }
}
