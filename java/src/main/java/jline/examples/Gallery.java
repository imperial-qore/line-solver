package jline.examples;

import jline.lang.*;
import jline.lang.constant.SchedStrategy;
import jline.lang.nodes.*;
import jline.lang.processes.MAP;
import jline.solvers.jmt.SolverJMT;

import javax.xml.parsers.ParserConfigurationException;
import java.util.List;
import java.util.ArrayList;
import jline.lang.distributions.*;

/**
 * Getting started examples
 */
public class Gallery {
    public static Network gallery_aphm1() {
        Network model = new Network("APH/M/1");
        Source source = new Source(model, "mySource");
        Queue queue = new Queue(model, "myQueue", SchedStrategy.FCFS);
        Sink sink = new Sink(model, "mySink");
        JobClass oclass = new OpenClass(model, "myClass");
        source.setArrival(oclass, APH.fitCentral(1,0.99,1.999));
        queue.setService(oclass, new Exp(2));
        model.link(Network.serialRouting(source,queue,sink));
        return model;
    }

    public static Network gallery_coxm1() {
        Network model = new Network("Cox/M/1");
        Source source = new Source(model, "mySource");
        Queue queue = new Queue(model, "myQueue", SchedStrategy.FCFS);
        Sink sink = new Sink(model, "mySink");
        JobClass oclass = new OpenClass(model, "myClass");
        source.setArrival(oclass, Coxian.fitCentral(1,0.99,1.999));
        queue.setService(oclass, new Exp(2));
        model.link(Network.serialRouting(source,queue,sink));
        return model;
    }

    public static Network gallery_cqn() {
        Network model = new Network("model");
        // TODO
        return model;
    }

    public static Network gallery_cqn_multiclass() {
        Network model = new Network("model");
        // TODO
        return model;
    }

    public static Network gallery_dm1() {
        Network model = new Network("Det/M/1");
        Source source = new Source(model, "mySource");
        Queue queue = new Queue(model, "myQueue", SchedStrategy.FCFS);
        Sink sink = new Sink(model, "mySink");
        JobClass oclass = new OpenClass(model, "myClass");
        source.setArrival(oclass, new Det(1));
        queue.setService(oclass, new Exp(2));
        model.link(Network.serialRouting(source,queue,sink));
        return model;
    }

    public static Network gallery_erldk() {
        return gallery_erldk(2);
    }

    public static Network gallery_erldk(int k) {
        Network model = new Network("Erl/Det/k");
        Source source = new Source(model, "mySource");
        Queue queue = new Queue(model, "myQueue", SchedStrategy.FCFS);
        Sink sink = new Sink(model, "mySink");
        JobClass oclass = new OpenClass(model, "myClass");
        source.setArrival(oclass, Erlang.fitMeanAndOrder(1.0,5));
        queue.setService(oclass, new Det(2.0/k));
        queue.setNumberOfServers(k);
        model.link(Network.serialRouting(source,queue,sink));
        return model;
    }

    public static Network gallery_erlerl1_reentrant() {
        Network model = new Network("model");
        // TODO
        return model;
    }

    public static Network gallery_erlm1() {
        Network model = new Network("Erl/M/1");
        Source source = new Source(model, "mySource");
        Queue queue = new Queue(model, "myQueue", SchedStrategy.FCFS);
        Sink sink = new Sink(model, "mySink");
        JobClass oclass = new OpenClass(model, "myClass");
        source.setArrival(oclass, Erlang.fitMeanAndOrder(1.0,5));
        queue.setService(oclass, new Exp(2));
        model.link(Network.serialRouting(source,queue,sink));
        return model;
    }

    public static Network gallery_erlm1_ps() {
        Network model = new Network("Erl/M/1");
        Source source = new Source(model, "mySource");
        Queue queue = new Queue(model, "myQueue", SchedStrategy.PS);
        Sink sink = new Sink(model, "mySink");
        JobClass oclass = new OpenClass(model, "myClass");
        source.setArrival(oclass, Erlang.fitMeanAndOrder(1.0,5));
        queue.setService(oclass, new Exp(2));
        model.link(Network.serialRouting(source,queue,sink));
        return model;
    }

    public static Network gallery_erlm1_reentrant() {
        Network model = new Network("model");
        // TODO
        return model;
    }

    public static Network gallery_gamm1() {
        Network model = new Network("Gam/M/1");
        Source source = new Source(model, "mySource");
        Queue queue = new Queue(model, "myQueue", SchedStrategy.FCFS);
        Sink sink = new Sink(model, "mySink");
        JobClass oclass = new OpenClass(model, "myClass");
        source.setArrival(oclass, Gamma.fitMeanAndSCV(1.0,1/5));
        queue.setService(oclass, new Exp(2));
        model.link(Network.serialRouting(source,queue,sink));
        return model;
    }

    public static Network gallery_hyperl1_feedback() {
        Network model = new Network("model");
        // TODO
        return model;
    }

    public static Network gallery_hyperl1_reentrant() {
        Network model = new Network("model");
        // TODO
        return model;
    }

    public static Network gallery_hyperlk() {
        Network model = new Network("model");
        // TODO
        return model;
    }

    public static Network gallery_hyphyp1_linear() {
        Network model = new Network("model");
        // TODO
        return model;
    }

    public static Network gallery_hyphyp1_reentrant() {
        Network model = new Network("model");
        // TODO
        return model;
    }

    public static Network gallery_hyphyp1_tandem() {
        Network model = new Network("model");
        // TODO
        return model;
    }

    public static Network gallery_hypm1() {
        Network model = new Network("H2/M/1");
        Source source = new Source(model, "mySource");
        Queue queue = new Queue(model, "myQueue", SchedStrategy.FCFS);
        Sink sink = new Sink(model, "mySink");
        JobClass oclass = new OpenClass(model, "myClass");
        source.setArrival(oclass, HyperExp.fitMeanAndSCV(1,64));
        queue.setService(oclass, new Exp(2));
        model.link(Network.serialRouting(source,queue,sink));
        return model;
    }

    public static Network gallery_hypm1_reentrant() {
        Network model = new Network("model");
        // TODO
        return model;
    }

    public static Network gallery_mapm1() {
        return gallery_mapm1(MAP.rand(2));
    }

    public static Network gallery_mapm1(MAP map) {
        Network model = new Network("MAP/M/1");
        Source source = new Source(model, "mySource");
        Queue queue = new Queue(model, "myQueue", SchedStrategy.FCFS);
        Sink sink = new Sink(model, "mySink");
        JobClass oclass = new OpenClass(model, "myClass");
        source.setArrival(oclass, map);
        queue.setService(oclass, new Exp(2));
        model.link(Network.serialRouting(source,queue,sink));
        return model;
    }

    public static Network gallery_mapmk() {
        return gallery_mapmk(MAP.rand(2), 2);
    }

    public static Network gallery_mapmk(MAP map) {
        return gallery_mapmk(map, 2);
    }

    public static Network gallery_mapmk(MAP map, int k) {
        Network model = new Network("MAP/M/k");
        Source source = new Source(model, "mySource");
        Queue queue = new Queue(model, "myQueue", SchedStrategy.FCFS);
        Sink sink = new Sink(model, "mySink");
        JobClass oclass = new OpenClass(model, "myClass");
        source.setArrival(oclass, map);
        queue.setService(oclass, new Exp(2));
        queue.setNumberOfServers(k);
        model.link(Network.serialRouting(source,queue,sink));
        return model;
    }
    public static Network gallery_mdk() {
        return gallery_mdk(2);
    }

    public static Network gallery_mdk(int k) {
        Network model = new Network("M/D/k");
        Source source = new Source(model, "mySource");
        Queue queue = new Queue(model, "myQueue", SchedStrategy.FCFS);
        Sink sink = new Sink(model, "mySink");
        JobClass oclass = new OpenClass(model, "myClass");
        source.setArrival(oclass, Exp.fitMean(1));
        queue.setService(oclass, new Det(2.0/k));
        queue.setNumberOfServers(k);
        model.link(Network.serialRouting(source,queue,sink));
        return model;
    }

    public static Network gallery_merl1() {
        Network model = new Network("M/E/1");
        Source source = new Source(model, "mySource");
        Queue queue = new Queue(model, "myQueue", SchedStrategy.FCFS);
        Sink sink = new Sink(model, "mySink");
        JobClass oclass = new OpenClass(model, "myClass");
        source.setArrival(oclass, new Exp(1));
        queue.setService(oclass, Erlang.fitMeanAndOrder(0.5,2));
        model.link(Network.serialRouting(source,queue,sink));
        return model;
    }

    public static Network gallery_merl1_linear() {
        Network model = new Network("model");
        // TODO
        return model;
    }

    public static Network gallery_merl1_reentrant() {
        Network model = new Network("model");
        // TODO
        return model;
    }

    public static Network gallery_merl1_tandem() {
        Network model = new Network("model");
        // TODO
        return model;
    }

    public static Network gallery_merlk() {
        return gallery_merlk(2);
    }

    public static Network gallery_merlk(int k) {
        Network model = new Network("M/E/1");
        Source source = new Source(model, "mySource");
        Queue queue = new Queue(model, "myQueue", SchedStrategy.FCFS);
        Sink sink = new Sink(model, "mySink");
        JobClass oclass = new OpenClass(model, "myClass");
        source.setArrival(oclass, new Exp(1));
        queue.setService(oclass, Erlang.fitMeanAndOrder(0.5,2));
        queue.setNumberOfServers(k);
        model.link(Network.serialRouting(source,queue,sink));
        return model;
    }

    public static Network gallery_mhyp1() {
        Network model = new Network("M/H/1");
        Source source = new Source(model, "mySource");
        Queue queue = new Queue(model, "myQueue", SchedStrategy.FCFS);
        Sink sink = new Sink(model, "mySink");
        JobClass oclass = new OpenClass(model, "myClass");
        source.setArrival(oclass, new Exp(1));
        queue.setService(oclass, HyperExp.fitMeanAndSCV(0.5,4));
        model.link(Network.serialRouting(source,queue,sink));
        return model;
    }

    public static Network gallery_mhyp1_linear() {
        Network model = new Network("model");
        // TODO
        return model;
    }

    public static Network gallery_mhyp1_reentrant() {
        Network model = new Network("model");
        // TODO
        return model;
    }

    public static Network gallery_mhyp1_tandem() {
        Network model = new Network("model");
        // TODO
        return model;
    }

    public static Network gallery_mhypk() {
        return gallery_mhypk(2);
    }

    public static Network gallery_mhypk(int k) {
        Network model = new Network("M/H/k");
        Source source = new Source(model, "mySource");
        Queue queue = new Queue(model, "myQueue", SchedStrategy.FCFS);
        Sink sink = new Sink(model, "mySink");
        JobClass oclass = new OpenClass(model, "myClass");
        source.setArrival(oclass, new Exp(1));
        queue.setService(oclass, HyperExp.fitMeanAndSCV(0.5,4));
        queue.setNumberOfServers(k);
        model.link(Network.serialRouting(source,queue,sink));
        return model;
    }

    public static Network gallery_mm1() {
        Network model = new Network("M/M/1");
        Source source = new Source(model, "mySource");
        Queue queue = new Queue(model, "myQueue", SchedStrategy.FCFS);
        Sink sink = new Sink(model, "mySink");
        JobClass oclass = new OpenClass(model, "myClass");
        source.setArrival(oclass, new Exp(1));
        queue.setService(oclass, new Exp(2));
        model.link(Network.serialRouting(source,queue,sink));
        return model;
    }

    public static Network gallery_mm1_feedback() {
        Network model = new Network("model");
        // TODO
        return model;
    }

    public static Network gallery_mm1_linear() {
        Network model = new Network("model");
        // TODO
        return model;
    }

    public static Network gallery_mm1_linear(Integer n) {
        return gallery_mm1_linear(n, 0.9);
    }

    public static Network gallery_mm1_linear(Integer n, Double Umax) {

        Network model = new Network("M/M/1-Linear");

        // Block 1: nodes
        List<Node> nodes = new ArrayList<>();
        nodes.add(new Source(model, "mySource"));
        for (int i = 1; i <= n; i++) {
            nodes.add(new Queue(model, "Queue" + i, SchedStrategy.FCFS));
        }
        nodes.add(new Sink(model, "mySink"));

        // Block 2: classes
        OpenClass oclass = new OpenClass(model, "myClass");
        ((Source)nodes.get(0)).setArrival(oclass, new Exp(1));

        double[] means = new double[n];
        int half = n / 2;
        for (int i = 0; i < half; i++) {
            means[i] = Umax;
        }
        if (n % 2 == 0) {
            for (int i = 0; i < half; i++) {
                means[n - 1 - i] = means[i];
            }
        } else {
            for (int i = 0; i < half; i++) {
                means[n - 1 - i] = means[i];
            }
            means[half] = Umax;
        }

        for (int i = 1; i <= n; i++) {
            ((Queue)nodes.get(i)).setService(oclass, Exp.fitMean(means[i - 1]));
        }

        // Block 3: topology
        model.link(model.serialRouting(nodes));

        return model;
    }

    public static Network gallery_mm1_multiclass() {
        Network model = new Network("model");
        // TODO
        return model;
    }

    public static Network gallery_mm1_prio() {
        Network model = new Network("M[2]/M[2]/1");
        Source source = new Source(model, "mySource");
        Queue queue = new Queue(model, "myQueue", SchedStrategy.HOL);
        Sink sink = new Sink(model, "mySink");
        JobClass oclass1 = new OpenClass(model, "myClass1", 1);
        source.setArrival(oclass1, new Exp(1));
        queue.setService(oclass1, new Exp(4));
        JobClass oclass2 = new OpenClass(model, "myClass2", 0);
        source.setArrival(oclass2, new Exp(0.5));
        queue.setService(oclass2, new Exp(4));
        RoutingMatrix P = model.initRoutingMatrix();
        P.set(oclass1, Network.serialRouting(source,queue,sink));
        P.set(oclass2, Network.serialRouting(source,queue,sink));
        model.link(P);
        return model;
    }

    public static Network gallery_mm1_ps() {
        Network model = new Network("M/M/1-PS");
        Source source = new Source(model, "mySource");
        Queue queue = new Queue(model, "myQueue", SchedStrategy.PS);
        Sink sink = new Sink(model, "mySink");
        JobClass oclass = new OpenClass(model, "myClass");
        source.setArrival(oclass, new Exp(1));
        queue.setService(oclass, new Exp(2));
        model.link(Network.serialRouting(source,queue,sink));
        return model;
    }

    public static Network gallery_mm1_ps_feedback() {
        Network model = new Network("model");
        // TODO
        return model;
    }

    public static Network gallery_mm1_ps_multiclass() {
        Network model = new Network("model");
        // TODO
        return model;
    }

    public static Network gallery_mm1_ps_reentrant() {
        Network model = new Network("model");
        // TODO
        return model;
    }

    public static Network gallery_mm1_reentrant() {
        Network model = new Network("model");
        // TODO
        return model;
    }

    public static Network gallery_mm1_tandem() {
        return gallery_mm1_linear(2, 0.9);
    }

    public static Network gallery_mm1_tandem(Double Umax) {
        return gallery_mm1_linear(2, Umax);
    }

    public static Network gallery_mm1_tandem_multiclass() {
        Network model = new Network("model");
        // TODO
        return model;
    }

    public static Network gallery_mmap1() {
        Network model = new Network("model");
        // TODO
        return model;
    }

    public static Network gallery_mmap1_multiclass() {
        Network model = new Network("model");
        // TODO
        return model;
    }

    public static Network gallery_mmapk() {
        Network model = new Network("model");
        // TODO
        return model;
    }

    public static Network gallery_mmk() {
        Network model = new Network("model");
        // TODO
        return model;
    }

    public static Network gallery_mpar1() {
        Network model = new Network("model");
        // TODO
        return model;
    }

    public static Network gallery_parm1() {
        Network model = new Network("model");
        // TODO
        return model;
    }

    public static Network gallery_repairmen() {
        Network model = new Network("model");
        // TODO
        return model;
    }

    public static Network gallery_replayerm1() {
        Network model = new Network("model");
        // TODO
        return model;
    }

    public static Network gallery_um1() {
        Network model = new Network("model");
        // TODO
        return model;
    }

    public static void main(String[] args) throws IllegalAccessException, ParserConfigurationException {
        Network model = gallery_mm1_prio();
        NetworkStruct sn = model.getStruct(false);
        new SolverJMT(model).getAvgTable().print();
        model.jsimgView();
        //new SolverMVA(gallery_mm1_tandem()).getAvgTable().tget("RespT", "Queue1","myClass").print();
    }
}
