package jline.opt;

import jline.lang.ClosedClass;
import jline.lang.Network;
import jline.lang.NetworkStruct;
import jline.lang.OpenClass;
import jline.lang.constant.SchedStrategy;
import jline.lang.nodes.Delay;
import jline.lang.nodes.Node;
import jline.lang.nodes.Queue;
import jline.lang.nodes.Sink;
import jline.lang.nodes.Source;
import jline.lang.nodes.Station;
import jline.lang.processes.Exp;
import jline.opt.variables.ClassPriority;
import jline.opt.variables.JobPopulation;
import jline.opt.variables.RoutingProbabilities;
import jline.opt.variables.ServerAllocation;
import jline.opt.variables.ServiceRate;
import org.junit.jupiter.api.Test;

import java.util.ArrayList;
import java.util.List;

import static org.junit.jupiter.api.Assertions.assertEquals;

/** Verifies each decision variable's decode + apply mutates a copied model. */
public class VariablesTest {

    private Network openModel() {
        Network model = new Network("mm1");
        Source src = new Source(model, "Src");
        Queue q = new Queue(model, "Server", SchedStrategy.FCFS);
        Sink snk = new Sink(model, "Snk");
        OpenClass oc = new OpenClass(model, "Jobs");
        src.setArrival(oc, new Exp(1.0));
        q.setService(oc, new Exp(2.0));
        model.link(Network.serialRouting(src, q, snk));
        return model;
    }

    private int stationIndexOf(Network m, String name) {
        List<Node> nodes = m.getNodes();
        NetworkStruct sn = m.getStruct();
        for (int i = 0; i < nodes.size(); i++) {
            if (nodes.get(i).getName().equals(name)) {
                return (int) sn.nodeToStation.get(i);
            }
        }
        return -1;
    }

    @Test
    public void testServerAllocationDecodeApply() {
        Network model = openModel();
        Station q = (Station) model.getNodes().get(1);
        ServerAllocation sa = new ServerAllocation(q, 1, 8);
        assertEquals(1, (int) (Integer) sa.decode(new double[]{0.0}));
        assertEquals(8, (int) (Integer) sa.decode(new double[]{1.0}));

        Network copy = model.copy();
        sa.apply(copy, 5);
        for (Node n : copy.getNodes()) {
            if (n.getName().equals("Server")) {
                assertEquals(5, ((Station) n).getNumberOfServers());
            }
        }
    }

    @Test
    public void testServiceRateApply() {
        Network model = openModel();
        Queue q = (Queue) model.getNodes().get(1);
        OpenClass oc = (OpenClass) model.getClasses().get(0);
        ServiceRate sr = new ServiceRate(q, oc, 0.5, 5.0);
        assertEquals(0.5, (Double) sr.decode(new double[]{0.0}), 1e-12);
        assertEquals(5.0, (Double) sr.decode(new double[]{1.0}), 1e-12);

        Network copy = model.copy();
        sr.apply(copy, 3.0);
        NetworkStruct sn = copy.getStruct();
        int st = stationIndexOf(copy, "Server");
        assertEquals(3.0, sn.rates.get(st, 0), 1e-9);
    }

    @Test
    public void testJobPopulationApply() {
        Network model = new Network("cqn");
        Delay delay = new Delay(model, "Think");
        Queue q = new Queue(model, "Q1", SchedStrategy.PS);
        ClosedClass cc = new ClosedClass(model, "C1", 5, delay);
        delay.setService(cc, new Exp(1.0));
        q.setService(cc, new Exp(2.0));
        model.link(Network.serialRouting(delay, q));

        JobPopulation jp = new JobPopulation(cc, 1, 50);
        assertEquals(20, (int) (Integer) jp.decode(new double[]{(20.0 - 1) / (50 - 1)}));

        Network copy = model.copy();
        jp.apply(copy, 20);
        NetworkStruct sn = copy.getStruct();
        assertEquals(20.0, sn.njobs.get(0), 1e-9);
    }

    @Test
    public void testClassPriorityApply() {
        Network model = new Network("cqn2");
        Delay delay = new Delay(model, "Think");
        Queue q = new Queue(model, "Q1", SchedStrategy.HOL);
        ClosedClass cc = new ClosedClass(model, "C1", 5, delay);
        delay.setService(cc, new Exp(1.0));
        q.setService(cc, new Exp(2.0));
        model.link(Network.serialRouting(delay, q));

        List<jline.lang.JobClass> classes = new ArrayList<jline.lang.JobClass>();
        classes.add(cc);
        ClassPriority cp = new ClassPriority(classes, "levels", 0, 5, "prio");
        Network copy = model.copy();
        cp.apply(copy, cp.decode(new double[]{1.0}));
        assertEquals(5, copy.getClasses().get(0).getPriority());
    }

    @Test
    public void testRoutingProbabilitiesDecode() {
        Network model = openModel();
        Source src = (Source) model.getNodes().get(0);
        Queue q = (Queue) model.getNodes().get(1);
        Sink snk = (Sink) model.getNodes().get(2);
        OpenClass oc = (OpenClass) model.getClasses().get(0);
        List<Node> targets = new ArrayList<Node>();
        targets.add(q);
        targets.add(snk);
        RoutingProbabilities rp = new RoutingProbabilities(oc, src, targets);
        double[] probs = (double[]) rp.decode(new double[]{0.743});
        assertEquals(1.0, probs[0] + probs[1], 1e-12);
        assertEquals(0.743, probs[0], 1e-12);
    }
}
