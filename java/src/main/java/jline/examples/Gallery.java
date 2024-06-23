package jline.examples;

import jline.lang.*;
import jline.lang.constant.SchedStrategy;
import jline.lang.distributions.Exp;
import jline.lang.nodes.*;
import jline.solvers.mva.SolverMVA;

import javax.xml.parsers.ParserConfigurationException;
import java.util.List;
import java.util.ArrayList;

/**
 * Getting started examples
 */
public class Gallery {
    public static Network gallery_aphm1() {
        Network model = new Network("model");
        // TODO
        return model;
    }

    public static Network gallery_coxm1() {
        Network model = new Network("model");
        // TODO
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
        Network model = new Network("model");
        // TODO
        return model;
    }

    public static Network gallery_erldk() {
        Network model = new Network("model");
        // TODO
        return model;
    }

    public static Network gallery_erlerl1_reentrant() {
        Network model = new Network("model");
        // TODO
        return model;
    }

    public static Network gallery_erlm1() {
        Network model = new Network("model");
        // TODO
        return model;
    }

    public static Network gallery_erlm1_ps() {
        Network model = new Network("model");
        // TODO
        return model;
    }

    public static Network gallery_erlm1_reentrant() {
        Network model = new Network("model");
        // TODO
        return model;
    }

    public static Network gallery_gamm1() {
        Network model = new Network("model");
        // TODO
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
        Network model = new Network("model");
        // TODO
        return model;
    }

    public static Network gallery_hypm1_reentrant() {
        Network model = new Network("model");
        // TODO
        return model;
    }

    public static Network gallery_mapm1() {
        Network model = new Network("model");
        // TODO
        return model;
    }

    public static Network gallery_mapmk() {
        Network model = new Network("model");
        // TODO
        return model;
    }

    public static Network gallery_mdk() {
        Network model = new Network("model");
        // TODO
        return model;
    }

    public static Network gallery_merl1() {
        Network model = new Network("model");
        // TODO
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
        Network model = new Network("model");
        // TODO
        return model;
    }

    public static Network gallery_mhyp1() {
        Network model = new Network("model");
        // TODO
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
        Network model = new Network("model");
        // TODO
        return model;
    }

    public static Network gallery_mm1() {
        return gallery_mm1_linear(1, 0.5);
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

        Network model = new Network("M/Hyp/1-Linear");

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
        Network model = new Network("model");
        // TODO
        return model;
    }

    public static Network gallery_mm1_ps() {
        Network model = new Network("model");
        // TODO
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
        new SolverMVA(gallery_mm1_tandem()).getAvgTable().tget("RespT", "Queue1","myClass").print();
    }
}
