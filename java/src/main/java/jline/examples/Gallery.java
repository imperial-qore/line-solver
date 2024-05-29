package jline.examples;

import jline.lang.*;
import jline.lang.constant.SchedStrategy;
import jline.lang.constant.SolverType;
import jline.lang.distributions.APH;
import jline.lang.distributions.Erlang;
import jline.lang.distributions.Exp;
import jline.lang.nodes.*;
import jline.solvers.NetworkAvgTable;
import jline.solvers.SolverOptions;
import jline.solvers.jmt.SolverJMT;
import jline.solvers.mva.SolverMVA;

import javax.xml.parsers.ParserConfigurationException;
import java.util.Arrays;
import java.util.Collections;
import java.util.List;
import java.util.ArrayList;
import java.util.List;

/**
 * Getting started examples
 */
public class Gallery {

    public static Network gallery_mm1() {
        return gallery_mm1_linear(1, 0.5);
    }

    public static Network gallery_mm1_tandem() {
        return gallery_mm1_linear(2, 0.9);
    }

    public static Network gallery_mm1_tandem(Double Umax) {
        return gallery_mm1_linear(2, Umax);
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

    public static void main(String[] args) throws IllegalAccessException, ParserConfigurationException {
        //new SolverJMT(gallery_mm1()).getAvgTable().print();
        new SolverMVA(gallery_mm1()).getAvgQLenTable().tget("Queue1","myClass").print();
    }
}
