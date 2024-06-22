package jline.solvers.jmt;

import jline.examples.*;
import jline.lang.Network;
import jline.lang.OpenClass;
import jline.lang.constant.SchedStrategy;
import jline.lang.constant.SolverType;
import jline.lang.constant.VerboseLevel;
import jline.lang.distributions.Exp;
import jline.lang.nodes.Queue;
import jline.lang.nodes.Sink;
import jline.lang.nodes.Source;
import jline.solvers.NetworkAvgTable;
import jline.solvers.SolverOptions;
import org.junit.jupiter.api.Test;

import javax.xml.parsers.ParserConfigurationException;
import java.util.List;

import static org.junit.jupiter.api.Assertions.assertEquals;

public class SolverJMTIOTest {

    @Test
    public void test_view() throws ParserConfigurationException {
        Network model = ClosedModel.example_closedModel_1();

        SolverOptions options = new SolverOptions(SolverType.JMT);
        options.seed = 23000;
        options.keep = true;
        SolverJMT solver = new SolverJMT(model, options);

        solver.runAnalyzer();
        //solver.jsimgView();
    }

    @Test
    public void test_getting_started_example_1() {
        Network model = new Network("M/M/1");
        Source source = new Source(model, "Source");
        Queue queue = new Queue(model, "Queue", SchedStrategy.FCFS);
        Sink sink = new Sink(model, "Sink");

        OpenClass oclass = new OpenClass(model, "Class1");
        source.setArrival(oclass, new Exp(1));
        queue.setService(oclass, new Exp(2));

        model.link(model.serialRouting(source, queue, sink));

        SolverOptions options = new SolverOptions(SolverType.JMT);
        options.seed = 23000;
        options.keep = true;
        options.verbose = VerboseLevel.STD;
        options.samples = 10000;
        SolverJMT solver = new SolverJMT(model, options);

        NetworkAvgTable avgTable = solver.getAvgTable();

        List<Double> QLen = avgTable.get(0);
        assertEquals(0, QLen.get(0), 1e-13);
        assertEquals(0.955501010809008, QLen.get(1), 1e-13);

        List<Double> Util = avgTable.get(1);
        assertEquals(0, Util.get(0), 1e-13);
        assertEquals(0.487360218010475, Util.get(1), 1e-13);


        List<Double> RespT = avgTable.get(2);
        assertEquals(0, RespT.get(0), 1e-13);
        assertEquals(0.954292928096683, RespT.get(1), 1e-13);


        List<Double> ResidT = avgTable.get(3);
        assertEquals(0, ResidT.get(0), 1e-13);
        assertEquals(0.954292928096683, ResidT.get(1), 1e-13);

        List<Double> ArvR = avgTable.get(4);
        List<Double> Tput = avgTable.get(5);
        assertEquals(0.998941736137035, Tput.get(0), 1e-13);
        assertEquals(0.99986838711003, Tput.get(1), 1e-13);
    }

    @Test
    public void test_example_forkJoin_1() throws ParserConfigurationException {
        Network model = ForkJoinModel.ex1_line();

        SolverOptions options = new SolverOptions(SolverType.JMT);
        options.seed = 23000;
//        options.keep = true;
        SolverJMT solver = new SolverJMT(model, options);

        NetworkAvgTable avgTable = solver.getAvgTable();
//        solver.jsimgView();

        // Expected values
        double[] expectedQLen = {0, 0.055452, 0.026193, 0.046605};
        double[] expectedUtil = {0, 0.051901, 0.025582, 0};
        double[] expectedRespT = {0, 1.0825, 0.50222, 0.452};
        double[] expectedResidT = {0, 0.54124, 0.25111, 0.452};
        double[] expectedTput = {0.051125, 0.050904, 0.051125, 0.050904};

        for (int i = 0; i < expectedQLen.length; i++) {
            assertEquals(expectedQLen[i], avgTable.get(0).get(i), 1e-3);
            assertEquals(expectedUtil[i], avgTable.get(1).get(i), 1e-3);
            assertEquals(expectedRespT[i], avgTable.get(2).get(i), 1e-3);
            assertEquals(expectedResidT[i], avgTable.get(3).get(i), 1e-3);
            assertEquals(expectedTput[i], avgTable.get(5).get(i), 1e-3);
        }
    }

    @Test void test_example_forkJoin_2() {
        Network model = ForkJoinModel.ex2_line();

        SolverOptions options = new SolverOptions(SolverType.JMT);
        options.seed = 23000;
//        options.keep = true;
        SolverJMT solver = new SolverJMT(model, options);

        NetworkAvgTable avgTable = solver.getAvgTable();

        // Expected values
        double[] expectedQLen = {0, 0, 1.3998, 0, 11.34, 4.1457, 20.023, 8.1158};
        double[] expectedUtil = {0, 0, 0.51269, 0, 0.69983, 0.2507, 0, 0};
        double[] expectedRespT = {0, 0, 3.0043, 0, 20.608, 8.5604, 19.275, 8.1805};
        double[] expectedResidT = {0, 0, 1.5021, 0, 10.304, 4.2802, 19.275, 8.1805};
        double[] expectedTput = {0.251, 0.24768, 0.52094, 0.50036, 0.51149, 0.50612, 0.2507, 0.24802};

        for (int i = 0; i < expectedQLen.length; i++) {
            assertEquals(expectedQLen[i], avgTable.get(0).get(i), 1e-3);
            assertEquals(expectedUtil[i], avgTable.get(1).get(i), 1e-3);
            assertEquals(expectedRespT[i], avgTable.get(2).get(i), 1e-3);
            assertEquals(expectedResidT[i], avgTable.get(3).get(i), 1e-3);
            assertEquals(expectedTput[i], avgTable.get(5).get(i), 1e-3);
        }
    }

    @Test void test_example_forkJoin_3() throws ParserConfigurationException {
        Network model = ForkJoinModel.ex3_line();

        SolverOptions options = new SolverOptions(SolverType.JMT);
        options.seed = 23000;
//        options.keep = true;
        SolverJMT solver = new SolverJMT(model, options);
//        solver.jsimgView();

        NetworkAvgTable avgTable = solver.getAvgTable();

        // Expected values
        double[] expectedQLen = {1.9934, 1.0125, 2.1943, 1.0062, 0, 0.51907, 1.4034, 0.85773, 2.3717, 1.0364};
        double[] expectedUtil = {1.9934, 1.0125, 0, 0, 0, 0, 0.50771, 0.25614, 0.66741, 0.25295};
        double[] expectedRespT = {3.9661, 3.9831, 2.2219, 0.99526, 0, 1.02, 2.694, 1.6624, 4.7882, 2.1242};
        double[] expectedResidT = {3.9661, 3.9831, 2.2219, 0.99526, 0, 1.02, 1.347, 0.83119, 2.3941, 1.0621};
        double[] expectedTput = {0.50379, 0.25574, 0.50404, 0.50666, 0, 0.25395, 0.5015, 0.51104, 0.50815, 0.5069};

        for (int i = 0; i < expectedQLen.length; i++) {
            assertEquals(expectedQLen[i], avgTable.get(0).get(i), 1e-2);
            assertEquals(expectedUtil[i], avgTable.get(1).get(i), 1e-2);
            assertEquals(expectedRespT[i], avgTable.get(2).get(i), 1e-2);
            assertEquals(expectedResidT[i], avgTable.get(3).get(i), 1e-2);
            assertEquals(expectedTput[i], avgTable.get(5).get(i), 1e-2);
        }
    }

    @Test
    public void test_example_mixedModel_1() throws ParserConfigurationException {
        Network model = MixedModel.ex1_line();

        SolverOptions options = new SolverOptions(SolverType.JMT);
        options.seed = 23000;
        options.keep = true;
        options.cutoff = 3;
        options.verbose = VerboseLevel.STD;
        SolverJMT solver = new SolverJMT(model, options);

        NetworkAvgTable avgTable = solver.getAvgTable();
//        solver.jsimgView();

        // Expected values
        double[] expectedQLen = {1.4528, 0.021961, 0.54717, 0.16391, 0, 0};
        double[] expectedUtil = {1.4528, 0.021961, 0.40216, 0.095781, 0, 0};
        double[] expectedRespT = {0.66644, 0.21581, 0.25983, 1.711, 0, 0};
        double[] expectedResidT = {0.66644, 0.21581, 0.25983, 1.711, 0, 0};
        double[] expectedTput = {2.1536, 0.10001, 2.1602, 0.10001, 0, 0.10001};

        for (int i = 0; i < expectedQLen.length; i++) {
            assertEquals(expectedQLen[i], avgTable.get(0).get(i), 1e-3);
            assertEquals(expectedUtil[i], avgTable.get(1).get(i), 1e-3);
            assertEquals(expectedRespT[i], avgTable.get(2).get(i), 1e-3);
            assertEquals(expectedResidT[i], avgTable.get(3).get(i), 1e-3);
            assertEquals(expectedTput[i], avgTable.get(5).get(i), 1e-3);
        }
    }
}
