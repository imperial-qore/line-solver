package jline.solvers.jmt;

import jline.examples.ClosedModel;
import jline.lang.Network;
import jline.lang.constant.SolverType;
import jline.lang.constant.VerboseLevel;
import jline.solvers.NetworkAvgTable;
import jline.solvers.SolverOptions;
import org.junit.jupiter.api.Test;

import javax.xml.parsers.ParserConfigurationException;
import java.io.IOException;
import java.util.List;

import static org.junit.jupiter.api.Assertions.assertEquals;

public class SolverJMTClosedExamplesTest {
    @Test
    public void test_example_closedModel_1() throws ParserConfigurationException, IOException {
        Network model = ClosedModel.ex1();

        SolverOptions options = new SolverOptions(SolverType.JMT);
        options.seed = 23000;
        options.keep = false;
        SolverJMT solver = new SolverJMT(model, options);
//        solver.jsimgView();

        NetworkAvgTable avgTable = solver.getAvgTable();

        List<Double> QLen = avgTable.get(0);
        assertEquals(2.19091543487542, QLen.get(0), 1e-13);
        assertEquals(7.73312520184885, QLen.get(1), 1e-13);

        List<Double> Util = avgTable.get(1);
        assertEquals(2.19091543487542, Util.get(0), 1e-13);
        assertEquals(0.999774485968553, Util.get(1), 1e-13);

        List<Double> RespT = avgTable.get(2);
        assertEquals(1.00806433579534, RespT.get(0), 1e-13);
        assertEquals(11.5733935011948, RespT.get(1), 1e-13);

        List<Double> ResidT = avgTable.get(3);
        assertEquals(1.00806433579534, ResidT.get(0), 1e-13);
        assertEquals(3.47201805035843, ResidT.get(1), 1e-13);

        List<Double> Tput = avgTable.get(4);
        assertEquals(2.21287923853873, Tput.get(0), 1e-13);
        assertEquals(0.680220567268009, Tput.get(1), 1e-13);
    }

    @Test
    public void test_example_closedModel_2() {
        Network model = ClosedModel.ex2_line();

        SolverOptions options = new SolverOptions(SolverType.JMT);
        options.seed = 23000;
        options.keep = false;
        options.verbose = VerboseLevel.STD;
        options.samples = 5000;
        SolverJMT solver = new SolverJMT(model, options);

        NetworkAvgTable avgTable = solver.getAvgTable();

        List<Double> QLen = avgTable.get(0);
        assertEquals(0.930076508980514, QLen.get(0), 1e-13);
        assertEquals(0.202015817235601, QLen.get(1), 1e-13);
        assertEquals(0.0746719839497427, QLen.get(2), 1e-13);
        assertEquals(2.80758266283344, QLen.get(3), 1e-13);

        List<Double> Util = avgTable.get(1);
        assertEquals(0.930076508980514, Util.get(0), 1e-13);
        assertEquals(0.202015817235601, Util.get(1), 1e-13);
        assertEquals(0.018631696765107, Util.get(2), 1e-13);
        assertEquals(0.953009799376815, Util.get(3), 1e-13);


        List<Double> RespT = avgTable.get(2);
        assertEquals(0.669928699090526, RespT.get(0), 1e-13);
        assertEquals(0.214316949764258, RespT.get(1), 1e-13);
        assertEquals(0.549908013640484, RespT.get(2), 1e-13);
        assertEquals(2.93018093134332, RespT.get(3), 1e-13);


        List<Double> ResidT = avgTable.get(3);
        assertEquals(0.39876708279198, ResidT.get(0), 1e-13);
        assertEquals(0.0867473368093427, ResidT.get(1), 1e-13);
        assertEquals(0.0327326198595526, ResidT.get(2), 1e-13);
        assertEquals(1.18602561506753, ResidT.get(3), 1e-13);

        List<Double> Tput = avgTable.get(4);
        assertEquals(1.40780162218347, Tput.get(0), 1e-13);
        assertEquals(0.950637005374731, Tput.get(1), 1e-13);
        assertEquals(0.137244652741601, Tput.get(2), 1e-13);
        assertEquals(0.95062989904158, Tput.get(3), 1e-13);
    }

    @Test
    public void test_example_closedModel_3() throws ParserConfigurationException, IOException {
        Network model = ClosedModel.ex3_line();

        SolverOptions options = new SolverOptions(SolverType.JMT);
        options.seed = 23000;
        options.verbose = VerboseLevel.STD;
        options.samples = 5000;
        options.keep = false;
        SolverJMT solver = new SolverJMT(model, options);

        NetworkAvgTable avgTable = solver.getAvgTable();

        List<Double> QLen = avgTable.get(0);
        assertEquals(1.13010048435194, QLen.get(0), 1e-13);
        assertEquals(0.243556890473459, QLen.get(1), 1e-13);
        assertEquals(0.738007949161062, QLen.get(2), 1e-13);
        assertEquals(0.0380384564942096, QLen.get(3), 1e-13);
        assertEquals(0.600989509876208, QLen.get(4), 1e-13);
        assertEquals(0.261992050838938, QLen.get(5), 1e-13);

        List<Double> Util = avgTable.get(1);
        assertEquals(1.13010048435194, Util.get(0), 1e-13);
        assertEquals(0.243556890473459, Util.get(1), 1e-13);
        assertEquals(0.738007949161062, Util.get(2), 1e-13);
        assertEquals(0.0142993151385417, Util.get(3), 1e-13);
        assertEquals(0.303434113583001, Util.get(4), 1e-13);
        assertEquals(0.117663081760052, Util.get(5), 1e-13);

        List<Double> RespT = avgTable.get(2);
        assertEquals(0.660075733334815, RespT.get(0), 1e-13);
        assertEquals(0.209990808303743, RespT.get(1), 1e-13);
        assertEquals(1.00236400023138, RespT.get(2), 1e-13);
        assertEquals(0.201574834645298, RespT.get(3), 1e-13);
        assertEquals(0.514282772377555, RespT.get(4), 1e-13);
        assertEquals(0.354291512693584, RespT.get(5), 1e-13);

        List<Double> ResidT = avgTable.get(3);
        assertEquals(0.392902222223104, ResidT.get(0), 1e-13);
        assertEquals(0.0849962795515152, ResidT.get(1), 1e-13);
        assertEquals(1.00236400023138, ResidT.get(2), 1e-13);
        assertEquals(0.0119985020622201, ResidT.get(3), 1e-13);
        assertEquals(0.208162074533773, ResidT.get(4), 1e-13);
        assertEquals(0.354291512693584, ResidT.get(5), 1e-13);

        List<Double> Tput = avgTable.get(4);
        assertEquals(1.70923082948479, Tput.get(0), 1e-13);
        assertEquals(1.16177133451427, Tput.get(1), 1e-13);
        assertEquals(0.740279575628389, Tput.get(2), 1e-13);
        assertEquals(0.166115602085003, Tput.get(3), 1e-13);
        assertEquals(1.16219867056842, Tput.get(4), 1e-13);
        assertEquals(0.740289840681653, Tput.get(5), 1e-13);

    }

    @Test
    public void test_example_closedModel_4() throws ParserConfigurationException {
        Network model = ClosedModel.ex4_line();
        SolverOptions options = new SolverOptions(SolverType.JMT);
        options.seed = 23000;
        options.keep = true;
        options.verbose = VerboseLevel.STD;
        SolverJMT solver = new SolverJMT(model, options);

        NetworkAvgTable avgTable = solver.getAvgTable();

        // Expected values
        double[] expectedQLen = {0.96237, 0.47101, 0.13316, 0.67278, 1.4139, 1.178, 0.56534, 0.94042};
        double[] expectedUtil = {0.96237, 0.47101, 0.13316, 0.67278, 0.33229, 0.33028, 0.033717, 0.23274};
        double[] expectedRespT = {0.99545, 0.99639, 0.10003, 0.99114, 1.4149, 2.4004, 0.56114, 1.356};
        double[] expectedResidT = {0.66363, 0.33213, 0.066687, 0.33038, 0.94326, 0.80015, 0.28057, 0.45199};
        double[] expectedTput = {0.97749, 0.48837, 1.3475, 0.68114, 0.97798, 0.48832, 1.0216, 0.68119};

        // Assert each value with a tolerance of 1e-5
        for (int i = 0; i < expectedQLen.length; i++) {
            assertEquals(expectedQLen[i], avgTable.get(0).get(i), 1e-4);
            assertEquals(expectedUtil[i], avgTable.get(1).get(i), 1e-4);
            assertEquals(expectedRespT[i], avgTable.get(2).get(i), 1e-4);
            assertEquals(expectedResidT[i], avgTable.get(3).get(i), 1e-4);
            assertEquals(expectedTput[i], avgTable.get(4).get(i), 1e-4);
        }

//        solver.jsimgView(options);
    }

    @Test
    public void test_example_closedModel_9() throws ParserConfigurationException {
        Network model = ClosedModel.ex9_line();
        SolverOptions options = new SolverOptions(SolverType.JMT);
        options.seed = 23000;
        options.keep = true;
        options.verbose = VerboseLevel.STD;
        SolverJMT solver = new SolverJMT(model, options);

        NetworkAvgTable avgTable = solver.getAvgTable();

        // Expected values
        double[] expectedQLen = {0.16748, 0.17233, 0.47695, 0.49262, 9.3458, 9.3393};
        double[] expectedUtil = {0.16748, 0.17233, 0.24266, 0.24903, 0.49444, 0.50553};
        double[] expectedRespT = {0.98104, 0.99347, 2.8777, 2.9578, 55.78, 55.766};
        double[] expectedResidT = {0.98104, 0.99347, 2.8777, 2.9578, 55.78, 55.766};
        double[] expectedTput = {0.16787, 0.16675, 0.16768, 0.16618, 0.16742, 0.16685};

        // Assert each value with a tolerance of 1e-5
        for (int i = 0; i < expectedQLen.length; i++) {
            assertEquals(expectedQLen[i], avgTable.get(0).get(i), 1e-2);
            assertEquals(expectedUtil[i], avgTable.get(1).get(i), 1e-2);
            assertEquals(expectedRespT[i], avgTable.get(2).get(i), 1e-2);
            assertEquals(expectedResidT[i], avgTable.get(3).get(i), 1e-2);
            assertEquals(expectedTput[i], avgTable.get(4).get(i), 1e-2);
        }

//        solver.jsimgView(options);
    }
}
