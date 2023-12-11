package jline.solvers.jmt;

import jline.examples.OpenModel;
import jline.lang.Network;
import jline.lang.constant.VerboseLevel;
import jline.solvers.NetworkAvgTable;
import jline.solvers.Solver;
import jline.solvers.SolverOptions;
import org.junit.jupiter.api.Test;

import javax.xml.parsers.ParserConfigurationException;
import java.util.List;

import static org.junit.jupiter.api.Assertions.assertEquals;

public class SolverJMTIOOpenExamplesTest {
    @Test
    public void test_example_openModel_1() throws ParserConfigurationException {
        Network model = OpenModel.ex1_line();

        SolverOptions options = Solver.defaultOptions();
        options.keep = false;
        options.verbose = VerboseLevel.STD;
        options.cutoff = 10;
        options.seed = 23000;
        options.iter_max = 200;
        SolverJMT solver = new SolverJMT(model, options);
        //solver.jsimgView();

        NetworkAvgTable avgTable = solver.getAvgTable();
        avgTable.printTable();
        // QLen values
        List<Double> QLen = avgTable.get(0);
        assertEquals(0.0210224504535244, QLen.get(0), 1e-13);
        assertEquals(0.0970769740346852, QLen.get(1), 1e-13);
        assertEquals(0.0, QLen.get(2), 1e-13);

        // Util values
        List<Double> Util = avgTable.get(1);
        assertEquals(0.0210224504535244, Util.get(0), 1e-13);
        assertEquals(0.0970769740346852, Util.get(1), 1e-13);
        assertEquals(0.0, Util.get(2), 1e-13);

        // RespT values
        List<Double> RespT = avgTable.get(2);
        assertEquals(0.2130525434758123, RespT.get(0), 1e-13);
        assertEquals(0.9963834920584885, RespT.get(1), 1e-13);
        assertEquals(0.0, RespT.get(2), 1e-13);

        // ResidT values
        List<Double> ResidT = avgTable.get(3);
        assertEquals(0.2130525434758124, ResidT.get(0), 1e-13);
        assertEquals(0.9963834920584888, ResidT.get(1), 1e-13);
        assertEquals(0.0, ResidT.get(2), 1e-13);

        // Tput values
        List<Double> ArvR = avgTable.get(4);
        List<Double> Tput = avgTable.get(5);
        assertEquals(0.0999978032520005, Tput.get(0), 1e-13);
        assertEquals(0.0999991505973458, Tput.get(1), 1e-13);
        assertEquals(0.1000000000000000, Tput.get(2), 1e-13);
    }

    @Test
    public void test_example_openModel_3() throws ParserConfigurationException {
        Network model = OpenModel.ex3_line();

        SolverOptions options = Solver.defaultOptions();
        options.keep = true;
        options.verbose = VerboseLevel.STD;
        options.cutoff = 7;
        options.seed = 23000;
        SolverJMT solver = new SolverJMT(model, options);
//        solver.jsimgView();
        NetworkAvgTable avgTable = solver.getAvgTable();
        // Expected values
        double[] expectedQLen = {0.0, 0.0, 0.0, 1.1982, 1.0325, 0.0, 0.0, 0.0, 0.81572};
        double[] expectedUtil = {0.0, 0.0, 0.0, 0.39921, 0.28192, 0.0, 0.0, 0.0, 0.44596};
        double[] expectedRespT = {0.0, 0.0, 0.0, 0.64707, 0.92558, 0.0, 0.0, 0.0, 0.27468};
        double[] expectedResidT = {0.0, 0.0, 0.0, 0.43138, 0.30853, 0.0, 0.0, 0.0, 0.27468};
        double[] expectedTput = {2.0204, 1.0009, 0.0, 2.008, 0.991, 0.0, 0.0, 0.0, 2.9856};

        for (int i = 0; i < expectedQLen.length; i++) {
            assertEquals(expectedQLen[i], avgTable.get(0).get(i), 1e-3);
            assertEquals(expectedUtil[i], avgTable.get(1).get(i), 1e-3);
            assertEquals(expectedRespT[i], avgTable.get(2).get(i), 1e-3);
            assertEquals(expectedResidT[i], avgTable.get(3).get(i), 1e-3);
            assertEquals(expectedTput[i], avgTable.get(5).get(i), 1e-3);
        }
    }
}
