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
//        solver.jsimgView();

        NetworkAvgTable avgTable = solver.getAvgTable();

        // QLen values
        List<Double> QLen = avgTable.get(0);
        assertEquals(0.0226789994053706, QLen.get(0), 1e-13);
        assertEquals(0.113319721462768, QLen.get(1), 1e-13);
        assertEquals(0.0, QLen.get(2), 1e-13);

        // Util values
        List<Double> Util = avgTable.get(1);
        assertEquals(0.0226789994053706, Util.get(0), 1e-13);
        assertEquals(0.102223795456197, Util.get(1), 1e-13);
        assertEquals(0.0, Util.get(2), 1e-13);

        // RespT values
        List<Double> RespT = avgTable.get(2);
        assertEquals(0.214832145125023, RespT.get(0), 1e-13);
        assertEquals(1.1075290123417, RespT.get(1), 1e-13);
        assertEquals(0.0, RespT.get(2), 1e-13);

        // ResidT values
        List<Double> ResidT = avgTable.get(3);
        assertEquals(0.214832145125023, ResidT.get(0), 1e-13);
        assertEquals(1.1075290123417, ResidT.get(1), 1e-13);
        assertEquals(0.0, ResidT.get(2), 1e-13);

        // Tput values
        List<Double> Tput = avgTable.get(4);
        assertEquals(0.0996073116118573, Tput.get(0), 1e-13);
        assertEquals(0.100011985776734, Tput.get(1), 1e-13);
        assertEquals(0.0996206521870407, Tput.get(2), 1e-13);
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
            assertEquals(expectedTput[i], avgTable.get(4).get(i), 1e-3);
        }
    }
}
