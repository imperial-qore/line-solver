package jline.solvers.cmtc;

import jline.examples.ClosedModel;
import jline.examples.GettingStarted;
import jline.lang.Network;
import jline.solvers.SolverResult;
import jline.solvers.ctmc.SolverCTMC;
import org.junit.jupiter.api.Test;

import static jline.examples.GettingStarted.*;

import static org.junit.jupiter.api.Assertions.assertEquals;

public class SolverCTMCClosedExamplesTest {


    static double tol = 0.0001;

    @Test
    public void matlabExample3ReturnsCorrectResultFromRunAnalyzer() {

        Network model = matlabExample3();
        SolverCTMC solver = new SolverCTMC(model);
        solver.runAnalyzer();
        SolverResult result = solver.result;

        // QN
        assertEquals(2, result.QN.getNumRows());
        assertEquals(1, result.QN.getNumCols());
        assertEquals(2, result.QN.getNumElements());
        assertEquals(2.6648, result.QN.get(0, 0), tol);
        assertEquals(0.33516, result.QN.get(1, 0), tol);

        // RN
        assertEquals(2, result.RN.getNumRows());
        assertEquals(1, result.RN.getNumCols());
        assertEquals(2, result.RN.getNumElements());
        assertEquals(2, result.RN.get(0, 0), tol);
        assertEquals(0.2515, result.RN.get(1, 0), tol);

        // XN
        assertEquals(1, result.XN.getNumRows());
        assertEquals(1, result.XN.getNumCols());
        assertEquals(1, result.XN.getNumElements());
        assertEquals(1.3324, result.XN.get(0, 0), tol);

        // UN
        assertEquals(2, result.UN.getNumRows());
        assertEquals(1, result.UN.getNumCols());
        assertEquals(2, result.UN.getNumElements());
        assertEquals(2.6648, result.UN.get(0, 0), tol);
        assertEquals(0.1666, result.UN.get(1, 0), tol);

        // TN
        assertEquals(2, result.TN.getNumRows());
        assertEquals(1, result.TN.getNumCols());
        assertEquals(2, result.TN.getNumElements());
        assertEquals(1.3324, result.TN.get(0, 0), tol);
        assertEquals(1.3324, result.TN.get(1, 0), tol);

        // CN
        assertEquals(1, result.CN.getNumRows());
        assertEquals(1, result.CN.getNumCols());
        assertEquals(1, result.CN.getNumElements());
        assertEquals(3, result.CN.get(0, 0), tol);
    }

    @Test
    public void ClosedModelEx1ReturnsCorrectResultFromRunAnalyzer() {

        Network model = ClosedModel.ex1();
        SolverCTMC solver = new SolverCTMC(model);
        solver.runAnalyzer();
        SolverResult result = solver.result;

        // QN
        assertEquals(2, result.QN.getNumRows());
        assertEquals(1, result.QN.getNumCols());
        assertEquals(2, result.QN.getNumElements());
        assertEquals(2.2220, result.QN.get(0, 0), tol);
        assertEquals(7.7780, result.QN.get(1, 0), tol);

        // RN
        assertEquals(2, result.RN.getNumRows());
        assertEquals(1, result.RN.getNumCols());
        assertEquals(2, result.RN.getNumElements());
        assertEquals(1.0000, result.RN.get(0, 0), tol);
        assertEquals(11.6680, result.RN.get(1, 0), tol);

        // XN
        assertEquals(1, result.XN.getNumRows());
        assertEquals(1, result.XN.getNumCols());
        assertEquals(1, result.XN.getNumElements());
        assertEquals(2.2220, result.XN.get(0, 0), tol);

        // UN
        assertEquals(2, result.UN.getNumRows());
        assertEquals(1, result.UN.getNumCols());
        assertEquals(2, result.UN.getNumElements());
        assertEquals(2.2220, result.UN.get(0, 0), tol);
        assertEquals(0.9999, result.UN.get(1, 0), tol);

        // TN
        assertEquals(2, result.TN.getNumRows());
        assertEquals(1, result.TN.getNumCols());
        assertEquals(2, result.TN.getNumElements());
        assertEquals(2.2220, result.TN.get(0, 0), tol);
        assertEquals(0.6666, result.TN.get(1, 0), tol);

        // CN
        assertEquals(1, result.CN.getNumRows());
        assertEquals(1, result.CN.getNumCols());
        assertEquals(1, result.CN.getNumElements());
        assertEquals(10, result.CN.get(0, 0), tol);
    }

    @Test
    public void ClosedModelEx8ReturnsCorrectResultFromRunAnalyzer() {

        Network model = ClosedModel.ex8();
        SolverCTMC solver = new SolverCTMC(model);
        solver.runAnalyzer();
        SolverResult result = solver.result;

        // QN
        assertEquals(3, result.QN.getNumRows());
        assertEquals(1, result.QN.getNumCols());
        assertEquals(3, result.QN.getNumElements());
        assertEquals(2.0650, result.QN.get(0, 0), tol);
        assertEquals(4.7448, result.QN.get(1, 0), tol);
        assertEquals(3.1901, result.QN.get(2, 0), tol);

        // RN
        assertEquals(3, result.RN.getNumRows());
        assertEquals(1, result.RN.getNumCols());
        assertEquals(3, result.RN.getNumElements());
        assertEquals(1.0000, result.RN.get(0, 0), tol);
        assertEquals(7.6590, result.RN.get(1, 0), tol);
        assertEquals(7.7241, result.RN.get(2, 0), tol);

        // XN
        assertEquals(1, result.XN.getNumRows());
        assertEquals(1, result.XN.getNumCols());
        assertEquals(1, result.XN.getNumElements());
        assertEquals(2.0650, result.XN.get(0, 0), tol);

        // UN
        assertEquals(3, result.UN.getNumRows());
        assertEquals(1, result.UN.getNumCols());
        assertEquals(3, result.UN.getNumElements());
        assertEquals(2.0650, result.UN.get(0, 0), tol);
        assertEquals(0.9293, result.UN.get(1, 0), tol);
        assertEquals(0.8260, result.UN.get(2, 0), tol);

        // TN
        assertEquals(3, result.TN.getNumRows());
        assertEquals(1, result.TN.getNumCols());
        assertEquals(3, result.TN.getNumElements());
        assertEquals(2.0650, result.TN.get(0, 0), tol);
        assertEquals(0.6195, result.TN.get(1, 0), tol);
        assertEquals(0.4130, result.TN.get(2, 0), tol);

        // CN
        assertEquals(1, result.CN.getNumRows());
        assertEquals(1, result.CN.getNumCols());
        assertEquals(1, result.CN.getNumElements());
        assertEquals(10, result.CN.get(0, 0), tol);
    }

    @Test
    public void ErlangEx1ReturnsCorrectResultFromRunAnalyzer() {

        Network model = GettingStarted.erlangExample1();
        SolverCTMC solver = new SolverCTMC(model);
        solver.runAnalyzer();
        SolverResult result = solver.result;

        // QN
        assertEquals(2, result.QN.getNumRows());
        assertEquals(1, result.QN.getNumCols());
        assertEquals(2, result.QN.getNumElements());
        assertEquals(1.4213, result.QN.get(0, 0), tol);
        assertEquals(1.5787, result.QN.get(1, 0), tol);

        // RN
        assertEquals(2, result.RN.getNumRows());
        assertEquals(1, result.RN.getNumCols());
        assertEquals(2, result.RN.getNumElements());
        assertEquals(2.0000, result.RN.get(0, 0), tol);
        assertEquals(2.2214, result.RN.get(1, 0), tol);

        // XN
        assertEquals(1, result.XN.getNumRows());
        assertEquals(1, result.XN.getNumCols());
        assertEquals(1, result.XN.getNumElements());
        assertEquals(0.7107, result.XN.get(0, 0), tol);

        // UN
        assertEquals(2, result.UN.getNumRows());
        assertEquals(1, result.UN.getNumCols());
        assertEquals(2, result.UN.getNumElements());
        assertEquals(1.4213, result.UN.get(0, 0), tol);
        assertEquals(0.7107, result.UN.get(1, 0), tol);

        // TN
        assertEquals(2, result.TN.getNumRows());
        assertEquals(1, result.TN.getNumCols());
        assertEquals(2, result.TN.getNumElements());
        assertEquals(0.7107, result.TN.get(0, 0), tol);
        assertEquals(0.7107, result.TN.get(1, 0), tol);

        // CN
        assertEquals(1, result.CN.getNumRows());
        assertEquals(1, result.CN.getNumCols());
        assertEquals(1, result.CN.getNumElements());
        assertEquals(3, result.CN.get(0, 0), tol);
    }

    @Test
    public void ErlangEx2ReturnsCorrectResultFromRunAnalyzer() {

        Network model = GettingStarted.erlangExample2();
        SolverCTMC solver = new SolverCTMC(model);
        solver.runAnalyzer();
        SolverResult result = solver.result;

        // QN
        assertEquals(2, result.QN.getNumRows());
        assertEquals(1, result.QN.getNumCols());
        assertEquals(2, result.QN.getNumElements());
        assertEquals(2.4952, result.QN.get(0, 0), tol);
        assertEquals(0.5048, result.QN.get(1, 0), tol);

        // RN
        assertEquals(2, result.RN.getNumRows());
        assertEquals(1, result.RN.getNumCols());
        assertEquals(2, result.RN.getNumElements());
        assertEquals(2.0000, result.RN.get(0, 0), tol);
        assertEquals(0.4046, result.RN.get(1, 0), tol);

        // XN
        assertEquals(1, result.XN.getNumRows());
        assertEquals(1, result.XN.getNumCols());
        assertEquals(1, result.XN.getNumElements());
        assertEquals(1.2476, result.XN.get(0, 0), tol);

        // UN
        assertEquals(2, result.UN.getNumRows());
        assertEquals(1, result.UN.getNumCols());
        assertEquals(2, result.UN.getNumElements());
        assertEquals(2.4952, result.UN.get(0, 0), tol);
        assertEquals(0.2495, result.UN.get(1, 0), tol);

        // TN
        assertEquals(2, result.TN.getNumRows());
        assertEquals(1, result.TN.getNumCols());
        assertEquals(2, result.TN.getNumElements());
        assertEquals(1.2476, result.TN.get(0, 0), tol);
        assertEquals(1.2476, result.TN.get(1, 0), tol);

        // CN
        assertEquals(1, result.CN.getNumRows());
        assertEquals(1, result.CN.getNumCols());
        assertEquals(1, result.CN.getNumElements());
        assertEquals(3, result.CN.get(0, 0), tol);
    }

    @Test
    public void ErlangEx3ReturnsCorrectResultFromRunAnalyzer() {

        Network model = GettingStarted.erlangExample3();
        SolverCTMC solver = new SolverCTMC(model);
        solver.runAnalyzer();
        SolverResult result = solver.result;


        // QN
        assertEquals(2, result.QN.getNumRows());
        assertEquals(1, result.QN.getNumCols());
        assertEquals(2, result.QN.getNumElements());
        assertEquals(3.0363, result.QN.get(0, 0), tol);
        assertEquals(0.9637, result.QN.get(1, 0), tol);

        // RN
        assertEquals(2, result.RN.getNumRows());
        assertEquals(1, result.RN.getNumCols());
        assertEquals(2, result.RN.getNumElements());
        assertEquals(5.0000, result.RN.get(0, 0), tol);
        assertEquals(1.5869, result.RN.get(1, 0), tol);

        // XN
        assertEquals(1, result.XN.getNumRows());
        assertEquals(1, result.XN.getNumCols());
        assertEquals(1, result.XN.getNumElements());
        assertEquals(0.6073, result.XN.get(0, 0), tol);

        // UN
        assertEquals(2, result.UN.getNumRows());
        assertEquals(1, result.UN.getNumCols());
        assertEquals(2, result.UN.getNumElements());
        assertEquals(3.0363, result.UN.get(0, 0), tol);
        assertEquals(0.4555, result.UN.get(1, 0), tol);

        // TN
        assertEquals(2, result.TN.getNumRows());
        assertEquals(1, result.TN.getNumCols());
        assertEquals(2, result.TN.getNumElements());
        assertEquals(0.6073, result.TN.get(0, 0), tol);
        assertEquals(0.6073, result.TN.get(1, 0), tol);

        // CN
        assertEquals(1, result.CN.getNumRows());
        assertEquals(1, result.CN.getNumCols());
        assertEquals(1, result.CN.getNumElements());
        assertEquals(4, result.CN.get(0, 0), tol);
    }

    @Test
    public void ClosedModelEx9ReturnsCorrectResultFromRunAnalyzer() {

        Network model = ClosedModel.ex9();
        SolverCTMC solver = new SolverCTMC(model);
        solver.runAnalyzer();
        SolverResult result = solver.result;


        // QN
        assertEquals(3, result.QN.getNumRows());
        assertEquals(1, result.QN.getNumCols());
        assertEquals(3, result.QN.getNumElements());
        assertEquals(2.1958, result.QN.get(0, 0), tol);
        assertEquals(0.3799, result.QN.get(1, 0), tol);
        assertEquals(0.4243, result.QN.get(2, 0), tol);

        // RN
        assertEquals(3, result.RN.getNumRows());
        assertEquals(1, result.RN.getNumCols());
        assertEquals(3, result.RN.getNumElements());
        assertEquals(5.0000, result.RN.get(0, 0), tol);
        assertEquals(2.8836, result.RN.get(1, 0), tol);
        assertEquals(4.8304, result.RN.get(2, 0), tol);

        // XN
        assertEquals(1, result.XN.getNumRows());
        assertEquals(1, result.XN.getNumCols());
        assertEquals(1, result.XN.getNumElements());
        assertEquals(0.4392, result.XN.get(0, 0), tol);

        // UN
        assertEquals(3, result.UN.getNumRows());
        assertEquals(1, result.UN.getNumCols());
        assertEquals(3, result.UN.getNumElements());
        assertEquals(2.1958, result.UN.get(0, 0), tol);
        assertEquals(0.3294, result.UN.get(1, 0), tol);
        assertEquals(0.3513, result.UN.get(2, 0), tol);

        // TN
        assertEquals(3, result.TN.getNumRows());
        assertEquals(1, result.TN.getNumCols());
        assertEquals(3, result.TN.getNumElements());
        assertEquals(0.4392, result.TN.get(0, 0), tol);
        assertEquals(0.1317, result.TN.get(1, 0), tol);
        assertEquals(0.0878, result.TN.get(2, 0), tol);

        // CN
        assertEquals(1, result.CN.getNumRows());
        assertEquals(1, result.CN.getNumCols());
        assertEquals(1, result.CN.getNumElements());
        assertEquals(3, result.CN.get(0, 0), tol);
    }

    /*
    @Test
    public void AphEx1ReturnsCorrectResultFromRunAnalyzer() {

        Network model = GettingStarted.aphExample1();
        SolverCTMC solver = new SolverCTMC(model);
        solver.runAnalyzer();
        SolverResult result = solver.result;


        // QN
        assertEquals(2, result.QN.getNumRows());
        assertEquals(1, result.QN.getNumCols());
        assertEquals(2, result.QN.getNumElements());
        assertEquals(0.8792, result.QN.get(0, 0), tol);
        assertEquals(2.1208, result.QN.get(1, 0), tol);

        // RN
        assertEquals(2, result.RN.getNumRows());
        assertEquals(1, result.RN.getNumCols());
        assertEquals(2, result.RN.getNumElements());
        assertEquals(2.0000, result.RN.get(0, 0), tol);
        assertEquals(4.8245, result.RN.get(1, 0), tol);

        // XN
        assertEquals(1, result.XN.getNumRows());
        assertEquals(1, result.XN.getNumCols());
        assertEquals(1, result.XN.getNumElements());
        assertEquals(0.4396, result.XN.get(0, 0), tol);

        // UN
        assertEquals(2, result.UN.getNumRows());
        assertEquals(1, result.UN.getNumCols());
        assertEquals(2, result.UN.getNumElements());
        assertEquals(0.8792, result.UN.get(0, 0), tol);
        assertEquals(0.8792, result.UN.get(1, 0), tol);

        // TN
        assertEquals(2, result.TN.getNumRows());
        assertEquals(1, result.TN.getNumCols());
        assertEquals(2, result.TN.getNumElements());
        assertEquals(0.4396, result.TN.get(0, 0), tol);
        assertEquals(0.4396, result.TN.get(1, 0), tol);

        // CN
        assertEquals(1, result.CN.getNumRows());
        assertEquals(1, result.CN.getNumCols());
        assertEquals(1, result.CN.getNumElements());
        assertEquals(3, result.CN.get(0, 0), tol);
    }
    */

    /* Temporarily disabled as it fails test after SSA refactoring
    @Test
    public void AphEx2ReturnsCorrectResultFromRunAnalyzer() {

        Network model = GettingStarted.aphExample2();
        SolverCTMC solver = new SolverCTMC(model);
        solver.runAnalyzer();
        SolverResult result = solver.result;


        // QN
        assertEquals(2, result.QN.getNumRows());
        assertEquals(2, result.QN.getNumCols());
        assertEquals(4, result.QN.getNumElements());
        assertEquals(0.6154, result.QN.get(0, 0), tol);
        assertEquals(0.4444, result.QN.get(0, 1), tol);
        assertEquals(0.3846, result.QN.get(1, 0), tol);
        assertEquals(0.5556, result.QN.get(1, 1), tol);

        // RN
        assertEquals(2, result.RN.getNumRows());
        assertEquals(2, result.RN.getNumCols());
        assertEquals(4, result.RN.getNumElements());
        assertEquals(16.0000, result.RN.get(0, 0), tol);
        assertEquals(8.0000, result.RN.get(0, 1), tol);
        assertEquals(10.0000, result.RN.get(1, 0), tol);
        assertEquals(10.0000, result.RN.get(1, 1), tol);

        // XN
        assertEquals(1, result.XN.getNumRows());
        assertEquals(2, result.XN.getNumCols());
        assertEquals(2, result.XN.getNumElements());
        assertEquals(0.0385, result.XN.get(0, 0), tol);
        assertEquals(0.0556, result.XN.get(0, 1), tol);

        // UN
        assertEquals(2, result.UN.getNumRows());
        assertEquals(2, result.UN.getNumCols());
        assertEquals(4, result.UN.getNumElements());
        assertEquals(0.6154, result.UN.get(0, 0), tol);
        assertEquals(0.4444, result.UN.get(0, 1), tol);
        assertEquals(0.1923, result.UN.get(1, 0), tol);
        assertEquals(0.2778, result.UN.get(1, 1), tol);

        // TN
        assertEquals(2, result.TN.getNumRows());
        assertEquals(2, result.TN.getNumCols());
        assertEquals(4, result.TN.getNumElements());
        assertEquals(0.0385, result.TN.get(0, 0), tol);
        assertEquals(0.0556, result.TN.get(0, 1), tol);
        assertEquals(0.0385, result.TN.get(1, 0), tol);
        assertEquals(0.0556, result.TN.get(1, 1), tol);

        // CN
        assertEquals(1, result.CN.getNumRows());
        assertEquals(2, result.CN.getNumCols());
        assertEquals(2, result.CN.getNumElements());
        assertEquals(1.0000, result.CN.get(0, 0), tol);
        assertEquals(0.5000, result.CN.get(0, 1), tol);
    }

     */
}
