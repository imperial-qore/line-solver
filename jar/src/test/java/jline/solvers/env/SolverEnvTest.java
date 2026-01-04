package jline.solvers.env;

import jline.solvers.SolverResult;
import jline.util.matrix.Matrix;
import org.junit.jupiter.api.Test;
import org.junit.jupiter.api.Disabled;
import org.junit.jupiter.api.Timeout;

import java.util.ArrayList;
import java.util.List;

import static org.junit.jupiter.api.Assertions.*;

import jline.solvers.env.SolverENV.SamplePathResult;

// @Disabled("SolverENV results have 16-28% value mismatches vs CTMC ground truth - algorithm needs revision")
public class SolverEnvTest {

    static double tol = 0.03;

    @Test
    public void exampleRandomEnvironment1ReturnsCorrectResult() {

        SolverENV envSolver = SolverEnvTestFixtures.renv_twostages_repairmen();
        envSolver.getAvg();
        SolverResult result = envSolver.result;

        // QN
        assertEquals(2, result.QN.getNumRows());
        assertEquals(1, result.QN.getNumCols());
        assertEquals(2, result.QN.getNumElements());
        assertEquals(0.5590, result.QN.value(), 0.5590 * tol);
        assertEquals(0.4410, result.QN.get(1, 0), 0.4410 * tol);

        // UN
        assertEquals(2, result.UN.getNumRows());
        assertEquals(1, result.UN.getNumCols());
        assertEquals(2, result.UN.getNumElements());
        assertEquals(0.5590, result.UN.value(), 0.5590 * tol);
        assertEquals(0.4410, result.UN.get(1, 0), 0.4410 * tol);

        // TN
        assertEquals(2, result.TN.getNumRows());
        assertEquals(1, result.TN.getNumCols());
        assertEquals(2, result.TN.getNumElements());
        assertEquals(0.6946, result.TN.value(), 0.6946 * tol);
        assertEquals(0.6842, result.TN.get(1, 0), 0.6842 * tol);
    }

    @Test
    @Disabled("1s, Test fails - value mismatch")
    public void exampleRandomEnvironment2ReturnsCorrectResult() {

        SolverENV envSolver = SolverEnvTestFixtures.renv_fourstages_repairmen();
        envSolver.getAvg();
        SolverResult result = envSolver.result;

        //  QN
        assertEquals(2, result.QN.getNumRows());
        assertEquals(1, result.QN.getNumCols());
        assertEquals(2, result.QN.getNumElements());
        assertEquals(0.4444, result.QN.value(), 0.4444 * tol);
        assertEquals(29.556, result.QN.get(1, 0), 29.556 * tol);

        //  UN
        assertEquals(2, result.UN.getNumRows());
        assertEquals(1, result.UN.getNumCols());
        assertEquals(2, result.UN.getNumElements());
        assertEquals(0.4444, result.UN.value(), 0.4444 * tol);
        assertEquals(1, result.UN.get(1, 0), 1 * tol);

        // TN
        assertEquals(2, result.TN.getNumRows());
        assertEquals(1, result.TN.getNumCols());
        assertEquals(2, result.TN.getNumElements());
        assertEquals(0.9784, result.TN.value(), 0.9784 * tol);
        assertEquals(1, result.TN.get(1, 0), 1 * tol);
    }

    @Test
    @Disabled("1s, Test fails - value mismatch")
    public void exampleRandomEnvironment3ReturnsCorrectResult() {

        SolverENV envSolver = SolverEnvTestFixtures.renv_threestages_repairmen();
        envSolver.getAvg();
        SolverResult result = envSolver.result;

        // QN
        assertEquals(2, result.QN.getNumRows());
        assertEquals(1, result.QN.getNumCols());
        assertEquals(2, result.QN.getNumElements());
        assertEquals(0.9446, result.QN.value(), 0.9446 * tol);
        assertEquals(1.0554, result.QN.get(1, 0), 1.0554 * tol);

        // UN
        assertEquals(2, result.UN.getNumRows());
        assertEquals(1, result.UN.getNumCols());
        assertEquals(2, result.UN.getNumElements());
        assertEquals(0.9446, result.UN.value(), 0.9446 * tol);
        assertEquals(0.8444, result.UN.get(1, 0), 0.8444 * tol);

        // TN
        assertEquals(2, result.TN.getNumRows());
        assertEquals(1, result.TN.getNumCols());
        assertEquals(2, result.TN.getNumElements());

        assertEquals(1.5542, result.TN.value(), 1.5542 * tol);
        assertEquals(1.5339, result.TN.get(1, 0), 1.5339 * tol);
    }

    @Test
    public void StateIndependentOneDelayStationTwoMMPPQueuesReturnsCorrectResult() {
        SolverENV envSolver = SolverEnvTestFixtures.example_randomEnvironment_4();
        envSolver.setStateDepMethod("stateindep");
        envSolver.getAvg();
        // envSolver.printAvgTable();
        SolverResult result = envSolver.result;

        // QN
        assertEquals(3, result.QN.getNumRows());
        assertEquals(1, result.QN.getNumCols());
        assertEquals(3, result.QN.getNumElements());
        assertEquals(41.7292, result.QN.value(), 41.7292* tol);
        assertEquals(56.2397, result.QN.get(1, 0), 56.2397 * tol);
        assertEquals(2.0233, result.QN.get(2, 0), 2.0233 * tol);

    }

    @Test
    public void D1SQ1TwoClassesOfMMPPSolvingByCTMCAndFluid() {
        SolverENV envSolverCTMC = SolverEnvTestFixtures.example_randomEnvironment_5();
        SolverResult CTMCresult = envSolverCTMC.runAnalyzerByCTMC();
        // System.out.println(CTMCresult.QN);
        SolverENV envSolver = SolverEnvTestFixtures.example_randomEnvironment_5();
        envSolver.getAvg();
        // envSolver.printAvgTable();
        SolverResult result = envSolver.result;

        // QN
        assertEquals(2, result.QN.getNumRows());
        assertEquals(2, result.QN.getNumCols());
        assertEquals(4, result.QN.getNumElements());
        assertEquals(CTMCresult.QN.value(), result.QN.value(), CTMCresult.QN.value() * tol);
        assertEquals(CTMCresult.QN.get(1, 0), result.QN.get(1, 0), CTMCresult.QN.get(1, 0) * tol);
        assertEquals(CTMCresult.QN.get(0, 1), result.QN.get(0, 1), CTMCresult.QN.get(0, 1) * tol);
        assertEquals(CTMCresult.QN.get(1, 1), result.QN.get(1, 1), CTMCresult.QN.get(1, 1) * tol);
    }

    @Test
    public void StateDependentD1SQ2stage2() {
        SolverENV envSolver = SolverEnvTestFixtures.example_randomEnvironment_6();
        envSolver.getAvg();
        // envSolver.printAvgTable();
        SolverResult result = envSolver.result;

        // QN
        assertEquals(3, result.QN.getNumRows());
        assertEquals(1, result.QN.getNumCols());
        assertEquals(3, result.QN.getNumElements());
        assertEquals(23.8671, result.QN.value(), 23.8671 * tol);
        assertEquals(74.9777, result.QN.get(1, 0), 74.9777 * tol);
        assertEquals(1.1685, result.QN.get(2, 0), 1.1685 * tol);
    }



    @Test
    @Disabled("38s, Test fails - value mismatch")
    public void StateDependentD1SQ3stage2() {
        SolverENV envSolver = SolverEnvTestFixtures.example_randomEnvironment_7();
        envSolver.setRef(0);
        envSolver.getAvg();
        // envSolver.printAvgTable();
        SolverResult result = envSolver.result;

        // QN
        assertEquals(4, result.QN.getNumRows());
        assertEquals(1, result.QN.getNumCols());
        assertEquals(4, result.QN.getNumElements());
        assertEquals(23.9966, result.QN.value(), 23.9966 * tol);
        assertEquals(73.5379, result.QN.get(1, 0), 73.5379 * tol);
        assertEquals(1.1808, result.QN.get(2, 0), 1.1808 * tol);
        assertEquals(1.3067, result.QN.get(3, 0), 1.3067 * tol);
    }

    @Test
    @Disabled("9s, Test fails - value mismatch")
    public void StateDependentD1SQ2stage2Parallel() {
        SolverENV envSolver = SolverEnvTestFixtures.example_randomEnvironment_8();
        envSolver.setRef(0);
        envSolver.getAvg();
        // envSolver.printAvgTable();
        SolverResult result = envSolver.result;

        // QN
        assertEquals(3, result.QN.getNumRows());
        assertEquals(1, result.QN.getNumCols());
        assertEquals(3, result.QN.getNumElements());
        assertEquals(47.999, result.QN.value(), 47.999 * tol);
        assertEquals(50.893, result.QN.get(1, 0), 50.893 * tol);
        assertEquals(1.164, result.QN.get(2, 0), 1.164 * tol);
    }

    @Test
    @Disabled("4s, Test fails - value mismatch")
    public void ErlangStateIndependentD1SQ2stage2() {
        SolverENV envSolver = SolverEnvTestFixtures.example_randomEnvironment_9();
        envSolver.getAvg();
        // envSolver.printAvgTable();
        SolverResult result = envSolver.result;

        // QN
        assertEquals(3, result.QN.getNumRows());
        assertEquals(1, result.QN.getNumCols());
        assertEquals(3, result.QN.getNumElements());
        assertEquals(42.0259, result.QN.value(), 42.0259 * tol);
        assertEquals(55.8969, result.QN.get(1, 0), 55.8969 * tol);
        assertEquals(2.0349, result.QN.get(2, 0), 2.0349 * tol);
    }


    @Test
    @Disabled("9s, Test fails - value mismatch")
    public void ErlangStateDependentD1SQ2stage2() {
        SolverENV envSolver = SolverEnvTestFixtures.example_randomEnvironment_9();
        envSolver.setStateDepMethod("statedep");
        envSolver.setRef(0);
        envSolver.getAvg();
        // envSolver.printAvgTable();
        SolverResult result = envSolver.result;

        // QN
        assertEquals(3, result.QN.getNumRows());
        assertEquals(1, result.QN.getNumCols());
        assertEquals(3, result.QN.getNumElements());
        assertEquals(42.2081, result.QN.value(), 42.2081 * tol);
        assertEquals(56.0814, result.QN.get(1, 0), 56.0814 * tol);
        assertEquals(2.0469, result.QN.get(2, 0), 2.0469 * tol);
    }

    @Test
    public void CTMCSolverMulticlassOneMMPPDelayTest() {
        SolverENV envSolverCTMC = SolverEnvTestFixtures.renv_twostages_repairmen1();
        SolverResult CTMCresult = envSolverCTMC.runAnalyzerByCTMC();
        // System.out.println(CTMCresult.QN);
        SolverENV envSolver = SolverEnvTestFixtures.renv_twostages_repairmen1();
        envSolver.setStateDepMethod("stateindep");
        envSolver.getAvg();
        // envSolver.printAvgTable();
        SolverResult result = envSolver.result;

        // QN
        assertEquals(1, result.QN.getNumRows());
        assertEquals(2, result.QN.getNumCols());
        assertEquals(2, result.QN.getNumElements());
        assertEquals(CTMCresult.QN.value(), result.QN.value(), CTMCresult.QN.value() * tol);
        assertEquals(CTMCresult.QN.get(0, 1), result.QN.get(0, 1), CTMCresult.QN.get(0, 1) * tol);
    }

    @Test
    public void CTMCSolverMulticlassOneExpOneMMPPTest() {
        SolverENV envSolverCTMC = SolverEnvTestFixtures.renv_twostages_repairmen2();
        SolverResult CTMCresult = envSolverCTMC.runAnalyzerByCTMC();
        // System.out.println(CTMCresult.QN);
        SolverENV envSolver = SolverEnvTestFixtures.renv_twostages_repairmen2();
        envSolver.setStateDepMethod("stateindep");
        envSolver.getAvg();
        // envSolver.printAvgTable();
        SolverResult result = envSolver.result;

        // QN
        assertEquals(2, result.QN.getNumRows());
        assertEquals(2, result.QN.getNumCols());
        assertEquals(4, result.QN.getNumElements());
        assertEquals(CTMCresult.QN.value(), result.QN.value(), CTMCresult.QN.value() * tol);
        assertEquals(CTMCresult.QN.get(0, 1), result.QN.get(0, 1), CTMCresult.QN.get(0, 1) * tol);
    }

    @Test
    @Disabled("33s, wrong final value")
    public void CTMCSolverMulticlassD1SQ2TwoClassesOfMMPPTest() {
        SolverENV envSolverCTMC = SolverEnvTestFixtures.renv_twostages_repairmen3();
        SolverResult CTMCresult = envSolverCTMC.runAnalyzerByCTMC();
        // System.out.println(CTMCresult.QN);
        SolverENV envSolver = SolverEnvTestFixtures.renv_twostages_repairmen3();
        envSolver.setStateDepMethod("stateindep");
        envSolver.getAvg();
        // envSolver.printAvgTable();
        SolverResult result = envSolver.result;

        // QN
        assertEquals(3, result.QN.getNumRows());
        assertEquals(2, result.QN.getNumCols());
        assertEquals(6, result.QN.getNumElements());
        assertEquals(CTMCresult.QN.value(), result.QN.value(), CTMCresult.QN.value() * tol);
        assertEquals(CTMCresult.QN.get(0, 1), result.QN.get(0, 1), CTMCresult.QN.get(0, 1) * tol);
        assertEquals(CTMCresult.QN.get(1, 0), result.QN.get(1, 0), CTMCresult.QN.get(1, 0) * tol);
        assertEquals(CTMCresult.QN.get(1, 1), result.QN.get(1, 1), CTMCresult.QN.get(1, 1) * tol);
        assertEquals(CTMCresult.QN.get(2, 0), result.QN.get(2, 0), CTMCresult.QN.get(2, 0) * tol);
        assertEquals(CTMCresult.QN.get(2, 1), result.QN.get(2, 1), CTMCresult.QN.get(2, 1) * tol);
    }

    @Test
    @Disabled("242s, incorrect value")
    public void CTMCSolverSingleClassErlangTest() {
        SolverENV envSolverCTMC = SolverEnvTestFixtures.renv_twostages_repairmen4();
        SolverResult CTMCresult = envSolverCTMC.runAnalyzerByCTMC();
        // System.out.println(CTMCresult.QN);
        SolverENV envSolver = SolverEnvTestFixtures.renv_twostages_repairmen4();
        envSolver.setStateDepMethod("stateindep");
        envSolver.getAvg();
        // envSolver.printAvgTable();
        SolverResult result = envSolver.result;

        // QN
        assertEquals(3, result.QN.getNumRows());
        assertEquals(1, result.QN.getNumCols());
        assertEquals(3, result.QN.getNumElements());
        assertEquals(CTMCresult.QN.value(), result.QN.value(), CTMCresult.QN.value() * tol);
        assertEquals(CTMCresult.QN.get(1, 0), result.QN.get(1, 0), CTMCresult.QN.get(1, 0) * tol);
        assertEquals(CTMCresult.QN.get(2, 0), result.QN.get(2, 0), CTMCresult.QN.get(2, 0) * tol);
    }

    @Test
    @Disabled("178s, incorrect value")
    public void CTMCSolverMultiClassErlangTest() {
        SolverENV envSolverCTMC = SolverEnvTestFixtures.renv_twostages_repairmen5();
        SolverResult CTMCresult = envSolverCTMC.runAnalyzerByCTMC();
        // System.out.println(CTMCresult.QN);
        SolverENV envSolver = SolverEnvTestFixtures.renv_twostages_repairmen5();
        envSolver.setStateDepMethod("stateindep");
        envSolver.getAvg();
        // envSolver.printAvgTable();
        SolverResult result = envSolver.result;

        // QN
        assertEquals(2, result.QN.getNumRows());
        assertEquals(2, result.QN.getNumCols());
        assertEquals(4, result.QN.getNumElements());
        assertEquals(CTMCresult.QN.value(), result.QN.value(), CTMCresult.QN.value() * tol);
        assertEquals(CTMCresult.QN.get(0, 1), result.QN.get(0, 1), CTMCresult.QN.get(0, 1) * tol);
        assertEquals(CTMCresult.QN.get(1, 0), result.QN.get(1, 0), CTMCresult.QN.get(1, 0) * tol);
        assertEquals(CTMCresult.QN.get(1, 1), result.QN.get(1, 1), CTMCresult.QN.get(1, 1) * tol);
    }

    @Test
    @Disabled("Stack overflow error")
    public void CTMCSolverSingleClassCoxianTest() {
        SolverENV envSolverCTMC = SolverEnvTestFixtures.renv_twostages_repairmen6();
        SolverResult CTMCresult = envSolverCTMC.runAnalyzerByCTMC();
        // System.out.println(CTMCresult.QN);
        SolverENV envSolver = SolverEnvTestFixtures.renv_twostages_repairmen6();
        envSolver.setStateDepMethod("stateindep");
        envSolver.getAvg();
        // envSolver.printAvgTable();
        SolverResult result = envSolver.result;

        // QN
        assertEquals(3, result.QN.getNumRows());
        assertEquals(1, result.QN.getNumCols());
        assertEquals(3, result.QN.getNumElements());
        assertEquals(CTMCresult.QN.value(), result.QN.value(), CTMCresult.QN.value() * tol);
        assertEquals(CTMCresult.QN.get(1, 0), result.QN.get(1, 0), CTMCresult.QN.get(1, 0) * tol);
        assertEquals(CTMCresult.QN.get(2, 0), result.QN.get(2, 0), CTMCresult.QN.get(2, 0) * tol);
    }

    /*@Test
    @Disabled("Exceeds 25 minutes runtime")
    public void CTMCSolverThreeClassMMPPTest() {
        SolverENV envSolverCTMC = SolverEnvTestFixtures.renv_twostages_repairmen7();
        SolverResult CTMCresult = envSolverCTMC.runAnalyzerByCTMC();
        // System.out.println(CTMCresult.QN);
        SolverENV envSolver = SolverEnvTestFixtures.renv_twostages_repairmen7();
        envSolver.setStateDepMethod("stateindep");
        envSolver.getAvg();
        // envSolver.printAvgTable();
        SolverResult result = envSolver.result;

        // QN
        assertEquals(3, result.QN.getNumRows());
        assertEquals(3, result.QN.getNumCols());
        assertEquals(9, result.QN.getNumElements());
        assertEquals(CTMCresult.QN.value(), result.QN.value(), CTMCresult.QN.value() * tol);
        assertEquals(CTMCresult.QN.get(0, 1), result.QN.get(0, 1), CTMCresult.QN.get(0, 1) * tol);
        assertEquals(CTMCresult.QN.get(1, 0), result.QN.get(1, 0), CTMCresult.QN.get(1, 0) * tol);
        assertEquals(CTMCresult.QN.get(1, 1), result.QN.get(1, 1), CTMCresult.QN.get(1, 1) * tol);
        assertEquals(CTMCresult.QN.get(2, 0), result.QN.get(2, 0), CTMCresult.QN.get(2, 0) * tol);
        assertEquals(CTMCresult.QN.get(2, 1), result.QN.get(2, 1), CTMCresult.QN.get(2, 1) * tol);
    }*/

    @Test
    public void getSamplePathTableReturnsCorrectStructure() {
        // Use the two-stage repairmen model
        SolverENV envSolver = SolverEnvTestFixtures.renv_twostages_repairmen();

        // Create a sample path: Stage1 for 5.0 time units, Stage2 for 10.0, back to Stage1 for 3.0
        List<Object[]> samplePath = new ArrayList<>();
        samplePath.add(new Object[]{"Stage1", 5.0});
        samplePath.add(new Object[]{"Stage2", 10.0});
        samplePath.add(new Object[]{"Stage1", 3.0});

        // Run the sample path analysis
        SamplePathResult result = envSolver.getSamplePathTable(samplePath);

        // Verify structure
        assertNotNull(result);
        assertNotNull(result.segments);
        assertEquals(3, result.segments.size());

        // Verify segment 0 (Stage1, duration 5.0)
        SamplePathResult.SamplePathSegment seg0 = result.segments.get(0);
        assertEquals(0, seg0.segmentIndex);
        assertEquals(0, seg0.stageIndex);
        assertEquals("Stage1", seg0.stageName);
        assertEquals(5.0, seg0.duration, 1e-10);
        assertNotNull(seg0.initialQ);
        assertNotNull(seg0.finalQ);
        assertNotNull(seg0.QNt);

        // Verify segment 1 (Stage2, duration 10.0)
        SamplePathResult.SamplePathSegment seg1 = result.segments.get(1);
        assertEquals(1, seg1.segmentIndex);
        assertEquals(1, seg1.stageIndex);
        assertEquals("Stage2", seg1.stageName);
        assertEquals(10.0, seg1.duration, 1e-10);

        // Verify segment 2 (Stage1, duration 3.0)
        SamplePathResult.SamplePathSegment seg2 = result.segments.get(2);
        assertEquals(2, seg2.segmentIndex);
        assertEquals(0, seg2.stageIndex);
        assertEquals("Stage1", seg2.stageName);
        assertEquals(3.0, seg2.duration, 1e-10);

        // Verify state continuity: final Q of segment N should equal initial Q of segment N+1
        int M = seg0.finalQ.getNumRows();
        int K = seg0.finalQ.getNumCols();
        for (int i = 0; i < M; i++) {
            for (int k = 0; k < K; k++) {
                assertEquals(seg0.finalQ.get(i, k), seg1.initialQ.get(i, k), 0.01,
                    "State continuity violated between segment 0 and 1 at (" + i + "," + k + ")");
                assertEquals(seg1.finalQ.get(i, k), seg2.initialQ.get(i, k), 0.01,
                    "State continuity violated between segment 1 and 2 at (" + i + "," + k + ")");
            }
        }
    }

    @Test
    public void getSamplePathTableWorksWithStageIndices() {
        // Use the two-stage repairmen model
        SolverENV envSolver = SolverEnvTestFixtures.renv_twostages_repairmen();

        // Create a sample path using 0-based indices instead of names
        List<Object[]> samplePath = new ArrayList<>();
        samplePath.add(new Object[]{0, 5.0});  // Stage1
        samplePath.add(new Object[]{1, 10.0}); // Stage2

        // Run the sample path analysis
        SamplePathResult result = envSolver.getSamplePathTable(samplePath);

        // Verify structure
        assertNotNull(result);
        assertEquals(2, result.segments.size());
        assertEquals("Stage1", result.segments.get(0).stageName);
        assertEquals("Stage2", result.segments.get(1).stageName);
    }

    @Test
    public void getSamplePathTableThrowsOnEmptyPath() {
        SolverENV envSolver = SolverEnvTestFixtures.renv_twostages_repairmen();

        List<Object[]> emptyPath = new ArrayList<>();

        assertThrows(IllegalArgumentException.class, () -> {
            envSolver.getSamplePathTable(emptyPath);
        });
    }

    @Test
    public void getSamplePathTableThrowsOnInvalidStageName() {
        SolverENV envSolver = SolverEnvTestFixtures.renv_twostages_repairmen();

        List<Object[]> samplePath = new ArrayList<>();
        samplePath.add(new Object[]{"NonExistentStage", 5.0});

        assertThrows(IllegalArgumentException.class, () -> {
            envSolver.getSamplePathTable(samplePath);
        });
    }

    @Test
    public void getSamplePathTableThrowsOnNegativeDuration() {
        SolverENV envSolver = SolverEnvTestFixtures.renv_twostages_repairmen();

        List<Object[]> samplePath = new ArrayList<>();
        samplePath.add(new Object[]{"Stage1", -5.0});

        assertThrows(IllegalArgumentException.class, () -> {
            envSolver.getSamplePathTable(samplePath);
        });
    }
}
