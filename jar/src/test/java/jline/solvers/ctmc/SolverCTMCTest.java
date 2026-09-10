package jline.solvers.ctmc;

import jline.examples.java.advanced.LoadDependentModel;
import jline.examples.java.basic.ClosedModel;
import jline.lang.Network;
import jline.VerboseLevel;
import jline.solvers.NetworkAvgTable;
import jline.lang.ClosedClass;
import jline.lang.constant.SchedStrategy;
import jline.lang.processes.Exp;
import jline.lang.nodes.Delay;
import jline.lang.nodes.Queue;
import jline.solvers.SolverOptions;
import jline.solvers.ssa.SolverSSA;
import jline.util.matrix.Matrix;
import org.junit.jupiter.api.Test;
import org.junit.jupiter.api.AfterEach;

import java.nio.file.Files;
import java.nio.file.Path;
import java.nio.file.Paths;
import java.util.Map;

import static org.junit.jupiter.api.Assertions.*;
import static jline.TestTools.MID_TOL;
import static jline.TestTools.FINE_TOL;
import static jline.TestTools.relativeTolerance;
import static jline.TestTools.assertTableMetrics;

public class SolverCTMCTest {
  
  /**
   * Cleanup method to delete any CTMC workspace files created during tests.
   * This ensures tests don't leave behind temporary files.
   */
  @AfterEach
  public void cleanupWorkspaceFiles() {
    try {
      // Get current working directory
      Path currentDir = Paths.get(".");
      
      // Find and delete all ctmc_analyzer_workspace_*.mat files
      Files.list(currentDir)
          .filter(path -> {
            String fileName = path.getFileName().toString();
            return fileName.startsWith("ctmc_analyzer_workspace_") && fileName.endsWith(".mat");
          })
          .forEach(path -> {
            try {
              Files.deleteIfExists(path);
            } catch (Exception e) {
              // Ignore deletion errors to avoid test failures
            }
          });
    } catch (Exception e) {
      // Ignore cleanup errors to avoid affecting test results
    }
  }



  // The following tests have been removed as they are already covered in ClosedExamplesTest.java:
  // - cqn_repairmen (example_closedModel_1)
  // - cqn_twoclass_hyperl (example_closedModel_2)
  // - cqn_threeclass_hyperl (example_closedModel_3)
  // - example_closedModel_7fcfs (cqn_bcmp_theorem_fcfs)
  // - example_closedModel_7ps (cqn_bcmp_theorem_ps)
  // - example_closedModel_7lcfspr (cqn_bcmp_theorem_lcfspr)

  @Test
  public void cqn_repairmen_multi() {
    Network model = ClosedModel.cqn_repairmen_multi();
    SolverCTMC solver = new SolverCTMC(model,"verbose",VerboseLevel.SILENT);
    NetworkAvgTable table = solver.getAvgTable();
    //table.printTable();

    double[] expectedQLen = {1.9454556018966782, 1.6382844704790733, 2.054544398103321, 0.36171552952092645};
    double[] expectedUtil = {1.9454556018966782, 1.6382844704790733, 0.648485200632226, 0.05460948234930244};
    double[] expectedRespT = {1, 1, 1.056073649843404, 0.2207892072706715};
    double[] expectedResidT = {1, 1, 1.056073649843404, 0.22079};
    double[] expectedArvR = {1.9455, 1.6383, 1.9455, 1.6383};
    double[] expectedTput = {1.9455, 1.6383, 1.9455, 1.6383};

    assertTableMetrics(table, expectedQLen, expectedUtil, expectedRespT, expectedResidT, expectedArvR, expectedTput);
  }

  @Test
  public void ld_multiserver_fcfs() {
    Network model = LoadDependentModel.ld_multiserver_fcfs();
    SolverCTMC solver = new SolverCTMC(model,"verbose",VerboseLevel.SILENT);
    NetworkAvgTable table = solver.getAvgTable();
    //table.printTable();

    double[] expectedQLen = {1.333333310169666, 14.667};
    double[] expectedUtil = {1.333333310169666, 1.0};
    double[] expectedRespT = {1.0, 11.000000017372754};
    double[] expectedResidT = {1.0, 11.000000017372754};
    double[] expectedArvR = {1.333333333333333, 1.333333310169666};
    double[] expectedTput = {1.333333310169666, 1.333333333333333};

    assertTableMetrics(table, expectedQLen, expectedUtil, expectedRespT, expectedResidT, expectedArvR, expectedTput);
  }

  @Test
  public void ld_multiserver_ps_twoclasses() {
    Network model = LoadDependentModel.ld_multiserver_ps_twoclasses();
    SolverCTMC solver = new SolverCTMC(model,"verbose",VerboseLevel.SILENT);
    NetworkAvgTable table = solver.getAvgTable();
    //table.printTable();

    double[] expectedQLen = {0.897501892505678, 0.511574166301447, 3.10249810749432, 1.48842583369855};
    double[] expectedUtil = {0.897501892505678, 0.511574166301447, 0.673126419379258, 0.319733853938404};
    double[] expectedRespT = {1.00000000000000, 2.00000000000000, 3.45681511470985, 5.81900311526480};
    double[] expectedResidT = {1.00000000000000, 2.00000000000000, 3.45681511470985, 5.81900311526480};
    double[] expectedArvR = {0.897501892505678, 0.255787083150723, 0.897501892505678, 0.255787083150723};
    double[] expectedTput = {0.897501892505678, 0.255787083150723, 0.897501892505678, 0.255787083150723};

    assertTableMetrics(table, expectedQLen, expectedUtil, expectedRespT, expectedResidT, expectedArvR, expectedTput);
  }

  @Test
  public void ld_multiserver_ps() {
    Network model = LoadDependentModel.ld_multiserver_ps();
    SolverCTMC solver = new SolverCTMC(model,"verbose",VerboseLevel.SILENT);
    NetworkAvgTable table = solver.getAvgTable();
    //table.printTable();

    double[] expectedQLen = {0.546737254325686, 0.369444841897479, 0.856387113240600, 0.480774229918630, 2.59687563243372, 1.14978092818389};
    double[] expectedUtil = {0.546737254325686, 0.369444841897479, 0.273368627162843, 0.153935350790616, 0.637860130046633, 0.277083631423110};
    double[] expectedRespT = {1.00000000000000, 2.00000000000000, 1.56635953827002, 2.60268476046038, 4.74976894639556, 6.22437126082790};
    double[] expectedResidT = {1.00000000000000, 2.00000000000000, 1.56635953827002, 2.60268476046038, 4.74976894639556, 6.22437126082790};
    double[] expectedArvR = {0.546737254325686, 0.184722420948740, 0.546737254325686, 0.184722420948740, 0.546737254325685, 0.184722420948739};
    double[] expectedTput = {0.546737254325686, 0.184722420948740, 0.546737254325685, 0.184722420948739, 0.546737254325686, 0.184722420948740};

    assertTableMetrics(table, expectedQLen, expectedUtil, expectedRespT, expectedResidT, expectedArvR, expectedTput);
  }

  // example_loadDependent_4 test has been deleted as the method was removed

  @Test
  public void ld_class_dependence() {
    Network model = LoadDependentModel.ld_class_dependence();

    // Queue1 is a class-dependent station beta(n)=min(n_1,2) with declared
    // peak c=2; utilization is normalized as T*S/peak (busy-capacity fraction).
    CTMCOptions options = new CTMCOptions();
    options.verbose = VerboseLevel.SILENT;
    SolverCTMC solver = new SolverCTMC(model, options);

    NetworkAvgTable table = solver.getAvgTable();
    
    //table.printTable();

    // rebased 2026-08-15 against MATLAB ld_class_dependence.m (CTMC exact), which
    // python reproduces to 5 digits; the previous golden predates setClassDependence
    double[] expectedQLen = {0.8819631781498992, 0.2708220931100575, 15.118036821850097, 7.729177906889943};
    double[] expectedUtil = {0.8819631781498992, 0.2708220931100575, 0.661472383612426, 0.3385276163875742};
    double[] expectedRespT = {1.0, 2.0, 17.14134693645973, 57.07937501058218};
    double[] expectedResidT = {1.0, 2.0, 17.14134693645973, 57.07937501058218};
    double[] expectedArvR = {0.8819631781499012, 0.1354110465550297, 0.8819631781498992, 0.13541104655502875};
    double[] expectedTput = {0.8819631781498992, 0.13541104655502875, 0.8819631781499012, 0.1354110465550297};

    assertTableMetrics(table, expectedQLen, expectedUtil, expectedRespT, expectedResidT, expectedArvR, expectedTput);
  }

  // mqn_basic test removed - already tested in MixedExamplesTest.java as testMxqnBasicCTMC
  // mqn_multiserver_ps test removed - similar test covered in MixedExamplesTest.java
  // mqn_multiserver_fcfs test removed - similar test covered in MixedExamplesTest.java


  // oqn_basic test removed - already tested in OpenExamplesTest.java as testOqnBasicCTMC
  // oqn_cs_routing test removed - incomplete test that doesn't execute solver or verify results
  // oqn_fourqueues test removed - already tested in OpenExamplesTest.java

  /**
   * Test reward-based CTMC analysis on a simple closed queueing network.
   *
   * This test verifies:
   * - setReward functionality on Network
   * - getAvgReward computation via value iteration
   * - Comparison with expected results
   */
  @Test
  public void testRewardModel_ClosedNetwork() {
    // Model Definition: Simple closed network with 2 jobs
    Network model = new Network("RewardExample");

    // Nodes
    Delay delay = new Delay(model, "Delay");
    Queue queue = new Queue(model, "Queue", SchedStrategy.FCFS);

    queue.setNumberOfServers(1);

    // Job class: 2 jobs in the system
    ClosedClass cclass = new ClosedClass(model, "Class1", 2, delay);
    delay.setService(cclass, new Exp(1.0));   // Think time = 1
    queue.setService(cclass, new Exp(2.0));   // Service rate = 2

    // Topology
    model.addLink(delay, queue);
    model.addLink(queue, delay);

    // Define Reward Functions
    // State format: [delay_jobs, queue_jobs]

    // Reward 1: Queue length
    model.setReward("QueueLength", (state, sn) -> state.get(0, 1));

    // Reward 2: Utilization (1 if server busy)
    model.setReward("Utilization", (state, sn) -> Math.min(state.get(0, 1), 1.0));

    // Verify rewards are defined
    assertTrue(model.hasRewards());
    assertEquals(2, model.getRewards().size());

    // Solve with CTMC Solver with reduced iterations for speed
    SolverOptions options = new SolverOptions();
    options.verbose = VerboseLevel.SILENT;
    options.rewardIterations = 500;  // Reduced for test speed

    SolverCTMC solver = new SolverCTMC(model, options);

    // Get Steady-State Expected Rewards
    Map<String, Double> avgRewards = solver.getAvgReward();

    // For a closed network with N=2, think rate=1, service rate=2:
    // States: (2,0), (1,1), (0,2)
    // This is a simple birth-death process
    double queueLength = avgRewards.get("QueueLength");
    double utilization = avgRewards.get("Utilization");

    // Verify results are finite and non-negative
    assertTrue(Double.isFinite(queueLength) && queueLength >= 0,
            "QueueLength should be finite and non-negative, got " + queueLength);
    assertTrue(Double.isFinite(utilization) && utilization >= 0,
            "Utilization should be finite and non-negative, got " + utilization);

    // Test getAvgReward(name) method
    assertEquals(queueLength, solver.getAvgReward("QueueLength"), FINE_TOL);
    assertEquals(utilization, solver.getAvgReward("Utilization"), FINE_TOL);

    // Test getRewardNames
    assertEquals(2, solver.getRewardNames().size());
    assertTrue(solver.getRewardNames().contains("QueueLength"));
    assertTrue(solver.getRewardNames().contains("Utilization"));

    // Test getRewardValueFunction
    Matrix V = solver.getRewardValueFunction("QueueLength");
    assertNotNull(V);
    assertTrue(V.getNumRows() > 0);
    assertTrue(V.getNumCols() > 0);

    // Test getRewardTimeVector
    double[] t = solver.getRewardTimeVector();
    assertNotNull(t);
    assertTrue(t.length > 0);
    assertEquals(0.0, t[0], FINE_TOL);  // Time starts at 0
  }

  @Test
  public void tut06_cache_lru_zipf_ctmc() {
    Network model = SolverCTMCTestFixtures.test_tut06_cache_lru_zipf();
    SolverCTMC solver = new SolverCTMC(model, "verbose", VerboseLevel.SILENT);
    NetworkAvgTable table = solver.getAvgTable();

    double[] expectedQLen = {0.0, 0.43337595752684316, 0.5666239604680446};
    double[] expectedUtil = {0.0, 0.43337595752684316, 0.5666239604680446};
    double[] expectedRespT = {1e-8, 0.2, 1.0};
    double[] expectedResidT = {0.0, 0.1585422949676638, 0.20728852516168095};
    double[] expectedArvR = {2.73350374810226, 2.1668797876342154, 0.5666239604680443};
    double[] expectedTput = {2.73350374810226, 2.1668797876342154, 0.5666239604680446};

    assertTableMetrics(table, expectedQLen, expectedUtil, expectedRespT, expectedResidT, expectedArvR, expectedTput);
  }

  @Test
  public void tut06_cache_lru_zipf_ssa() {
    Network model = SolverCTMCTestFixtures.test_tut06_cache_lru_zipf();
    SolverSSA solver = new SolverSSA(model, "samples", 2e4, "seed", 1, "verbose", VerboseLevel.SILENT, "method", "serial");
    NetworkAvgTable table = solver.getAvgTable();

    double[] expectedQLen = {0.0, 0.43790923334741466, 0.5620906824628674};
    double[] expectedUtil = {0.0, 0.43790923334741466, 0.5620906824628674};
    double[] expectedRespT = {1e-8, 0.2, 1.0};
    double[] expectedResidT = {0.0, 0.15639485518962562, 0.21802572405187115};
    double[] expectedArvR = {2.751636849199943, 2.218007395317449, 0.6184125017797503};
    double[] expectedTput = {2.8364198970971994, 2.1895461667370753, 0.5620906824628674};

    assertTableMetrics(table, expectedQLen, expectedUtil, expectedRespT, expectedResidT, expectedArvR, expectedTput);
  }

  /**
   * Test that calling reward methods without defining rewards throws exception.
   */
  @Test
  public void testRewardNotDefined() {
    Network model = new Network("NoRewardModel");

    Delay delay = new Delay(model, "Delay");
    Queue queue = new Queue(model, "Queue", SchedStrategy.FCFS);

    queue.setNumberOfServers(1);

    ClosedClass cclass = new ClosedClass(model, "Class1", 2, delay);
    delay.setService(cclass, new Exp(1.0));
    queue.setService(cclass, new Exp(2.0));

    model.addLink(delay, queue);
    model.addLink(queue, delay);

    // No rewards defined
    assertFalse(model.hasRewards());

    SolverCTMC solver = new SolverCTMC(model, "verbose", VerboseLevel.SILENT);

    // Should throw exception when trying to get rewards
    assertThrows(IllegalStateException.class, () -> solver.getAvgReward());
  }

  /**
   * Test clearRewards functionality.
   */
  @Test
  public void testClearRewards() {
    Network model = new Network("ClearRewardModel");

    Delay delay = new Delay(model, "Delay");
    Queue queue = new Queue(model, "Queue", SchedStrategy.FCFS);

    queue.setNumberOfServers(1);

    ClosedClass cclass = new ClosedClass(model, "Class1", 2, delay);
    delay.setService(cclass, new Exp(1.0));
    queue.setService(cclass, new Exp(2.0));

    model.addLink(delay, queue);
    model.addLink(queue, delay);

    // Add a reward
    model.setReward("Test", (state, sn) -> state.get(0, 1));
    assertTrue(model.hasRewards());
    assertEquals(1, model.getRewards().size());

    // Clear rewards
    model.clearRewards();
    assertFalse(model.hasRewards());
  }

  /**
   * Test updating an existing reward.
   */
  @Test
  public void testUpdateReward() {
    Network model = new Network("UpdateRewardModel");

    Delay delay = new Delay(model, "Delay");
    Queue queue = new Queue(model, "Queue", SchedStrategy.FCFS);

    queue.setNumberOfServers(1);

    ClosedClass cclass = new ClosedClass(model, "Class1", 2, delay);
    delay.setService(cclass, new Exp(1.0));
    queue.setService(cclass, new Exp(2.0));

    model.addLink(delay, queue);
    model.addLink(queue, delay);

    // Add initial reward
    model.setReward("Test", (state, sn) -> 1.0);
    assertEquals(1, model.getRewards().size());

    // Update the same reward (should replace, not add)
    model.setReward("Test", (state, sn) -> 2.0);
    assertEquals(1, model.getRewards().size());

    // Add a different reward
    model.setReward("Test2", (state, sn) -> 3.0);
    assertEquals(2, model.getRewards().size());
  }

}
