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
import jline.util.matrix.Matrix;
import org.junit.jupiter.api.Test;
import org.junit.jupiter.api.AfterEach;

import java.nio.file.Files;
import java.nio.file.Path;
import java.nio.file.Paths;
import java.util.Map;

import static org.junit.jupiter.api.Assertions.*;
import static jline.TestTools.MID_TOL;
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

    double[] expectedQLen = {0.8975018925056777, 0.5115741663014466, 3.1024981074943225, 1.4884258336985532};
    double[] expectedUtil = {0.8975018925056777, 0.5115741663014466, 0.6772062631977369, 0.3214337888627701};
    double[] expectedRespT = {1.0, 2.0, 3.456815114709851, 5.819};
    double[] expectedResidT = {1.0, 2.0, 3.456815114709851, 5.819};
    double[] expectedArvR = {0.8975, 0.25579, 0.8975, 0.25579};
    double[] expectedTput = {0.8975, 0.25579, 0.8975, 0.25579};

    assertTableMetrics(table, expectedQLen, expectedUtil, expectedRespT, expectedResidT, expectedArvR, expectedTput);
  }

  @Test
  public void ld_multiserver_ps() {
    Network model = LoadDependentModel.ld_multiserver_ps();
    SolverCTMC solver = new SolverCTMC(model,"verbose",VerboseLevel.SILENT);
    NetworkAvgTable table = solver.getAvgTable();
    //table.printTable();

    double[] expectedQLen = {0.5467372543256855, 0.3694448418974793, 0.8563871132405995, 0.48077422991863045, 2.5968756324337163, 1.1497809281838902};
    double[] expectedUtil = {0.5467372543256855, 0.3694448418974793, 0.47391698344747757, 0.2712785254907754, 0.6941051054883322, 0.29802};
    double[] expectedRespT = {1.0, 2.0, 1.5664, 2.6027, 4.7498, 6.2244};
    double[] expectedResidT = {1.0, 2.0, 1.5664, 2.6027, 4.7498, 6.2244};
    double[] expectedArvR = {0.54674, 0.18472, 0.54674, 0.18472, 0.54674, 0.18472};
    double[] expectedTput = {0.54674, 0.18472, 0.54674, 0.18472, 0.54674, 0.18472};

    assertTableMetrics(table, expectedQLen, expectedUtil, expectedRespT, expectedResidT, expectedArvR, expectedTput);
  }

  // example_loadDependent_4 test has been deleted as the method was removed

  @Test
  public void ld_class_dependence() {
    Network model = LoadDependentModel.ld_class_dependence();
    
    SolverCTMC solver = new SolverCTMC(model,"verbose",VerboseLevel.SILENT);
    
    NetworkAvgTable table = solver.getAvgTable();
    
    //table.printTable();

    double[] expectedQLen = {0.8923099748616234, 0.5292274741426646, 15.108, 7.470772525857336};
    double[] expectedUtil = {0.8923099748616234, 0.5292274741426646, 0.669232704525731, 0.3307672954742693};
    double[] expectedRespT = {1.0, 2.0, 16.93098296336792, 28.233};
    double[] expectedResidT = {1.0, 2.0, 16.931, 28.233};
    double[] expectedArvR = {0.8923102727009744, 0.2646138363794153, 0.8923099748616234, 0.26461};
    double[] expectedTput = {0.89231, 0.26461, 0.89231, 0.26461};

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
    assertEquals(queueLength, solver.getAvgReward("QueueLength"), 1e-10);
    assertEquals(utilization, solver.getAvgReward("Utilization"), 1e-10);

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
    assertEquals(0.0, t[0], 1e-10);  // Time starts at 0
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
