package jline.solvers.des;

import java.util.List;

import jline.GlobalConstants;
import jline.VerboseLevel;
import jline.lang.Network;
import jline.lang.constant.SolverType;
import jline.solvers.NetworkAvgTable;
import jline.solvers.SolverOptions;
import jline.solvers.jmt.SolverJMT;
import jline.solvers.mva.SolverMVA;
import org.junit.jupiter.api.BeforeEach;

import static org.junit.jupiter.api.Assertions.assertNotNull;
import static org.junit.jupiter.api.Assertions.assertTrue;

/**
 * Base class for SolverDES test files.
 *
 * Provides shared infrastructure:
 * - Constants (BASE_SEED, REL_ERROR_TOL)
 * - Helper methods for creating DESOptions
 * - Validation methods for comparing simulation results
 *
 * All DES test files (SolverDESTestCore, SolverDESTestScheduling, etc.) extend this class
 * to share validation logic and avoid duplication.
 */
public abstract class SolverDESTestFixtures {

	/** Base seed for reproducibility across all tests */
	protected static final int BASE_SEED = 23000;

	/** Relative error tolerance for simulation results (5% for stochastic comparisons) */
	protected static final double REL_ERROR_TOL = 5e-2;

	/**
	 * Sets up verbose level to SILENT before each test.
	 * This ensures priority messages and other verbose output don't clutter test output.
	 */
	@BeforeEach
	public void setUpVerboseLevel() {
		GlobalConstants.setVerbose(VerboseLevel.SILENT);
	}

	/**
	 * Creates a DESOptions instance with standard test configuration.
	 * Uses default samples (200,000) from DESOptions.
	 *
	 * @return configured DESOptions for testing
	 */
	protected static DESOptions createDefaultTestOptions() {
		DESOptions options = new DESOptions();
		options.verbose = VerboseLevel.SILENT;
		options.seed = BASE_SEED;
		// Uses default samples from DESOptions (200,000)
		return options;
	}

	/**
	 * Creates a DESOptions instance with custom sample count.
	 * Use only when test requires different sample size than default.
	 *
	 * @param samples number of samples
	 * @return configured DESOptions for testing
	 */
	protected static DESOptions createTestOptions(int samples) {
		DESOptions options = createDefaultTestOptions();
		options.samples = samples;
		return options;
	}

	/**
	 * Validates DES results against exact MVA results using relative error.
	 *
	 * @param testName name of the test (for error messages)
	 * @param desModel network model to test with DES
	 */
	protected void validateDESAgainstMVA(String testName, Network desModel) {
		// Get exact results from MVA solver
		Network mvaModel = (Network) desModel.copy();
		SolverOptions mvaOptions = new SolverOptions(SolverType.MVA);
		mvaOptions.verbose = VerboseLevel.SILENT;
		SolverMVA solverMVA = new SolverMVA(mvaModel, mvaOptions);
		NetworkAvgTable mvaResults = solverMVA.getAvgTable();

		// Run DES simulation
		DESOptions desOptions = createDefaultTestOptions();
		SolverDES solverDES = new SolverDES(desModel, desOptions);
		NetworkAvgTable desResults = solverDES.getAvgTable();

		// Compare results against MVA
		List<Double> mvaQLen = mvaResults.getQLen();
		List<Double> mvaUtil = mvaResults.getUtil();
		List<Double> mvaTput = mvaResults.getTput();
		List<Double> mvaRespT = mvaResults.getRespT();

		List<Double> desQLen = desResults.getQLen();
		List<Double> desUtil = desResults.getUtil();
		List<Double> desTput = desResults.getTput();
		List<Double> desRespT = desResults.getRespT();

		// Validate each metric
		for (int i = 0; i < mvaQLen.size(); i++) {
			double mvaQ = mvaQLen.get(i);
			double desQ = desQLen.get(i);
			if (mvaQ > 1e-6) {
				double relErr = Math.abs(desQ - mvaQ) / Math.abs(mvaQ);
				assertTrue(relErr <= REL_ERROR_TOL,
						testName + ": Queue " + i + " QLen relative error " + String.format("%.2e", relErr)
								+ " exceeds tolerance " + String.format("%.2e", REL_ERROR_TOL));
			}

			double mvaU = mvaUtil.get(i);
			double desU = desUtil.get(i);
			if (mvaU > 1e-6) {
				double relErr = Math.abs(desU - mvaU) / Math.abs(mvaU);
				assertTrue(relErr <= REL_ERROR_TOL,
						testName + ": Queue " + i + " Util relative error " + String.format("%.2e", relErr)
								+ " exceeds tolerance " + String.format("%.2e", REL_ERROR_TOL));
			}

			double mvaT = mvaTput.get(i);
			double desT = desTput.get(i);
			if (mvaT > 1e-6) {
				double relErr = Math.abs(desT - mvaT) / Math.abs(mvaT);
				assertTrue(relErr <= REL_ERROR_TOL,
						testName + ": Queue " + i + " Tput relative error " + String.format("%.2e", relErr)
								+ " exceeds tolerance " + String.format("%.2e", REL_ERROR_TOL));
			}

			double mvaR = mvaRespT.get(i);
			double desR = desRespT.get(i);
			if (mvaR > 1e-6) {
				double relErr = Math.abs(desR - mvaR) / Math.abs(mvaR);
				assertTrue(relErr <= REL_ERROR_TOL,
						testName + ": Queue " + i + " RespT relative error " + String.format("%.2e", relErr)
								+ " exceeds tolerance " + String.format("%.2e", REL_ERROR_TOL));
			}
		}
	}

	/**
	 * Validates DES results against static reference values using relative error.
	 * Reference values were pre-computed using SolverJMT.
	 * Uses default samples from DESOptions.
	 *
	 * @param testName name of the test (for error messages)
	 * @param desModel network model to test with DES
	 * @param refQLen reference queue length values
	 * @param refUtil reference utilization values
	 * @param refTput reference throughput values
	 * @param refRespT reference response time values
	 */
	protected void validateDESAgainstStaticRef(String testName, Network desModel, double[] refQLen,
			double[] refUtil, double[] refTput, double[] refRespT) {
		validateDESAgainstStaticRef(testName, desModel, refQLen, refUtil, refTput, refRespT,
				createDefaultTestOptions());
	}

	/**
	 * Validates DES results against static reference values using relative error.
	 * Reference values were pre-computed using SolverJMT.
	 *
	 * @param testName name of the test (for error messages)
	 * @param desModel network model to test with DES
	 * @param refQLen reference queue length values
	 * @param refUtil reference utilization values
	 * @param refTput reference throughput values
	 * @param refRespT reference response time values
	 * @param samples number of simulation samples
	 */
	protected void validateDESAgainstStaticRef(String testName, Network desModel, double[] refQLen,
			double[] refUtil, double[] refTput, double[] refRespT, int samples) {
		validateDESAgainstStaticRef(testName, desModel, refQLen, refUtil, refTput, refRespT,
				createTestOptions(samples));
	}

	/**
	 * Validates DES results against static reference values using relative error.
	 * Reference values were pre-computed using SolverJMT.
	 *
	 * @param testName name of the test (for error messages)
	 * @param desModel network model to test with DES
	 * @param refQLen reference queue length values
	 * @param refUtil reference utilization values
	 * @param refTput reference throughput values
	 * @param refRespT reference response time values
	 * @param desOptions DES solver options
	 */
	protected void validateDESAgainstStaticRef(String testName, Network desModel, double[] refQLen,
			double[] refUtil, double[] refTput, double[] refRespT, DESOptions desOptions) {
		// Run DES simulation
		SolverDES solverDES = new SolverDES(desModel, desOptions);
		NetworkAvgTable desResults = solverDES.getAvgTable();

		List<Double> desQLen = desResults.getQLen();
		List<Double> desUtil = desResults.getUtil();
		List<Double> desTput = desResults.getTput();
		List<Double> desRespT = desResults.getRespT();

		// Validate each metric
		for (int i = 0; i < refQLen.length; i++) {
			double refQ = refQLen[i];
			double desQ = desQLen.get(i);
			if (refQ > 1e-6) {
				double relErr = Math.abs(desQ - refQ) / Math.abs(refQ);
				assertTrue(relErr <= REL_ERROR_TOL, testName + ": Queue " + i + " QLen relative error "
						+ String.format("%.2e", relErr) + " exceeds tolerance "
						+ String.format("%.2e", REL_ERROR_TOL));
			}

			double refU = refUtil[i];
			double desU = desUtil.get(i);
			if (refU > 1e-6) {
				double relErr = Math.abs(desU - refU) / Math.abs(refU);
				assertTrue(relErr <= REL_ERROR_TOL, testName + ": Queue " + i + " Util relative error "
						+ String.format("%.2e", relErr) + " exceeds tolerance "
						+ String.format("%.2e", REL_ERROR_TOL));
			}

			double refT = refTput[i];
			double desT = desTput.get(i);
			if (refT > 1e-6) {
				double relErr = Math.abs(desT - refT) / Math.abs(refT);
				assertTrue(relErr <= REL_ERROR_TOL, testName + ": Queue " + i + " Tput relative error "
						+ String.format("%.2e", relErr) + " exceeds tolerance "
						+ String.format("%.2e", REL_ERROR_TOL));
			}

			double refR = refRespT[i];
			double desR = desRespT.get(i);
			if (refR > 1e-6) {
				double relErr = Math.abs(desR - refR) / Math.abs(refR);
				assertTrue(relErr <= REL_ERROR_TOL, testName + ": Queue " + i + " RespT relative error "
						+ String.format("%.2e", relErr) + " exceeds tolerance "
						+ String.format("%.2e", REL_ERROR_TOL));
			}
		}
	}

	/**
	 * Validates two DES simulations produce identical results.
	 * Both models are run with the same seed, so results should match exactly.
	 * Uses very tight tolerance (1e-10) since deterministic match is expected.
	 *
	 * @param testName name of the test (for error messages)
	 * @param model1 first network model
	 * @param model2 second network model
	 */
	protected void validateDESAgainstDES(String testName, Network model1, Network model2) {
		DESOptions options = createDefaultTestOptions();

		SolverDES solver1 = new SolverDES(model1, options);
		NetworkAvgTable results1 = solver1.getAvgTable();

		SolverDES solver2 = new SolverDES(model2, options);
		NetworkAvgTable results2 = solver2.getAvgTable();

		List<Double> qlen1 = results1.getQLen();
		List<Double> qlen2 = results2.getQLen();

		// Results should be identical because logic and seed are identical
		for (int i = 0; i < qlen1.size(); i++) {
			double v1 = qlen1.get(i);
			double v2 = qlen2.get(i);
			if (v1 > 1e-6) {
				double relErr = Math.abs(v1 - v2) / Math.abs(v1);
				// Use very tight tolerance because it should be deterministic match
				assertTrue(relErr <= 1e-10,
						testName + ": Queue " + i + " QLen mismatch. Val1=" + v1 + ", Val2=" + v2);
			}
		}
	}

	/**
	 * Validates DES tardiness metrics against JMT tardiness metrics using relative error.
	 * Compares both station-level tardiness (Tard) and system-level tardiness (SysTard).
	 * Uses default samples from DESOptions.
	 *
	 * @param testName name of the test (for error messages)
	 * @param desModel network model to test with DES (must have soft deadlines)
	 */
	protected void validateDESTardinessAgainstJMT(String testName, Network desModel) {
		validateDESTardinessAgainstJMT(testName, desModel, createDefaultTestOptions());
	}

	/**
	 * Validates DES tardiness metrics against JMT tardiness metrics using relative error.
	 * Compares both station-level tardiness (Tard) and system-level tardiness (SysTard).
	 *
	 * @param testName name of the test (for error messages)
	 * @param desModel network model to test with DES (must have soft deadlines)
	 * @param samples number of simulation samples
	 */
	protected void validateDESTardinessAgainstJMT(String testName, Network desModel, int samples) {
		validateDESTardinessAgainstJMT(testName, desModel, createTestOptions(samples));
	}

	/**
	 * Validates DES tardiness metrics against JMT tardiness metrics using relative error.
	 * Compares both station-level tardiness (Tard) and system-level tardiness (SysTard).
	 *
	 * @param testName name of the test (for error messages)
	 * @param desModel network model to test with DES (must have soft deadlines)
	 * @param desOptions DES solver options
	 */
	protected void validateDESTardinessAgainstJMT(String testName, Network desModel, DESOptions desOptions) {
		// Get JMT tardiness results
		Network jmtModel = (Network) desModel.copy();
		SolverOptions jmtOptions = new SolverOptions(SolverType.JMT);
		jmtOptions.verbose = VerboseLevel.SILENT;
		jmtOptions.seed = desOptions.seed;
		jmtOptions.samples = desOptions.samples;
		SolverJMT solverJMT = new SolverJMT(jmtModel, jmtOptions);
		NetworkAvgTable jmtDeadlineTable = solverJMT.getDeadlineTable();

		// Get DES tardiness results
		SolverDES solverDES = new SolverDES(desModel, desOptions);
		NetworkAvgTable desDeadlineTable = solverDES.getDeadlineTable();

		// Validate that both solvers returned tardiness data
		assertNotNull(jmtDeadlineTable, testName + ": JMT did not return tardiness data");
		assertNotNull(desDeadlineTable, testName + ": DES did not return tardiness data");

		// Compare results
		List<Double> jmtRespT = jmtDeadlineTable.getRespT();
		List<Double> jmtTard = jmtDeadlineTable.getTard();
		List<Double> jmtSysTard = jmtDeadlineTable.getSysTard();

		List<Double> desRespT = desDeadlineTable.getRespT();
		List<Double> desTard = desDeadlineTable.getTard();
		List<Double> desSysTard = desDeadlineTable.getSysTard();

		// Validate response time (sanity check - should match since both are simulations)
		for (int i = 0; i < jmtRespT.size(); i++) {
			double jmtR = jmtRespT.get(i);
			double desR = desRespT.get(i);
			if (jmtR > 1e-6) {
				double relErr = Math.abs(desR - jmtR) / Math.abs(jmtR);
				assertTrue(relErr <= REL_ERROR_TOL,
						testName + ": Entry " + i + " RespT relative error " + String.format("%.2e", relErr)
								+ " exceeds tolerance " + String.format("%.2e", REL_ERROR_TOL)
								+ " (JMT=" + jmtR + ", DES=" + desR + ")");
			}
		}

		// Validate station-level tardiness (Tard)
		for (int i = 0; i < jmtTard.size(); i++) {
			double jmtT = jmtTard.get(i);
			double desT = desTard.get(i);
			if (jmtT > 1e-6) {
				double relErr = Math.abs(desT - jmtT) / Math.abs(jmtT);
				assertTrue(relErr <= REL_ERROR_TOL,
						testName + ": Entry " + i + " Tard relative error " + String.format("%.2e", relErr)
								+ " exceeds tolerance " + String.format("%.2e", REL_ERROR_TOL)
								+ " (JMT=" + jmtT + ", DES=" + desT + ")");
			}
		}

		// Validate system-level tardiness (SysTard)
		for (int i = 0; i < jmtSysTard.size(); i++) {
			double jmtST = jmtSysTard.get(i);
			double desST = desSysTard.get(i);
			if (jmtST > 1e-6) {
				double relErr = Math.abs(desST - jmtST) / Math.abs(jmtST);
				assertTrue(relErr <= REL_ERROR_TOL,
						testName + ": Entry " + i + " SysTard relative error " + String.format("%.2e", relErr)
								+ " exceeds tolerance " + String.format("%.2e", REL_ERROR_TOL)
								+ " (JMT=" + jmtST + ", DES=" + desST + ")");
			}
		}
	}
}
