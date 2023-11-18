// Copyright (c) 2012-2022, Imperial College London
// All rights reserved.

package jline.solvers;

import jline.lang.Model;
import jline.lang.Network;
import jline.lang.constant.SolverType;

import javax.xml.parsers.ParserConfigurationException;
import java.util.Random;

// Abstract class for model solution algorithms and tools
public abstract class Solver {

  public String name; // Solver name
  public SolverOptions options; // Data structure with solver options
  public SolverResult result; // Last result
  public boolean enableChecks;
  public Random random;

  protected Solver(String name, SolverOptions options) {
    this.name = name;
    this.setOptions(options);
    this.result = new SolverResult();
    this.enableChecks = true;
    this.random = new Random(options.seed);
  }

  protected Solver(String name) {
    this(name, defaultOptions());
  }

  // Generic method to run the solver
  protected abstract void runAnalyzer() throws IllegalAccessException, ParserConfigurationException;

  protected void setChecks(boolean bool) {
    enableChecks = bool;
  }


  // Check if the model has been solved
  protected boolean hasResults() {

    if (result.QN == null) {
      return false;
    } else {
      return !result.QN.isEmpty();
    }
  }

  // Dispose previously stored results
  public void resetResults() {
    this.result = new SolverResult();
  }

  // Set a new options data structure
  protected void setOptions(SolverOptions options) {
    this.options = options;
  }

  // Assign a new seed to the random number generator
  public void resetRandomGeneratorSeed(int seed) {
    this.random = new Random(seed);
  }

  // List valid fields for options data structure
  protected static void listValidOptions() {
    // TODO: implementation - note return type should likely not be void
    throw new RuntimeException("listValidOptions() has not yet been implemented in JLINE.");
  }

  public String getName() { return name; }

  // Return default options
  public static SolverOptions defaultOptions() {
    return new SolverOptions(null);
    // The line below (iter_max = 100) appears in LINE but is very out of place here
    // Removed for now but leaving commented in case it serves an important purpose
    // options.iter_max = 100;
  }

  // Parse option parameters into options data structure
  protected static SolverOptions parseOptions() {
    // TODO: implementation - note arguments should likely not be void
    throw new RuntimeException("parseOptions() has not yet been implemented in JLINE.");
  }

  // Returns a solver configured to run the chosen method
  protected static Solver load(String chosenMethod, Model model) {
    // TODO: implementation - note further arguments may be needed
    throw new RuntimeException("load() has not yet been implemented in JLINE.");
  }

  // NOTE: the following LINE methods have not been migrated to JLINE
  // a) getModel() - since model has moved to NetworkSolver/EnsembleSolver level
  // b) getStruct() - part of Network, not Solver
  // c) getName() - name is public and therefore no need for getter
  // d) getResults() - results is public and therefore no need for getter
  // e) getOptions() - options is public and therefore no need for getter
  // f) reset() - entirely duplicative with resetResults()
  // g) isAvailable() - always returns true, provides no functionality
  // h) isJavaAvailable() - unnecessary in JLINE
  // i) 4 x ODESolver() - accessible via options
  // j) isValidOption() - all options are always available as part of SolverOptions class
  // k) getDefaultOptions() - duplicative with static method defaultOptions()
  // l) supports() - static method at specific Solver level rather than abstract within Solver class


}
