package jline.opt.results;

import java.util.LinkedHashMap;
import java.util.Map;

/**
 * Result from solving one subproblem in a decomposition workflow. Mirrors
 * native-Python {@code line_solver.opt.results.SubProblemResult}.
 */
public class SubProblemResult {

    public String name = "";
    public OptimizationResult result = new OptimizationResult();
    public final Map<String, Object> variablesFixed = new LinkedHashMap<String, Object>();

    public SubProblemResult() {
    }

    public SubProblemResult(String name, OptimizationResult result) {
        this.name = name;
        this.result = result;
    }
}
