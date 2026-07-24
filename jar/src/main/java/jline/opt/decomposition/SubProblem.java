package jline.opt.decomposition;

import jline.opt.variables.DecisionVariable;

import java.util.ArrayList;
import java.util.LinkedHashMap;
import java.util.List;
import java.util.Map;

/**
 * A subproblem in a decomposition: a subset of variables to optimize while
 * others are held fixed. Mirrors native-Python
 * {@code line_solver.opt.decomposition.SubProblem}.
 */
public class SubProblem {

    public String name;
    public String variableType;
    public List<DecisionVariable> variables;
    public Map<String, Object> fixedValues = new LinkedHashMap<String, Object>();

    public SubProblem(String name, String variableType, List<DecisionVariable> variables) {
        this.name = name;
        this.variableType = variableType;
        this.variables = variables;
    }

    public List<String> getVariableNames() {
        List<String> names = new ArrayList<String>();
        for (DecisionVariable v : variables) {
            names.add(v.getName());
        }
        return names;
    }
}
