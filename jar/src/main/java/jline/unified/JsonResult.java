package jline.unified;

import java.util.ArrayList;
import java.util.HashMap;
import java.util.List;
import java.util.Map;

/**
 * Data class for JSON-formatted test results.
 */
public class JsonResult {
    public final String modelName;
    public final String language;
    public String status;
    public final Map<String, Object> timing;
    public final Map<String, Map<String, Object>> solverResults;
    public final List<String> errors;

    public JsonResult(String modelName) {
        this.modelName = modelName;
        this.language = "java";
        this.status = "passed";
        this.timing = new HashMap<String, Object>();
        this.timing.put("total_ms", 0L);
        this.timing.put("solvers", new HashMap<String, Long>());
        this.solverResults = new HashMap<String, Map<String, Object>>();
        this.errors = new ArrayList<String>();
    }
}
