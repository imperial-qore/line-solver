package jline.opt.results;

import java.util.LinkedHashMap;
import java.util.Map;

/**
 * Analytic performance sensitivities for a product-form model, mirroring the
 * nested dict returned by native-Python {@code compute_model_sensitivities}:
 * metric kind ('RespT'|'QLen'|'Tput'|'Util') -&gt; metric key -&gt; parameter
 * key -&gt; d(metric)/d(parameter).
 *
 * <p>Keys are canonical strings so the gradient path can match a decision
 * variable's parameter key against a sensitivity entry: a metric key is
 * {@code station} (for 'Util') or {@code station||jobclass}; a parameter key is
 * {@code rate||station||jobclass}. Use {@link #metricKey} / {@link #paramKey}
 * to build them consistently.</p>
 */
public class SensitivityData {

    public static final String SEP = "||";

    // kind -> metricKey -> paramKey -> value
    public final Map<String, Map<String, Map<String, Double>>> data =
            new LinkedHashMap<String, Map<String, Map<String, Double>>>();

    public static String metricKey(String station, String jobclass) {
        return station + SEP + jobclass;
    }

    public static String paramKey(String station, String jobclass) {
        return "rate" + SEP + station + SEP + jobclass;
    }

    public void add(String kind, String metricKey, String paramKey, double value) {
        Map<String, Map<String, Double>> byMetric = data.get(kind);
        if (byMetric == null) {
            byMetric = new LinkedHashMap<String, Map<String, Double>>();
            data.put(kind, byMetric);
        }
        Map<String, Double> byParam = byMetric.get(metricKey);
        if (byParam == null) {
            byParam = new LinkedHashMap<String, Double>();
            byMetric.put(metricKey, byParam);
        }
        byParam.put(paramKey, byParam.containsKey(paramKey)
                ? byParam.get(paramKey) + value : value);
    }

    public Map<String, Map<String, Double>> forKind(String kind) {
        return data.get(kind);
    }

    public boolean isEmpty() {
        return data.isEmpty();
    }
}
