package jline.lib.smc;

import static jline.io.InputOutput.line_warning;

import java.util.Collections;
import java.util.HashMap;
import java.util.List;
import java.util.Map;

/**
 * ParseOptPara - Utility class to parse optional parameters in QBD and related solvers.
 */
public final class ParseOptPara {
    private ParseOptPara() {}

    public static Map<String, Object> parse(Map<String, Object> defaults, Object[] optParams) {
        return parse(defaults, optParams, Collections.<String, String>emptyMap(),
                Collections.<String, List<String>>emptyMap());
    }

    public static Map<String, Object> parse(
            Map<String, Object> defaults,
            Object[] optParams,
            Map<String, String> validTypes,
            Map<String, List<String>> validValues) {

        Map<String, Object> result = new HashMap<String, Object>(defaults);

        int i = 0;
        while (i < optParams.length - 1) {
            Object pn = optParams[i];
            Object paramValue = optParams[i + 1];

            String paramName = (pn instanceof String) ? (String) pn : null;
            if (paramName == null) {
                line_warning("ParseOptPara", "Parameter name at position %d is not a string, ignoring", i);
                i += 2;
                continue;
            }

            if (!defaults.containsKey(paramName)) {
                line_warning("ParseOptPara", "Property name '%s' not recognized and ignored", paramName);
                i += 2;
                continue;
            }

            String expectedType = validTypes.get(paramName);
            if (expectedType != null) {
                boolean isValidType;
                if ("numeric".equals(expectedType)) {
                    isValidType = paramValue instanceof Number;
                } else if ("char".equals(expectedType)) {
                    isValidType = paramValue instanceof String;
                } else if ("logical".equals(expectedType)) {
                    isValidType = paramValue instanceof Boolean;
                } else {
                    isValidType = true;
                }
                if (!isValidType) {
                    line_warning("ParseOptPara", "Property value '%s' of '%s' has an incorrect type and is ignored",
                            paramValue, paramName);
                    i += 2;
                    continue;
                }
            }

            List<String> allowedValues = validValues.get(paramName);
            if (allowedValues != null && paramValue instanceof String) {
                if (!allowedValues.contains(paramValue)) {
                    line_warning("ParseOptPara", "Property value '%s' of '%s' not allowed and ignored",
                            paramValue, paramName);
                    i += 2;
                    continue;
                }
            }

            result.put(paramName, paramValue);
            i += 2;
        }

        if (optParams.length % 2 != 0) {
            Object lastParam = optParams[optParams.length - 1];
            line_warning("ParseOptPara",
                    "An odd number of optional parameters detected, last parameter '%s' ignored", lastParam);
        }

        return result;
    }

    public static Map<String, Object> parseSimple(Map<String, Object> defaults, Object... optParams) {
        return parse(defaults, optParams);
    }
}
