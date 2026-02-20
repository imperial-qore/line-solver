package jline.lib.smc

import jline.io.line_warning

/**
 * ParseOptPara - Utility function to parse optional parameters
 * 
 * This is a Kotlin port of the MATLAB ParseOptPara function for parsing
 * optional parameters in QBD and related solvers.
 */
class ParseOptPara {
    companion object {
        /**
         * Parse optional parameters from a map of parameter names to values
         * 
         * @param defaults Default parameter values as a map
         * @param optParams Optional parameters provided by user
         * @param validTypes Map of parameter names to their expected types
         * @param validValues Map of parameter names to lists of valid values (for string params)
         * @return Updated parameters map with parsed values
         */
        fun parse(
            defaults: MutableMap<String, Any?>,
            optParams: Array<out Any?>,
            validTypes: Map<String, String> = emptyMap(),
            validValues: Map<String, List<String>> = emptyMap()
        ): Map<String, Any?> {
            
            val result = defaults.toMutableMap()
            
            // Process pairs of parameter name and value
            var i = 0
            while (i < optParams.size - 1) {
                val paramName = optParams[i] as? String
                val paramValue = optParams[i + 1]
                
                if (paramName == null) {
                    line_warning("ParseOptPara", "Parameter name at position %d is not a string, ignoring", i)
                    i += 2
                    continue
                }

                if (!defaults.containsKey(paramName)) {
                    line_warning("ParseOptPara", "Property name '%s' not recognized and ignored", paramName)
                    i += 2
                    continue
                }
                
                // Check type if validation is specified
                val expectedType = validTypes[paramName]
                if (expectedType != null) {
                    val isValidType = when (expectedType) {
                        "numeric" -> paramValue is Number
                        "char" -> paramValue is String
                        "logical" -> paramValue is Boolean
                        else -> true
                    }
                    
                    if (!isValidType) {
                        line_warning("ParseOptPara", "Property value '%s' of '%s' has an incorrect type and is ignored", paramValue, paramName)
                        i += 2
                        continue
                    }
                }

                // Check valid values if specified
                val allowedValues = validValues[paramName]
                if (allowedValues != null && paramValue is String) {
                    if (!allowedValues.contains(paramValue)) {
                        line_warning("ParseOptPara", "Property value '%s' of '%s' not allowed and ignored", paramValue, paramName)
                        i += 2
                        continue
                    }
                }
                
                // Set the parameter value
                result[paramName] = paramValue
                i += 2
            }
            
            // Check for odd number of parameters
            if (optParams.size % 2 != 0) {
                val lastParam = optParams.last()
                line_warning("ParseOptPara", "An odd number of optional parameters detected, last parameter '%s' ignored", lastParam)
            }
            
            return result
        }
        
        /**
         * Simplified version for numeric and string parameters only
         */
        fun parseSimple(
            defaults: MutableMap<String, Any?>,
            vararg optParams: Any?
        ): Map<String, Any?> {
            return parse(defaults, optParams)
        }
    }
}