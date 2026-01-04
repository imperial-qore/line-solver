/**
 * @file NetworkStruct Modification Options
 *
 * Defines options for direct modification of NetworkStruct parameters via SN API methods.
 * These options allow users to control modification behavior for optimization scenarios.
 *
 * @since LINE 3.0
 */
package jline.api.sn

/**
 * Modification mode for SN setter methods.
 *
 * Controls whether modifications are made in-place on the input NetworkStruct
 * or on a copy that is returned.
 */
enum class ModifyMode {
    /**
     * Modify the input NetworkStruct in place and return it.
     * This is the fastest option, suitable for optimization loops.
     */
    IN_PLACE,

    /**
     * Create a deep copy of the NetworkStruct, modify the copy, and return it.
     * The original NetworkStruct remains unchanged.
     * Useful when you need to preserve the original state.
     */
    COPY
}

/**
 * Validation level for SN setter methods.
 *
 * Controls the amount of validation performed when modifying NetworkStruct fields.
 * Lower validation levels provide better performance at the cost of safety.
 */
enum class ValidationLevel {
    /**
     * Full validation of all constraints and consistency checks.
     * - Index bounds checking
     * - NaN/Inf detection
     * - Routing matrix stochasticity (row sums = 1.0)
     * - Population non-negative for closed classes
     * - nServers > 0
     * - Node type verification for fork/join operations
     */
    FULL,

    /**
     * Minimal validation - only basic bounds checking.
     * - Index bounds checking
     * - NaN detection
     */
    MINIMAL,

    /**
     * Skip all validation for maximum performance.
     * Use only when you are certain the input is valid.
     * Suitable for tight optimization loops with pre-validated data.
     */
    NONE
}

/**
 * Exception thrown when validation fails during NetworkStruct modification.
 *
 * @param message Description of the validation error
 * @param errors List of specific validation errors found
 */
class SnValidationException(
    message: String,
    val errors: List<String> = emptyList()
) : RuntimeException(message)
