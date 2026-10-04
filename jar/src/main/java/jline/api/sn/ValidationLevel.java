/**
 * @file NetworkStruct Validation Level
 *
 * @since LINE 3.0
 */
package jline.api.sn;

/**
 * Validation level for SN setter methods.
 *
 * Controls the amount of validation performed when modifying NetworkStruct fields.
 * Lower validation levels provide better performance at the cost of safety.
 */
public enum ValidationLevel {
    /**
     * Full validation of all constraints and consistency checks.
     */
    FULL,

    /**
     * Minimal validation - only basic bounds checking.
     */
    MINIMAL,

    /**
     * Skip all validation for maximum performance.
     */
    NONE
}
