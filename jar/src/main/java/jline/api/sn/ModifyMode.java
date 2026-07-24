/**
 * @file NetworkStruct Modification Options
 *
 * Defines options for direct modification of NetworkStruct parameters via SN API methods.
 *
 * @since LINE 3.0
 */
package jline.api.sn;

/**
 * Modification mode for SN setter methods.
 *
 * Controls whether modifications are made in-place on the input NetworkStruct
 * or on a copy that is returned.
 */
public enum ModifyMode {
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
