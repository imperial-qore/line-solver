package jline.opt.objectives;

/** Shared helpers for objectives/constraints reading decoded variable values. */
final class ObjectiveUtil {

    private ObjectiveUtil() {
    }

    /**
     * Numeric contribution of a decoded variable value: a scalar returns its
     * value; a vector (routing probabilities / priority levels) returns its
     * sum; anything else contributes 0. Mirrors native-Python's isinstance
     * checks for scalar vs np.ndarray.
     */
    static double numericValue(Object value) {
        if (value == null) {
            return 0.0;
        }
        if (value instanceof Number) {
            return ((Number) value).doubleValue();
        }
        if (value instanceof double[]) {
            double s = 0.0;
            for (double v : (double[]) value) {
                s += v;
            }
            return s;
        }
        if (value instanceof int[]) {
            double s = 0.0;
            for (int v : (int[]) value) {
                s += v;
            }
            return s;
        }
        return 0.0;
    }

    /** True only for a scalar Number (matches Python isinstance(value,(int,float))). */
    static boolean isScalar(Object value) {
        return value instanceof Number;
    }
}
