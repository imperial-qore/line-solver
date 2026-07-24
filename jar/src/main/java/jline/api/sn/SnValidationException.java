/**
 * @file NetworkStruct Validation Exception
 *
 * @since LINE 3.0
 */
package jline.api.sn;

import java.util.Collections;
import java.util.List;

/**
 * Exception thrown when validation fails during NetworkStruct modification.
 */
public class SnValidationException extends RuntimeException {
    private static final long serialVersionUID = 1L;
    private final List<String> errors;

    public SnValidationException(String message) {
        super(message);
        this.errors = Collections.emptyList();
    }

    public SnValidationException(String message, List<String> errors) {
        super(message);
        this.errors = errors == null ? Collections.<String>emptyList() : errors;
    }

    public List<String> getErrors() {
        return errors;
    }
}
