package jline.util;

import java.io.Serializable;
import java.util.function.Function;

/**
 * Interface used for routing functions
 */
public interface SerializableFunction<T, U> extends Function<T, U>, Serializable { }
