/*
 * Copyright (c) 2012-2026, QORE Lab, Imperial College London
 * All rights reserved.
 */

package jline.util;

import java.io.Serializable;
import java.util.function.Function;

/**
 * Interface used for routing functions
 */
public interface SerializableFunction<T, U> extends Function<T, U>, Serializable {
}
