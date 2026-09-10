/*
 * Copyright (c) 2012-2026, QORE Lab, Imperial College London
 * All rights reserved.
 */

package jline.util.matrix;

/**
 * A (row, column, value) triple returned by {@link BaseMatrix#nonZeroIterator()}.
 *
 * <p>Iterators reuse a single instance for all entries to avoid allocation;
 * callers must copy the fields if they need to retain an entry beyond the
 * next call to {@code next()}.</p>
 */
public final class MatrixEntry {

    /** Row index of the entry. */
    public int row;

    /** Column index of the entry. */
    public int col;

    /** Value of the entry. */
    public double value;
}
