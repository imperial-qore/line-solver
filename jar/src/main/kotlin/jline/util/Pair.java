/*
 * Copyright (c) 2012-2026, QORE Lab, Imperial College London
 * All rights reserved.
 */

package jline.util;

import java.io.Serializable;

/**
 * A generic pair container that holds two objects of potentially different types.
 * 
 * <p>This utility class provides a simple way to group two related values together.
 * It implements {@link Comparable} for ordering pairs and {@link Serializable} for
 * persistence. The comparison is based on the left element if it's comparable,
 * otherwise on the right element.</p>
 * 
 * <p>Common use cases:
 * <ul>
 *   <li>Returning multiple values from methods</li>
 *   <li>Storing coordinate pairs or ranges</li>
 *   <li>Key-value associations where Map is not suitable</li>
 *   <li>Temporary grouping of related objects</li>
 * </ul>
 * </p>
 * 
 * @param <T> the type of the left element
 * @param <U> the type of the right element
 * @since 1.0
 */
public class Pair<T, U> implements Comparable<Pair<T, U>>, Serializable {
    private T left;
    private U right;

    /**
     * Constructs a new pair with the specified left and right elements.
     * 
     * @param left  the left element of the pair
     * @param right the right element of the pair
     */
    public Pair(T left, U right) {
        this.left = left;
        this.right = right;
    }

    /**
     * Copy constructor that creates a new pair with the same elements as another pair.
     * 
     * @param other the pair to copy from
     */
    public Pair(Pair<T, U> other) {
        this.left = other.left;
        this.right = other.right;
    }

    @SuppressWarnings("unchecked")
    @Override
    public int compareTo(Pair<T, U> other) {
        if (this.left instanceof Comparable) {
            return ((Comparable) this.left).compareTo(other.getLeft());
        } else if (this.right instanceof Comparable) {
            return ((Comparable) this.right).compareTo(other.getRight());
        }

        return 0;
    }

    public T getLeft() {
        return this.left;
    }

    public void setLeft(T left) {
        this.left = left;
    }

    public U getRight() {
        return this.right;
    }

    public void setRight(U right) {
        this.right = right;
    }
}
