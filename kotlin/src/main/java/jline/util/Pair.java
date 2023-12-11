package jline.util;

import java.io.Serializable;

/**
 * A pair of objects of different classes.
 */
public class Pair<T, U> implements Comparable<Pair<T,U>>, Serializable {
    private T left;
    private U right;

    public Pair(T left, U right) {
        this.left = left;
        this.right = right;
    }

    // Copy constructor
    public Pair(Pair<T, U> other) {
        this.left = other.left;
        this.right = other.right;
    }

    public T getLeft() {
        return this.left;
    }

    public U getRight() {
        return this.right;
    }

    public void setLeft(T left) {
        this.left = left;
    }

    public void setRight(U right) {
        this.right = right;
    }
    @SuppressWarnings("unchecked")
    @Override
    public int compareTo(Pair<T,U> other) {
        if (this.left instanceof Comparable) {
            return ((Comparable)this.left).compareTo(other.getLeft());
        } else if (this.right instanceof Comparable) {
            return ((Comparable)this.right).compareTo(other.getRight());
        }

        return 0;
    }
}
