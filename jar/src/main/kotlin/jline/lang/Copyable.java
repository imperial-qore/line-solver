package jline.lang;

import java.io.ByteArrayInputStream;
import java.io.ByteArrayOutputStream;
import java.io.IOException;
import java.io.ObjectInputStream;
import java.io.ObjectOutputStream;
import java.io.Serializable;

/**
 * Copyable interface allows to perform deep-copy of objects via the copy() method.
 * Classes implementing this interface must also implement Serializable.
 * 
 * Copyright (c) 2012-2026, Imperial College London
 * All rights reserved.
 */
public interface Copyable extends Serializable {
    
    /**
     * Creates a deep copy of this object using serialization.
     * 
     * @return A deep copy of this object
     * @throws RuntimeException if the copy operation fails
     */
    @SuppressWarnings("unchecked")
    default <T extends Copyable> T copy() {
        try {
            ByteArrayOutputStream bos = new ByteArrayOutputStream();
            ObjectOutputStream out = new ObjectOutputStream(bos);
            out.writeObject(this);
            out.close();
            
            ByteArrayInputStream bis = new ByteArrayInputStream(bos.toByteArray());
            ObjectInputStream in = new ObjectInputStream(bis);
            T copy = (T) in.readObject();
            in.close();
            
            return copy;
        } catch (IOException | ClassNotFoundException e) {
            throw new RuntimeException("Failed to create deep copy of " + this.getClass().getSimpleName(), e);
        }
    }
}