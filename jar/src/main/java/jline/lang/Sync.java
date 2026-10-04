/*
 * Copyright (c) 2012-2026, QORE Lab, Imperial College London
 * All rights reserved.
 */

package jline.lang;

import java.io.Serializable;
import java.util.HashMap;
import java.util.Map;

/**
 * A declaration of a synchronization on a NetworkEvent
 */
public class Sync implements Serializable {

    //Index starts from 0
    public Map<Integer, Event> active;
    public Map<Integer, Event> passive;

    public Sync() {
        active = new HashMap<Integer, Event>();
        passive = new HashMap<Integer, Event>();
    }
}
