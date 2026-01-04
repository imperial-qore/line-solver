/*
 * Copyright (c) 2012-2026, QORE Lab, Imperial College London
 * All rights reserved.
 */

package jline.lang;

import java.io.Serializable;
import java.util.ArrayList;
import java.util.List;

/**
 * Class representing a global synchronization event with active and passive participants
 */
public class GlobalSync implements Serializable {
    
    private List<ModeEvent> active;
    private List<ModeEvent> passive;
    
    public GlobalSync() {
        this.active = new ArrayList<ModeEvent>();
        this.passive = new ArrayList<ModeEvent>();
    }
    
    public GlobalSync(List<ModeEvent> active, List<ModeEvent> passive) {
        this.active = active != null ? active : new ArrayList<ModeEvent>();
        this.passive = passive != null ? passive : new ArrayList<ModeEvent>();
    }
    
    // Getters and setters
    public List<ModeEvent> getActive() {
        return active;
    }
    
    public void setActive(List<ModeEvent> active) {
        this.active = active;
    }
    
    public List<ModeEvent> getPassive() {
        return passive;
    }
    
    public void setPassive(List<ModeEvent> passive) {
        this.passive = passive;
    }
    
    // Helper methods
    public void addActive(ModeEvent event) {
        if (active == null) {
            active = new ArrayList<ModeEvent>();
        }
        active.add(event);
    }
    
    public void addPassive(ModeEvent event) {
        if (passive == null) {
            passive = new ArrayList<ModeEvent>();
        }
        passive.add(event);
    }
}