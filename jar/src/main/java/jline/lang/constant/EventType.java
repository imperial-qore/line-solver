/*
 * Copyright (c) 2012-2026, QORE Lab, Imperial College London
 * All rights reserved.
 */

package jline.lang.constant;

/**
 * Constants for specifying events
 */
public enum EventType {
    INIT,
    LOCAL,
    ARV,
    DEP,
    PHASE,
    READ,
    STAGE,
    ENABLE,
    FIRE,
    PRE,
    POST,
    RENEGE,
    RETRY,
    /**
     * The server of a polling station advances its switchover timer.
     * <p>Unlike PHASE, which carries only the internal transitions of a
     * phase-type and leaves the absorption to DEP, this event carries both:
     * the completion of a switchover moves no job and so has no departure to
     * attach the absorption to. It is therefore emitted also for a
     * single-phase switchover, where it consists of the absorption alone.</p>
     */
    SWITCH,
    /**
     * The server of a station breaks down (goes from up to down).
     * <p>Only the server-status local variable changes: jobs in service are
     * not lost and, service being memoryless in the supported case, they
     * resume when the server is repaired. The passive half is LOCAL, so no
     * job moves.</p>
     */
    FAILURE,
    /**
     * The server of a station is repaired (goes from down to up).
     * <p>No START is emitted here on purpose: the supported breakdown model
     * resumes the job the server was holding rather than restarting it
     * (downServiceRates degrades the server, it does not evict).</p>
     */
    REPAIR,
    /**
     * A job begins, or resumes, holding a server at a station.
     * <p>Instantaneous tag on the arc of the ARV or DEP transition that causes
     * it, never the active half of an sn.sync entry: it carries no clock, adds
     * no state and leaves every numerical result unchanged.</p>
     */
    START,
    /**
     * A job holding a server is pushed back into the buffer.
     * <p>Spelled in full because PRE already names the Petri-net pre-arc.
     * Like START, an instantaneous tag on an existing arc.</p>
     */
    PREEMPT
}
