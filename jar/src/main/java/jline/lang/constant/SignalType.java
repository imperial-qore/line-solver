/*
 * Copyright (c) 2012-2026, QORE Lab, Imperial College London
 * All rights reserved.
 */

package jline.lang.constant;

import java.io.Serializable;

/**
 * Enumeration of signal types for signal classes in queueing networks.
 */
public enum SignalType implements Serializable {
    NEGATIVE(1),
    REPLY(0),
    CATASTROPHE(2);

    private final int id;

    SignalType(int id) {
        this.id = id;
    }

    public int getID() {
        return id;
    }

    public static SignalType fromID(int id) {
        for (SignalType type : values()) {
            if (type.id == id) {
                return type;
            }
        }
        return null;
    }

    public static String toText(SignalType type) {
        if (type == null) {
            return null;
        }
        switch (type) {
            case NEGATIVE:
                return "negative";
            case REPLY:
                return "reply";
            case CATASTROPHE:
                return "catastrophe";
            default:
                return type.name();
        }
    }

    public static SignalType fromText(String text) {
        if (text == null || text.isEmpty()) {
            return null;
        }
        switch (text.toLowerCase()) {
            case "negative":
                return NEGATIVE;
            case "reply":
                return REPLY;
            case "catastrophe":
                return CATASTROPHE;
            default:
                return null;
        }
    }
}
