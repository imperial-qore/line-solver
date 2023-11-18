package jline.lang.constant;

import java.io.Serializable;

/**
 *  Constants for specifying timing strategies at Petri net transitions
 */
public enum TimingStrategy implements Serializable {
    TIMED,
    IMMEDIATE
}