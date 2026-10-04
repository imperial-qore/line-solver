/*
 * Copyright (c) 2012-2026, QORE Lab, Imperial College London
 * All rights reserved.
 */

package jline.lang.constant;

/**
 * Constants for defining activity precedences in LayeredNetwork models
 */
public class ActivityPrecedenceType {
    public static final int ID_POST_AND = 12;
    public static final int ID_POST_CACHE = 15;
    public static final int ID_POST_LOOP = 14;
    public static final int ID_POST_OR = 13;
    public static final int ID_POST_SEQ = 11;
    public static final int ID_PRE_AND = 2;
    public static final int ID_PRE_OR = 3;
    public static final int ID_PRE_SEQ = 1;
    public static final String POST_AND = "post-AND";
    public static final String POST_CACHE = "post-CACHE";
    public static final String POST_LOOP = "post-LOOP";
    public static final String POST_OR = "post-OR";
    public static final String POST_SEQ = "post";
    public static final String PRE_AND = "pre-AND";
    public static final String PRE_OR = "pre-OR";
    public static final String PRE_SEQ = "pre";

    /**
     * The FeatureSet entry naming a precedence, as ActivityPrecedence stores it.
     *
     * The argument is the String an ActivityPrecedence carries in preType /
     * postType, not an instance of this constant holder: this class has no
     * instances, so the previous signature could only ever be passed null, and
     * every branch compared a String against it and fell through to the throw.
     * The names returned are the REGISTRY names (ActivityPrecedence_*), not the
     * class name with "Type" in it, which matched no registry entry either.
     *
     * POST_LOOP has no registry entry -- LINE models a loop by its pseudo-task,
     * not by a declared capability -- so it returns the empty string, the same
     * "no gated capability" convention the SchedStrategy and RoutingStrategy
     * helpers use for their internal markers.
     *
     * @param precedenceType the precedence type string
     * @return the registry name, or "" when the type gates nothing
     */
    public static String toFeature(String precedenceType) {
        if (PRE_SEQ.equals(precedenceType)) {
            return "ActivityPrecedence_PRE_SEQ";
        }
        if (PRE_AND.equals(precedenceType)) {
            return "ActivityPrecedence_PRE_AND";
        }
        if (PRE_OR.equals(precedenceType)) {
            return "ActivityPrecedence_PRE_OR";
        }
        if (POST_SEQ.equals(precedenceType)) {
            return "ActivityPrecedence_POST_SEQ";
        }
        if (POST_AND.equals(precedenceType)) {
            return "ActivityPrecedence_POST_AND";
        }
        if (POST_OR.equals(precedenceType)) {
            return "ActivityPrecedence_POST_OR";
        }
        if (POST_CACHE.equals(precedenceType)) {
            return "ActivityPrecedence_POST_CACHE";
        }
        if (POST_LOOP.equals(precedenceType)) {
            return "";
        }
        throw new RuntimeException("Unrecognized precedence type: " + precedenceType);
    }

    /**
     * The numeric id of a precedence, as ActivityPrecedence stores its type.
     *
     * Takes the type String for the same reason toFeature does: this class has
     * no instances, so the previous signature could only be passed null and
     * every branch fell through to the throw.
     *
     * @param precedenceType the precedence type string
     * @return the ID_ constant for that type
     */
    public static int toId(String precedenceType) {
        if (PRE_SEQ.equals(precedenceType)) {
            return ID_PRE_SEQ;
        }
        if (PRE_AND.equals(precedenceType)) {
            return ID_PRE_AND;
        }
        if (PRE_OR.equals(precedenceType)) {
            return ID_PRE_OR;
        }
        if (POST_SEQ.equals(precedenceType)) {
            return ID_POST_SEQ;
        }
        if (POST_AND.equals(precedenceType)) {
            return ID_POST_AND;
        }
        if (POST_OR.equals(precedenceType)) {
            return ID_POST_OR;
        }
        if (POST_LOOP.equals(precedenceType)) {
            return ID_POST_LOOP;
        }
        if (POST_CACHE.equals(precedenceType)) {
            return ID_POST_CACHE;
        }
        throw new RuntimeException("Unrecognized precedence type");
    }
}
