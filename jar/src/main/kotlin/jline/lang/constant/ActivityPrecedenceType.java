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

    public static String toFeature(ActivityPrecedenceType precedenceType) {
        if (PRE_SEQ.equals(precedenceType)) {
            return "ActivityPrecedenceType_PRE_SEQ";
        }
        if (PRE_AND.equals(precedenceType)) {
            return "ActivityPrecedenceType_PRE_AND";
        }
        if (PRE_OR.equals(precedenceType)) {
            return "ActivityPrecedenceType_PRE_OR";
        }
        if (POST_SEQ.equals(precedenceType)) {
            return "ActivityPrecedenceType_POST_SEQ";
        }
        if (POST_AND.equals(precedenceType)) {
            return "ActivityPrecedenceType_POST_AND";
        }
        if (POST_OR.equals(precedenceType)) {
            return "ActivityPrecedenceType_POST_OR";
        }
        if (POST_LOOP.equals(precedenceType)) {
            return "ActivityPrecedenceType_POST_LOOP";
        }
        if (POST_CACHE.equals(precedenceType)) {
            return "ActivityPrecedenceType_POST_CACHE";
        }
        throw new RuntimeException("Unrecognized precedence type");
    }

    public static int toId(ActivityPrecedenceType precedenceType) {
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
