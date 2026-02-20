/*
 * Copyright (c) 2012-2026, QORE Lab, Imperial College London
 * All rights reserved.
 */

package jline.lang.layered;

import jline.lang.constant.ActivityPrecedenceType;
import jline.util.matrix.Matrix;

import java.util.ArrayList;
import java.util.List;

/**
 * A class modeling precedence relationships among activities
 */
public class ActivityPrecedence {


    protected List<String> preActs;
    protected List<String> postActs;
    protected String preType;
    protected String postType;
    protected Matrix preParams;
    protected Matrix postParams;

    /**
     * Helper method to convert Activity objects to their names.
     *
     * @param activities the list of activities to convert
     * @return the list of activity names
     */
    private static List<String> convertActivitiesToNames(List<Activity> activities) {
        List<String> names = new ArrayList<>();
        for (Activity activity : activities) {
            names.add(activity.getName());
        }
        return names;
    }

    /**
     * Constructs an ActivityPrecedence with the specified parameters.
     *
     * @param preActs    the list of preceding activity names
     * @param postActs   the list of following activity names
     * @param preType    the type of the precedence relationship before the activity
     * @param postType   the type of the precedence relationship after the activity
     * @param preParams  the parameters for the preceding activities
     * @param postParams the parameters for the following activities
     */
    public ActivityPrecedence(List<String> preActs, List<String> postActs, String preType, String postType, Matrix preParams, Matrix postParams) {
        this.preActs = preActs;
        this.postActs = postActs;
        this.preType = preType;
        this.postType = postType;
        this.preParams = preParams;
        this.postParams = postParams;
    }

    /**
     * Constructs an ActivityPrecedence with the specified parameters, without postParams.
     *
     * @param preActs   the list of preceding activity names
     * @param postActs  the list of following activity names
     * @param preType   the type of the precedence relationship before the activity
     * @param postType  the type of the precedence relationship after the activity
     * @param preParams the parameters for the preceding activities
     */
    public ActivityPrecedence(List<String> preActs, List<String> postActs, String preType, String postType, Matrix preParams) {
        this.preActs = preActs;
        this.postActs = postActs;
        this.preType = preType;
        this.postType = postType;
        this.preParams = preParams;
        this.postParams = null;
    }

    /**
     * Constructs an ActivityPrecedence with the specified parameters, without preParams and postParams.
     *
     * @param preActs  the list of preceding activity names
     * @param postActs the list of following activity names
     * @param preType  the type of the precedence relationship before the activity
     * @param postType the type of the precedence relationship after the activity
     */
    public ActivityPrecedence(List<String> preActs, List<String> postActs, String preType, String postType) {
        this.preActs = preActs;
        this.postActs = postActs;
        this.preType = preType;
        this.postType = postType;
        this.preParams = null;
        this.postParams = null;
    }

    /**
     * Constructs an ActivityPrecedence with the specified parameters, assuming POST_SEQ as postType.
     *
     * @param preActs  the list of preceding activity names
     * @param postActs the list of following activity names
     * @param preType  the type of the precedence relationship before the activity
     */
    public ActivityPrecedence(List<String> preActs, List<String> postActs, String preType) {
        this.preActs = preActs;
        this.postActs = postActs;
        this.preType = preType;
        this.postType = ActivityPrecedenceType.POST_SEQ;
        this.preParams = null;
        this.postParams = null;
    }

    /**
     * Constructs an ActivityPrecedence with the specified parameters, assuming PRE_SEQ as preType and POST_SEQ as postType.
     *
     * @param preActs  the list of preceding activity names
     * @param postActs the list of following activity names
     */
    public ActivityPrecedence(List<String> preActs, List<String> postActs) {
        this.preActs = preActs;
        this.postActs = postActs;
        this.preType = ActivityPrecedenceType.PRE_SEQ;
        this.postType = ActivityPrecedenceType.POST_SEQ;
        this.preParams = null;
        this.postParams = null;
    }

    // Factory methods for Activity objects
    
    /**
     * Creates an ActivityPrecedence with Activity objects.
     *
     * @param preActs    the list of preceding activities
     * @param postActs   the list of following activities
     * @param preType    the type of the precedence relationship before the activity
     * @param postType   the type of the precedence relationship after the activity
     * @param preParams  the parameters for the preceding activities
     * @param postParams the parameters for the following activities
     * @return an ActivityPrecedence object
     */
    public static ActivityPrecedence fromActivities(List<Activity> preActs, List<Activity> postActs, String preType, String postType, Matrix preParams, Matrix postParams) {
        return new ActivityPrecedence(convertActivitiesToNames(preActs), convertActivitiesToNames(postActs), preType, postType, preParams, postParams);
    }

    /**
     * Creates an ActivityPrecedence with Activity objects, without postParams.
     *
     * @param preActs   the list of preceding activities
     * @param postActs  the list of following activities
     * @param preType   the type of the precedence relationship before the activity
     * @param postType  the type of the precedence relationship after the activity
     * @param preParams the parameters for the preceding activities
     * @return an ActivityPrecedence object
     */
    public static ActivityPrecedence fromActivities(List<Activity> preActs, List<Activity> postActs, String preType, String postType, Matrix preParams) {
        return new ActivityPrecedence(convertActivitiesToNames(preActs), convertActivitiesToNames(postActs), preType, postType, preParams);
    }

    /**
     * Creates an ActivityPrecedence with Activity objects, without preParams and postParams.
     *
     * @param preActs  the list of preceding activities
     * @param postActs the list of following activities
     * @param preType  the type of the precedence relationship before the activity
     * @param postType the type of the precedence relationship after the activity
     * @return an ActivityPrecedence object
     */
    public static ActivityPrecedence fromActivities(List<Activity> preActs, List<Activity> postActs, String preType, String postType) {
        return new ActivityPrecedence(convertActivitiesToNames(preActs), convertActivitiesToNames(postActs), preType, postType);
    }

    /**
     * Creates an ActivityPrecedence with Activity objects, assuming POST_SEQ as postType.
     *
     * @param preActs  the list of preceding activities
     * @param postActs the list of following activities
     * @param preType  the type of the precedence relationship before the activity
     * @return an ActivityPrecedence object
     */
    public static ActivityPrecedence fromActivities(List<Activity> preActs, List<Activity> postActs, String preType) {
        return new ActivityPrecedence(convertActivitiesToNames(preActs), convertActivitiesToNames(postActs), preType);
    }

    /**
     * Creates an ActivityPrecedence with Activity objects, assuming PRE_SEQ as preType and POST_SEQ as postType.
     *
     * @param preActs  the list of preceding activities
     * @param postActs the list of following activities
     * @return an ActivityPrecedence object
     */
    public static ActivityPrecedence fromActivities(List<Activity> preActs, List<Activity> postActs) {
        return new ActivityPrecedence(convertActivitiesToNames(preActs), convertActivitiesToNames(postActs));
    }






    /**
     * Creates an ActivityPrecedence object representing an AND-fork relationship.
     *
     * @param preAct   the preceding activity
     * @param postActs the list of following activities
     * @return an ActivityPrecedence object
     */
    public static ActivityPrecedence AndFork(Activity preAct, List<Activity> postActs) {
        Matrix fanout = new Matrix(1, postActs.size(), postActs.size());
        fanout.fill(1.0);
        return ActivityPrecedence.AndFork(preAct, postActs, fanout);
    }

    /**
     * Creates an ActivityPrecedence object representing an AND-fork relationship.
     *
     * @param preAct   the preceding activity
     * @param postActs the list of following activities
     * @return an ActivityPrecedence object
     */
    public static ActivityPrecedence AndFork(String preAct, List<String> postActs) {
        Matrix fanout = new Matrix(1, postActs.size(), postActs.size());
        fanout.fill(1.0);
        return ActivityPrecedence.AndFork(preAct, postActs, fanout);
    }

    /**
     * Creates an ActivityPrecedence object representing an AND-fork relationship with a specified fanout matrix.
     *
     * @param preAct   the preceding activity
     * @param postActs the list of following activities
     * @param fanout   the fanout matrix
     * @return an ActivityPrecedence object
     */
    public static ActivityPrecedence AndFork(String preAct, List<String> postActs, Matrix fanout) {
        List<String> preActs = new ArrayList<>();
        preActs.add(preAct);
        return new ActivityPrecedence(preActs, postActs, ActivityPrecedenceType.PRE_SEQ, ActivityPrecedenceType.POST_AND, null, fanout);
    }

    /**
     * Creates an ActivityPrecedence object representing an AND-fork relationship with a specified fanout matrix.
     *
     * @param preAct   the preceding activity
     * @param postActs the list of following activities
     * @param fanout   the fanout matrix
     * @return an ActivityPrecedence object
     */
    public static ActivityPrecedence AndFork(Activity preAct, List<Activity> postActs, Matrix fanout) {
        List<Activity> preActs = new ArrayList<>();
        preActs.add(preAct);
        return fromActivities(preActs, postActs, ActivityPrecedenceType.PRE_SEQ, ActivityPrecedenceType.POST_AND, null, fanout);
    }

    /**
     * Creates an ActivityPrecedence object representing an AND-join relationship.
     *
     * @param preActs the list of preceding activities
     * @param postAct the following activity
     * @return an ActivityPrecedence object
     */
    public static ActivityPrecedence AndJoin(List<String> preActs, String postAct) {
        Matrix quorum = new Matrix(0, 0); // Empty matrix to match MATLAB behavior
        List<String> postActs = new ArrayList<String>();
        postActs.add(postAct);
        return new ActivityPrecedence(preActs, postActs, ActivityPrecedenceType.PRE_AND, ActivityPrecedenceType.POST_SEQ, quorum, null);
    }

    /**
     * Creates an ActivityPrecedence object representing an AND-join relationship.
     *
     * @param preActs the list of preceding activities
     * @param postAct the following activity
     * @return an ActivityPrecedence object
     */
    public static ActivityPrecedence AndJoin(List<Activity> preActs, Activity postAct) {
        Matrix quorum = Matrix.ones(1, preActs.size());
        List<Activity> postActs = new ArrayList<Activity>();
        postActs.add(postAct);
        return fromActivities(preActs, postActs, ActivityPrecedenceType.PRE_AND, ActivityPrecedenceType.POST_SEQ, quorum, null);
    }

    /**
     * Creates an ActivityPrecedence object representing an AND-join relationship with a specified quorum matrix.
     *
     * @param preActs the list of preceding activities
     * @param postAct the following activity
     * @param quorum  the quorum matrix
     * @return an ActivityPrecedence object
     */
    public static ActivityPrecedence AndJoin(List<Activity> preActs, Activity postAct, Matrix quorum) {
        List<Activity> postActs = new ArrayList<Activity>();
        postActs.add(postAct);
        return fromActivities(preActs, postActs, ActivityPrecedenceType.PRE_AND, ActivityPrecedenceType.POST_SEQ, quorum, null);
    }

    /**
     * Creates an ActivityPrecedence object representing an AND-join relationship with a specified quorum matrix.
     *
     * @param preActs the list of preceding activities
     * @param postAct the following activity
     * @param quorum  the quorum matrix
     * @return an ActivityPrecedence object
     */
    public static ActivityPrecedence AndJoin(List<String> preActs, String postAct, Matrix quorum) {
        List<String> postActs = new ArrayList<String>();
        postActs.add(postAct);
        return new ActivityPrecedence(preActs, postActs, ActivityPrecedenceType.PRE_AND, ActivityPrecedenceType.POST_SEQ, quorum, null);
    }

    /**
     * Creates an ActivityPrecedence object representing a cache access relationship.
     *
     * @param preAct   the preceding activity
     * @param postActs the list of following activities
     * @return an ActivityPrecedence object
     */
    public static ActivityPrecedence CacheAccess(Activity preAct, List<Activity> postActs) {
        List<Activity> preActs = new ArrayList<>();
        preActs.add(preAct);
        return fromActivities(preActs, postActs, ActivityPrecedenceType.PRE_SEQ, ActivityPrecedenceType.POST_CACHE);
    }

    /**
     * Creates an ActivityPrecedence object representing a cache access relationship.
     *
     * @param preAct   the preceding activity
     * @param postActs the list of following activities
     * @return an ActivityPrecedence object
     */
    public static ActivityPrecedence CacheAccess(String preAct, List<String> postActs) {
        List<String> preActs = new ArrayList<>();
        preActs.add(preAct);
        return new ActivityPrecedence(preActs, postActs, ActivityPrecedenceType.PRE_SEQ, ActivityPrecedenceType.POST_CACHE);
    }

    /**
     * Creates an ActivityPrecedence object representing a loop relationship.
     *
     * @param preAct  the preceding activity
     * @param loopAct the activity that is looped
     * @param endAct  the activity after the loop completes
     * @param nloops  the number of loops
     * @return an ActivityPrecedence object
     */
    public static ActivityPrecedence Loop(Activity preAct, Activity loopAct, Activity endAct, double nloops) {
        List<Activity> preActs = new ArrayList<>();
        preActs.add(preAct);
        List<Activity> postActs = new ArrayList<>();
        postActs.add(loopAct);
        postActs.add(endAct);
        return fromActivities(preActs, postActs, ActivityPrecedenceType.PRE_SEQ, ActivityPrecedenceType.POST_LOOP, null, Matrix.singleton(nloops));
    }

    /**
     * Creates an ActivityPrecedence object representing a loop relationship.
     *
     * @param preAct  the preceding activity
     * @param loopAct the activity that is looped
     * @param endAct  the activity after the loop completes
     * @param nloops  the number of loops
     * @return an ActivityPrecedence object
     */
    public static ActivityPrecedence Loop(String preAct, String loopAct, String endAct, double nloops) {
        List<String> preActs = new ArrayList<>();
        preActs.add(preAct);
        List<String> postActs = new ArrayList<>();
        postActs.add(loopAct);
        postActs.add(endAct);
        return new ActivityPrecedence(preActs, postActs, ActivityPrecedenceType.PRE_SEQ, ActivityPrecedenceType.POST_LOOP, null, Matrix.singleton(nloops));
    }

    /**
     * Creates an ActivityPrecedence object representing a loop relationship.
     *
     * @param preAct   the preceding activity
     * @param postActs a list with loopAct and endAct (see 4-parameter Loop constructor)
     * @param nloops   the number of loops
     * @return an ActivityPrecedence object
     */
    public static ActivityPrecedence Loop(Activity preAct, List<Activity> postActs, double nloops) {
        List<Activity> preActs = new ArrayList<>();
        preActs.add(preAct);
        return fromActivities(preActs, postActs, ActivityPrecedenceType.PRE_SEQ, ActivityPrecedenceType.POST_LOOP, null, Matrix.singleton(nloops));
    }

    /**
     * Creates an ActivityPrecedence object representing a loop relationship.
     *
     * @param preAct   the preceding activity
     * @param postActs a list with loopAct and endAct (see 4-parameter Loop constructor)
     * @param nloops   the number of loops
     * @return an ActivityPrecedence object
     */
    public static ActivityPrecedence Loop(String preAct, List<String> postActs, double nloops) {
        List<String> preActs = new ArrayList<>();
        preActs.add(preAct);
        return new ActivityPrecedence(preActs, postActs, ActivityPrecedenceType.PRE_SEQ, ActivityPrecedenceType.POST_LOOP, null, Matrix.singleton(nloops));
    }

    public static ActivityPrecedence Loop(String preAct, List<String> postActs, String endAct, Matrix nloops) {
        List<String> postActsNew = new ArrayList<>(postActs);
        postActsNew.add(endAct);
        return Loop(preAct, postActsNew, nloops);
    }

    public static ActivityPrecedence Loop(Activity preAct, List<Activity> postActs, Activity endAct, Matrix nloops) {
        List<Activity> postActsNew = new ArrayList<>(postActs);
        postActsNew.add(endAct);
        return Loop(preAct, postActsNew, nloops);
    }

    /**
     * Creates an ActivityPrecedence object representing a loop relationship.
     *
     * @param preAct   the preceding activity
     * @param postActs a list with loopAct and endAct (see 4-parameter Loop constructor)
     * @param nloops   the number of loops
     * @return an ActivityPrecedence object
     */
    public static ActivityPrecedence Loop(Activity preAct, List<Activity> postActs, Matrix nloops) {
        List<Activity> preActs = new ArrayList<>();
        preActs.add(preAct);
        return fromActivities(preActs, postActs, ActivityPrecedenceType.PRE_SEQ, ActivityPrecedenceType.POST_LOOP, null, nloops);
    }

    /**
     * Creates an ActivityPrecedence object representing a loop relationship.
     *
     * @param preAct   the preceding activity
     * @param postActs a list with loopAct and endAct (see 4-parameter Loop constructor)
     * @param nloops   the number of loops
     * @return an ActivityPrecedence object
     */
    public static ActivityPrecedence Loop(String preAct, List<String> postActs, Matrix nloops) {
        List<String> preActs = new ArrayList<>();
        preActs.add(preAct);
        return new ActivityPrecedence(preActs, postActs, ActivityPrecedenceType.PRE_SEQ, ActivityPrecedenceType.POST_LOOP, null, nloops);
    }

    /**
     * Creates an ActivityPrecedence object representing a loop relationship.
     *
     * @param preAct  the preceding activity
     * @param loopAct the activity that is looped
     * @param endAct  the activity after the loop completes
     * @param nloops  the number of loops
     * @return an ActivityPrecedence object
     */
    public static ActivityPrecedence Loop(Activity preAct, Activity loopAct, Activity endAct, Matrix nloops) {
        List<Activity> preActs = new ArrayList<>();
        preActs.add(preAct);
        List<Activity> postActs = new ArrayList<>();
        postActs.add(loopAct);
        postActs.add(endAct);
        return fromActivities(preActs, postActs, ActivityPrecedenceType.PRE_SEQ, ActivityPrecedenceType.POST_LOOP, null, nloops);
    }

    /**
     * Creates an ActivityPrecedence object representing a loop relationship.
     *
     * @param preAct  the preceding activity
     * @param loopAct the activity that is looped
     * @param endAct  the activity after the loop completes
     * @param nloops  the number of loops
     * @return an ActivityPrecedence object
     */
    public static ActivityPrecedence Loop(String preAct, String loopAct, String endAct, Matrix nloops) {
        List<String> preActs = new ArrayList<>();
        preActs.add(preAct);
        List<String> postActs = new ArrayList<>();
        postActs.add(loopAct);
        postActs.add(endAct);
        return new ActivityPrecedence(preActs, postActs, ActivityPrecedenceType.PRE_SEQ, ActivityPrecedenceType.POST_LOOP, null, nloops);
    }

    /**
     * Creates an ActivityPrecedence object representing an OR-fork relationship with a specified probability matrix.
     *
     * @param preAct   the preceding activity
     * @param postActs the list of following activities
     * @param probs    the probability matrix
     * @return an ActivityPrecedence object
     */
    public static ActivityPrecedence OrFork(Activity preAct, List<Activity> postActs, Matrix probs) {
        List<Activity> preActs = new ArrayList<>();
        preActs.add(preAct);
        return fromActivities(preActs, postActs, ActivityPrecedenceType.PRE_SEQ, ActivityPrecedenceType.POST_OR, null, probs);
    }

    /**
     * Creates an ActivityPrecedence object representing an OR-fork relationship with a specified probability matrix.
     *
     * @param preAct   the preceding activity
     * @param postActs the list of following activities
     * @param probs    the probability matrix
     * @return an ActivityPrecedence object
     */
    public static ActivityPrecedence OrFork(String preAct, List<String> postActs, Matrix probs) {
        List<String> preActs = new ArrayList<>();
        preActs.add(preAct);
        return new ActivityPrecedence(preActs, postActs, ActivityPrecedenceType.PRE_SEQ, ActivityPrecedenceType.POST_OR, null, probs);
    }

    /**
     * Creates an ActivityPrecedence object representing an OR-join relationship.
     *
     * @param preActs the list of preceding activities
     * @param postAct the following activity
     * @return an ActivityPrecedence object
     */
    public static ActivityPrecedence OrJoin(List<String> preActs, String postAct) {
        List<String> postActs = new ArrayList<String>();
        postActs.add(postAct);
        return new ActivityPrecedence(preActs, postActs, ActivityPrecedenceType.PRE_OR, ActivityPrecedenceType.POST_SEQ);
    }

    /**
     * Creates an ActivityPrecedence object representing an OR-join relationship.
     *
     * @param preActs the list of preceding activities
     * @param postAct the following activity
     * @return an ActivityPrecedence object
     */
    public static ActivityPrecedence OrJoin(List<Activity> preActs, Activity postAct) {
        List<Activity> postActs = new ArrayList<Activity>();
        postActs.add(postAct);
        return fromActivities(preActs, postActs, ActivityPrecedenceType.PRE_OR, ActivityPrecedenceType.POST_SEQ);
    }

    /**
     * Creates an array of ActivityPrecedence objects representing a sequence.
     *
     * @param varargin the list of activity names
     * @return an array of ActivityPrecedence objects
     */
    public static ActivityPrecedence[] Serial(List<?> varargin) {
        int length = varargin.size();
        if (length < 2) {
            throw new IllegalArgumentException("Serial precedence requires at least 2 activities");
        }

        ActivityPrecedence[] ap = new ActivityPrecedence[length - 1];

        if (varargin.get(0) instanceof String) {
            for (int i = 0; i < length - 1; i++) {
                ap[i] = ActivityPrecedence.Serial((String) varargin.get(i), (String) varargin.get(i + 1));
            }
        } else if (varargin.get(0) instanceof Activity) {
            for (int i = 0; i < length - 1; i++) {
                ap[i] = ActivityPrecedence.Serial((Activity) varargin.get(i), (Activity) varargin.get(i + 1));
            }
        }

        return ap;
    }

    /**
     * Creates an ActivityPrecedence object representing a sequence relationship.
     *
     * @param preAct  the preceding activity
     * @param postAct the following activity
     * @return an ActivityPrecedence object
     */
    public static ActivityPrecedence Serial(String preAct, String postAct) {
        List<String> preActs = new ArrayList<String>();
        preActs.add(preAct);
        List<String> postActs = new ArrayList<String>();
        postActs.add(postAct);
        return new ActivityPrecedence(preActs, postActs);
    }

    /**
     * Creates an ActivityPrecedence object representing a sequence relationship.
     *
     * @param preAct  the preceding activity
     * @param postAct the following activity
     * @return an ActivityPrecedence object
     */
    public static ActivityPrecedence Serial(Activity preAct, Activity postAct) {
        List<Activity> preActs = new ArrayList<Activity>();
        preActs.add(preAct);
        List<Activity> postActs = new ArrayList<Activity>();
        postActs.add(postAct);
        return fromActivities(preActs, postActs);
    }

    /**
     * Creates an array of ActivityPrecedence objects representing a sequence from varargs.
     *
     * @param activities the activities to form a serial sequence
     * @return an array of ActivityPrecedence objects
     */
    public static ActivityPrecedence[] Serial(Activity... activities) {
        List<Activity> activityList = new ArrayList<>();
        for (Activity activity : activities) {
            activityList.add(activity);
        }
        return Serial(activityList);
    }

    /**
     * Creates an array of ActivityPrecedence objects representing a sequence from varargs.
     *
     * @param activities the activity names to form a serial sequence
     * @return an array of ActivityPrecedence objects
     */
    public static ActivityPrecedence[] Serial(String... activities) {
        List<String> activityList = new ArrayList<>();
        for (String activity : activities) {
            activityList.add(activity);
        }
        return Serial(activityList);
    }

    /**
     * Creates an ActivityPrecedence object representing an XOR relationship.
     *
     * @param preAct   the preceding activity
     * @param postActs the list of following activities
     * @param probs    the probability matrix
     * @return an ActivityPrecedence object
     */
    public static ActivityPrecedence Xor(Activity preAct, List<Activity> postActs, Matrix probs) {

        return ActivityPrecedence.OrFork(preAct, postActs, probs);
    }

    /**
     * Creates an ActivityPrecedence object representing an XOR relationship.
     *
     * @param preAct   the preceding activity
     * @param postActs the list of following activities
     * @param probs    the probability matrix
     * @return an ActivityPrecedence object
     */
    public static ActivityPrecedence Xor(String preAct, List<String> postActs, Matrix probs) {

        return ActivityPrecedence.OrFork(preAct, postActs, probs);
    }

    // Getter methods

    /**
     * Returns the list of preceding activity names.
     *
     * @return the list of preceding activity names
     */
    public List<String> getPreActs() {
        return preActs;
    }

    /**
     * Returns the list of following activity names.
     *
     * @return the list of following activity names
     */
    public List<String> getPostActs() {
        return postActs;
    }

    /**
     * Returns the type of the precedence relationship before the activity.
     *
     * @return the pre-type string
     */
    public String getPreType() {
        return preType;
    }

    /**
     * Returns the type of the precedence relationship after the activity.
     *
     * @return the post-type string
     */
    public String getPostType() {
        return postType;
    }

    /**
     * Returns the parameters for the preceding activities.
     *
     * @return the pre-parameters matrix
     */
    public Matrix getPreParams() {
        return preParams;
    }

    /**
     * Returns the parameters for the following activities.
     *
     * @return the post-parameters matrix
     */
    public Matrix getPostParams() {
        return postParams;
    }

    /**
     * Retrieves the precedence ID based on the precedence type string.
     *
     * @param precedence the precedence type string
     * @return the precedence ID
     */
    public static int getPrecedenceId(String precedence) {
        int typeId = -1;

        switch (precedence) {
            case ActivityPrecedenceType.PRE_SEQ:
                typeId = ActivityPrecedenceType.ID_PRE_SEQ;
                break;
            case ActivityPrecedenceType.PRE_AND:
                typeId = ActivityPrecedenceType.ID_PRE_AND;
                break;
            case ActivityPrecedenceType.PRE_OR:
                typeId = ActivityPrecedenceType.ID_PRE_OR;
                break;
            case ActivityPrecedenceType.POST_AND:
                typeId = ActivityPrecedenceType.ID_POST_AND;
                break;
            case ActivityPrecedenceType.POST_SEQ:
                typeId = ActivityPrecedenceType.ID_POST_SEQ;
                break;
            case ActivityPrecedenceType.POST_CACHE:
                typeId = ActivityPrecedenceType.ID_POST_CACHE;
                break;
            case ActivityPrecedenceType.POST_OR:
                typeId = ActivityPrecedenceType.ID_POST_OR;
                break;
            case ActivityPrecedenceType.POST_LOOP:
                typeId = ActivityPrecedenceType.ID_POST_LOOP;
                break;
        }

        return typeId;
    }

}
