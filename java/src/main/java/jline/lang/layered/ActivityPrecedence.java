package jline.lang.layered;

import jline.lang.constant.ActivityPrecedenceType;
import jline.util.Matrix;

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

    // TODO: this should use activity objects rather than strings
    public ActivityPrecedence(List<String> preActs, List<String> postActs, String preType, String postType, Matrix preParams, Matrix postParams) {
        this.preActs = preActs;
        this.postActs = postActs;
        this.preType = preType;
        this.postType = postType;
        this.preParams = preParams;
        this.postParams = postParams;
    }

    public ActivityPrecedence(List<String> preActs, List<String> postActs, String preType, String postType, Matrix preParams) {
        this.preActs = preActs;
        this.postActs = postActs;
        this.preType = preType;
        this.postType = postType;
        this.preParams = preParams;
        this.postParams = null;
    }

    public ActivityPrecedence(List<String> preActs, List<String> postActs, String preType, String postType) {
        this.preActs = preActs;
        this.postActs = postActs;
        this.preType = preType;
        this.postType = postType;
        this.preParams = null;
        this.postParams = null;
    }

    public ActivityPrecedence(List<String> preActs, List<String> postActs, String preType) {
        this.preActs = preActs;
        this.postActs = postActs;
        this.preType = preType;
        this.postType = ActivityPrecedenceType.POST_SEQ;
        this.preParams = null;
        this.postParams = null;
    }

    public ActivityPrecedence(List<String > preActs, List<String> postActs){
        this.preActs = preActs;
        this.postActs = postActs;
        this.preType = ActivityPrecedenceType.PRE_SEQ;
        this.postType = ActivityPrecedenceType.POST_SEQ;
        this.preParams = null;
        this.postParams = null;
    }

    public ActivityPrecedence(boolean isActivity, List<Activity> preActs, List<Activity> postActs, String preType, String postType, Matrix preParams, Matrix postParams) {
        List<String> preActsName = new ArrayList<>();
        List<String> postActsName = new ArrayList<>();
        for(Activity preAct:preActs){
            preActsName.add(preAct.getName());
        }
        for(Activity postAct: postActs){
            postActsName.add(postAct.getName());
        }
        this.preActs = preActsName;
        this.postActs = postActsName;
        this.preType = preType;
        this.postType = postType;
        this.preParams = preParams;
        this.postParams = postParams;
    }

    public ActivityPrecedence(boolean isActivity, List<Activity> preActs, List<Activity> postActs,String preType, String postType, Matrix preParams) {
        List<String> preActsName = new ArrayList<>();
        List<String> postActsName = new ArrayList<>();
        for(Activity preAct:preActs){
            preActsName.add(preAct.getName());
        }
        for(Activity postAct: postActs){
            postActsName.add(postAct.getName());
        }
        this.preActs = preActsName;
        this.postActs = postActsName;
        this.preType = preType;
        this.postType = postType;
        this.preParams = preParams;
        this.postParams = null;
    }

    public ActivityPrecedence(boolean isActivity, List<Activity> preActs, List<Activity> postActs,String preType, String postType) {
        List<String> preActsName = new ArrayList<>();
        List<String> postActsName = new ArrayList<>();
        for(Activity preAct:preActs){
            preActsName.add(preAct.getName());
        }
        for(Activity postAct: postActs){
            postActsName.add(postAct.getName());
        }
        this.preActs = preActsName;
        this.postActs = postActsName;
        this.preType = preType;
        this.postType = postType;
        this.preParams = null;
        this.postParams = null;
    }

    public ActivityPrecedence(boolean isActivity, List<Activity> preActs, List<Activity> postActs,String preType) {
        List<String> preActsName = new ArrayList<>();
        List<String> postActsName = new ArrayList<>();
        for(Activity preAct:preActs){
            preActsName.add(preAct.getName());
        }
        for(Activity postAct: postActs){
            postActsName.add(postAct.getName());
        }
        this.preActs = preActsName;
        this.postActs = postActsName;
        this.preType = preType;
        this.postType = ActivityPrecedenceType.POST_SEQ;
        this.preParams = null;
        this.postParams = null;
    }

    public ActivityPrecedence(boolean isActivity, List<Activity> preActs, List<Activity> postActs) {
        List<String> preActsName = new ArrayList<>();
        List<String> postActsName = new ArrayList<>();
        for(Activity preAct:preActs){
            preActsName.add(preAct.getName());
        }
        for(Activity postAct: postActs){
            postActsName.add(postAct.getName());
        }
        this.preActs = preActsName;
        this.postActs = postActsName;
        this.preType = ActivityPrecedenceType.PRE_SEQ;
        this.postType = ActivityPrecedenceType.POST_SEQ;
        this.preParams = null;
        this.postParams = null;
    }



    public static ActivityPrecedence[] Serial(List<String> varargin){
        int length = varargin.size();

        ActivityPrecedence[] ap = new ActivityPrecedence[length];//TODO

        for(int i = 0; i<length;i++){
            ap[i] = ActivityPrecedence.Sequence(varargin.get(i),varargin.get(i+1));
        }

        return ap;

    }

    public static int getPrecedenceId(String precedence){
        int typeId = -1;

        switch (precedence){
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

    public static ActivityPrecedence Sequence(String preAct, String postAct) {
        List<String> preActs = new ArrayList<String>();
        preActs.add(preAct);
        List<String> postActs = new ArrayList<String>();
        postActs.add(postAct);
        return new ActivityPrecedence(preActs, postActs);
    }

    public static ActivityPrecedence AndJoin(List<String> preActs, String postAct, Matrix quorum) {
        List<String> postActs = new ArrayList<String>();
        postActs.add(postAct);
        return new ActivityPrecedence(preActs, postActs,ActivityPrecedenceType.PRE_AND,ActivityPrecedenceType.POST_SEQ,quorum,null);
    }

    public static ActivityPrecedence OrJoin(List<String> preActs, String postAct) {
        List<String> postActs = new ArrayList<String>();
        postActs.add(postAct);
        return new ActivityPrecedence(preActs,postActs,ActivityPrecedenceType.PRE_OR,ActivityPrecedenceType.POST_SEQ);
    }

    public static ActivityPrecedence AndFork(String preAct, List<String> postActs, Matrix fanout) {
        List<String> preActs = new ArrayList<>();
        preActs.add(preAct);
        return new ActivityPrecedence(preActs,postActs,ActivityPrecedenceType.PRE_SEQ,ActivityPrecedenceType.POST_AND, null, fanout);
    }

    public static ActivityPrecedence Xor(String preAct, List<String> postActs, Matrix probs) {

        return ActivityPrecedence.OrFork(preAct,postActs,probs);
    }

    public static ActivityPrecedence OrFork(String preAct, List<String> postActs, Matrix probs) {
        List<String> preActs = new ArrayList<>();
        preActs.add(preAct);
        return new ActivityPrecedence(preActs,postActs,ActivityPrecedenceType.PRE_SEQ,ActivityPrecedenceType.POST_OR,null,probs);
    }

    public static ActivityPrecedence Loop(String preAct, List<String> postActs, Matrix counts){
        List<String> preActs = new ArrayList<>();
        preActs.add(preAct);
        return new ActivityPrecedence(preActs,postActs,ActivityPrecedenceType.PRE_SEQ,ActivityPrecedenceType.POST_LOOP,null,counts);
    }

    public static ActivityPrecedence CacheAccess(String preAct, List<String> postActs){
        List<String> preActs = new ArrayList<>();
        preActs.add(preAct);
        return new ActivityPrecedence(preActs,postActs,ActivityPrecedenceType.PRE_SEQ,ActivityPrecedenceType.POST_CACHE);
    }
}
