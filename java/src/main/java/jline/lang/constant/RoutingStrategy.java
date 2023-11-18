package jline.lang.constant;

/**
 *  Constants for specifying routing strategies at stations
 */
public enum RoutingStrategy {
    RAND,
    PROB,
    RROBIN,
    WRROBIN,
    JSQ,
    DISABLED,
    FIRING,
    KCHOICES;

    public static String toFeature(RoutingStrategy routing){
        switch (routing){
            case RAND:
                return "RoutingStrategy_RAND";
            case PROB:
                return "RoutingStrategy_PROB";
            case RROBIN:
                return "RoutingStrategy_RROBIN";
            case WRROBIN:
                return "RoutingStrategy_WRROBIN";
            case KCHOICES:
                return "RoutingStrategy_KCHOICES";
            default:
                return "";
        }
    }
}
