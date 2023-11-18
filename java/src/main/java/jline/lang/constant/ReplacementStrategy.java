package jline.lang.constant;

/**
 *  Constants for specifying a cache replacement strategy
 */
public enum ReplacementStrategy {
    RR,
    FIFO,
    SFIFO,
    LRU;

    public static int toId(ReplacementStrategy r){
        switch(r){
            case RR:
                return 0;
            case FIFO:
                return 1;
            case SFIFO:
                return 2;
            case LRU:
                return 3;
            default:
                throw new RuntimeException("Unrecognized replacement strategy");
        }
    }

    public static String toFeature(ReplacementStrategy r){
        switch(r){
            case RR:
                return "ReplacementStrategy_RR";
            case FIFO:
                return "ReplacementStrategy_FIFO";
            case SFIFO:
                return "ReplacementStrategy_SFIFO";
            case LRU:
                return "ReplacementStrategy_LRU";
            default:
                throw new RuntimeException("Unrecognized replacement strategy");
        }
    }
}
