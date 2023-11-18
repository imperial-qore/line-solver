package jline.lang.constant;

import java.io.Serializable;

/**
 *  Constants for specifying drop strategies at stations
 */
public enum DropStrategy implements Serializable {
    WaitingQueue,
    Drop,
    BlockingAfterService;

    public static String toText(DropStrategy strategy) {
        switch (strategy) {
            case WaitingQueue:
                return "waiting queue";
            case Drop:
                return "drop";
            case BlockingAfterService:
                return "blocking after service";
            default:
                return strategy.name();
        }
    }
}
