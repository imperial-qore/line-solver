package jline.api.sn;

import jline.lang.NetworkStruct;

public final class SnHasMultiServer {
    private SnHasMultiServer() {}

    public static boolean snHasMultiServer(NetworkStruct sn) {
        // Infinite servers are delays, not multiserver queues: counting them
        // made every model with a Delay read as multiserver.
        if (sn.nservers == null) {
            return false;
        }
        for (int i = 0; i < sn.nservers.getNumRows(); i++) {
            for (int j = 0; j < sn.nservers.getNumCols(); j++) {
                double c = sn.nservers.get(i, j);
                if (!Double.isInfinite(c) && c > 1) {
                    return true;
                }
            }
        }
        return false;
    }
}
