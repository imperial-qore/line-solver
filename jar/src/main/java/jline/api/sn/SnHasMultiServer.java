package jline.api.sn;

import jline.lang.NetworkStruct;

public final class SnHasMultiServer {
    private SnHasMultiServer() {}

    public static boolean snHasMultiServer(NetworkStruct sn) {
        return sn.nservers.elementMax() > 1;
    }
}
