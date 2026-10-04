package jline.api.sn;

import jline.lang.NetworkStruct;

public final class SnHasMultiChain {
    private SnHasMultiChain() {}

    public static boolean snHasMultiChain(NetworkStruct sn) {
        return sn.nchains > 1;
    }
}
