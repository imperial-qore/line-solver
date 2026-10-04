package jline.api.sn;

import jline.lang.NetworkStruct;

public final class SnHasSingleChain {
    private SnHasSingleChain() {}

    public static boolean snHasSingleChain(NetworkStruct sn) {
        return sn.nchains == 1;
    }
}
