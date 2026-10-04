package jline.api.sn;

import jline.lang.NetworkStruct;

public final class SnHasSingleClass {
    private SnHasSingleClass() {}

    public static boolean snHasSingleClass(NetworkStruct sn) {
        return sn.nclasses == 1;
    }
}
