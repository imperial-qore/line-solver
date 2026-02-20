/*
 * Copyright (c) 2012-2026, QORE Lab, Imperial College London
 * All rights reserved.
 */

package jline.lang.sections;

import java.io.Serializable;

/**
 * A section that forwards jobs without introducing delays in a Log node
 */
public class LogTunnel extends ServiceTunnel implements Serializable {
    public LogTunnel(String name) {
        super(name);

        this.numberOfServers = 1;
    }

    public LogTunnel() {
        this("LogTunnel");
    }
}
