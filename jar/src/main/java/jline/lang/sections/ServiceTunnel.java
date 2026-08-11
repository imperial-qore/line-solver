/*
 * Copyright (c) 2012-2026, QORE Lab, Imperial College London
 * All rights reserved.
 */

package jline.lang.sections;

import java.io.Serializable;

/**
 * A section that forwards jobs without introducing delays in a service station
 */
public class ServiceTunnel extends ServiceSection implements Serializable {
    public ServiceTunnel(String name) {
        super(name);

        this.numberOfServers = 1;
    }

    public ServiceTunnel() {
        this("ServiceTunnel");
    }
}
