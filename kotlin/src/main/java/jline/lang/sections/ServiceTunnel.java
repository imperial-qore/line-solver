package jline.lang.sections;

import java.io.Serializable;

/**
 * A section that forwards jobs without introducing delays
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
