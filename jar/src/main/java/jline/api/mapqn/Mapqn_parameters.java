/**
 * @file MAPQN Parameters Base Class
 *
 * @since LINE 3.0
 */
package jline.api.mapqn;

import jline.lang.NetworkStruct;

/**
 * Base class for MAPQN model parameters.
 */
public abstract class Mapqn_parameters {
    public abstract int getM();
    public abstract int getN();

    public void validate() {
        if (getM() <= 0) throw new IllegalArgumentException("M must be positive");
        if (getN() <= 0) throw new IllegalArgumentException("N must be positive");
    }

    public static Mapqn_parameters fromNetworkStruct(NetworkStruct networkStruct) {
        return Mapqn_parameters_factory.createFromNetworkStruct(networkStruct);
    }
}
