package jline.lang.distributions;

import java.util.Arrays;
import java.util.Collections;

/**
 * This class is a subclass of the general Coxian distribution class, but simplifies the usage for
 * the special case of a two-phase Coxian distribution.
 */
@SuppressWarnings("unchecked")
public class Cox2 extends Coxian {
    public Cox2(double mu1, double mu2, double phi1) {
        super(Arrays.asList(mu1, mu2), Collections.singletonList(phi1));
    }
}