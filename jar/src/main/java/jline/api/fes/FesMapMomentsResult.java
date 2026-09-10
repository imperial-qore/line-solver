/**
 * @file Descriptors of an inter-departure MAP
 *
 * @since LINE 3.0
 */
package jline.api.fes;

/**
 * Moments and index of dispersion of an inter-departure MAP.
 *
 * These are the four descriptors a MAP(2) is fitted against in Section 5.2.2 of Casale,
 * Mi, Cherkasova and Smirni, IEEE Trans. Soft. Eng. 37(5), 2011.
 *
 * Copyright (c) 2012-2026, Imperial College London
 * All rights reserved.
 */
public final class FesMapMomentsResult {
    /** Mean inter-departure time. */
    public final double e1;
    /** Second moment of the inter-departure times. */
    public final double e2;
    /** Third moment of the inter-departure times. */
    public final double e3;
    /** Joint moment of consecutive inter-departure times. */
    public final double e11;
    /** Asymptotic index of dispersion. */
    public final double idc;

    public FesMapMomentsResult(double e1, double e2, double e3, double e11, double idc) {
        this.e1 = e1;
        this.e2 = e2;
        this.e3 = e3;
        this.e11 = e11;
        this.idc = idc;
    }
}
