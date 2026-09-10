/*
 * Copyright (c) 2012-2026, QORE Lab, Imperial College London
 * All rights reserved.
 */
package jline.lib.ode;

import java.util.Arrays;

/**
 * Result returned by RKF45.integrate.
 */
public final class RKF45Result {
    public final double[] y;
    public final double[] yp;
    public final double t;
    public final int flag;
    public final double relerr;

    public RKF45Result(double[] y, double[] yp, double t, int flag, double relerr) {
        this.y = y;
        this.yp = yp;
        this.t = t;
        this.flag = flag;
        this.relerr = relerr;
    }

    public double[] getY() { return y; }
    public double[] getYp() { return yp; }
    public double getT() { return t; }
    public int getFlag() { return flag; }
    public double getRelerr() { return relerr; }

    @Override
    public boolean equals(Object o) {
        if (this == o) return true;
        if (!(o instanceof RKF45Result)) return false;
        RKF45Result other = (RKF45Result) o;
        return Double.compare(t, other.t) == 0
                && flag == other.flag
                && Double.compare(relerr, other.relerr) == 0
                && Arrays.equals(y, other.y)
                && Arrays.equals(yp, other.yp);
    }

    @Override
    public int hashCode() {
        int h = Arrays.hashCode(y);
        h = 31 * h + Arrays.hashCode(yp);
        h = 31 * h + Double.valueOf(t).hashCode();
        h = 31 * h + flag;
        h = 31 * h + Double.valueOf(relerr).hashCode();
        return h;
    }

    @Override
    public String toString() {
        return "RKF45Result(y=" + Arrays.toString(y) + ", yp=" + Arrays.toString(yp)
                + ", t=" + t + ", flag=" + flag + ", relerr=" + relerr + ")";
    }
}
