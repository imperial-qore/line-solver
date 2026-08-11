/*
 * Copyright (c) 2012-2026, Imperial College London
 * All rights reserved.
 */
package jline.lib.fjcodes;

import java.util.Objects;

import jline.util.matrix.Matrix;

/**
 * Service representation for FJ_codes.
 */
public final class FJService {
    public final double mu;          // Service rate
    public final Matrix ST;          // PH sub-generator
    private final Matrix St;          // PH exit rate vector (column)
    public final Matrix tau_st;      // PH initial probability (column)
    private final int SerChoice;      // Service type: 1=Exponential, other=PH

    public FJService(double mu, Matrix ST, Matrix St, Matrix tau_st, int SerChoice) {
        this.mu = mu;
        this.ST = ST;
        this.St = St;
        this.tau_st = tau_st;
        this.SerChoice = SerChoice;
    }

    public double getMu() { return mu; }
    public Matrix getST() { return ST; }
    public Matrix getSt() { return St; }
    public Matrix getTau_st() { return tau_st; }
    public int getSerChoice() { return SerChoice; }

    @Override
    public boolean equals(Object o) {
        if (this == o) return true;
        if (!(o instanceof FJService)) return false;
        FJService that = (FJService) o;
        return Double.compare(that.mu, mu) == 0
                && SerChoice == that.SerChoice
                && Objects.equals(ST, that.ST)
                && Objects.equals(St, that.St)
                && Objects.equals(tau_st, that.tau_st);
    }

    @Override
    public int hashCode() {
        return Objects.hash(mu, ST, St, tau_st, SerChoice);
    }

    @Override
    public String toString() {
        return "FJService(mu=" + mu + ", ST=" + ST + ", St=" + St
                + ", tau_st=" + tau_st + ", SerChoice=" + SerChoice + ")";
    }
}
