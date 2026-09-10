/*
 * Copyright (c) 2012-2026, QORE Lab, Imperial College London
 * All rights reserved.
 */
/**
 * Discrete-time (slotted) single-queue formulas.
 *
 * <p>Every function here observes the system on a lattice of unit slots rather
 * than on the continuous time axis, under the late-arrival rule with the
 * departure resolved before the arrival (Daduna's LA and D/A rules). The
 * continuous-time counterparts live in {@link jline.api.qsys}.</p>
 */
package jline.api.dqsys;
