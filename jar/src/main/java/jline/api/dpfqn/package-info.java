/*
 * Copyright (c) 2012-2026, QORE Lab, Imperial College London
 * All rights reserved.
 */
/**
 * Discrete-time (slotted) product-form queueing network algorithms.
 *
 * <p>Counterparts of {@link jline.api.pfqn} for models that live on a lattice
 * of unit slots rather than on the continuous time axis. In each slot a busy
 * server completes with probability p and an arrival occurs with probability b,
 * both recorded at the end of the slot with the departure resolved before the
 * arrival (Daduna's LA rule and D/A rule).</p>
 *
 * <p>Reference: H. Daduna, Queueing Networks with Discrete Time Scale, LNCS
 * 2046, Springer, 2001.</p>
 */
package jline.api.dpfqn;
