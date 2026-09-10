/**
 * Package containing Discrete Phase-Type (DPH) distribution functions.
 * Ported from BUTools V2.0.
 *
 * DPH distributions are the discrete-time analogues of PH distributions.
 * They are characterized by an initial probability vector (alpha) and
 * a substochastic transition probability matrix (A).
 *
 * This package also includes Matrix-Geometric (MG) distributions,
 * which are the discrete-time analogues of Matrix-Exponential (ME)
 * distributions.
 *
 * Functions include:
 * - CheckMGRepresentation: Validate matrix-geometric representation
 * - CheckDPHRepresentation: Validate discrete phase-type representation
 * - MomentsFromMG: Compute moments from MG distribution
 * - MomentsFromDPH: Compute moments from DPH distribution
 * - CdfFromMG: Compute CDF from MG distribution
 * - CdfFromDPH: Compute CDF from DPH distribution
 * - PmfFromMG: Compute PMF from MG distribution
 * - PmfFromDPH: Compute PMF from DPH distribution
 * - MGFromMoments: Create MG distribution from moments
 * - DPH2From3Moments: Create DPH(2) from 3 moments
 * - DPH3From5Moments: Create DPH(3) from 5 moments
 * - CanonicalFromDPH2: Convert DPH(2) to canonical form
 * - CanonicalFromDPH3: Convert DPH(3) to canonical form
 * - SamplesFromDPH: Generate random samples from DPH distribution
 * - RandomDPH: Generate random DPH distribution
 * - DPHFromMG: Transform MG to DPH representation
 * - AcyclicDPHFromMG: Transform MG to acyclic DPH
 */
package jline.lib.butools.dph;
