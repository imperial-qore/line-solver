/**
 * Procedural algorithms for solving stochastic models and analyzing queueing systems.
 *
 * <p>This package provides low-level algorithmic implementations for various analysis methods
 * used in queueing theory, stochastic modeling, and performance evaluation. The algorithms
 * are organized into specialized subpackages based on the modeling domain and analysis technique.
 *
 * <h2>Subpackages</h2>
 * <ul>
 * <li><strong>{@link jline.api.cache}</strong> - Cache modeling algorithms (LRU, TTL, hit/miss probabilities)</li>
 * <li><strong>{@link jline.api.lossn}</strong> - Loss network algorithms (Erlang formulas, blocking probabilities)</li>
 * <li><strong>{@link jline.api.lsn}</strong> - Layered stochastic network (LayeredNetworkStruct) data structure manipulation</li>
 * <li><strong>{@link jline.api.mam}</strong> - Matrix Analytic Methods (MAP, APH, MMAP processes)</li>
 * <li><strong>{@link jline.api.mapqn}</strong> - MAP-based queueing network analysis algorithms</li>
 * <li><strong>{@link jline.api.mc}</strong> - Markov Chain algorithms (CTMC, DTMC analysis)</li>
 * <li><strong>{@link jline.api.measures}</strong> - Distance and similarity metrics for performance analysis</li>
 * <li><strong>{@link jline.api.npfqn}</strong> - Non-Product Form Queueing Network algorithms</li>
 * <li><strong>{@link jline.api.pfqn}</strong> - Product Form Queueing Network algorithms (MVA, convolution, normalizing constant, load-dependent)</li>
 * <li><strong>{@link jline.api.polling}</strong> - Polling system algorithms (gated, exhaustive, limited service)</li>
 * <li><strong>{@link jline.api.qsys}</strong> - Single queue system algorithms (M/M/1, M/G/1, GI/GI/1, etc.)</li>
 * <li><strong>{@link jline.api.rl}</strong> - Reinforcement learning and optimization algorithms</li>
 * <li><strong>{@link jline.api.sn}</strong> - Stochastic network data structure (NetworkStruct) manipulation</li>
 * <li><strong>{@link jline.api.trace}</strong> - Trace analysis algorithms (moments, statistics)</li>
 * <li><strong>{@link jline.api.wf}</strong> - Workflow and job scheduling analysis algorithms</li>
 * </ul>
 *
 * <p><strong>Note:</strong> These are low-level procedural algorithms primarily implemented in Kotlin.
 * For high-level object-oriented modeling, use the {@link jline.lang} package. For complete solvers,
 * use the {@link jline.solvers} package.
 *
 * @see jline.lang
 * @see jline.solvers
 * @since LINE 2.0
 */
package jline.api;