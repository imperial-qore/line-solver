package jline.lang.nodeparam;

import jline.lang.JobClass;
import jline.lang.constant.PollingType;
import jline.lang.processes.Distribution;

import java.util.HashMap;
import java.util.Map;

/**
 * Parameter container for queue nodes in queueing networks.
 *
 * <p>This class extends {@link ServiceNodeParam} to add queue-specific parameters,
 * particularly for polling scheduling strategies. Queue nodes may use polling
 * disciplines (gated, exhaustive, or k-limited) with associated switchover times
 * between job classes.</p>
 *
 * @see jline.lang.nodes.Queue
 * @see PollingType
 * @since 1.0
 */
public class QueueNodeParam extends ServiceNodeParam {

    /**
     * Polling type for queues with polling scheduling strategy.
     * <p>Specifies the polling discipline: GATED, EXHAUSTIVE, or KLIMITED.</p>
     */
    public PollingType pollingType;

    /**
     * K parameter for K-LIMITED polling.
     * <p>Specifies the maximum number of jobs to serve per polling cycle.</p>
     */
    public Integer pollingPar;

    /**
     * Switchover time distributions by job class for polling servers.
     * <p>Defines the time distribution for switching between job classes.</p>
     */
    public Map<JobClass, Distribution> switchoverTime;

    /**
     * Setup (cold-start) time distributions by job class.
     * <p>Populated from {@code Queue.setDelayOff}, as in MATLAB
     * refreshLocalVars.m. Always paired with {@link #delayoffTime}.</p>
     */
    public Map<JobClass, Distribution> setupTime;

    /**
     * Delay-off (teardown) time distributions by job class.
     * <p>Populated from {@code Queue.setDelayOff}, as in MATLAB
     * refreshLocalVars.m. Always paired with {@link #setupTime}.</p>
     */
    public Map<JobClass, Distribution> delayoffTime;

    /**
     * Memoized derived description of the polling controller (server position,
     * switchover phase-type, visit budget) for a POLLING queue.
     * <p>Built on demand by {@code State.Polling.info}, which is read once per
     * state per synchronization while the CTMC generator is built; rebuilding it
     * per call is orders of magnitude slower than the analysis itself.</p>
     */
    public jline.lang.state.Polling.Info pollinfo;

    /**
     * Pass-and-swap (PAS) class compatibility/swap graph (nclasses x nclasses).
     */
    public jline.util.matrix.Matrix swapGraph;

    /**
     * Pass-and-swap (PAS) total service rate function mu(c) of the ordered state.
     */
    public jline.util.SerializableFunction<jline.util.matrix.Matrix, Double> svcRateFun;

    /**
     * Constructs an empty QueueNodeParam with initialized parameter maps.
     */
    public QueueNodeParam() {
        super();
        switchoverTime = new HashMap<JobClass, Distribution>();
        setupTime = new HashMap<JobClass, Distribution>();
        delayoffTime = new HashMap<JobClass, Distribution>();
    }

    /**
     * Checks if this queue parameter container is empty.
     *
     * @return true if no parameters are specified, false otherwise
     */
    @Override
    public boolean isEmpty() {
        return super.isEmpty() &&
                pollingType == null &&
                pollingPar == null &&
                (switchoverTime == null || switchoverTime.isEmpty());
    }
}
