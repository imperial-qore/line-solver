/*
 * Copyright (c) 2012-2026, QORE Lab, Imperial College London
 * All rights reserved.
 */

package jline.lang.reward;

import jline.lang.JobClass;
import jline.lang.NetworkStruct;
import jline.lang.nodes.Node;
import jline.util.matrix.Matrix;

/**
 * Factory for common reward function templates.
 *
 * The templates return a {@link RewardDescriptor}, which is itself a
 * {@link RewardFunction} and can therefore be passed straight to
 * {@code Network.setReward}. Unlike a plain lambda, a descriptor records what the
 * reward measures, which is what makes the reward serializable to JSON.
 *
 * <pre>
 * model.setReward("QLen_Q1", Reward.queueLength(queue1));
 * model.setReward("QLen_Q1_C1", Reward.queueLength(queue1, class1));
 * model.setReward("Util_Q1", Reward.utilization(queue1));
 * model.setReward("Block_Q1", Reward.blocking(queue1));
 * </pre>
 *
 * The reward functions operate on a row of the aggregate state space, laid out as
 * [n_{1,1}, ..., n_{1,K}, n_{2,1}, ..., n_{M,K}], matching {@link RewardFunction}.
 *
 * @see RewardDescriptor
 */
public class Reward {

    private Reward() {
    }

    /**
     * Resolves the station index of a node through the network structure.
     */
    private static int stationIndexOf(Node node, NetworkStruct sn) {
        int nodeIdx = node.getNodeIndex();
        // nodeToStation is read through the linear accessor, as everywhere else in the
        // codebase: it is built as a row vector in some paths and a column vector in
        // others, so a (row, col) read would be orientation-dependent.
        if (sn.nodeToStation == null || nodeIdx < 0 || nodeIdx >= sn.nodeToStation.getNumElements()) {
            throw new IllegalArgumentException(
                    "Reward refers to node \"" + node.getName() + "\", which is not a station in this model.");
        }
        int stationIdx = (int) sn.nodeToStation.get(nodeIdx);
        if (stationIdx < 0) {
            throw new IllegalArgumentException(
                    "Reward refers to node \"" + node.getName() + "\", which is not a station in this model.");
        }
        return stationIdx;
    }

    /**
     * Resolves the class index of a job class through the owning network.
     */
    private static int classIndexOf(Node node, JobClass jobClass) {
        int classIdx = node.getModel().getJobClassIndex(jobClass);
        if (classIdx < 0) {
            throw new IllegalArgumentException(
                    "Reward refers to class \"" + jobClass.getName() + "\", which is not defined in this model.");
        }
        return classIdx;
    }

    /**
     * Reads the number of jobs of a given class at a given station from a state row.
     */
    private static double jobsAt(Matrix state, int stationIdx, int classIdx, int nclasses) {
        return state.get(0, stationIdx * nclasses + classIdx);
    }

    /**
     * Sums the jobs of all classes at a given station in a state row.
     */
    private static double totalJobsAt(Matrix state, int stationIdx, int nclasses) {
        double total = 0;
        for (int r = 0; r < nclasses; r++) {
            total += jobsAt(state, stationIdx, r, nclasses);
        }
        return total;
    }

    /**
     * Queue length reward: total number of jobs at a node, over all classes.
     *
     * @param node the station
     * @return a serializable reward descriptor
     */
    public static RewardDescriptor queueLength(final Node node) {
        RewardFunction fn = new RewardFunction() {
            private static final long serialVersionUID = 1L;

            @Override
            public double compute(Matrix state, NetworkStruct sn) {
                return totalJobsAt(state, stationIndexOf(node, sn), sn.nclasses);
            }
        };
        return new RewardDescriptor(RewardDescriptor.Kind.QLen, node, null, fn);
    }

    /**
     * Queue length reward: number of jobs of a given class at a node.
     *
     * @param node     the station
     * @param jobClass the job class
     * @return a serializable reward descriptor
     */
    public static RewardDescriptor queueLength(final Node node, final JobClass jobClass) {
        if (jobClass == null) {
            return queueLength(node);
        }
        RewardFunction fn = new RewardFunction() {
            private static final long serialVersionUID = 1L;

            @Override
            public double compute(Matrix state, NetworkStruct sn) {
                return jobsAt(state, stationIndexOf(node, sn), classIndexOf(node, jobClass), sn.nclasses);
            }
        };
        return new RewardDescriptor(RewardDescriptor.Kind.QLen, node, jobClass, fn);
    }

    /**
     * Utilization reward: min(jobs, nservers) at a node, over all classes.
     *
     * @param node the station
     * @return a serializable reward descriptor
     */
    public static RewardDescriptor utilization(final Node node) {
        RewardFunction fn = new RewardFunction() {
            private static final long serialVersionUID = 1L;

            @Override
            public double compute(Matrix state, NetworkStruct sn) {
                int stationIdx = stationIndexOf(node, sn);
                double jobs = totalJobsAt(state, stationIdx, sn.nclasses);
                return Math.min(jobs, sn.nservers.get(stationIdx));
            }
        };
        return new RewardDescriptor(RewardDescriptor.Kind.Util, node, null, fn);
    }

    /**
     * Utilization reward: min(class jobs, nservers) at a node.
     *
     * @param node     the station
     * @param jobClass the job class
     * @return a serializable reward descriptor
     */
    public static RewardDescriptor utilization(final Node node, final JobClass jobClass) {
        if (jobClass == null) {
            return utilization(node);
        }
        RewardFunction fn = new RewardFunction() {
            private static final long serialVersionUID = 1L;

            @Override
            public double compute(Matrix state, NetworkStruct sn) {
                int stationIdx = stationIndexOf(node, sn);
                double jobs = jobsAt(state, stationIdx, classIndexOf(node, jobClass), sn.nclasses);
                return Math.min(jobs, sn.nservers.get(stationIdx));
            }
        };
        return new RewardDescriptor(RewardDescriptor.Kind.Util, node, jobClass, fn);
    }

    /**
     * Blocking reward: 1 when the node holds at least as many jobs as its capacity, 0 otherwise.
     *
     * @param node the station
     * @return a serializable reward descriptor
     */
    public static RewardDescriptor blocking(final Node node) {
        RewardFunction fn = new RewardFunction() {
            private static final long serialVersionUID = 1L;

            @Override
            public double compute(Matrix state, NetworkStruct sn) {
                int stationIdx = stationIndexOf(node, sn);
                double jobs = totalJobsAt(state, stationIdx, sn.nclasses);
                return (jobs >= sn.cap.get(stationIdx)) ? 1.0 : 0.0;
            }
        };
        return new RewardDescriptor(RewardDescriptor.Kind.Blocking, node, null, fn);
    }

    /**
     * Wraps an arbitrary user reward function.
     *
     * The result is marked as {@link RewardDescriptor.Kind#Custom}, which is deliberately
     * not serializable: an arbitrary function cannot be reproduced from JSON, so the
     * writer warns and omits it rather than emitting a wrong reward.
     *
     * @param userFn the user reward function
     * @return a reward descriptor of kind Custom
     */
    public static RewardDescriptor custom(RewardFunction userFn) {
        return new RewardDescriptor(RewardDescriptor.Kind.Custom, null, null, userFn);
    }
}
