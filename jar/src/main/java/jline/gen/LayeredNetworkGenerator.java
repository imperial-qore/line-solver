/*
 * Copyright (c) 2012-2026, QORE Lab, Imperial College London
 * All rights reserved.
 */

package jline.gen;

import jline.GlobalConstants;
import jline.lang.constant.SchedStrategy;
import jline.lang.layered.Activity;
import jline.lang.layered.Entry;
import jline.lang.layered.LayeredNetwork;
import jline.lang.layered.Processor;
import jline.lang.layered.Task;
import jline.lang.processes.Distribution;
import jline.lang.processes.Exp;
import jline.lang.processes.Immediate;

import java.util.ArrayList;
import java.util.List;
import java.util.Random;

/**
 * A generator object that generates layered queueing network models
 * based on user specification. Characteristics of generated
 * models can be configured via the generator's properties.
 */
public class LayeredNetworkGenerator {

    private double[] populationRange;
    private double[] thinkTimeRange;
    private double taskInfProbability;
    private double procInfProbability;
    private double[] taskMultiRange;
    private double[] procMultiRange;
    private double[] hostDemandRange;
    private double[] synchCallRange;

    private List<Activity> cActivities;
    private List<Entry> cEntries;
    private List<Task> cTasks;
    private List<Processor> cProcessors;
    private List<Activity> activities;
    private List<Entry> entries;
    private List<Task> tasks;
    private List<Processor> processors;
    private int[] numTasksPerLevel;
    private int[] numTasksPerProcessor;

    private Random random = new Random();

    /**
     * Constructor with default settings
     */
    public LayeredNetworkGenerator() {
        this(new double[]{1, 1}, new double[]{1, 1}, 0, 0, new double[]{1, 1},
             new double[]{1, 1}, new double[]{1, 1}, new double[]{1, 1});
    }

    /**
     * Reseeds this generator's internal random source for reproducible generation.
     *
     * @param seed the seed value
     * @return this generator (for chaining)
     */
    public LayeredNetworkGenerator setSeed(long seed) {
        this.random = new Random(seed);
        return this;
    }

    /**
     * Constructor with custom settings
     *
     * @param populationRange    range of reference-task populations
     * @param thinkTimeRange     range of client think times
     * @param taskInfProbability probability that a task is infinite-server
     * @param procInfProbability probability that a processor is infinite-server
     * @param taskMultiRange     range of task multiplicities
     * @param procMultiRange     range of processor multiplicities
     * @param hostDemandRange    range of activity host demands
     * @param synchCallRange     range of synchronous call means
     */
    public LayeredNetworkGenerator(double[] populationRange, double[] thinkTimeRange,
                                   double taskInfProbability, double procInfProbability,
                                   double[] taskMultiRange, double[] procMultiRange,
                                   double[] hostDemandRange, double[] synchCallRange) {
        setPopulationRange(populationRange);
        setThinkTimeRange(thinkTimeRange);
        setTaskInfProbability(taskInfProbability);
        setProcInfProbability(procInfProbability);
        setTaskMultiRange(taskMultiRange);
        setProcMultiRange(procMultiRange);
        setHostDemandRange(hostDemandRange);
        setSynchCallRange(synchCallRange);
    }

    /**
     * Main function to call. Returns a generated LQN model according to
     * specified properties of the LayeredNetworkGenerator object
     *
     * @param numClients    number of client (reference) tasks
     * @param numLevels     number of task layers
     * @param numTasks      number of non-client tasks, distributed over the levels
     * @param numProcessors number of processors hosting the non-client tasks
     * @return the generated layered network
     */
    public LayeredNetwork generate(int numClients, int numLevels, int numTasks, int numProcessors) {
        validateArgs(numClients, numLevels, numTasks, numProcessors);
        LayeredNetwork model = new LayeredNetwork("lnw");
        createClients(model, numClients);
        createTasks(model, numTasks);
        createProcessors(model, numProcessors);
        assignTasks(numLevels, numTasks, numProcessors);
        connectClientsToTasks(numClients);
        connectTasksToTasks(numLevels);
        connectTasksToProcessors(numProcessors);
        return model;
    }

    // Property setters with validation

    public void setPopulationRange(double[] range) {
        this.populationRange = checkRange(range, "Population range", true);
    }

    public void setThinkTimeRange(double[] range) {
        this.thinkTimeRange = checkRange(range, "Think time range", false);
    }

    public void setTaskInfProbability(double probability) {
        this.taskInfProbability = checkProbability(probability, "Task infinite probability");
    }

    public void setProcInfProbability(double probability) {
        this.procInfProbability = checkProbability(probability, "Processor infinite probability");
    }

    public void setTaskMultiRange(double[] range) {
        this.taskMultiRange = checkRange(range, "Task multiplicity range", true);
    }

    public void setProcMultiRange(double[] range) {
        this.procMultiRange = checkRange(range, "Processor multiplicity range", true);
    }

    public void setHostDemandRange(double[] range) {
        this.hostDemandRange = checkRange(range, "Host demand range", false);
    }

    public void setSynchCallRange(double[] range) {
        this.synchCallRange = checkRange(range, "Synchronous call range", true);
    }

    // Getters

    public double[] getPopulationRange() { return populationRange; }
    public double[] getThinkTimeRange() { return thinkTimeRange; }
    public double getTaskInfProbability() { return taskInfProbability; }
    public double getProcInfProbability() { return procInfProbability; }
    public double[] getTaskMultiRange() { return taskMultiRange; }
    public double[] getProcMultiRange() { return procMultiRange; }
    public double[] getHostDemandRange() { return hostDemandRange; }
    public double[] getSynchCallRange() { return synchCallRange; }

    // Private methods

    private static double[] checkRange(double[] range, String name, boolean strictlyPositive) {
        if (range == null || range.length != 2) {
            throw new IllegalArgumentException(name + " must have two elements");
        }
        boolean lowerOk = strictlyPositive ? range[0] > 0 : range[0] >= 0;
        if (!(lowerOk && range[0] <= range[1])) {
            throw new IllegalArgumentException(name + " is not valid");
        }
        return new double[]{range[0], range[1]};
    }

    private static double checkProbability(double value, String name) {
        if (!(value >= 0 && value <= 1)) {
            throw new IllegalArgumentException(name + " is not valid");
        }
        return value;
    }

    /**
     * Validates that parameter values for the network are sound
     */
    private void validateArgs(int numClients, int numLevels, int numTasks, int numProcessors) {
        if (numClients < 1) {
            throw new IllegalArgumentException("Number of clients is less than one");
        } else if (numLevels < 1) {
            throw new IllegalArgumentException("Number of levels is less than one");
        } else if (numTasks < 1) {
            throw new IllegalArgumentException("Number of tasks is less than one");
        } else if (numProcessors < 1) {
            throw new IllegalArgumentException("Number of processors is less than one");
        } else if (numLevels > numTasks) {
            throw new IllegalArgumentException("Number of levels is greater than that of tasks");
        } else if (numProcessors > numTasks) {
            throw new IllegalArgumentException("Number of processors is greater than that of tasks");
        }
    }

    /**
     * Creates the clients in the layered network
     */
    private void createClients(LayeredNetwork model, int numClients) {
        cActivities = new ArrayList<Activity>();
        cEntries = new ArrayList<Entry>();
        cTasks = new ArrayList<Task>();
        cProcessors = new ArrayList<Processor>();

        for (int c = 0; c < numClients; c++) {
            int population = sampleIntegerValue(populationRange);
            double thinkTime = sampleRealValue(thinkTimeRange);

            Activity activity = new Activity(model, "c_activity_" + (c + 1), asDistribution(thinkTime));
            Entry entry = new Entry(model, "c_entry_" + (c + 1));
            Task task = new Task(model, "c_task_" + (c + 1), population, SchedStrategy.REF);
            Processor processor = new Processor(model, "c_processor_" + (c + 1),
                                                Integer.MAX_VALUE, SchedStrategy.INF);

            activity.on(task).boundTo(entry);
            entry.on(task);
            task.on(processor);

            cActivities.add(activity);
            cEntries.add(entry);
            cTasks.add(task);
            cProcessors.add(processor);
        }
    }

    /**
     * Creates the tasks in the layered network
     */
    private void createTasks(LayeredNetwork model, int numTasks) {
        activities = new ArrayList<Activity>();
        entries = new ArrayList<Entry>();
        tasks = new ArrayList<Task>();

        for (int t = 0; t < numTasks; t++) {
            int multiplicity;
            SchedStrategy scheduling;
            if (chooseBooleanValue(taskInfProbability)) {
                multiplicity = Integer.MAX_VALUE;
                scheduling = SchedStrategy.INF;
            } else {
                multiplicity = sampleIntegerValue(taskMultiRange);
                scheduling = SchedStrategy.FCFS;
            }
            double hostDemand = sampleRealValue(hostDemandRange);

            Activity activity = new Activity(model, "activity_" + (t + 1), asDistribution(hostDemand));
            Entry entry = new Entry(model, "entry_" + (t + 1));
            Task task = new Task(model, "task_" + (t + 1), multiplicity, scheduling);

            activity.on(task).boundTo(entry).repliesTo(entry);
            entry.on(task);

            activities.add(activity);
            entries.add(entry);
            tasks.add(task);
        }
    }

    /**
     * Creates the processors in the layered network
     */
    private void createProcessors(LayeredNetwork model, int numProcessors) {
        processors = new ArrayList<Processor>();

        for (int i = 0; i < numProcessors; i++) {
            int multiplicity;
            SchedStrategy scheduling;
            if (chooseBooleanValue(procInfProbability)) {
                multiplicity = Integer.MAX_VALUE;
                scheduling = SchedStrategy.INF;
            } else {
                multiplicity = sampleIntegerValue(procMultiRange);
                scheduling = SchedStrategy.PS;
            }
            processors.add(new Processor(model, "processor_" + (i + 1), multiplicity, scheduling));
        }
    }

    /**
     * Assigns the tasks to different levels and processors
     */
    private void assignTasks(int numLevels, int numTasks, int numProcessors) {
        numTasksPerLevel = makeIntegerVector(numLevels, numTasks);
        numTasksPerProcessor = makeIntegerVector(numProcessors, numTasks);
    }

    /**
     * Connects the clients to the first-level tasks
     */
    private void connectClientsToTasks(int numClients) {
        boolean[] clientConnected = new boolean[numClients];

        for (int t = 0; t < numTasksPerLevel[0]; t++) {
            double synchCall = sampleRealValue(synchCallRange);

            int c = sampleIntegerValue(new double[]{1, numClients}) - 1;
            cActivities.get(c).synchCall(entries.get(t), synchCall);
            clientConnected[c] = true;
        }

        for (int c = 0; c < numClients; c++) {
            if (clientConnected[c]) {
                continue;
            }

            double synchCall = sampleRealValue(synchCallRange);

            int t = sampleIntegerValue(new double[]{1, numTasksPerLevel[0]}) - 1;
            cActivities.get(c).synchCall(entries.get(t), synchCall);
            clientConnected[c] = true;
        }
    }

    /**
     * Connects the tasks between adjacent levels
     */
    private void connectTasksToTasks(int numLevels) {
        int numTasks = numTasksPerLevel[0];
        for (int l = 1; l < numLevels; l++) {
            for (int t2 = numTasks; t2 < numTasks + numTasksPerLevel[l]; t2++) {
                double synchCall = sampleRealValue(synchCallRange);

                // Caller drawn among the tasks of the previous level
                int t1 = sampleIntegerValue(new double[]{numTasks - numTasksPerLevel[l - 1] + 1,
                                                        numTasks}) - 1;
                activities.get(t1).synchCall(entries.get(t2), synchCall);
            }
            numTasks += numTasksPerLevel[l];
        }
    }

    /**
     * Connects the tasks to the processors
     */
    private void connectTasksToProcessors(int numProcessors) {
        int numTasks = 0;
        for (int p = 0; p < numProcessors; p++) {
            for (int t = numTasks; t < numTasks + numTasksPerProcessor[p]; t++) {
                tasks.get(t).on(processors.get(p));
            }
            numTasks += numTasksPerProcessor[p];
        }
    }

    /**
     * Coerces a scalar mean into a distribution, as MATLAB's setHostDemand does
     */
    private static Distribution asDistribution(double mean) {
        if (mean <= GlobalConstants.FineTol) {
            return new Immediate();
        }
        return new Exp(1.0 / mean);
    }

    /**
     * Samples an integer value from a given range
     *
     * @param range two-element range, whose bounds are rounded inwards
     * @return an integer in [ceil(range[0]), floor(range[1])]
     */
    public int sampleIntegerValue(double[] range) {
        int lower = (int) Math.ceil(range[0]);
        int upper = (int) Math.floor(range[1]);
        return lower + random.nextInt(upper - lower + 1);
    }

    /**
     * Samples a real value from a given range
     *
     * @param range two-element range
     * @return a real value in [range[0], range[1]]
     */
    public double sampleRealValue(double[] range) {
        return range[0] + (range[1] - range[0]) * random.nextDouble();
    }

    /**
     * Chooses a boolean value for a given probability
     *
     * @param probability probability of returning true
     * @return true with the given probability
     */
    public boolean chooseBooleanValue(double probability) {
        return random.nextDouble() < probability;
    }

    /**
     * Makes a vector of positive integers with given length and sum
     *
     * @param length number of elements
     * @param sum    required sum of the elements
     * @return a vector of length positive integers summing to sum
     */
    public int[] makeIntegerVector(int length, int sum) {
        int[] vector = new int[length];
        for (int i = 0; i < length; i++) {
            vector[i] = 1;
        }
        for (int s = 0; s < sum - length; s++) {
            vector[random.nextInt(length)]++;
        }
        return vector;
    }
}
