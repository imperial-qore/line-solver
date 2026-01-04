/*
 * Copyright (c) 2012-2026, QORE Lab, Imperial College London
 * All rights reserved.
 */

package jline.lang.constant;

import jline.lang.Element;
import jline.lang.JobClass;
import jline.lang.nodes.Queue;

import java.io.Serializable;
import java.util.ArrayList;
import java.util.List;

/**
 * Represents a type of server within a heterogeneous multiserver queue.
 * <p>
 * A server type defines a group of identical servers with:
 * <ul>
 * <li>A unique name identifying this server type</li>
 * <li>A count of servers of this type</li>
 * <li>A list of job classes that are compatible with (can be served by) this type</li>
 * </ul>
 * <p>
 * Server types enable modeling of heterogeneous multiserver queues where different
 * servers may have different service rates and serve different subsets of job classes.
 * <p>
 * Example usage:
 * <pre>
 * ServerType fastServer = new ServerType("Fast", 2);
 * fastServer.addCompatibleClass(classA);
 * fastServer.addCompatibleClass(classB);
 * queue.addServerType(fastServer);
 * queue.setService(classA, fastServer, new Exp(2.0));
 * </pre>
 */
public class ServerType extends Element implements Serializable {
    private static final long serialVersionUID = 1L;

    /** Unique identifier for this server type within the queue */
    private int id;

    /** Number of servers of this type */
    private int numOfServers;

    /** List of job classes that can be served by this server type */
    private List<JobClass> compatibleClasses;

    /** The queue this server type belongs to (set when added to queue) */
    private Queue parentQueue;

    /**
     * Creates a new server type with the specified name and number of servers.
     * By default, no job classes are compatible until explicitly added.
     *
     * @param name the name identifying this server type (e.g., "Fast", "Slow")
     * @param numOfServers the number of servers of this type (must be >= 1)
     * @throws IllegalArgumentException if numOfServers is less than 1
     */
    public ServerType(String name, int numOfServers) {
        super(name);
        if (numOfServers < 1) {
            throw new IllegalArgumentException("Number of servers must be at least 1");
        }
        this.numOfServers = numOfServers;
        this.compatibleClasses = new ArrayList<JobClass>();
        this.id = -1; // Will be set when added to a queue
    }

    /**
     * Creates a new server type with the specified name, number of servers, and compatible classes.
     *
     * @param name the name identifying this server type
     * @param numOfServers the number of servers of this type (must be >= 1)
     * @param compatibleClasses initial list of compatible job classes
     * @throws IllegalArgumentException if numOfServers is less than 1
     */
    public ServerType(String name, int numOfServers, List<JobClass> compatibleClasses) {
        super(name);
        if (numOfServers < 1) {
            throw new IllegalArgumentException("Number of servers must be at least 1");
        }
        this.numOfServers = numOfServers;
        this.compatibleClasses = new ArrayList<JobClass>(compatibleClasses);
        this.id = -1;
    }

    /**
     * Gets the unique identifier of this server type within its queue.
     *
     * @return the server type ID, or -1 if not yet added to a queue
     */
    public int getId() {
        return id;
    }

    /**
     * Sets the unique identifier of this server type.
     * This is typically called by the Queue when the server type is added.
     *
     * @param id the ID to assign
     */
    public void setId(int id) {
        this.id = id;
    }

    /**
     * Gets the number of servers of this type.
     *
     * @return the number of servers
     */
    public int getNumOfServers() {
        return numOfServers;
    }

    /**
     * Sets the number of servers of this type.
     *
     * @param numOfServers the number of servers (must be >= 1)
     * @throws IllegalArgumentException if numOfServers is less than 1
     */
    public void setNumOfServers(int numOfServers) {
        if (numOfServers < 1) {
            throw new IllegalArgumentException("Number of servers must be at least 1");
        }
        this.numOfServers = numOfServers;
    }

    /**
     * Gets the list of job classes compatible with this server type.
     *
     * @return a new list containing the compatible job classes
     */
    public List<JobClass> getCompatibleClasses() {
        return new ArrayList<JobClass>(compatibleClasses);
    }

    /**
     * Adds a job class to the list of classes that can be served by this server type.
     *
     * @param jobClass the job class to add as compatible
     * @throws IllegalArgumentException if jobClass is null
     */
    public void addCompatibleClass(JobClass jobClass) {
        if (jobClass == null) {
            throw new IllegalArgumentException("Job class cannot be null");
        }
        if (!compatibleClasses.contains(jobClass)) {
            compatibleClasses.add(jobClass);
        }
    }

    /**
     * Removes a job class from the list of compatible classes.
     *
     * @param jobClass the job class to remove
     * @return true if the class was removed, false if it was not in the list
     */
    public boolean removeCompatibleClass(JobClass jobClass) {
        return compatibleClasses.remove(jobClass);
    }

    /**
     * Sets the list of compatible job classes, replacing any existing list.
     *
     * @param classes the list of job classes to set as compatible
     * @throws IllegalArgumentException if classes is null
     */
    public void setCompatibleClasses(List<JobClass> classes) {
        if (classes == null) {
            throw new IllegalArgumentException("Classes list cannot be null");
        }
        this.compatibleClasses = new ArrayList<JobClass>(classes);
    }

    /**
     * Checks if a job class is compatible with this server type.
     *
     * @param jobClass the job class to check
     * @return true if the class can be served by this server type
     */
    public boolean isCompatible(JobClass jobClass) {
        return compatibleClasses.contains(jobClass);
    }

    /**
     * Gets the number of compatible job classes.
     *
     * @return the count of compatible classes
     */
    public int getNumCompatibleClasses() {
        return compatibleClasses.size();
    }

    /**
     * Checks if this server type has any compatible classes defined.
     *
     * @return true if at least one compatible class is defined
     */
    public boolean hasCompatibleClasses() {
        return !compatibleClasses.isEmpty();
    }

    /**
     * Gets the parent queue this server type belongs to.
     *
     * @return the parent queue, or null if not yet added to a queue
     */
    public Queue getParentQueue() {
        return parentQueue;
    }

    /**
     * Sets the parent queue for this server type.
     * This is typically called by the Queue when the server type is added.
     *
     * @param queue the parent queue
     */
    public void setParentQueue(Queue queue) {
        this.parentQueue = queue;
    }

    @Override
    public String toString() {
        return "ServerType{" +
                "name='" + getName() + '\'' +
                ", id=" + id +
                ", numOfServers=" + numOfServers +
                ", compatibleClasses=" + compatibleClasses.size() +
                '}';
    }
}
