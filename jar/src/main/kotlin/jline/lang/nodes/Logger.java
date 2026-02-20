/*
 * Copyright (c) 2012-2026, QORE Lab, Imperial College London
 * All rights reserved.
 */

package jline.lang.nodes;

import jline.lang.Network;
import jline.lang.sections.LogTunnel;

import java.io.File;
import java.io.Serializable;

/**
 * A node that logs passage of jobs
 */
public class Logger extends Node implements Serializable {

    public String fileName;
    public String filePath;
    protected Network model;
    boolean wantStartTime;
    boolean wantLoggerName;
    boolean wantTimestamp;
    boolean wantJobID;
    boolean wantJobClass;
    boolean wantTimeSameClass;
    boolean wantTimeAnyClass;

    public Logger(Network model, String name, String logfile) {
        super("Logger");
        this.setName(name);
        this.setModel(model);
        model.addNode(this);
        this.model = model;

        // Set up LogTunnel as the server section
        this.server = new LogTunnel();
        
        // Initialize input/output sections (required by JMT)
        // These are already initialized in parent Node constructor

        // Split logfile into path and name
        File file = new File(logfile);
        this.fileName = file.getName();
        this.filePath = file.getParent();
        if (this.filePath == null) {
            this.filePath = "";
        }
        if (!this.filePath.endsWith(File.separator) && !this.filePath.isEmpty()) {
            this.filePath += File.separator;
        }
    }

    public String getFileName() {
        return fileName;
    }

    public String getFilePath() {
        return filePath;
    }

    public boolean getJobClass() {
        return wantJobClass;
    }

    public void setJobClass(boolean value) {
        wantJobClass = value;
    }

    public boolean getJobID() {
        return wantJobID;
    }

    public void setJobID(boolean value) {
        wantJobID = value;
    }

    public boolean getLoggerName() {
        return wantLoggerName;
    }

    public void setLoggerName(boolean value) {
        wantLoggerName = value;
    }

    @Override
    public Network getModel() {
        return this.model;
    }

    public boolean getStartTime() {
        return wantStartTime;
    }

    public void setStartTime(boolean value) {
        wantStartTime = value;
    }

    public boolean getTimeAnyClass() {
        return wantTimeAnyClass;
    }

    public void setTimeAnyClass(boolean value) {
        wantTimeAnyClass = value;
    }

    public boolean getTimeSameClass() {
        return wantTimeSameClass;
    }

    public void setTimeSameClass(boolean value) {
        wantTimeSameClass = value;
    }

    public boolean getTimestamp() {
        return wantTimestamp;
    }

    public void setTimestamp(boolean value) {
        wantTimestamp = value;
    }

}
