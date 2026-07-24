/*
 * Copyright (c) 2012-2026, QORE Lab, Imperial College London
 * All rights reserved.
 */

package jline.io.tikz;

/**
 * Configuration options for UML sequence diagram visualization of LayeredNetwork models.
 */
public class SequenceDiagramOptions {

    /** Horizontal spacing between lifelines (in cm) */
    private double lifelineSpacing = 3.0;

    /** Vertical spacing between messages (in cm) */
    private double messageSpacing = 1.0;

    /** Whether to show host demand as notes on activation boxes */
    private boolean showHostDemand = true;

    /** Whether to show call multiplicity on messages */
    private boolean showCallMeans = true;

    /** Whether to group tasks by processor in swim lanes */
    private boolean showProcessorFrames = true;

    /** Whether to show entry names on messages */
    private boolean showEntryNames = true;

    /** Whether to underline object names (UML style) */
    private boolean underlineObjectNames = true;

    /** Border padding around the diagram (in pt) */
    private int borderPadding = 10;

    /** Whether to use dashed arrows for async calls */
    private boolean asyncDashed = true;

    /** Whether to show reply messages explicitly */
    private boolean showReplies = true;

    /** Scale factor for the diagram (1.0 = normal size) */
    private double scale = 1.0;

    public SequenceDiagramOptions() {
    }

    // Getters and fluent setters

    public double getLifelineSpacing() {
        return lifelineSpacing;
    }

    public SequenceDiagramOptions setLifelineSpacing(double lifelineSpacing) {
        this.lifelineSpacing = lifelineSpacing;
        return this;
    }

    public double getMessageSpacing() {
        return messageSpacing;
    }

    public SequenceDiagramOptions setMessageSpacing(double messageSpacing) {
        this.messageSpacing = messageSpacing;
        return this;
    }

    public boolean isShowHostDemand() {
        return showHostDemand;
    }

    public SequenceDiagramOptions setShowHostDemand(boolean showHostDemand) {
        this.showHostDemand = showHostDemand;
        return this;
    }

    public boolean isShowCallMeans() {
        return showCallMeans;
    }

    public SequenceDiagramOptions setShowCallMeans(boolean showCallMeans) {
        this.showCallMeans = showCallMeans;
        return this;
    }

    public boolean isShowProcessorFrames() {
        return showProcessorFrames;
    }

    public SequenceDiagramOptions setShowProcessorFrames(boolean showProcessorFrames) {
        this.showProcessorFrames = showProcessorFrames;
        return this;
    }

    public boolean isShowEntryNames() {
        return showEntryNames;
    }

    public SequenceDiagramOptions setShowEntryNames(boolean showEntryNames) {
        this.showEntryNames = showEntryNames;
        return this;
    }

    public boolean isUnderlineObjectNames() {
        return underlineObjectNames;
    }

    public SequenceDiagramOptions setUnderlineObjectNames(boolean underlineObjectNames) {
        this.underlineObjectNames = underlineObjectNames;
        return this;
    }

    public int getBorderPadding() {
        return borderPadding;
    }

    public SequenceDiagramOptions setBorderPadding(int borderPadding) {
        this.borderPadding = borderPadding;
        return this;
    }

    public boolean isAsyncDashed() {
        return asyncDashed;
    }

    public SequenceDiagramOptions setAsyncDashed(boolean asyncDashed) {
        this.asyncDashed = asyncDashed;
        return this;
    }

    public boolean isShowReplies() {
        return showReplies;
    }

    public SequenceDiagramOptions setShowReplies(boolean showReplies) {
        this.showReplies = showReplies;
        return this;
    }

    public double getScale() {
        return scale;
    }

    public SequenceDiagramOptions setScale(double scale) {
        this.scale = scale;
        return this;
    }
}
