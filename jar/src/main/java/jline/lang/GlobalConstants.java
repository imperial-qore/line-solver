/*
 * Copyright (c) 2012-2026, QORE Lab, Imperial College London
 * All rights reserved.
 */

package jline.lang;

/**
 * Global constants and settings for LINE solver.
 * Maintains singleton state for configuration flags.
 */
public class GlobalConstants {

    // Library attribution flag
    private static boolean libraryAttributionShown = false;

    /**
     * Check if library attribution message has been shown.
     *
     * @return true if attribution has been shown, false otherwise
     */
    public static boolean isLibraryAttributionShown() {
        return libraryAttributionShown;
    }

    /**
     * Set the library attribution shown flag.
     *
     * @param value true if attribution message has been shown
     */
    public static void setLibraryAttributionShown(boolean value) {
        libraryAttributionShown = value;
    }
}
