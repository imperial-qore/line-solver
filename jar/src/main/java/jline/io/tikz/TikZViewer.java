/*
 * Copyright (c) 2012-2026, QORE Lab, Imperial College London
 * All rights reserved.
 */

package jline.io.tikz;

import java.awt.Desktop;
import java.io.File;
import java.io.IOException;

/**
 * Displays PDF files generated from TikZ diagrams using external viewers.
 */
public class TikZViewer {

    /**
     * Displays a PDF file using an external viewer.
     *
     * @param pdfFile The PDF file to display
     */
    public static void displayPDF(File pdfFile) {
        displayPDF(pdfFile, null);
    }

    /**
     * Displays a PDF file using an external viewer.
     *
     * @param pdfFile The PDF file to display
     * @param title   Optional title for logging
     */
    public static void displayPDF(File pdfFile, String title) {
        if (!pdfFile.exists()) {
            System.err.println("PDF file does not exist: " + pdfFile.getAbsolutePath());
            return;
        }

        // Prefer standalone PDF viewers over xdg-open (which may open in browser)
        String[] standaloneViewers = {"evince", "okular", "mupdf", "zathura", "qpdfview", "acroread"};
        for (String viewer : standaloneViewers) {
            if (isCommandAvailable(viewer)) {
                try {
                    ProcessBuilder pb = new ProcessBuilder(viewer, pdfFile.getAbsolutePath());
                    pb.start();
                    if (title != null) {
                        System.out.println("Opened diagram for: " + title);
                    }
                    return;
                } catch (IOException e) {
                    // Try next viewer
                }
            }
        }

        // Fallback to xdg-open or Desktop.open
        if (Desktop.isDesktopSupported()) {
            Desktop desktop = Desktop.getDesktop();
            if (desktop.isSupported(Desktop.Action.OPEN)) {
                try {
                    desktop.open(pdfFile);
                    if (title != null) {
                        System.out.println("Opened diagram for: " + title);
                    }
                    return;
                } catch (IOException e) {
                    System.err.println("Failed to open PDF with default viewer: " + e.getMessage());
                }
            }
        }

        // If all else fails, print the path
        System.out.println("Could not open PDF viewer. File saved at: " + pdfFile.getAbsolutePath());
    }

    /**
     * Checks if a command is available on the system.
     */
    private static boolean isCommandAvailable(String command) {
        try {
            ProcessBuilder pb = new ProcessBuilder("which", command);
            Process p = pb.start();
            int exitCode = p.waitFor();
            return exitCode == 0;
        } catch (Exception e) {
            return false;
        }
    }

    /**
     * Checks if a PDF viewer is available on the system.
     */
    public static boolean isPdfViewerAvailable() {
        String[] viewers = {"evince", "okular", "mupdf", "zathura", "qpdfview", "acroread", "xdg-open"};
        for (String viewer : viewers) {
            if (isCommandAvailable(viewer)) {
                return true;
            }
        }
        return Desktop.isDesktopSupported() && Desktop.getDesktop().isSupported(Desktop.Action.OPEN);
    }
}
