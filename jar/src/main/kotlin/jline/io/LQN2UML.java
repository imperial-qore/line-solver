/*
 * Copyright (c) 2012-2026, QORE Lab, Imperial College London
 * All rights reserved.
 */
package jline.io;

import jline.io.tikz.SequenceDiagramExporter;
import jline.io.tikz.SequenceDiagramOptions;
import jline.lang.layered.LayeredNetwork;

import java.io.File;
import java.io.IOException;

/**
 * Converts a LayeredNetwork (LQN) model to a UML sequence diagram.
 *
 * <p>Default behavior exports to TEX (LaTeX/TikZ). Use options to specify other formats.
 *
 * <p>Example usage:
 * <pre>
 * // Export to TEX (default)
 * File tex = LQN2UML.export(model);
 *
 * // Export to PDF
 * LQN2UML.export(model, new LQN2UML.Options().setFormat(Format.PDF).setOutputPath("/tmp/diagram.pdf"));
 *
 * // Export to PNG
 * LQN2UML.export(model, new LQN2UML.Options().setFormat(Format.PNG).setOutputPath("/tmp/diagram.png"));
 * </pre>
 */
public class LQN2UML {

    /**
     * Output format for the sequence diagram.
     */
    public enum Format {
        TEX,    // LaTeX/TikZ source file (default)
        PDF,    // PDF file
        PNG     // PNG image
    }

    /**
     * Options for sequence diagram export.
     */
    public static class Options {
        private Format format = Format.TEX;
        private String outputPath = null;
        private int dpi = 150;
        private SequenceDiagramOptions diagramOptions = new SequenceDiagramOptions();

        public Options() {}

        public Format getFormat() {
            return format;
        }

        public Options setFormat(Format format) {
            this.format = format;
            return this;
        }

        public String getOutputPath() {
            return outputPath;
        }

        public Options setOutputPath(String outputPath) {
            this.outputPath = outputPath;
            return this;
        }

        public int getDpi() {
            return dpi;
        }

        public Options setDpi(int dpi) {
            this.dpi = dpi;
            return this;
        }

        public SequenceDiagramOptions getDiagramOptions() {
            return diagramOptions;
        }

        public Options setDiagramOptions(SequenceDiagramOptions diagramOptions) {
            this.diagramOptions = diagramOptions;
            return this;
        }

        // Convenience setters for common diagram options

        public Options setScale(double scale) {
            this.diagramOptions.setScale(scale);
            return this;
        }

        public Options setShowProcessorFrames(boolean show) {
            this.diagramOptions.setShowProcessorFrames(show);
            return this;
        }

        public Options setShowEntryNames(boolean show) {
            this.diagramOptions.setShowEntryNames(show);
            return this;
        }

        public Options setShowCallMeans(boolean show) {
            this.diagramOptions.setShowCallMeans(show);
            return this;
        }

        public Options setShowReplies(boolean show) {
            this.diagramOptions.setShowReplies(show);
            return this;
        }
    }

    /**
     * Exports a LayeredNetwork model to a UML sequence diagram (TEX format by default).
     *
     * @param model The LayeredNetwork model to convert
     * @return The generated TEX file
     * @throws IOException If export fails
     */
    public static File export(LayeredNetwork model) throws IOException {
        return export(model, new Options());
    }

    /**
     * Exports a LayeredNetwork model to a UML sequence diagram with custom options.
     *
     * @param model The LayeredNetwork model to convert
     * @param options Export options (format, output path, etc.)
     * @return The generated file (TEX, PDF, or PNG)
     * @throws IOException If export fails
     */
    public static File export(LayeredNetwork model, Options options) throws IOException {
        SequenceDiagramExporter exporter = new SequenceDiagramExporter(model, options.getDiagramOptions());

        switch (options.getFormat()) {
            case PDF:
                File pdfFile = exporter.exportToPDF();
                if (options.getOutputPath() != null) {
                    File destFile = new File(options.getOutputPath());
                    java.nio.file.Files.copy(pdfFile.toPath(), destFile.toPath(),
                            java.nio.file.StandardCopyOption.REPLACE_EXISTING);
                    return destFile;
                }
                return pdfFile;

            case PNG:
                String pngPath = options.getOutputPath();
                if (pngPath == null) {
                    pngPath = System.getProperty("java.io.tmpdir") + File.separator +
                              model.getName() + "_sequence.png";
                }
                exporter.exportToPNG(pngPath, options.getDpi());
                return new File(pngPath);

            case TEX:
            default:
                String texPath = options.getOutputPath();
                if (texPath == null) {
                    texPath = System.getProperty("java.io.tmpdir") + File.separator +
                              model.getName() + "_sequence.tex";
                }
                exporter.exportToFile(texPath);
                return new File(texPath);
        }
    }

    /**
     * Generates the TikZ/LaTeX source code for the sequence diagram.
     *
     * @param model The LayeredNetwork model to convert
     * @return The TikZ source code as a string
     */
    public static String toTikZ(LayeredNetwork model) {
        return toTikZ(model, new SequenceDiagramOptions());
    }

    /**
     * Generates the TikZ/LaTeX source code with custom diagram options.
     *
     * @param model The LayeredNetwork model to convert
     * @param diagramOptions Options for the diagram layout
     * @return The TikZ source code as a string
     */
    public static String toTikZ(LayeredNetwork model, SequenceDiagramOptions diagramOptions) {
        return new SequenceDiagramExporter(model, diagramOptions).generateSequenceDiagram();
    }

    /**
     * Displays the sequence diagram in a PDF viewer.
     *
     * @param model The LayeredNetwork model to display
     */
    public static void view(LayeredNetwork model) {
        view(model, new SequenceDiagramOptions());
    }

    /**
     * Displays the sequence diagram in a PDF viewer with custom options.
     *
     * @param model The LayeredNetwork model to display
     * @param diagramOptions Options for the diagram layout
     */
    public static void view(LayeredNetwork model, SequenceDiagramOptions diagramOptions) {
        new SequenceDiagramExporter(model, diagramOptions).display();
    }
}
