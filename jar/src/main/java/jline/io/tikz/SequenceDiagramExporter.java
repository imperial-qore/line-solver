/*
 * Copyright (c) 2012-2026, QORE Lab, Imperial College London
 * All rights reserved.
 */

package jline.io.tikz;

import jline.lang.layered.*;

import java.io.*;
import java.nio.charset.StandardCharsets;
import java.nio.file.Files;
import java.nio.file.Path;
import java.util.*;

/**
 * Exports a LayeredNetwork model as a UML sequence diagram using TikZ/pgf-umlsd.
 * Generates LaTeX code that can be compiled to PDF.
 */
public class SequenceDiagramExporter {

    private final LayeredNetwork model;
    private final SequenceDiagramOptions options;
    private final SequenceDiagramLayoutEngine layoutEngine;
    private final SequenceDiagramTraverser traverser;

    /**
     * Creates a sequence diagram exporter with default options.
     */
    public SequenceDiagramExporter(LayeredNetwork model) {
        this(model, new SequenceDiagramOptions());
    }

    /**
     * Creates a sequence diagram exporter with custom options.
     */
    public SequenceDiagramExporter(LayeredNetwork model, SequenceDiagramOptions options) {
        this.model = model;
        this.options = options;
        this.layoutEngine = new SequenceDiagramLayoutEngine(model, options);
        this.traverser = new SequenceDiagramTraverser(model, layoutEngine);
    }

    /**
     * Generates the complete LaTeX document with UML sequence diagram.
     *
     * @return Complete LaTeX document as a string
     */
    public String generateSequenceDiagram() {
        // Compute layout first (must be done before preamble generation)
        layoutEngine.computeLayout();

        // Extract interactions and fragments
        List<SequenceDiagramTraverser.Interaction> interactions = traverser.extractInteractions();
        List<SequenceDiagramTraverser.Fragment> fragments = traverser.extractFragments();

        StringBuilder sb = new StringBuilder();

        // Generate preamble (uses layout info for thread color definitions)
        sb.append(generatePreamble());
        sb.append("\n");

        // Document body
        sb.append("\\begin{document}\n");
        sb.append("\\begin{sequencediagram}\n");
        sb.append("\n");

        // Generate lifelines (grouped by processor)
        sb.append(generateLifelines());
        sb.append("\n");

        // Generate combined fragments first (they need to wrap interactions)
        // For now, we'll output fragments as comments/annotations
        if (!fragments.isEmpty()) {
            sb.append("% Combined Fragments\n");
            for (SequenceDiagramTraverser.Fragment frag : fragments) {
                sb.append(generateFragment(frag));
            }
            sb.append("\n");
        }

        // Generate interactions (messages)
        sb.append("% Interactions\n");
        sb.append(generateInteractions(interactions));
        sb.append("\n");

        sb.append("\\end{sequencediagram}\n");
        sb.append("\\end{document}\n");

        return sb.toString();
    }

    /**
     * Generates the LaTeX preamble with required packages and settings.
     */
    private String generatePreamble() {
        StringBuilder sb = new StringBuilder();

        sb.append("\\documentclass[tikz,border=").append(options.getBorderPadding()).append("pt]{standalone}\n");
        sb.append("\\usepackage{pgf-umlsd}\n");
        sb.append("\\usepgflibrary{arrows}\n");
        sb.append("\n");

        // Custom styling for processor frames
        sb.append("% Custom styles for processor grouping\n");
        sb.append("\\tikzset{\n");
        sb.append("  processor frame/.style={draw=blue!50, fill=blue!5, rounded corners=3pt, thick},\n");
        sb.append("  processor label/.style={font=\\bfseries\\small, fill=blue!20, rounded corners=2pt}\n");
        sb.append("}\n");
        sb.append("\n");

        // Define thread colors for all tasks (required by pgf-umlsd)
        sb.append("% Thread instance colors\n");
        for (Map.Entry<Host, List<Task>> entry : layoutEngine.getProcessorTaskGroups().entrySet()) {
            for (Task task : entry.getValue()) {
                String threadId = sanitizeId(task.getName());
                sb.append("\\tikzset{instcolor").append(threadId).append("/.style={}}\n");
            }
        }
        sb.append("\n");

        // Scale if needed
        if (options.getScale() != 1.0) {
            sb.append("\\tikzset{every picture/.style={scale=")
              .append(String.format("%.2f", options.getScale()))
              .append("}}\n");
        }

        return sb.toString();
    }

    /**
     * Generates lifelines grouped by processor.
     */
    private String generateLifelines() {
        StringBuilder sb = new StringBuilder();
        Map<Host, List<Task>> groups = layoutEngine.getProcessorTaskGroups();

        if (options.isShowProcessorFrames()) {
            // Generate processor frames with nested tasks
            for (Map.Entry<Host, List<Task>> entry : groups.entrySet()) {
                Host processor = entry.getKey();
                List<Task> tasks = entry.getValue();

                if (tasks.isEmpty()) continue;

                sb.append("% Processor: ").append(sanitizeName(processor.getName())).append("\n");
                sb.append("\\begin{sdblock}{").append(escapeLatex(processor.getName())).append("}{}\n");

                for (Task task : tasks) {
                    String threadId = sanitizeId(task.getName());
                    String threadLabel = formatTaskLabel(task);
                    sb.append("  \\newthread{").append(threadId).append("}{")
                      .append(threadLabel).append("}\n");
                }

                sb.append("\\end{sdblock}\n");
                sb.append("\n");
            }
        } else {
            // Generate flat lifelines without processor grouping
            for (Map.Entry<Host, List<Task>> entry : groups.entrySet()) {
                for (Task task : entry.getValue()) {
                    String threadId = sanitizeId(task.getName());
                    String threadLabel = formatTaskLabel(task);
                    sb.append("\\newthread{").append(threadId).append("}{")
                      .append(threadLabel).append("}\n");
                }
            }
        }

        return sb.toString();
    }

    /**
     * Formats a task label for the lifeline.
     */
    private String formatTaskLabel(Task task) {
        StringBuilder label = new StringBuilder();

        if (options.isUnderlineObjectNames()) {
            label.append("\\underline{").append(escapeLatex(task.getName())).append("}");
        } else {
            label.append(escapeLatex(task.getName()));
        }

        label.append(":Task");

        return label.toString();
    }

    /**
     * Generates interaction messages.
     */
    private String generateInteractions(List<SequenceDiagramTraverser.Interaction> interactions) {
        StringBuilder sb = new StringBuilder();

        // Group synchronous calls with their potential replies
        Map<String, List<SequenceDiagramTraverser.Interaction>> callsBySource = new LinkedHashMap<>();
        for (SequenceDiagramTraverser.Interaction interaction : interactions) {
            String key = sanitizeId(interaction.sourceTask.getName());
            if (!callsBySource.containsKey(key)) {
                callsBySource.put(key, new ArrayList<>());
            }
            callsBySource.get(key).add(interaction);
        }

        // Generate calls
        for (SequenceDiagramTraverser.Interaction interaction : interactions) {
            String sourceId = sanitizeId(interaction.sourceTask.getName());
            String targetId = sanitizeId(interaction.targetTask.getName());

            // Skip self-calls for now (they need special handling)
            if (sourceId.equals(targetId)) {
                sb.append("% Self-call from ").append(sourceId).append("\n");
                sb.append("\\begin{callself}{").append(sourceId).append("}{")
                  .append(formatMessageLabel(interaction)).append("}\n");
                sb.append("\\end{callself}\n");
                continue;
            }

            if (interaction.type == SequenceDiagramTraverser.Interaction.Type.SYNCH_CALL) {
                // Synchronous call with reply
                sb.append("\\begin{call}{").append(sourceId).append("}{")
                  .append(formatMessageLabel(interaction)).append("}{")
                  .append(targetId).append("}{");

                if (options.isShowReplies()) {
                    sb.append("reply");
                }
                sb.append("}\n");
                sb.append("\\end{call}\n");
            } else if (interaction.type == SequenceDiagramTraverser.Interaction.Type.ASYNCH_CALL) {
                // Asynchronous message (no blocking wait)
                if (options.isAsyncDashed()) {
                    sb.append("\\mess[1]{").append(sourceId).append("}{")
                      .append(formatMessageLabel(interaction)).append("}{")
                      .append(targetId).append("}\n");
                } else {
                    sb.append("\\mess{").append(sourceId).append("}{")
                      .append(formatMessageLabel(interaction)).append("}{")
                      .append(targetId).append("}\n");
                }
            }
        }

        return sb.toString();
    }

    /**
     * Formats a message label for an interaction.
     */
    private String formatMessageLabel(SequenceDiagramTraverser.Interaction interaction) {
        StringBuilder label = new StringBuilder();

        if (options.isShowEntryNames() && interaction.targetEntry != null) {
            label.append(escapeLatex(interaction.targetEntry.getName()));
        } else {
            label.append("call");
        }

        if (options.isShowCallMeans() && interaction.callMean > 1.0) {
            label.append(" [").append(String.format("%.1f", interaction.callMean)).append("]");
        }

        return label.toString();
    }

    /**
     * Generates a combined fragment (loop, par, alt).
     */
    private String generateFragment(SequenceDiagramTraverser.Fragment fragment) {
        StringBuilder sb = new StringBuilder();

        switch (fragment.type) {
            case LOOP:
                sb.append("\\begin{sdloop}{").append(escapeLatex(fragment.label)).append("}\n");
                sb.append("  % Loop body activities at levels ")
                  .append(fragment.startLevel).append("-").append(fragment.endLevel).append("\n");
                sb.append("\\end{sdloop}\n");
                break;

            case PAR:
                sb.append("\\begin{sdblock}{").append(escapeLatex(fragment.label)).append("}{}\n");
                sb.append("  % Parallel activities at levels ")
                  .append(fragment.startLevel).append("-").append(fragment.endLevel).append("\n");
                sb.append("\\end{sdblock}\n");
                break;

            case ALT:
                sb.append("\\begin{sdblock}{").append(escapeLatex(fragment.label)).append("}{}\n");
                for (SequenceDiagramTraverser.Fragment operand : fragment.operands) {
                    sb.append("  % Alternative: ").append(operand.label).append("\n");
                }
                sb.append("\\end{sdblock}\n");
                break;

            default:
                sb.append("% Unknown fragment type\n");
        }

        return sb.toString();
    }

    /**
     * Sanitizes a name for use as a TikZ identifier.
     */
    private String sanitizeId(String name) {
        return name.replaceAll("[^a-zA-Z0-9]", "");
    }

    /**
     * Sanitizes a name for comments.
     */
    private String sanitizeName(String name) {
        return name.replaceAll("[%\\\\]", "");
    }

    /**
     * Escapes special LaTeX characters.
     */
    private String escapeLatex(String text) {
        if (text == null) return "";
        return text
                .replace("\\", "\\textbackslash{}")
                .replace("_", "\\_")
                .replace("&", "\\&")
                .replace("%", "\\%")
                .replace("$", "\\$")
                .replace("#", "\\#")
                .replace("{", "\\{")
                .replace("}", "\\}")
                .replace("~", "\\textasciitilde{}")
                .replace("^", "\\textasciicircum{}");
    }

    /**
     * Exports the sequence diagram to a .tex file.
     *
     * @param filePath Path to the output .tex file
     */
    public void exportToFile(String filePath) throws IOException {
        String tikzCode = generateSequenceDiagram();
        Path path = new File(filePath).toPath();
        Files.write(path, tikzCode.getBytes(StandardCharsets.UTF_8));
    }

    /**
     * Compiles the sequence diagram to PDF using pdflatex.
     *
     * @return The generated PDF file
     * @throws IOException If file operations fail
     * @throws RuntimeException If pdflatex fails
     */
    public File exportToPDF() throws IOException {
        if (!TikZExporter.isPdfLatexAvailable()) {
            throw new RuntimeException("pdflatex is not available. Please install TeX Live or MiKTeX.");
        }

        // Create temporary directory
        Path tempDir = Files.createTempDirectory("line-seqdiag-");
        Path texFile = tempDir.resolve("sequence.tex");
        Path pdfFile = tempDir.resolve("sequence.pdf");

        // Write .tex file
        String tikzCode = generateSequenceDiagram();
        Files.write(texFile, tikzCode.getBytes(StandardCharsets.UTF_8));

        // Compile with pdflatex
        ProcessBuilder pb = new ProcessBuilder(
                "pdflatex",
                "-interaction=nonstopmode",
                "-output-directory=" + tempDir.toString(),
                texFile.toString()
        );
        pb.redirectErrorStream(true);
        pb.directory(tempDir.toFile());

        Process process = pb.start();

        // Read output
        BufferedReader reader = new BufferedReader(new InputStreamReader(process.getInputStream()));
        StringBuilder output = new StringBuilder();
        String line;
        while ((line = reader.readLine()) != null) {
            output.append(line).append("\n");
        }

        try {
            process.waitFor();
            // Only fail if PDF was not generated - pdflatex may return non-zero on warnings
            if (!Files.exists(pdfFile)) {
                throw new RuntimeException("pdflatex compilation failed - no PDF generated:\n" + output.toString());
            }
        } catch (InterruptedException e) {
            Thread.currentThread().interrupt();
            throw new RuntimeException("pdflatex was interrupted", e);
        }

        return pdfFile.toFile();
    }

    /**
     * Exports the sequence diagram to a PNG file.
     *
     * @param filePath Path to the output .png file
     * @param dpi Resolution in dots per inch
     */
    public void exportToPNG(String filePath, int dpi) throws IOException {
        // First generate PDF
        File pdfFile = exportToPDF();

        // Convert PDF to PNG using pdftoppm
        String pngPathWithoutExt = filePath.replaceAll("\\.png$", "");
        ProcessBuilder pb = new ProcessBuilder(
                "pdftoppm",
                "-png",
                "-r", String.valueOf(dpi),
                "-singlefile",
                pdfFile.getAbsolutePath(),
                pngPathWithoutExt
        );
        pb.redirectErrorStream(true);

        Process process = pb.start();
        BufferedReader reader = new BufferedReader(new InputStreamReader(process.getInputStream()));
        StringBuilder output = new StringBuilder();
        String line;
        while ((line = reader.readLine()) != null) {
            output.append(line).append("\n");
        }

        try {
            int exitCode = process.waitFor();
            if (exitCode != 0) {
                throw new RuntimeException("pdftoppm conversion failed:\n" + output.toString());
            }
        } catch (InterruptedException e) {
            Thread.currentThread().interrupt();
            throw new RuntimeException("pdftoppm was interrupted", e);
        }

        System.out.println("PNG exported to: " + filePath);
    }

    /**
     * Exports the sequence diagram to a PNG file with default 150 DPI.
     */
    public void exportToPNG(String filePath) throws IOException {
        exportToPNG(filePath, 150);
    }

    /**
     * Displays the sequence diagram in a PDF viewer.
     */
    public void display() {
        if (!TikZExporter.isPdfLatexAvailable()) {
            String texPath = "sequence-diagram.tex";
            try {
                exportToFile(texPath);
                System.out.println("pdflatex not found. TikZ code saved to: " + texPath);
                System.out.println("Compile manually with: pdflatex " + texPath);
            } catch (IOException e) {
                System.err.println("Failed to save TikZ file: " + e.getMessage());
            }
            return;
        }

        try {
            File pdf = exportToPDF();
            TikZViewer.displayPDF(pdf, model.getName());
        } catch (Exception e) {
            System.err.println("Failed to generate visualization: " + e.getMessage());
            try {
                String texPath = "sequence-diagram.tex";
                exportToFile(texPath);
                System.out.println("TikZ code saved to: " + texPath);
            } catch (IOException ex) {
                System.err.println("Failed to save TikZ file: " + ex.getMessage());
            }
        }
    }

    /**
     * Gets the options used by this exporter.
     */
    public SequenceDiagramOptions getOptions() {
        return options;
    }
}
