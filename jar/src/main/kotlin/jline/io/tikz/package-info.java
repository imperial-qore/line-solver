/*
 * Copyright (c) 2012-2026, QORE Lab, Imperial College London
 * All rights reserved.
 */

/**
 * TikZ visualization package for LINE queueing networks.
 *
 * <p>This package provides functionality to export queueing network models
 * as TikZ/LaTeX diagrams. The generated diagrams can be compiled to PDF
 * using pdflatex and displayed in a viewer.</p>
 *
 * <h2>Main Classes</h2>
 * <ul>
 *   <li>{@link jline.io.tikz.TikZExporter} - Main entry point for generating TikZ code</li>
 *   <li>{@link jline.io.tikz.TikZNodeRenderer} - Renders individual node types</li>
 *   <li>{@link jline.io.tikz.TikZLayoutEngine} - Computes automatic node positions</li>
 *   <li>{@link jline.io.tikz.TikZViewer} - Displays generated PDF diagrams</li>
 *   <li>{@link jline.io.tikz.TikZOptions} - Configuration options</li>
 * </ul>
 *
 * <h2>Usage Example</h2>
 * <pre>{@code
 * Network model = new Network("MyModel");
 * // ... add nodes and connections ...
 *
 * // Generate and display TikZ diagram
 * model.tikzView();
 *
 * // Or generate TikZ code
 * String tikz = model.toTikZ();
 *
 * // Or export to PDF file
 * File pdf = model.exportTikZ("mymodel");
 * }</pre>
 *
 * <h2>Node Shapes</h2>
 * <table>
 *   <tr><th>Node Type</th><th>Shape</th><th>Color</th></tr>
 *   <tr><td>Queue</td><td>Rectangle + server circle</td><td>Blue</td></tr>
 *   <tr><td>Delay</td><td>Dashed ellipse</td><td>Green</td></tr>
 *   <tr><td>Source</td><td>Circle</td><td>Yellow</td></tr>
 *   <tr><td>Sink</td><td>Circle</td><td>Red</td></tr>
 *   <tr><td>Fork</td><td>Diamond</td><td>Orange</td></tr>
 *   <tr><td>Join</td><td>Diamond</td><td>Purple</td></tr>
 *   <tr><td>Router</td><td>Hexagon</td><td>Cyan</td></tr>
 *   <tr><td>ClassSwitch</td><td>Trapezium</td><td>Pink</td></tr>
 *   <tr><td>Cache</td><td>Stacked rectangle</td><td>Gray</td></tr>
 * </table>
 *
 * @see jline.lang.Network#toTikZ()
 * @see jline.lang.Network#tikzView()
 */
package jline.io.tikz;
