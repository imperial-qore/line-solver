/*
 * Copyright (c) 2012-2026, QORE Lab, Imperial College London
 * All rights reserved.
 */

package jline.io;

import jline.lang.Network;
import jline.lang.JobClass;
import jline.lang.layered.*;
import jline.lang.nodes.*;

import edu.uci.ics.jung.graph.DirectedSparseGraph;
import edu.uci.ics.jung.graph.Graph;
import edu.uci.ics.jung.algorithms.layout.Layout;
import edu.uci.ics.jung.algorithms.layout.StaticLayout;
import edu.uci.ics.jung.visualization.VisualizationViewer;
import edu.uci.ics.jung.visualization.control.DefaultModalGraphMouse;
import edu.uci.ics.jung.visualization.control.ModalGraphMouse;
import edu.uci.ics.jung.visualization.decorators.EdgeShape;
import edu.uci.ics.jung.visualization.renderers.Renderer;

import com.google.common.base.Function;

import javax.swing.Icon;
import javax.swing.JButton;
import javax.swing.JDialog;
import javax.swing.JFileChooser;
import javax.swing.JFrame;
import javax.swing.JLabel;
import javax.swing.JOptionPane;
import javax.swing.JPanel;
import javax.swing.JScrollPane;
import javax.swing.JSlider;
import javax.swing.JTextArea;
import javax.swing.JTextField;
import javax.swing.JToolBar;
import javax.swing.BorderFactory;
import javax.swing.ButtonGroup;
import javax.swing.JToggleButton;
import javax.swing.SwingUtilities;
import javax.swing.UIManager;
import javax.swing.event.ChangeEvent;
import javax.swing.event.ChangeListener;
import javax.swing.BoxLayout;
import javax.swing.filechooser.FileNameExtensionFilter;
import javax.swing.JComboBox;
import javax.swing.JSpinner;
import javax.swing.SpinnerNumberModel;
import javax.swing.JSeparator;
import javax.swing.JTabbedPane;
import javax.swing.JTable;
import javax.swing.SwingWorker;
import javax.swing.table.DefaultTableModel;
import javax.swing.table.DefaultTableCellRenderer;

import com.formdev.flatlaf.FlatLightLaf;
import javax.imageio.ImageIO;
import java.awt.Component;
import java.awt.RenderingHints;
import java.awt.geom.AffineTransform;
import java.awt.geom.GeneralPath;
import java.awt.Graphics;
import java.awt.Graphics2D;
import java.awt.image.BufferedImage;
import java.io.File;
import java.io.IOException;
import java.net.URL;
import java.awt.BasicStroke;
import java.awt.BorderLayout;
import java.awt.Color;
import java.awt.Dimension;
import java.awt.Font;
import java.awt.Paint;
import java.awt.Polygon;
import java.awt.Shape;
import java.awt.Stroke;
import java.awt.event.ActionEvent;
import java.awt.event.ActionListener;
import java.awt.event.ComponentAdapter;
import java.awt.event.ComponentEvent;
import java.awt.event.MouseAdapter;
import java.awt.event.MouseEvent;
import java.awt.geom.Ellipse2D;
import java.awt.geom.Point2D;
import java.awt.geom.Rectangle2D;
import java.util.ArrayList;
import java.util.Collections;
import java.util.Comparator;
import java.util.HashMap;
import java.util.HashSet;
import java.util.Iterator;
import java.util.List;
import java.util.Map;
import java.util.Set;

import jline.lang.constant.SchedStrategy;
import jline.solvers.SolverResult;
import jline.solvers.LayeredNetworkAvgTable;
import jline.solvers.ln.SolverLN;
import jline.util.matrix.Matrix;

/**
 * Unified visualizer for both Network and LayeredNetwork models using the JUNG library.
 *
 * <p>This visualizer accepts either a Network or LayeredNetwork object and provides
 * appropriate visual styling for each model type:
 *
 * <h3>Network models (queueing networks):</h3>
 * <ul>
 * <li>White shapes (except Sink which is black)</li>
 * <li>9 node types: Source, Sink, Queue, Delay, Fork, Join, ClassSwitch, Cache, Router</li>
 * <li>Sugiyama-style left-to-right hierarchical layout</li>
 * <li>Solid black edges</li>
 * </ul>
 *
 * <h3>LayeredNetwork models (LQN):</h3>
 * <ul>
 * <li>Red shapes (yellow when selected)</li>
 * <li>4 node types: Host (triangle), Task (parallelogram), Entry (rectangle), Activity (circle)</li>
 * <li>4-level hierarchical layout (Host→Task→Entry→Activity)</li>
 * <li>6 edge styles: solid (hierarchy), dashed (sync), dot-dashed (async), dotted (forwarding)</li>
 * </ul>
 */
public class ModelVisualizer {

    // Static initializer to set up FlatLaf modern Look and Feel
    static {
        try {
            FlatLightLaf.setup();
            Font defaultFont = new Font("Segoe UI", Font.PLAIN, 13);
            Font smallFont = new Font("Segoe UI", Font.PLAIN, 11);
            Font titleFont = new Font("Segoe UI", Font.BOLD, 14);

            UIManager.put("defaultFont", defaultFont);
            UIManager.put("Button.font", defaultFont);
            UIManager.put("Label.font", defaultFont);
            UIManager.put("TextField.font", defaultFont);
            UIManager.put("TextArea.font", new Font("JetBrains Mono", Font.PLAIN, 12));
            UIManager.put("TitledBorder.font", titleFont);
            UIManager.put("ToolTip.font", smallFont);
            UIManager.put("Slider.font", smallFont);
            UIManager.put("OptionPane.messageFont", defaultFont);
            UIManager.put("OptionPane.buttonFont", defaultFont);

            UIManager.put("Button.arc", 8);
            UIManager.put("Component.arc", 8);
            UIManager.put("TextComponent.arc", 8);
            UIManager.put("ScrollBar.width", 12);
            UIManager.put("TitlePane.unifiedBackground", true);

        } catch (Exception e) {
            System.err.println("Warning: Could not initialize FlatLaf theme: " + e.getMessage());
            try {
                UIManager.setLookAndFeel(UIManager.getSystemLookAndFeelClassName());
            } catch (Exception ex) {
                // Ignore
            }
        }
    }

    /** Horizontal spacing between layers (in pixels) for Network models */
    private static final double LAYER_SPACING = 120.0;

    /** Vertical spacing between nodes within the same layer (in pixels) for Network models */
    private static final double NODE_SPACING = 80.0;

    /**
     * Unified node type enumeration covering both Network and LayeredNetwork elements.
     */
    public enum NodeType {
        // Network types (9)
        SOURCE,      // Circle (white)
        SINK,        // Circle (black)
        QUEUE,       // Rectangle
        DELAY,       // Circle
        FORK,        // Diamond
        JOIN,        // Diamond
        CLASSSWITCH, // Hexagon
        CACHE,       // Rectangle
        ROUTER,      // Octagon
        // LayeredNetwork types (4)
        HOST,        // Triangle (pyramid)
        TASK,        // Parallelogram
        ENTRY,       // Rectangle
        ACTIVITY,    // Circle
        // Default
        OTHER        // Default circle
    }

    /**
     * Unified edge type enumeration covering both Network and LayeredNetwork relationships.
     */
    public enum EdgeType {
        ROUTING,        // Network: solid black
        PARENT_CHILD,   // LQN: solid (Host→Task, Task→Entry)
        BOUND_TO,       // LQN: solid (Entry→Activity binding)
        PRECEDENCE,     // LQN: solid (Activity→Activity precedence)
        SYNCH_CALL,     // LQN: dashed (Activity→Entry sync call)
        ASYNCH_CALL,    // LQN: dot-dashed (Activity→Entry async call)
        FORWARDING      // LQN: dotted (Entry→Entry forwarding)
    }

    /**
     * Unified wrapper class for graph vertices representing model elements.
     */
    public static class ModelVertex {
        private final String name;
        private final NodeType type;
        private final int index;
        private final Object element;  // Node or LayeredNetworkElement

        public ModelVertex(String name, NodeType type, int index, Object element) {
            this.name = name;
            this.type = type;
            this.index = index;
            this.element = element;
        }

        public String getName() { return name; }
        public NodeType getType() { return type; }
        public int getIndex() { return index; }
        public Object getElement() { return element; }

        /**
         * Returns the short name for this vertex based on its type.
         */
        public String getShortName() {
            String prefix;
            switch (type) {
                // Network types
                case SOURCE: prefix = "Src"; break;
                case SINK: prefix = "Snk"; break;
                case QUEUE: prefix = "Q"; break;
                case DELAY: prefix = "D"; break;
                case FORK: prefix = "F"; break;
                case JOIN: prefix = "J"; break;
                case CLASSSWITCH: prefix = "CS"; break;
                case CACHE: prefix = "C"; break;
                case ROUTER: prefix = "R"; break;
                // LQN types
                case HOST: prefix = "P"; break;
                case TASK: prefix = "T"; break;
                case ENTRY: prefix = "E"; break;
                case ACTIVITY: prefix = "A"; break;
                default: prefix = "N"; break;
            }
            return prefix + index;
        }

        @Override
        public String toString() { return name; }

        @Override
        public boolean equals(Object o) {
            if (this == o) return true;
            if (o == null || getClass() != o.getClass()) return false;
            ModelVertex that = (ModelVertex) o;
            return name.equals(that.name) && type == that.type;
        }

        @Override
        public int hashCode() {
            return 31 * name.hashCode() + type.hashCode();
        }
    }

    /**
     * Unified wrapper class for graph edges representing relationships.
     */
    public static class ModelEdge {
        private final String label;
        private final EdgeType type;
        private static int edgeCounter = 0;
        private final int id;

        public ModelEdge(String label, EdgeType type) {
            this.label = label;
            this.type = type;
            this.id = edgeCounter++;
        }

        public String getLabel() { return label; }
        public EdgeType getType() { return type; }

        @Override
        public String toString() { return label; }

        @Override
        public boolean equals(Object o) {
            if (this == o) return true;
            if (o == null || getClass() != o.getClass()) return false;
            ModelEdge that = (ModelEdge) o;
            return id == that.id;
        }

        @Override
        public int hashCode() { return id; }
    }

    // Model references (only one will be non-null)
    private final Network networkModel;
    private final LayeredNetwork layeredModel;
    private final boolean isLayeredNetwork;

    private Graph<ModelVertex, ModelEdge> graph;
    private Map<String, ModelVertex> vertexMap;
    private List<List<ModelVertex>> layers;  // For Network Sugiyama layout

    /**
     * Creates a visualizer for a Network model.
     *
     * @param model the Network to visualize
     */
    public ModelVisualizer(Network model) {
        this.networkModel = model;
        this.layeredModel = null;
        this.isLayeredNetwork = false;
        this.vertexMap = new HashMap<String, ModelVertex>();
        this.layers = new ArrayList<List<ModelVertex>>();
    }

    /**
     * Creates a visualizer for a LayeredNetwork model.
     *
     * @param model the LayeredNetwork to visualize
     */
    public ModelVisualizer(LayeredNetwork model) {
        this.networkModel = null;
        this.layeredModel = model;
        this.isLayeredNetwork = true;
        this.vertexMap = new HashMap<String, ModelVertex>();
        this.layers = new ArrayList<List<ModelVertex>>();
    }

    /**
     * Builds a JUNG graph from the underlying model.
     *
     * @return the constructed JUNG graph
     */
    public Graph<ModelVertex, ModelEdge> buildGraph() {
        if (isLayeredNetwork) {
            return buildLayeredNetworkGraph();
        } else {
            return buildNetworkGraph();
        }
    }

    // ========================================================================
    // NETWORK GRAPH BUILDING
    // ========================================================================

    /**
     * Determines the node type for a Network node.
     */
    private NodeType getNetworkNodeType(Node node) {
        if (node instanceof Source) return NodeType.SOURCE;
        if (node instanceof Sink) return NodeType.SINK;
        if (node instanceof Queue) return NodeType.QUEUE;
        if (node instanceof Delay) return NodeType.DELAY;
        if (node instanceof Fork) return NodeType.FORK;
        if (node instanceof Join) return NodeType.JOIN;
        if (node instanceof ClassSwitch) return NodeType.CLASSSWITCH;
        if (node instanceof Cache) return NodeType.CACHE;
        if (node instanceof Router) return NodeType.ROUTER;
        return NodeType.OTHER;
    }

    /**
     * Builds a JUNG graph from the Network model.
     */
    private Graph<ModelVertex, ModelEdge> buildNetworkGraph() {
        graph = new DirectedSparseGraph<ModelVertex, ModelEdge>();
        vertexMap.clear();

        List<Node> nodes = networkModel.getNodes();

        // Add all nodes as vertices
        for (int i = 0; i < nodes.size(); i++) {
            Node node = nodes.get(i);
            NodeType type = getNetworkNodeType(node);
            ModelVertex vertex = new ModelVertex(node.getName(), type, i, node);
            graph.addVertex(vertex);
            vertexMap.put("N:" + node.getName(), vertex);
        }

        // Add edges from connection matrix
        Matrix connectionMatrix = networkModel.getConnectionMatrix();
        if (connectionMatrix != null) {
            int numNodes = nodes.size();
            for (int i = 0; i < numNodes; i++) {
                for (int j = 0; j < numNodes; j++) {
                    if (connectionMatrix.get(i, j) > 0) {
                        ModelVertex source = vertexMap.get("N:" + nodes.get(i).getName());
                        ModelVertex dest = vertexMap.get("N:" + nodes.get(j).getName());
                        if (source != null && dest != null) {
                            graph.addEdge(new ModelEdge("", EdgeType.ROUTING), source, dest);
                        }
                    }
                }
            }
        }

        return graph;
    }

    // ========================================================================
    // LAYERED NETWORK GRAPH BUILDING
    // ========================================================================

    /**
     * Builds a JUNG graph from the LayeredNetwork model.
     */
    private Graph<ModelVertex, ModelEdge> buildLayeredNetworkGraph() {
        graph = new DirectedSparseGraph<ModelVertex, ModelEdge>();
        vertexMap.clear();

        // Build a map from activity name to Activity object for lookups
        Map<String, Activity> activityByName = new HashMap<String, Activity>();
        for (Map.Entry<Integer, Activity> entry : layeredModel.getActivities().entrySet()) {
            activityByName.put(entry.getValue().getName(), entry.getValue());
        }

        // 1. Add all Hosts as vertices
        for (Map.Entry<Integer, Host> entry : layeredModel.getHosts().entrySet()) {
            int idx = entry.getKey();
            Host host = entry.getValue();
            ModelVertex vertex = new ModelVertex(host.getName(), NodeType.HOST, idx, host);
            graph.addVertex(vertex);
            vertexMap.put("H:" + host.getName(), vertex);
        }

        // 2. Add all Tasks as vertices with edges to parent Hosts
        for (Map.Entry<Integer, Task> entry : layeredModel.getTasks().entrySet()) {
            int idx = entry.getKey();
            Task task = entry.getValue();
            ModelVertex vertex = new ModelVertex(task.getName(), NodeType.TASK, idx, task);
            graph.addVertex(vertex);
            vertexMap.put("T:" + task.getName(), vertex);

            if (task.getParent() != null) {
                ModelVertex parentVertex = vertexMap.get("H:" + task.getParent().getName());
                if (parentVertex != null) {
                    graph.addEdge(new ModelEdge("", EdgeType.PARENT_CHILD), parentVertex, vertex);
                }
            }
        }

        // 3. Add all Entries as vertices with edges to parent Tasks
        for (Map.Entry<Integer, Entry> entry : layeredModel.getEntries().entrySet()) {
            int idx = entry.getKey();
            Entry lqnEntry = entry.getValue();
            ModelVertex vertex = new ModelVertex(lqnEntry.getName(), NodeType.ENTRY, idx, lqnEntry);
            graph.addVertex(vertex);
            vertexMap.put("E:" + lqnEntry.getName(), vertex);

            if (lqnEntry.getParent() != null) {
                ModelVertex parentVertex = vertexMap.get("T:" + lqnEntry.getParent().getName());
                if (parentVertex != null) {
                    graph.addEdge(new ModelEdge("", EdgeType.PARENT_CHILD), parentVertex, vertex);
                }
            }
        }

        // 3b. Add forwarding edges (Entry → Entry)
        for (Map.Entry<Integer, Entry> entry : layeredModel.getEntries().entrySet()) {
            Entry lqnEntry = entry.getValue();
            ModelVertex sourceVertex = vertexMap.get("E:" + lqnEntry.getName());
            if (sourceVertex == null) continue;

            for (Map.Entry<Integer, String> fwdDest : lqnEntry.getForwardingDests().entrySet()) {
                String destEntryName = fwdDest.getValue();
                ModelVertex destVertex = vertexMap.get("E:" + destEntryName);
                if (destVertex != null) {
                    graph.addEdge(new ModelEdge("fwd", EdgeType.FORWARDING), sourceVertex, destVertex);
                }
            }
        }

        // 4. Add all Activities as vertices
        for (Map.Entry<Integer, Activity> entry : layeredModel.getActivities().entrySet()) {
            int idx = entry.getKey();
            Activity activity = entry.getValue();
            ModelVertex vertex = new ModelVertex(activity.getName(), NodeType.ACTIVITY, idx, activity);
            graph.addVertex(vertex);
            vertexMap.put("A:" + activity.getName(), vertex);

            // Edge from Entry to Activity (bound-to relationship)
            if (activity.getBoundToEntry() != null && !activity.getBoundToEntry().isEmpty()) {
                ModelVertex entryVertex = vertexMap.get("E:" + activity.getBoundToEntry());
                if (entryVertex != null) {
                    graph.addEdge(new ModelEdge("bound", EdgeType.BOUND_TO), entryVertex, vertex);
                }
            }

            // Synch call edges (Activity → Entry)
            for (Map.Entry<Integer, String> syncCall : activity.getSyncCallDests().entrySet()) {
                String destEntryName = syncCall.getValue();
                ModelVertex destVertex = vertexMap.get("E:" + destEntryName);
                if (destVertex != null) {
                    graph.addEdge(new ModelEdge("synch", EdgeType.SYNCH_CALL), vertex, destVertex);
                }
            }

            // Asynch call edges (Activity → Entry)
            for (Map.Entry<Integer, String> asyncCall : activity.getAsyncCallDests().entrySet()) {
                String destEntryName = asyncCall.getValue();
                ModelVertex destVertex = vertexMap.get("E:" + destEntryName);
                if (destVertex != null) {
                    graph.addEdge(new ModelEdge("asynch", EdgeType.ASYNCH_CALL), vertex, destVertex);
                }
            }
        }

        // 5. Add Activity→Activity precedence edges from Task precedences
        for (Map.Entry<Integer, Task> entry : layeredModel.getTasks().entrySet()) {
            Task task = entry.getValue();
            for (ActivityPrecedence prec : task.getPrecedences()) {
                for (String preActName : prec.getPreActs()) {
                    ModelVertex preVertex = vertexMap.get("A:" + preActName);
                    if (preVertex == null) continue;

                    for (String postActName : prec.getPostActs()) {
                        ModelVertex postVertex = vertexMap.get("A:" + postActName);
                        if (postVertex == null) continue;

                        graph.addEdge(new ModelEdge("", EdgeType.PRECEDENCE), preVertex, postVertex);
                    }
                }
            }
        }

        return graph;
    }

    // ========================================================================
    // SHOW METHODS
    // ========================================================================

    /**
     * Displays the graph in a JFrame window with default settings.
     */
    public void show() {
        show(getDefaultTitle(), 800, 600);
    }

    /**
     * Displays the graph with a custom title.
     *
     * @param title the window title
     */
    public void show(String title) {
        show(title, 800, 600);
    }

    /**
     * Displays the graph with custom title and dimensions.
     *
     * @param title  the window title
     * @param width  the window width
     * @param height the window height
     */
    public void show(String title, int width, int height) {
        show(title, width, height, true);
    }

    /**
     * Displays the graph with custom title, dimensions, and toolbar visibility.
     *
     * @param title        the window title
     * @param width        the window width
     * @param height       the window height
     * @param showToolbars if true, show toolbars; if false, show only the graph image
     */
    public void show(String title, int width, int height, boolean showToolbars) {
        if (graph == null) {
            buildGraph();
        }

        Layout<ModelVertex, ModelEdge> layout = createHierarchicalLayout(width, height);
        showWithLayout(layout, title, width, height, showToolbars);
    }

    private String getDefaultTitle() {
        if (isLayeredNetwork) {
            return "LayeredNetwork: " + layeredModel.getName();
        } else {
            return "Network: " + networkModel.getName();
        }
    }

    /**
     * Displays the graph using the specified layout with optional toolbars.
     */
    public void showWithLayout(Layout<ModelVertex, ModelEdge> layout, String title,
                               final int width, final int height, boolean showToolbars) {
        layout.setSize(new Dimension(width - 50, height - 50));

        final VisualizationViewer<ModelVertex, ModelEdge> vv =
            new VisualizationViewer<ModelVertex, ModelEdge>(layout);
        vv.setPreferredSize(new Dimension(width, height));

        // Set vertex label position
        vv.getRenderer().getVertexLabelRenderer().setPosition(
            isLayeredNetwork ? Renderer.VertexLabel.Position.N : Renderer.VertexLabel.Position.S);

        // Track selected vertex for highlighting (LQN only)
        final ModelVertex[] selectedVertex = {null};

        // Color theme: 0=Color, 1=Greyscale, 2=Black/White
        final int[] colorTheme = {0};

        // Apply vertex transformers with selection-aware paint and color theme
        vv.getRenderContext().setVertexFillPaintTransformer(getVertexPaintTransformer(selectedVertex, colorTheme));
        vv.getRenderContext().setVertexShapeTransformer(getVertexShapeTransformer());
        vv.getRenderContext().setVertexDrawPaintTransformer(getVertexOutlineTransformer());

        // Label display mode: 0=Full, 1=Short, 2=None
        final int[] labelMode = {1};  // Start with Short names

        // Vertex label transformer
        final Function<ModelVertex, String> labelTransformer = new Function<ModelVertex, String>() {
            @Override
            public String apply(ModelVertex v) {
                switch (labelMode[0]) {
                    case 1: return v.getShortName();
                    case 2: return "";
                    default:
                        String name = v.getName();
                        if (name.length() > 20) return v.getShortName();
                        return name;
                }
            }
        };
        vv.getRenderContext().setVertexLabelTransformer(labelTransformer);

        // Custom vertex label renderer with background
        vv.getRenderContext().setVertexLabelRenderer(new edu.uci.ics.jung.visualization.renderers.DefaultVertexLabelRenderer(Color.BLACK) {
            @Override
            public <V> Component getVertexLabelRendererComponent(javax.swing.JComponent vvComponent, Object value,
                    Font font, boolean isSelected, V vertex) {
                Component c = super.getVertexLabelRendererComponent(vvComponent, value, font, isSelected, vertex);
                if (c instanceof JLabel) {
                    JLabel label = (JLabel) c;
                    String text = value != null ? value.toString() : "";
                    if (text.isEmpty()) {
                        label.setOpaque(false);
                        label.setBorder(null);
                        label.setPreferredSize(new Dimension(0, 0));
                    } else {
                        label.setOpaque(true);
                        label.setBackground(new Color(255, 255, 255, 191));
                        label.setBorder(BorderFactory.createCompoundBorder(
                            BorderFactory.createLineBorder(new Color(180, 180, 180), 1),
                            BorderFactory.createEmptyBorder(1, 3, 1, 3)
                        ));
                        label.setFont(new Font("Segoe UI", Font.PLAIN, 11));
                        label.setPreferredSize(null);
                    }
                }
                return c;
            }
        });

        // Apply edge transformers
        vv.getRenderContext().setEdgeDrawPaintTransformer(getEdgePaintTransformer());
        vv.getRenderContext().setEdgeStrokeTransformer(getEdgeStrokeTransformer());
        vv.getRenderContext().setEdgeShapeTransformer(EdgeShape.line(graph));

        // Add mouse interaction
        edu.uci.ics.jung.visualization.control.PluggableGraphMouse graphMouse =
            new edu.uci.ics.jung.visualization.control.PluggableGraphMouse();
        graphMouse.add(new edu.uci.ics.jung.visualization.control.PickingGraphMousePlugin<ModelVertex, ModelEdge>() {
            @Override
            public void mousePressed(MouseEvent e) {
                down = e.getPoint();
                ModelVertex v = vertex;
                if (v == null && vv.getPickSupport() != null) {
                    v = vv.getPickSupport().getVertex(vv.getGraphLayout(), e.getX(), e.getY());
                }
                if (v != null) {
                    super.mousePressed(e);
                }
            }
            @Override
            public void mouseDragged(MouseEvent e) {
                if (vertex != null) {
                    super.mouseDragged(e);
                }
            }
        });
        graphMouse.add(new edu.uci.ics.jung.visualization.control.ScalingGraphMousePlugin(
            new edu.uci.ics.jung.visualization.control.CrossoverScalingControl(), 0, 1.1f, 0.9f));
        vv.setGraphMouse(graphMouse);

        // Add background panning
        final Point2D[] dragStart = {null};
        final Point2D[] viewStart = {null};
        vv.addMouseListener(new MouseAdapter() {
            @Override
            public void mousePressed(MouseEvent e) {
                ModelVertex v = vv.getPickSupport() != null ?
                    vv.getPickSupport().getVertex(vv.getGraphLayout(), e.getX(), e.getY()) : null;
                if (v == null) {
                    dragStart[0] = e.getPoint();
                    edu.uci.ics.jung.visualization.transform.MutableTransformer viewTransformer =
                        vv.getRenderContext().getMultiLayerTransformer()
                            .getTransformer(edu.uci.ics.jung.visualization.Layer.VIEW);
                    viewStart[0] = new Point2D.Double(viewTransformer.getTranslateX(), viewTransformer.getTranslateY());
                }
            }
            @Override
            public void mouseReleased(MouseEvent e) {
                dragStart[0] = null;
                viewStart[0] = null;
            }
        });

        vv.addMouseMotionListener(new java.awt.event.MouseMotionAdapter() {
            @Override
            public void mouseDragged(MouseEvent e) {
                if (dragStart[0] != null && viewStart[0] != null) {
                    double dx = e.getX() - dragStart[0].getX();
                    double dy = e.getY() - dragStart[0].getY();
                    edu.uci.ics.jung.visualization.transform.MutableTransformer viewTransformer =
                        vv.getRenderContext().getMultiLayerTransformer()
                            .getTransformer(edu.uci.ics.jung.visualization.Layer.VIEW);
                    viewTransformer.setTranslate(viewStart[0].getX() + dx, viewStart[0].getY() + dy);
                    vv.repaint();
                }
            }
        });

        // Create toolbar
        JToolBar toolBar = new JToolBar();
        toolBar.setFloatable(false);
        toolBar.setBorder(BorderFactory.createCompoundBorder(
            BorderFactory.createMatteBorder(0, 0, 1, 0, new Color(200, 200, 200)),
            BorderFactory.createEmptyBorder(4, 8, 4, 8)
        ));
        toolBar.setBackground(new Color(250, 250, 252));

        final int ICON_SIZE = 27;

        // Zoom slider
        final JSlider zoomSlider = new JSlider(JSlider.VERTICAL, -100, 100, 0);

        // Track spacing factors
        final double[] hSpaceFactor = {1.0};
        final double[] vSpaceFactor = {1.0};

        // Redraw button
        JButton redrawButton = new JButton(createRefreshIcon(ICON_SIZE, new Color(230, 126, 34)));
        redrawButton.setToolTipText("Recalculate layout to fit current view");
        redrawButton.addActionListener(new ActionListener() {
            @Override
            public void actionPerformed(ActionEvent e) {
                Dimension size = vv.getSize();
                int newWidth = Math.max(400, size.width);
                int newHeight = Math.max(300, size.height);
                Layout<ModelVertex, ModelEdge> newLayout = createHierarchicalLayout(newWidth, newHeight, vSpaceFactor[0]);
                newLayout.setSize(new Dimension(newWidth - 50, newHeight - 50));
                vv.getRenderContext().getMultiLayerTransformer()
                    .getTransformer(edu.uci.ics.jung.visualization.Layer.VIEW).setToIdentity();
                vv.getRenderContext().getMultiLayerTransformer()
                    .getTransformer(edu.uci.ics.jung.visualization.Layer.LAYOUT).setToIdentity();
                vv.setGraphLayout(newLayout);
                zoomSlider.setValue(0);
                vv.repaint();
            }
        });
        toolBar.add(redrawButton);

        toolBar.addSeparator();

        // Load background image button (LQN only feature but available for both)
        JButton loadBgButton = new JButton(createMapIcon(ICON_SIZE, new Color(192, 57, 43)));
        loadBgButton.setToolTipText("Load background map/image");
        loadBgButton.addActionListener(new ActionListener() {
            @Override
            public void actionPerformed(ActionEvent e) {
                String[] options = {"From File", "From URL", "Clear Background"};
                int choice = JOptionPane.showOptionDialog(vv, "Load background image:",
                    "Background Image", JOptionPane.DEFAULT_OPTION, JOptionPane.QUESTION_MESSAGE,
                    null, options, options[0]);

                if (choice == 0) {
                    JFileChooser fileChooser = new JFileChooser();
                    fileChooser.setFileFilter(new FileNameExtensionFilter(
                        "Image files", "jpg", "jpeg", "png", "gif", "bmp"));
                    if (fileChooser.showOpenDialog(vv) == JFileChooser.APPROVE_OPTION) {
                        try {
                            BufferedImage img = ImageIO.read(fileChooser.getSelectedFile());
                            setBackgroundImage(vv, img);
                        } catch (IOException ex) {
                            JOptionPane.showMessageDialog(vv, "Error loading image: " + ex.getMessage(),
                                "Error", JOptionPane.ERROR_MESSAGE);
                        }
                    }
                } else if (choice == 1) {
                    String urlStr = JOptionPane.showInputDialog(vv, "Enter image URL:");
                    if (urlStr != null && !urlStr.isEmpty()) {
                        try {
                            BufferedImage img = ImageIO.read(new URL(urlStr));
                            setBackgroundImage(vv, img);
                        } catch (IOException ex) {
                            JOptionPane.showMessageDialog(vv, "Error loading image: " + ex.getMessage(),
                                "Error", JOptionPane.ERROR_MESSAGE);
                        }
                    }
                } else if (choice == 2) {
                    clearBackgroundImage(vv);
                }
            }
        });
        toolBar.add(loadBgButton);

        toolBar.addSeparator();

        // Play button for solver (LQN only)
        if (isLayeredNetwork) {
            final JTabbedPane statusPane = new JTabbedPane();
            statusPane.setFont(new Font("Segoe UI", Font.PLAIN, 12));

            final JTextArea consoleArea = new JTextArea();
            consoleArea.setEditable(false);
            consoleArea.setFont(new Font("JetBrains Mono", Font.PLAIN, 12));
            consoleArea.setBackground(new Color(30, 30, 30));
            consoleArea.setForeground(new Color(212, 212, 212));
            consoleArea.setBorder(BorderFactory.createEmptyBorder(8, 8, 8, 8));
            statusPane.addTab("Console", new JScrollPane(consoleArea));

            final JPanel resultsPanel = new JPanel(new BorderLayout());
            resultsPanel.setBackground(new Color(250, 250, 252));
            statusPane.addTab("Results", resultsPanel);

            final JButton playButton = new JButton(createPlayIcon(ICON_SIZE, new Color(46, 204, 113)));
            playButton.setToolTipText("Run LINE solver on this model");
            playButton.addActionListener(new ActionListener() {
                @Override
                public void actionPerformed(ActionEvent e) {
                    runSolver(layeredModel, playButton, consoleArea, resultsPanel, statusPane);
                }
            });
            toolBar.add(playButton);

            toolBar.addSeparator();
        }

        // Label mode dropdown
        final String[] labelModeNames = {"Full", "Short", "None"};
        final JComboBox<String> labelCombo = new JComboBox<String>(labelModeNames);
        labelCombo.setSelectedIndex(1);
        labelCombo.setToolTipText("Select label display mode");
        labelCombo.setFont(new Font("Segoe UI", Font.PLAIN, 11));
        labelCombo.setMaximumSize(new Dimension(80, 28));
        labelCombo.addActionListener(new ActionListener() {
            @Override
            public void actionPerformed(ActionEvent e) {
                labelMode[0] = labelCombo.getSelectedIndex();
                vv.repaint();
            }
        });
        toolBar.add(labelCombo);

        toolBar.addSeparator();

        // Recenter button
        JButton recenterButton = new JButton(createCenterIcon(ICON_SIZE, new Color(155, 89, 182)));
        recenterButton.setToolTipText("Fit and center graph in view");
        recenterButton.addActionListener(new ActionListener() {
            @Override
            public void actionPerformed(ActionEvent e) {
                fitGraphToView(vv, graph, zoomSlider);
            }
        });
        toolBar.add(recenterButton);

        // Export button
        JButton exportButton = new JButton(createExportIcon(ICON_SIZE, new Color(52, 73, 94)));
        exportButton.setToolTipText("Export graph as PNG or TikZ");
        exportButton.addActionListener(new ActionListener() {
            @Override
            public void actionPerformed(ActionEvent e) {
                String[] options = {"PNG Image", "TikZ Code"};
                int choice = JOptionPane.showOptionDialog(vv, "Export format:",
                    "Export Graph", JOptionPane.DEFAULT_OPTION, JOptionPane.QUESTION_MESSAGE,
                    null, options, options[0]);

                if (choice < 0) return;

                JFileChooser fileChooser = new JFileChooser();
                String extension = choice == 0 ? "png" : "tex";
                fileChooser.setSelectedFile(new java.io.File("graph." + extension));

                if (fileChooser.showSaveDialog(vv) == JFileChooser.APPROVE_OPTION) {
                    java.io.File file = fileChooser.getSelectedFile();
                    try {
                        if (choice == 0) {
                            exportToPNG(vv, file);
                        } else {
                            exportToTikZ(vv, graph, file);
                        }
                        JOptionPane.showMessageDialog(vv, "Exported to: " + file.getAbsolutePath(),
                            "Export Successful", JOptionPane.INFORMATION_MESSAGE);
                    } catch (Exception ex) {
                        JOptionPane.showMessageDialog(vv, "Export failed: " + ex.getMessage(),
                            "Error", JOptionPane.ERROR_MESSAGE);
                    }
                }
            }
        });
        toolBar.add(exportButton);

        toolBar.addSeparator();

        // Color theme selector dropdown
        final String[] themeNames = {"Color", "Multicolor", "Greyscale", "Black/White"};
        final JComboBox<String> themeCombo = new JComboBox<String>(themeNames);
        themeCombo.setFont(new Font("Segoe UI", Font.PLAIN, 11));
        themeCombo.setToolTipText("Select color theme");
        themeCombo.addActionListener(new ActionListener() {
            @Override
            public void actionPerformed(ActionEvent e) {
                colorTheme[0] = themeCombo.getSelectedIndex();
                vv.repaint();
            }
        });
        toolBar.add(themeCombo);

        // Create second toolbar for sliders
        JToolBar sliderToolBar = new JToolBar();
        sliderToolBar.setFloatable(false);
        sliderToolBar.setBorder(BorderFactory.createCompoundBorder(
            BorderFactory.createMatteBorder(0, 0, 1, 0, new Color(200, 200, 200)),
            BorderFactory.createEmptyBorder(2, 8, 2, 8)
        ));
        sliderToolBar.setBackground(new Color(250, 250, 252));

        JPanel slidersPanel = new JPanel(new java.awt.GridLayout(1, 3, 12, 0));
        slidersPanel.setBackground(new Color(250, 250, 252));
        slidersPanel.setBorder(BorderFactory.createEmptyBorder(0, 0, 0, 0));

        // Zoom slider panel
        JPanel zoomPanel = new JPanel(new BorderLayout(4, 0));
        zoomPanel.setBackground(new Color(250, 250, 252));
        zoomPanel.setBorder(BorderFactory.createCompoundBorder(
            BorderFactory.createLineBorder(new Color(41, 128, 185), 1),
            BorderFactory.createEmptyBorder(2, 6, 2, 6)
        ));
        JLabel zoomLabel = new JLabel("Zoom");
        zoomLabel.setFont(new Font("Segoe UI", Font.BOLD, 10));
        zoomLabel.setForeground(new Color(41, 128, 185));
        zoomPanel.add(zoomLabel, BorderLayout.WEST);

        zoomSlider.setOrientation(JSlider.HORIZONTAL);
        zoomSlider.setMajorTickSpacing(50);
        zoomSlider.setMinorTickSpacing(25);
        zoomSlider.setPaintTicks(true);
        zoomSlider.setPaintLabels(false);
        zoomSlider.setBackground(new Color(250, 250, 252));
        zoomSlider.putClientProperty("JSlider.isFilled", Boolean.TRUE);
        zoomSlider.putClientProperty("Slider.thumbSize", new Dimension(12, 12));
        zoomSlider.putClientProperty("Slider.trackWidth", 4);

        zoomSlider.addChangeListener(new ChangeListener() {
            @Override
            public void stateChanged(ChangeEvent e) {
                double value = zoomSlider.getValue();
                double scale = Math.pow(2, value / 50.0);
                vv.getRenderContext().getMultiLayerTransformer()
                    .getTransformer(edu.uci.ics.jung.visualization.Layer.VIEW).setToIdentity();
                vv.getRenderContext().getMultiLayerTransformer()
                    .getTransformer(edu.uci.ics.jung.visualization.Layer.VIEW)
                    .scale(scale, scale, vv.getCenter());
                vv.repaint();
            }
        });
        zoomPanel.add(zoomSlider, BorderLayout.CENTER);
        slidersPanel.add(zoomPanel);

        // HSpace slider panel
        JPanel hSpacePanel = new JPanel(new BorderLayout(4, 0));
        hSpacePanel.setBackground(new Color(250, 250, 252));
        hSpacePanel.setBorder(BorderFactory.createCompoundBorder(
            BorderFactory.createLineBorder(new Color(39, 174, 96), 1),
            BorderFactory.createEmptyBorder(2, 6, 2, 6)
        ));
        JLabel hSpaceLabel = new JLabel("HSpace");
        hSpaceLabel.setFont(new Font("Segoe UI", Font.BOLD, 10));
        hSpaceLabel.setForeground(new Color(39, 174, 96));
        hSpacePanel.add(hSpaceLabel, BorderLayout.WEST);

        final JSlider hSpaceSlider = new JSlider(JSlider.HORIZONTAL, 50, 150, 100);
        hSpaceSlider.setMajorTickSpacing(25);
        hSpaceSlider.setMinorTickSpacing(5);
        hSpaceSlider.setPaintTicks(true);
        hSpaceSlider.setPaintLabels(false);
        hSpaceSlider.setBackground(new Color(250, 250, 252));
        hSpaceSlider.putClientProperty("JSlider.isFilled", Boolean.TRUE);
        hSpaceSlider.putClientProperty("Slider.thumbSize", new Dimension(12, 12));
        hSpaceSlider.putClientProperty("Slider.trackWidth", 4);

        hSpaceSlider.addChangeListener(new ChangeListener() {
            @Override
            public void stateChanged(ChangeEvent e) {
                if (!hSpaceSlider.getValueIsAdjusting()) {
                    hSpaceFactor[0] = hSpaceSlider.getValue() / 100.0;
                    Dimension size = vv.getSize();
                    int scaledWidth = (int) (Math.max(400, size.width) * hSpaceFactor[0]);
                    int baseHeight = Math.max(300, size.height);
                    Layout<ModelVertex, ModelEdge> newLayout = createHierarchicalLayout(scaledWidth, baseHeight, vSpaceFactor[0]);
                    newLayout.setSize(new Dimension(scaledWidth - 50, baseHeight - 50));
                    vv.setGraphLayout(newLayout);
                    vv.repaint();
                }
            }
        });
        hSpacePanel.add(hSpaceSlider, BorderLayout.CENTER);
        slidersPanel.add(hSpacePanel);

        // VSpace slider panel
        JPanel vSpacePanel = new JPanel(new BorderLayout(4, 0));
        vSpacePanel.setBackground(new Color(250, 250, 252));
        vSpacePanel.setBorder(BorderFactory.createCompoundBorder(
            BorderFactory.createLineBorder(new Color(230, 126, 34), 1),
            BorderFactory.createEmptyBorder(2, 6, 2, 6)
        ));
        JLabel vSpaceLabel = new JLabel("VSpace");
        vSpaceLabel.setFont(new Font("Segoe UI", Font.BOLD, 10));
        vSpaceLabel.setForeground(new Color(230, 126, 34));
        vSpacePanel.add(vSpaceLabel, BorderLayout.WEST);

        final JSlider vSpaceSlider = new JSlider(JSlider.HORIZONTAL, 50, 150, 100);
        vSpaceSlider.setMajorTickSpacing(25);
        vSpaceSlider.setMinorTickSpacing(5);
        vSpaceSlider.setPaintTicks(true);
        vSpaceSlider.setPaintLabels(false);
        vSpaceSlider.setBackground(new Color(250, 250, 252));
        vSpaceSlider.putClientProperty("JSlider.isFilled", Boolean.TRUE);
        vSpaceSlider.putClientProperty("Slider.thumbSize", new Dimension(12, 12));
        vSpaceSlider.putClientProperty("Slider.trackWidth", 4);

        vSpaceSlider.addChangeListener(new ChangeListener() {
            @Override
            public void stateChanged(ChangeEvent e) {
                if (!vSpaceSlider.getValueIsAdjusting()) {
                    vSpaceFactor[0] = vSpaceSlider.getValue() / 100.0;
                    Dimension size = vv.getSize();
                    int scaledWidth = (int) (Math.max(400, size.width) * hSpaceFactor[0]);
                    int baseHeight = Math.max(300, size.height);
                    Layout<ModelVertex, ModelEdge> newLayout = createHierarchicalLayout(scaledWidth, baseHeight, vSpaceFactor[0]);
                    newLayout.setSize(new Dimension(scaledWidth - 50, baseHeight - 50));
                    vv.setGraphLayout(newLayout);
                    vv.repaint();
                }
            }
        });
        vSpacePanel.add(vSpaceSlider, BorderLayout.CENTER);
        slidersPanel.add(vSpacePanel);

        sliderToolBar.add(slidersPanel);

        // Create panel to hold both toolbars
        JPanel toolBarPanel = new JPanel(new BorderLayout());
        toolBarPanel.add(toolBar, BorderLayout.NORTH);
        toolBarPanel.add(sliderToolBar, BorderLayout.SOUTH);

        // Auto-fit graph at startup
        fitGraphToView(vv, graph, zoomSlider);
        fitGraphToView(vv, graph, zoomSlider);

        // Create info panel for node details
        final JPanel infoPanel = new JPanel(new BorderLayout());
        infoPanel.setBorder(BorderFactory.createCompoundBorder(
            BorderFactory.createMatteBorder(0, 1, 0, 0, new Color(200, 200, 200)),
            BorderFactory.createEmptyBorder(8, 8, 8, 8)
        ));
        infoPanel.setPreferredSize(new Dimension(320, 0));
        infoPanel.setBackground(new Color(250, 250, 252));

        final JPanel propsContent = new JPanel();
        propsContent.setLayout(new BoxLayout(propsContent, BoxLayout.Y_AXIS));
        propsContent.setBackground(new Color(250, 250, 252));
        JScrollPane infoScroll = new JScrollPane(propsContent);
        infoScroll.setBorder(BorderFactory.createEmptyBorder());
        infoScroll.getVerticalScrollBar().setUnitIncrement(16);
        infoPanel.add(infoScroll, BorderLayout.CENTER);

        final JLabel headerLabel = new JLabel(" Node Details");
        headerLabel.setFont(new Font("Segoe UI", Font.BOLD, 14));
        headerLabel.setForeground(new Color(60, 60, 60));

        JButton closeInfoButton = new JButton("\u2715");
        closeInfoButton.setToolTipText("Close info panel");
        closeInfoButton.setFont(new Font("Segoe UI", Font.PLAIN, 14));
        closeInfoButton.setBorderPainted(false);
        closeInfoButton.setContentAreaFilled(false);
        closeInfoButton.setFocusPainted(false);
        closeInfoButton.addActionListener(new ActionListener() {
            @Override
            public void actionPerformed(ActionEvent e) {
                selectedVertex[0] = null;
                vv.repaint();
                infoPanel.setVisible(false);
                infoPanel.getParent().revalidate();
            }
        });
        JPanel infoHeader = new JPanel(new BorderLayout());
        infoHeader.setBackground(new Color(250, 250, 252));
        infoHeader.add(headerLabel, BorderLayout.CENTER);
        infoHeader.add(closeInfoButton, BorderLayout.EAST);
        infoHeader.setBorder(BorderFactory.createEmptyBorder(0, 0, 8, 0));
        infoPanel.add(infoHeader, BorderLayout.NORTH);
        infoPanel.setVisible(false);

        // Add double-click listener for node properties
        vv.addMouseListener(new MouseAdapter() {
            @Override
            public void mouseClicked(MouseEvent e) {
                if (e.getClickCount() == 2) {
                    Point2D p = e.getPoint();
                    edu.uci.ics.jung.algorithms.layout.GraphElementAccessor<ModelVertex, ModelEdge> pickSupport =
                        vv.getPickSupport();
                    if (pickSupport != null) {
                        ModelVertex vertex = pickSupport.getVertex(vv.getGraphLayout(), p.getX(), p.getY());
                        if (vertex != null) {
                            selectedVertex[0] = vertex;
                            vv.repaint();

                            propsContent.removeAll();
                            headerLabel.setText(" " + vertex.getName());
                            buildPropertiesPanel(propsContent, vertex, graph, vv);
                            propsContent.revalidate();
                            propsContent.repaint();
                            infoPanel.setVisible(true);
                            infoPanel.getParent().revalidate();
                        }
                    }
                }
            }
        });

        // Create main content panel
        JPanel graphPanel = new JPanel(new BorderLayout());
        graphPanel.add(vv, BorderLayout.CENTER);
        if (showToolbars) {
            graphPanel.add(infoPanel, BorderLayout.EAST);
        }

        // Create and display JFrame
        JFrame frame = new JFrame(title);
        frame.setDefaultCloseOperation(JFrame.DISPOSE_ON_CLOSE);
        frame.setLayout(new BorderLayout());
        if (showToolbars) {
            frame.add(toolBarPanel, BorderLayout.NORTH);
        }
        frame.add(graphPanel, BorderLayout.CENTER);
        frame.pack();
        frame.setLocationRelativeTo(null);
        frame.setVisible(true);
    }

    // ========================================================================
    // LAYOUT ALGORITHMS
    // ========================================================================

    /**
     * Creates a hierarchical layout appropriate for the model type.
     */
    public Layout<ModelVertex, ModelEdge> createHierarchicalLayout(int width, int height) {
        return createHierarchicalLayout(width, height, 1.0);
    }

    /**
     * Creates a hierarchical layout with vertical spacing factor.
     */
    public Layout<ModelVertex, ModelEdge> createHierarchicalLayout(int width, int height, double vSpaceFactor) {
        if (isLayeredNetwork) {
            return createLQNHierarchicalLayout(width, height, vSpaceFactor);
        } else {
            return createNetworkSugiyamaLayout(width, height, vSpaceFactor);
        }
    }

    // ========== NETWORK SUGIYAMA LAYOUT ==========

    /**
     * Creates a Sugiyama-style hierarchical layout for Network models.
     */
    private Layout<ModelVertex, ModelEdge> createNetworkSugiyamaLayout(int width, int height, double vSpaceFactor) {
        StaticLayout<ModelVertex, ModelEdge> layout = new StaticLayout<ModelVertex, ModelEdge>(graph);

        // Build adjacency information
        Map<ModelVertex, List<ModelVertex>> successors = buildSuccessorMap();
        Map<ModelVertex, List<ModelVertex>> predecessors = buildPredecessorMap(successors);

        // Step 1: Assign nodes to layers using topological sort
        assignNetworkLayers(predecessors, successors);

        // Step 2: Order nodes within layers to minimize crossings
        minimizeCrossings(successors, predecessors);

        // Step 3: Assign coordinates
        assignNetworkCoordinates(layout, width, height, vSpaceFactor);

        return layout;
    }

    private Map<ModelVertex, List<ModelVertex>> buildSuccessorMap() {
        Map<ModelVertex, List<ModelVertex>> successors = new HashMap<ModelVertex, List<ModelVertex>>();
        for (ModelVertex v : graph.getVertices()) {
            successors.put(v, new ArrayList<ModelVertex>());
        }

        for (ModelEdge edge : graph.getEdges()) {
            ModelVertex from = graph.getSource(edge);
            ModelVertex to = graph.getDest(edge);
            if (from != null && to != null) {
                successors.get(from).add(to);
            }
        }

        return successors;
    }

    private Map<ModelVertex, List<ModelVertex>> buildPredecessorMap(Map<ModelVertex, List<ModelVertex>> successors) {
        Map<ModelVertex, List<ModelVertex>> predecessors = new HashMap<ModelVertex, List<ModelVertex>>();
        for (ModelVertex v : graph.getVertices()) {
            predecessors.put(v, new ArrayList<ModelVertex>());
        }

        for (Map.Entry<ModelVertex, List<ModelVertex>> entry : successors.entrySet()) {
            ModelVertex from = entry.getKey();
            for (ModelVertex to : entry.getValue()) {
                predecessors.get(to).add(from);
            }
        }

        return predecessors;
    }

    private void assignNetworkLayers(Map<ModelVertex, List<ModelVertex>> predecessors,
                                      Map<ModelVertex, List<ModelVertex>> successors) {
        layers.clear();
        Set<ModelVertex> assigned = new HashSet<ModelVertex>();
        List<ModelVertex> allVertices = new ArrayList<ModelVertex>(graph.getVertices());

        // Find all sources and nodes with no predecessors
        List<ModelVertex> layer0 = new ArrayList<ModelVertex>();
        for (ModelVertex v : allVertices) {
            if (v.getType() == NodeType.SOURCE || predecessors.get(v).isEmpty()) {
                layer0.add(v);
                assigned.add(v);
            }
        }

        if (layer0.isEmpty() && !allVertices.isEmpty()) {
            ModelVertex first = allVertices.get(0);
            layer0.add(first);
            assigned.add(first);
        }

        layers.add(layer0);

        int currentLayerIdx = 0;
        while (assigned.size() < allVertices.size()) {
            List<ModelVertex> nextLayer = new ArrayList<ModelVertex>();

            for (ModelVertex v : allVertices) {
                if (!assigned.contains(v)) {
                    List<ModelVertex> preds = predecessors.get(v);
                    boolean allPredsAssigned = true;
                    for (ModelVertex pred : preds) {
                        if (!assigned.contains(pred)) {
                            allPredsAssigned = false;
                            break;
                        }
                    }
                    if (allPredsAssigned) {
                        nextLayer.add(v);
                    }
                }
            }

            if (nextLayer.isEmpty()) {
                for (ModelVertex v : allVertices) {
                    if (!assigned.contains(v)) {
                        nextLayer.add(v);
                        break;
                    }
                }
            }

            if (!nextLayer.isEmpty()) {
                for (ModelVertex v : nextLayer) {
                    assigned.add(v);
                }
                layers.add(nextLayer);
            }

            currentLayerIdx++;
            if (currentLayerIdx > allVertices.size()) {
                break;
            }
        }

        // Move sinks to the last layer
        if (layers.size() > 1) {
            List<ModelVertex> lastLayer = layers.get(layers.size() - 1);
            for (int i = 0; i < layers.size() - 1; i++) {
                List<ModelVertex> layer = layers.get(i);
                Iterator<ModelVertex> it = layer.iterator();
                while (it.hasNext()) {
                    ModelVertex v = it.next();
                    if (v.getType() == NodeType.SINK) {
                        it.remove();
                        if (!lastLayer.contains(v)) {
                            lastLayer.add(v);
                        }
                    }
                }
            }
            List<List<ModelVertex>> nonEmptyLayers = new ArrayList<List<ModelVertex>>();
            for (List<ModelVertex> layer : layers) {
                if (!layer.isEmpty()) {
                    nonEmptyLayers.add(layer);
                }
            }
            layers = nonEmptyLayers;
        }
    }

    private void minimizeCrossings(Map<ModelVertex, List<ModelVertex>> successors,
                                    Map<ModelVertex, List<ModelVertex>> predecessors) {
        for (int pass = 0; pass < 4; pass++) {
            for (int i = 1; i < layers.size(); i++) {
                reorderLayerByBarycenter(layers.get(i), predecessors, layers.get(i - 1));
            }
            for (int i = layers.size() - 2; i >= 0; i--) {
                reorderLayerByBarycenter(layers.get(i), successors, layers.get(i + 1));
            }
        }
    }

    private void reorderLayerByBarycenter(List<ModelVertex> layer,
                                           Map<ModelVertex, List<ModelVertex>> connections,
                                           List<ModelVertex> referenceLayer) {
        final Map<ModelVertex, Integer> refPositions = new HashMap<ModelVertex, Integer>();
        for (int i = 0; i < referenceLayer.size(); i++) {
            refPositions.put(referenceLayer.get(i), i);
        }

        final Map<ModelVertex, Double> barycenters = new HashMap<ModelVertex, Double>();
        for (ModelVertex v : layer) {
            List<ModelVertex> connected = connections.get(v);
            if (connected == null || connected.isEmpty()) {
                barycenters.put(v, Double.MAX_VALUE);
            } else {
                double sum = 0;
                int count = 0;
                for (ModelVertex conn : connected) {
                    Integer pos = refPositions.get(conn);
                    if (pos != null) {
                        sum += pos;
                        count++;
                    }
                }
                barycenters.put(v, count > 0 ? sum / count : Double.MAX_VALUE);
            }
        }

        Collections.sort(layer, new Comparator<ModelVertex>() {
            @Override
            public int compare(ModelVertex a, ModelVertex b) {
                return Double.compare(barycenters.get(a), barycenters.get(b));
            }
        });
    }

    private void assignNetworkCoordinates(StaticLayout<ModelVertex, ModelEdge> layout, int width, int height, double vSpaceFactor) {
        int margin = 60;

        double effectiveLayerSpacing = layers.size() > 1
            ? (width - 2 * margin) / (double)(layers.size() - 1)
            : LAYER_SPACING;
        effectiveLayerSpacing = Math.min(effectiveLayerSpacing, LAYER_SPACING);

        for (int layerIdx = 0; layerIdx < layers.size(); layerIdx++) {
            List<ModelVertex> layer = layers.get(layerIdx);
            double x = margin + layerIdx * effectiveLayerSpacing;

            double effectiveNodeSpacing = layer.size() > 1
                ? Math.min((height - 2 * margin) / (double)(layer.size() - 1), NODE_SPACING * vSpaceFactor)
                : NODE_SPACING * vSpaceFactor;
            double totalHeight = (layer.size() - 1) * effectiveNodeSpacing;
            double startY = (height - totalHeight) / 2.0;

            for (int nodeIdx = 0; nodeIdx < layer.size(); nodeIdx++) {
                ModelVertex v = layer.get(nodeIdx);
                double y = startY + nodeIdx * effectiveNodeSpacing;
                layout.setLocation(v, new Point2D.Double(x, y));
            }
        }
    }

    // ========== LQN HIERARCHICAL LAYOUT ==========

    /**
     * Creates a hierarchical layout for LayeredNetwork models.
     */
    private Layout<ModelVertex, ModelEdge> createLQNHierarchicalLayout(int width, int height, double vSpaceFactor) {
        StaticLayout<ModelVertex, ModelEdge> layout = new StaticLayout<ModelVertex, ModelEdge>(graph);

        // Collect vertices by type
        List<ModelVertex> hosts = new ArrayList<ModelVertex>();
        List<ModelVertex> tasks = new ArrayList<ModelVertex>();
        List<ModelVertex> entries = new ArrayList<ModelVertex>();
        List<ModelVertex> activities = new ArrayList<ModelVertex>();

        for (ModelVertex v : graph.getVertices()) {
            switch (v.getType()) {
                case HOST: hosts.add(v); break;
                case TASK: tasks.add(v); break;
                case ENTRY: entries.add(v); break;
                case ACTIVITY: activities.add(v); break;
            }
        }

        // Build parent mapping for ordering
        final Map<ModelVertex, ModelVertex> parentMap = new HashMap<ModelVertex, ModelVertex>();
        final Map<ModelVertex, ModelVertex> boundToMap = new HashMap<ModelVertex, ModelVertex>();

        for (ModelEdge edge : graph.getEdges()) {
            if (edge.getType() == EdgeType.PARENT_CHILD) {
                ModelVertex parent = graph.getSource(edge);
                ModelVertex child = graph.getDest(edge);
                parentMap.put(child, parent);
            } else if (edge.getType() == EdgeType.BOUND_TO) {
                ModelVertex entry = graph.getSource(edge);
                ModelVertex activity = graph.getDest(edge);
                boundToMap.put(activity, entry);
            }
        }

        // Step 1: Assign sublayers to activities based on precedence depth
        Map<ModelVertex, Integer> activitySublayer = assignActivitySublayers(activities);
        int maxSublayer = 0;
        for (Integer sublayer : activitySublayer.values()) {
            maxSublayer = Math.max(maxSublayer, sublayer);
        }

        // Step 2: Order hosts
        Collections.sort(hosts, new Comparator<ModelVertex>() {
            @Override
            public int compare(ModelVertex a, ModelVertex b) {
                return a.getName().compareTo(b.getName());
            }
        });

        // Assign initial X positions to hosts
        final Map<ModelVertex, Double> xPositions = new HashMap<ModelVertex, Double>();
        double hostSpacing = hosts.size() > 0 ? (width - 100.0) / hosts.size() : width - 100;
        for (int i = 0; i < hosts.size(); i++) {
            xPositions.put(hosts.get(i), 50 + hostSpacing * (i + 0.5));
        }

        // Step 3: Order tasks by parent host position
        orderByBarycenterLQN(tasks, parentMap, xPositions);
        double taskSpacing = tasks.size() > 0 ? (width - 100.0) / tasks.size() : width - 100;
        for (int i = 0; i < tasks.size(); i++) {
            xPositions.put(tasks.get(i), 50 + taskSpacing * (i + 0.5));
        }

        // Step 4: Order entries by parent task position
        orderByBarycenterLQN(entries, parentMap, xPositions);
        double entrySpacing = entries.size() > 0 ? (width - 100.0) / entries.size() : width - 100;
        for (int i = 0; i < entries.size(); i++) {
            xPositions.put(entries.get(i), 50 + entrySpacing * (i + 0.5));
        }

        // Step 5: Group activities by sublayer and order each sublayer
        List<List<ModelVertex>> activityLayers = new ArrayList<List<ModelVertex>>();
        for (int i = 0; i <= maxSublayer; i++) {
            activityLayers.add(new ArrayList<ModelVertex>());
        }
        for (ModelVertex act : activities) {
            int sublayer = activitySublayer.get(act);
            activityLayers.get(sublayer).add(act);
        }

        // Order first activity sublayer by bound entry position
        if (!activityLayers.isEmpty() && !activityLayers.get(0).isEmpty()) {
            orderByBarycenterLQN(activityLayers.get(0), boundToMap, xPositions);
            double actSpacing = activityLayers.get(0).size() > 0
                ? (width - 100.0) / activityLayers.get(0).size() : width - 100;
            for (int i = 0; i < activityLayers.get(0).size(); i++) {
                xPositions.put(activityLayers.get(0).get(i), 50 + actSpacing * (i + 0.5));
            }
        }

        // Order subsequent sublayers by precedence connections
        for (int layer = 1; layer <= maxSublayer; layer++) {
            List<ModelVertex> currentLayer = activityLayers.get(layer);
            if (currentLayer.isEmpty()) continue;

            final Map<ModelVertex, List<ModelVertex>> predecessors = new HashMap<ModelVertex, List<ModelVertex>>();
            for (ModelVertex act : currentLayer) {
                predecessors.put(act, new ArrayList<ModelVertex>());
            }

            for (ModelEdge edge : graph.getEdges()) {
                if (edge.getType() == EdgeType.PRECEDENCE) {
                    ModelVertex source = graph.getSource(edge);
                    ModelVertex dest = graph.getDest(edge);
                    if (predecessors.containsKey(dest) && xPositions.containsKey(source)) {
                        predecessors.get(dest).add(source);
                    }
                }
            }

            Collections.sort(currentLayer, new Comparator<ModelVertex>() {
                @Override
                public int compare(ModelVertex a, ModelVertex b) {
                    double bcA = computeBarycenter(predecessors.get(a), xPositions);
                    double bcB = computeBarycenter(predecessors.get(b), xPositions);
                    return Double.compare(bcA, bcB);
                }
            });

            double actSpacing = currentLayer.size() > 0
                ? (width - 100.0) / currentLayer.size() : width - 100;
            for (int i = 0; i < currentLayer.size(); i++) {
                xPositions.put(currentLayer.get(i), 50 + actSpacing * (i + 0.5));
            }
        }

        // Step 6: Calculate Y positions
        int baseY = 50;
        int baseLayerHeight = Math.max(80, (height - 100) / (4 + maxSublayer));
        int layerHeight = (int)(baseLayerHeight * vSpaceFactor);

        int hostY = baseY;
        int taskY = baseY + layerHeight;
        int entryY = baseY + 2 * layerHeight;
        int activityBaseY = baseY + 3 * layerHeight;

        // Assign final positions
        for (ModelVertex v : hosts) {
            layout.setLocation(v, new Point2D.Double(xPositions.get(v), hostY));
        }
        for (ModelVertex v : tasks) {
            layout.setLocation(v, new Point2D.Double(xPositions.get(v), taskY));
        }
        for (ModelVertex v : entries) {
            layout.setLocation(v, new Point2D.Double(xPositions.get(v), entryY));
        }
        for (ModelVertex v : activities) {
            int sublayer = activitySublayer.get(v);
            int actY = activityBaseY + sublayer * (layerHeight / 2);
            layout.setLocation(v, new Point2D.Double(xPositions.get(v), actY));
        }

        return layout;
    }

    private Map<ModelVertex, Integer> assignActivitySublayers(List<ModelVertex> activities) {
        Map<ModelVertex, Integer> sublayers = new HashMap<ModelVertex, Integer>();
        Map<ModelVertex, Set<ModelVertex>> predecessors = new HashMap<ModelVertex, Set<ModelVertex>>();

        for (ModelVertex act : activities) {
            sublayers.put(act, 0);
            predecessors.put(act, new HashSet<ModelVertex>());
        }

        Set<ModelVertex> activitySet = new HashSet<ModelVertex>(activities);
        for (ModelEdge edge : graph.getEdges()) {
            if (edge.getType() == EdgeType.PRECEDENCE) {
                ModelVertex source = graph.getSource(edge);
                ModelVertex dest = graph.getDest(edge);
                if (activitySet.contains(source) && activitySet.contains(dest)) {
                    predecessors.get(dest).add(source);
                }
            }
        }

        boolean changed = true;
        while (changed) {
            changed = false;
            for (ModelVertex act : activities) {
                int maxPredSublayer = -1;
                for (ModelVertex pred : predecessors.get(act)) {
                    maxPredSublayer = Math.max(maxPredSublayer, sublayers.get(pred));
                }
                int newSublayer = maxPredSublayer + 1;
                if (newSublayer > sublayers.get(act)) {
                    sublayers.put(act, newSublayer);
                    changed = true;
                }
            }
        }

        return sublayers;
    }

    private void orderByBarycenterLQN(List<ModelVertex> vertices,
                                       final Map<ModelVertex, ModelVertex> parentMap,
                                       final Map<ModelVertex, Double> xPositions) {
        Collections.sort(vertices, new Comparator<ModelVertex>() {
            @Override
            public int compare(ModelVertex a, ModelVertex b) {
                ModelVertex parentA = parentMap.get(a);
                ModelVertex parentB = parentMap.get(b);
                double posA = parentA != null && xPositions.containsKey(parentA)
                    ? xPositions.get(parentA) : Double.MAX_VALUE;
                double posB = parentB != null && xPositions.containsKey(parentB)
                    ? xPositions.get(parentB) : Double.MAX_VALUE;
                int result = Double.compare(posA, posB);
                if (result == 0) {
                    return a.getName().compareTo(b.getName());
                }
                return result;
            }
        });
    }

    private double computeBarycenter(List<ModelVertex> vertices, Map<ModelVertex, Double> xPositions) {
        if (vertices == null || vertices.isEmpty()) {
            return Double.MAX_VALUE;
        }
        double sum = 0;
        int count = 0;
        for (ModelVertex v : vertices) {
            if (xPositions.containsKey(v)) {
                sum += xPositions.get(v);
                count++;
            }
        }
        return count > 0 ? sum / count : Double.MAX_VALUE;
    }

    // ========================================================================
    // VISUAL TRANSFORMERS
    // ========================================================================

    /**
     * Returns a transformer for vertex fill colors based on color theme.
     * Theme 0=Color, 1=Multicolor, 2=Greyscale, 3=Black/White
     */
    public Function<ModelVertex, Paint> getVertexPaintTransformer(final ModelVertex[] selectedVertex, final int[] colorTheme) {
        return new Function<ModelVertex, Paint>() {
            @Override
            public Paint apply(ModelVertex v) {
                int theme = colorTheme[0];

                // Black/White theme: all white (sink black)
                if (theme == 3) {
                    if (v.getType() == NodeType.SINK) {
                        return Color.BLACK;
                    }
                    return Color.WHITE;
                }

                // Greyscale theme
                if (theme == 2) {
                    if (v == selectedVertex[0]) {
                        return new Color(200, 200, 200);  // Light grey for selected
                    }
                    if (isLayeredNetwork) {
                        switch (v.getType()) {
                            case HOST: return new Color(80, 80, 80);      // Dark grey
                            case TASK: return new Color(120, 120, 120);   // Medium-dark grey
                            case ENTRY: return new Color(160, 160, 160);  // Medium grey
                            case ACTIVITY: return new Color(200, 200, 200); // Light grey
                            default: return new Color(180, 180, 180);
                        }
                    } else {
                        if (v.getType() == NodeType.SINK) {
                            return Color.BLACK;
                        }
                        switch (v.getType()) {
                            case SOURCE: return new Color(240, 240, 240);   // Very light grey
                            case QUEUE: return new Color(180, 180, 180);    // Medium grey
                            case DELAY: return new Color(200, 200, 200);    // Light grey
                            case FORK:
                            case JOIN: return new Color(140, 140, 140);     // Medium-dark grey
                            case CLASSSWITCH: return new Color(160, 160, 160);
                            case CACHE: return new Color(170, 170, 170);
                            case ROUTER: return new Color(150, 150, 150);
                            default: return new Color(190, 190, 190);
                        }
                    }
                }

                // Multicolor theme: different color per element type
                if (theme == 1) {
                    if (v == selectedVertex[0]) {
                        return new Color(255, 220, 50);  // Yellow for selected
                    }
                    if (isLayeredNetwork) {
                        switch (v.getType()) {
                            case HOST: return Color.GRAY;
                            case TASK: return Color.RED;
                            case ENTRY: return new Color(255, 255, 0);      // Yellow
                            case ACTIVITY: return Color.WHITE;
                            default: return new Color(200, 200, 200);
                        }
                    } else {
                        // Network: same as Color theme
                        switch (v.getType()) {
                            case SOURCE: return new Color(144, 238, 144);   // Light green
                            case SINK: return Color.BLACK;
                            case QUEUE: return new Color(135, 206, 250);    // Light blue
                            case DELAY: return new Color(255, 218, 185);    // Peach
                            case FORK:
                            case JOIN: return new Color(255, 182, 193);     // Light pink
                            case CLASSSWITCH: return new Color(221, 160, 221); // Plum
                            case CACHE: return new Color(176, 224, 230);    // Powder blue
                            case ROUTER: return new Color(240, 230, 140);   // Khaki
                            default: return Color.WHITE;
                        }
                    }
                }

                // Color theme (default, theme == 0)
                if (isLayeredNetwork) {
                    // LQN: red fill, yellow when selected
                    if (v == selectedVertex[0]) {
                        return new Color(255, 220, 50);  // Yellow for selected
                    }
                    return new Color(220, 60, 60);  // Red for all LQN elements
                } else {
                    // Network: colored by node type
                    switch (v.getType()) {
                        case SOURCE: return Color.WHITE;
                        case SINK: return Color.BLACK;
                        case QUEUE: return Color.WHITE;
                        case DELAY: return new Color(255, 218, 185);    // Peach
                        case FORK:
                        case JOIN:
                        case ROUTER: return Color.GRAY;
                        case CLASSSWITCH: return new Color(221, 160, 221); // Plum
                        case CACHE: return new Color(176, 224, 230);    // Powder blue
                        default: return Color.WHITE;
                    }
                }
            }
        };
    }

    /**
     * Returns a transformer for vertex outline colors.
     */
    public Function<ModelVertex, Paint> getVertexOutlineTransformer() {
        return new Function<ModelVertex, Paint>() {
            @Override
            public Paint apply(ModelVertex v) {
                return Color.BLACK;
            }
        };
    }

    /**
     * Returns a transformer for vertex shapes based on node type.
     */
    public Function<ModelVertex, Shape> getVertexShapeTransformer() {
        return new Function<ModelVertex, Shape>() {
            @Override
            public Shape apply(ModelVertex v) {
                switch (v.getType()) {
                    // Network types
                    case SOURCE:
                        return new Ellipse2D.Double(-6, -6, 12, 12);
                    case SINK:
                        return new Ellipse2D.Double(-6, -6, 12, 12);
                    case QUEUE:
                        // TikZ-style queue: buffer rectangle + server circle
                        GeneralPath queue = new GeneralPath();
                        int bufferWidth = 24;
                        int bufferHeight = 16;
                        int serverRadius = 8;
                        // Buffer starts at left, server circle at right
                        int bufferLeft = -bufferWidth - serverRadius;
                        int bufferTop = -bufferHeight / 2;
                        // Draw buffer outline
                        queue.moveTo(bufferLeft, bufferTop);
                        queue.lineTo(bufferLeft + bufferWidth, bufferTop);
                        queue.lineTo(bufferLeft + bufferWidth, bufferTop + bufferHeight);
                        queue.lineTo(bufferLeft, bufferTop + bufferHeight);
                        queue.closePath();
                        // Draw server circle (right side)
                        queue.append(new Ellipse2D.Double(-serverRadius, -serverRadius, serverRadius * 2, serverRadius * 2), false);
                        return queue;
                    case DELAY:
                        return new Ellipse2D.Double(-12, -12, 24, 24);
                    case FORK:
                    case JOIN:
                        Polygon diamond = new Polygon();
                        diamond.addPoint(0, -7);
                        diamond.addPoint(7, 0);
                        diamond.addPoint(0, 7);
                        diamond.addPoint(-7, 0);
                        return diamond;
                    case CLASSSWITCH:
                        Polygon hexagon = new Polygon();
                        for (int i = 0; i < 6; i++) {
                            double angle = Math.PI / 3 * i - Math.PI / 2;
                            hexagon.addPoint((int) (12 * Math.cos(angle)), (int) (12 * Math.sin(angle)));
                        }
                        return hexagon;
                    case CACHE:
                        return new Rectangle2D.Double(-14, -12, 28, 24);
                    case ROUTER:
                        Polygon octagon = new Polygon();
                        for (int i = 0; i < 8; i++) {
                            double angle = Math.PI / 4 * i - Math.PI / 8;
                            octagon.addPoint((int) (6 * Math.cos(angle)), (int) (6 * Math.sin(angle)));
                        }
                        return octagon;

                    // LQN types
                    case HOST:
                        Polygon pyramid = new Polygon();
                        pyramid.addPoint(-8, -8);
                        pyramid.addPoint(8, -8);
                        pyramid.addPoint(0, 8);
                        return pyramid;
                    case TASK:
                        Polygon parallelogram = new Polygon();
                        parallelogram.addPoint(-5, -6);
                        parallelogram.addPoint(10, -6);
                        parallelogram.addPoint(5, 6);
                        parallelogram.addPoint(-10, 6);
                        return parallelogram;
                    case ENTRY:
                        return new Rectangle2D.Double(-9, -6, 18, 12);
                    case ACTIVITY:
                        return new Ellipse2D.Double(-8, -8, 16, 16);

                    default:
                        return new Ellipse2D.Double(-12, -12, 24, 24);
                }
            }
        };
    }

    /**
     * Returns a transformer for edge colors.
     */
    public Function<ModelEdge, Paint> getEdgePaintTransformer() {
        return new Function<ModelEdge, Paint>() {
            @Override
            public Paint apply(ModelEdge e) {
                return Color.BLACK;
            }
        };
    }

    /**
     * Returns a transformer for edge strokes based on edge type.
     */
    public Function<ModelEdge, Stroke> getEdgeStrokeTransformer() {
        return new Function<ModelEdge, Stroke>() {
            @Override
            public Stroke apply(ModelEdge e) {
                switch (e.getType()) {
                    case ROUTING:
                    case PARENT_CHILD:
                    case BOUND_TO:
                    case PRECEDENCE:
                        return new BasicStroke(1.5f);
                    case SYNCH_CALL:
                        return new BasicStroke(1.0f, BasicStroke.CAP_BUTT, BasicStroke.JOIN_MITER,
                                               10.0f, new float[]{8.0f, 4.0f}, 0.0f);
                    case ASYNCH_CALL:
                        return new BasicStroke(1.0f, BasicStroke.CAP_BUTT, BasicStroke.JOIN_MITER,
                                               10.0f, new float[]{8.0f, 4.0f, 2.0f, 4.0f}, 0.0f);
                    case FORWARDING:
                        return new BasicStroke(1.0f, BasicStroke.CAP_ROUND, BasicStroke.JOIN_ROUND,
                                               10.0f, new float[]{2.0f, 4.0f}, 0.0f);
                    default:
                        return new BasicStroke(1.0f);
                }
            }
        };
    }

    /**
     * Returns the underlying JUNG graph.
     */
    public Graph<ModelVertex, ModelEdge> getGraph() {
        return graph;
    }

    // ========================================================================
    // HELPER METHODS
    // ========================================================================

    // Store background images per visualization viewer
    private static Map<VisualizationViewer<?, ?>, BufferedImage> backgroundImages =
        new HashMap<VisualizationViewer<?, ?>, BufferedImage>();
    private static Map<VisualizationViewer<?, ?>, ComponentAdapter> resizeListeners =
        new HashMap<VisualizationViewer<?, ?>, ComponentAdapter>();

    /**
     * Sets a background image for the visualization viewer.
     */
    public static void setBackgroundImage(final VisualizationViewer<?, ?> vv, BufferedImage img) {
        backgroundImages.put(vv, img);

        if (!resizeListeners.containsKey(vv)) {
            vv.addPreRenderPaintable(new VisualizationViewer.Paintable() {
                @Override
                public void paint(Graphics g) {
                    BufferedImage bgImg = backgroundImages.get(vv);
                    if (bgImg != null) {
                        Graphics2D g2d = (Graphics2D) g;
                        g2d.setRenderingHint(RenderingHints.KEY_INTERPOLATION,
                            RenderingHints.VALUE_INTERPOLATION_BILINEAR);
                        Dimension size = vv.getSize();
                        g2d.drawImage(bgImg, 0, 0, size.width, size.height, null);
                    }
                }
                @Override
                public boolean useTransform() { return false; }
            });

            ComponentAdapter resizeListener = new ComponentAdapter() {
                @Override
                public void componentResized(ComponentEvent e) {
                    vv.repaint();
                }
            };
            vv.addComponentListener(resizeListener);
            resizeListeners.put(vv, resizeListener);
        }

        vv.repaint();
    }

    /**
     * Clears the background image from the visualization viewer.
     */
    public static void clearBackgroundImage(VisualizationViewer<?, ?> vv) {
        backgroundImages.remove(vv);
        vv.repaint();
    }

    /**
     * Fits the graph to the current view.
     */
    private void fitGraphToView(VisualizationViewer<ModelVertex, ModelEdge> vv,
            Graph<ModelVertex, ModelEdge> graph, JSlider zoomSlider) {
        vv.getRenderContext().getMultiLayerTransformer()
            .getTransformer(edu.uci.ics.jung.visualization.Layer.VIEW).setToIdentity();
        vv.getRenderContext().getMultiLayerTransformer()
            .getTransformer(edu.uci.ics.jung.visualization.Layer.LAYOUT).setToIdentity();

        Layout<ModelVertex, ModelEdge> layout = vv.getGraphLayout();
        double minX = Double.MAX_VALUE, minY = Double.MAX_VALUE;
        double maxX = Double.MIN_VALUE, maxY = Double.MIN_VALUE;
        for (ModelVertex v : graph.getVertices()) {
            Point2D p = layout.apply(v);
            if (p != null) {
                minX = Math.min(minX, p.getX());
                minY = Math.min(minY, p.getY());
                maxX = Math.max(maxX, p.getX());
                maxY = Math.max(maxY, p.getY());
            }
        }

        if (minX < Double.MAX_VALUE) {
            double margin = 20;
            minX -= margin;
            minY -= margin;
            maxX += margin;
            maxY += margin;

            double graphWidth = maxX - minX;
            double graphHeight = maxY - minY;
            double graphCenterX = (minX + maxX) / 2;

            Dimension viewSize = vv.getSize();
            double viewCenterX = viewSize.width / 2.0;

            double scaleX = viewSize.width / graphWidth;
            double scaleY = viewSize.height / graphHeight;
            double scale = Math.min(scaleX, scaleY);
            scale = Math.min(scale, 2.0);
            scale = Math.max(scale, 0.25);

            edu.uci.ics.jung.visualization.transform.MutableTransformer viewTransformer =
                vv.getRenderContext().getMultiLayerTransformer()
                    .getTransformer(edu.uci.ics.jung.visualization.Layer.VIEW);
            viewTransformer.scale(scale, scale, new Point2D.Double(0, 0));
            double offsetX = viewCenterX - graphCenterX * scale;
            double topMargin = 15;
            double offsetY = topMargin - minY * scale;
            viewTransformer.translate(offsetX / scale, offsetY / scale);

            int sliderValue = (int)(Math.log(scale) / Math.log(2) * 50);
            sliderValue = Math.max(-100, Math.min(100, sliderValue));
            zoomSlider.setValue(sliderValue);
        }
        vv.repaint();
    }

    // ========================================================================
    // PROPERTIES PANEL
    // ========================================================================

    /**
     * Builds the properties panel for a vertex.
     */
    private void buildPropertiesPanel(JPanel panel, ModelVertex vertex,
            Graph<ModelVertex, ModelEdge> graph, VisualizationViewer<ModelVertex, ModelEdge> vv) {

        addSectionHeader(panel, vertex.getType().toString());
        addReadOnlyField(panel, "Name", vertex.getName());
        addReadOnlyField(panel, "Index", String.valueOf(vertex.getIndex()));

        Object element = vertex.getElement();
        if (element != null) {
            if (isLayeredNetwork) {
                if (element instanceof Host) {
                    buildHostProperties(panel, (Host) element);
                } else if (element instanceof Task) {
                    buildTaskProperties(panel, (Task) element);
                } else if (element instanceof Entry) {
                    buildEntryProperties(panel, (Entry) element);
                } else if (element instanceof Activity) {
                    buildActivityProperties(panel, (Activity) element);
                }
            } else {
                if (element instanceof Node) {
                    buildNetworkNodeProperties(panel, (Node) element);
                }
            }
        }

        // Add connections section
        addConnectionsSection(panel, vertex, graph);
    }

    private void buildNetworkNodeProperties(JPanel panel, Node node) {
        addSectionHeader(panel, "Configuration");

        if (node instanceof Station) {
            Station station = (Station) node;
            addReadOnlyField(panel, "Servers", String.valueOf(station.getNumberOfServers()));

            if (node instanceof Queue) {
                Queue queue = (Queue) node;
                SchedStrategy sched = queue.getSchedStrategy();
                if (sched != null) {
                    addReadOnlyField(panel, "Scheduling", sched.toString());
                }
            }
        }

        if (node instanceof Source) {
            addSectionHeader(panel, "Arrival Rates");
            List<JobClass> classes = networkModel.getClasses();
            for (JobClass jc : classes) {
                try {
                    Source source = (Source) node;
                    Object dist = source.getArrivalDistribution(jc);
                    if (dist != null) {
                        addReadOnlyField(panel, jc.getName(), dist.toString());
                    }
                } catch (Exception e) {
                    // Ignore
                }
            }
        }
    }

    private void buildHostProperties(JPanel panel, final Host host) {
        addSectionHeader(panel, "Configuration");
        addReadOnlyField(panel, "Multiplicity", String.valueOf(host.getMultiplicity()));
        addReadOnlyField(panel, "Replication", String.valueOf(host.getReplication()));
        addReadOnlyField(panel, "Scheduling", host.getScheduling().toString());
        addReadOnlyField(panel, "Quantum", String.format("%.4f", host.getQuantum()));
        addReadOnlyField(panel, "Speed Factor", String.format("%.4f", host.getSpeedFactor()));
    }

    private void buildTaskProperties(JPanel panel, final Task task) {
        addSectionHeader(panel, "Configuration");
        addReadOnlyField(panel, "Multiplicity", String.valueOf(task.getMultiplicity()));
        addReadOnlyField(panel, "Replication", String.valueOf(task.getReplication()));
        addReadOnlyField(panel, "Scheduling", task.getScheduling().toString());
        addReadOnlyField(panel, "Think Time Mean", String.format("%.4f", task.getThinkTimeMean()));

        if (task.getParent() != null) {
            addSectionHeader(panel, "Parent");
            addReadOnlyField(panel, "Processor", task.getParent().getName());
        }
    }

    private void buildEntryProperties(JPanel panel, final Entry entry) {
        addSectionHeader(panel, "Configuration");

        if (entry.getArrival() != null) {
            addReadOnlyField(panel, "Arrival Dist", entry.getArrival().getClass().getSimpleName());
            addReadOnlyField(panel, "Arrival Mean", String.format("%.4f", entry.getArrival().getMean()));
        } else {
            addReadOnlyField(panel, "Arrival", "None (closed)");
        }

        if (entry.getParent() != null) {
            addSectionHeader(panel, "Parent");
            addReadOnlyField(panel, "Task", entry.getParent().getName());
        }
    }

    private void buildActivityProperties(JPanel panel, final Activity activity) {
        addSectionHeader(panel, "Host Demand");
        addReadOnlyField(panel, "Mean", String.format("%.4f", activity.getHostDemandMean()));
        addReadOnlyField(panel, "SCV", String.format("%.4f", activity.getHostDemandSCV()));

        if (activity.getBoundToEntry() != null && !activity.getBoundToEntry().isEmpty()) {
            addReadOnlyField(panel, "Bound To", activity.getBoundToEntry());
        }

        if (activity.getParent() != null) {
            addSectionHeader(panel, "Parent");
            addReadOnlyField(panel, "Task", activity.getParent().getName());
        }
    }

    private void addConnectionsSection(JPanel panel, ModelVertex vertex, Graph<ModelVertex, ModelEdge> graph) {
        List<String> incoming = new ArrayList<String>();
        List<String> outgoing = new ArrayList<String>();

        for (ModelEdge edge : graph.getEdges()) {
            if (graph.getDest(edge).equals(vertex)) {
                ModelVertex source = graph.getSource(edge);
                incoming.add(source.getName() + " (" + edge.getType() + ")");
            }
            if (graph.getSource(edge).equals(vertex)) {
                ModelVertex dest = graph.getDest(edge);
                outgoing.add(dest.getName() + " (" + edge.getType() + ")");
            }
        }

        if (!incoming.isEmpty() || !outgoing.isEmpty()) {
            addSectionHeader(panel, "Connections");

            for (String conn : incoming) {
                addReadOnlyField(panel, "\u2190 From", conn);
            }

            for (String conn : outgoing) {
                addReadOnlyField(panel, "\u2192 To", conn);
            }
        }
    }

    private static void addSectionHeader(JPanel panel, String title) {
        panel.add(javax.swing.Box.createVerticalStrut(12));
        JLabel header = new JLabel(title);
        header.setFont(new Font("Segoe UI", Font.BOLD, 12));
        header.setForeground(new Color(80, 80, 80));
        header.setAlignmentX(Component.LEFT_ALIGNMENT);
        header.setBorder(BorderFactory.createEmptyBorder(0, 0, 4, 0));
        panel.add(header);
        JSeparator sep = new JSeparator();
        sep.setMaximumSize(new Dimension(Integer.MAX_VALUE, 1));
        panel.add(sep);
        panel.add(javax.swing.Box.createVerticalStrut(4));
    }

    private static void addReadOnlyField(JPanel panel, String label, String value) {
        JPanel row = new JPanel(new BorderLayout(8, 0));
        row.setBackground(new Color(250, 250, 252));
        row.setMaximumSize(new Dimension(Integer.MAX_VALUE, 28));
        row.setAlignmentX(Component.LEFT_ALIGNMENT);

        JLabel lbl = new JLabel(label + ":");
        lbl.setFont(new Font("Segoe UI", Font.PLAIN, 12));
        lbl.setForeground(new Color(100, 100, 100));
        lbl.setPreferredSize(new Dimension(100, 24));
        row.add(lbl, BorderLayout.WEST);

        JLabel val = new JLabel(value);
        val.setFont(new Font("Segoe UI", Font.PLAIN, 12));
        val.setForeground(new Color(40, 40, 40));
        row.add(val, BorderLayout.CENTER);

        panel.add(row);
    }

    // ========================================================================
    // SOLVER EXECUTION (LQN ONLY)
    // ========================================================================

    private static void runSolver(final LayeredNetwork model, final JButton playButton,
            final JTextArea consoleArea, final JPanel resultsPanel, final JTabbedPane statusPane) {

        final Icon playIcon = playButton.getIcon();
        playButton.setIcon(createStopIcon(27, new Color(231, 76, 60)));
        playButton.setToolTipText("Solver running...");
        playButton.setEnabled(false);

        consoleArea.setText("");
        statusPane.setSelectedIndex(0);

        appendToConsole(consoleArea, "=== LINE Solver ===\n");
        appendToConsole(consoleArea, "Model: " + model.getName() + "\n");
        appendToConsole(consoleArea, "Starting solver...\n\n");

        SwingWorker<LayeredNetworkAvgTable, String> worker = new SwingWorker<LayeredNetworkAvgTable, String>() {
            @Override
            protected LayeredNetworkAvgTable doInBackground() throws Exception {
                publish("Creating SolverLN instance...\n");
                SolverLN solver = new SolverLN(model);
                publish("Running analysis...\n");
                long startTime = System.currentTimeMillis();
                solver.runAnalyzer();
                long elapsed = System.currentTimeMillis() - startTime;
                publish("Analysis completed in " + elapsed + " ms\n\n");
                return (LayeredNetworkAvgTable) solver.getAvgTable();
            }

            @Override
            protected void process(java.util.List<String> chunks) {
                for (String chunk : chunks) {
                    appendToConsole(consoleArea, chunk);
                }
            }

            @Override
            protected void done() {
                try {
                    LayeredNetworkAvgTable avgTable = get();
                    if (avgTable != null) {
                        appendToConsole(consoleArea, "Solver completed successfully.\n");
                        statusPane.setSelectedIndex(1);
                        appendToConsole(consoleArea, "\nResults ready. See 'Results' tab.\n");
                    } else {
                        appendToConsole(consoleArea, "\nERROR: Solver returned null result.\n");
                    }
                } catch (Exception e) {
                    appendToConsole(consoleArea, "\nERROR: " + e.getMessage() + "\n");
                } finally {
                    playButton.setIcon(playIcon);
                    playButton.setToolTipText("Run LINE solver on this model");
                    playButton.setEnabled(true);
                }
            }
        };

        worker.execute();
    }

    private static void appendToConsole(final JTextArea consoleArea, final String text) {
        SwingUtilities.invokeLater(new Runnable() {
            @Override
            public void run() {
                consoleArea.append(text);
                consoleArea.setCaretPosition(consoleArea.getDocument().getLength());
            }
        });
    }

    // ========================================================================
    // EXPORT METHODS
    // ========================================================================

    private void exportToPNG(VisualizationViewer<ModelVertex, ModelEdge> vv, java.io.File file) throws Exception {
        Dimension size = vv.getSize();
        BufferedImage image = new BufferedImage(size.width, size.height, BufferedImage.TYPE_INT_ARGB);
        Graphics2D g2 = image.createGraphics();
        g2.setRenderingHint(RenderingHints.KEY_ANTIALIASING, RenderingHints.VALUE_ANTIALIAS_ON);
        vv.paint(g2);
        g2.dispose();
        ImageIO.write(image, "PNG", file);
    }

    private void exportToTikZ(VisualizationViewer<ModelVertex, ModelEdge> vv,
            Graph<ModelVertex, ModelEdge> graph, java.io.File file) throws Exception {
        Layout<ModelVertex, ModelEdge> layout = vv.getGraphLayout();

        StringBuilder tikz = new StringBuilder();
        tikz.append("% TikZ export from ModelVisualizer\n");
        tikz.append("\\begin{tikzpicture}[scale=0.02,\n");

        if (isLayeredNetwork) {
            tikz.append("  host/.style={draw, fill=red!70, regular polygon, regular polygon sides=3, minimum size=16pt},\n");
            tikz.append("  task/.style={draw, fill=red!70, trapezium, minimum size=12pt},\n");
            tikz.append("  entry/.style={draw, fill=red!70, rectangle, minimum size=10pt},\n");
            tikz.append("  activity/.style={draw, fill=red!70, circle, minimum size=10pt},\n");
            tikz.append("  synch/.style={->, dashed},\n");
            tikz.append("  asynch/.style={->, dashdotted},\n");
            tikz.append("  forward/.style={->, dotted},\n");
            tikz.append("  parent/.style={->}]\n\n");
        } else {
            tikz.append("  source/.style={draw, fill=white, circle, minimum size=12pt},\n");
            tikz.append("  sink/.style={draw, fill=black, circle, minimum size=12pt},\n");
            tikz.append("  queue/.style={draw, fill=white, rectangle, minimum size=12pt},\n");
            tikz.append("  delay/.style={draw, fill=white, circle, minimum size=12pt},\n");
            tikz.append("  fork/.style={draw, fill=white, diamond, minimum size=12pt},\n");
            tikz.append("  join/.style={draw, fill=white, diamond, minimum size=12pt},\n");
            tikz.append("  routing/.style={->}]\n\n");
        }

        tikz.append("% Nodes\n");
        for (ModelVertex v : graph.getVertices()) {
            Point2D p = layout.apply(v);
            if (p != null) {
                String style = v.getType().toString().toLowerCase();
                tikz.append(String.format("\\node[%s] (%s) at (%.1f, %.1f) {};\n",
                    style, "n" + v.getIndex(), p.getX(), -p.getY()));
            }
        }

        tikz.append("\n% Edges\n");
        for (ModelEdge edge : graph.getEdges()) {
            ModelVertex source = graph.getSource(edge);
            ModelVertex dest = graph.getDest(edge);
            String style = "";
            switch (edge.getType()) {
                case SYNCH_CALL: style = "synch"; break;
                case ASYNCH_CALL: style = "asynch"; break;
                case FORWARDING: style = "forward"; break;
                case ROUTING: style = "routing"; break;
                default: style = "parent"; break;
            }
            tikz.append(String.format("\\draw[%s] (n%d) -- (n%d);\n",
                style, source.getIndex(), dest.getIndex()));
        }

        tikz.append("\n\\end{tikzpicture}\n");

        java.io.FileWriter writer = new java.io.FileWriter(file);
        writer.write(tikz.toString());
        writer.close();
    }

    // ========================================================================
    // ICON CREATION METHODS
    // ========================================================================

    private static void configureHighQualityRendering(Graphics2D g2d) {
        g2d.setRenderingHint(RenderingHints.KEY_ANTIALIASING, RenderingHints.VALUE_ANTIALIAS_ON);
        g2d.setRenderingHint(RenderingHints.KEY_RENDERING, RenderingHints.VALUE_RENDER_QUALITY);
        g2d.setRenderingHint(RenderingHints.KEY_STROKE_CONTROL, RenderingHints.VALUE_STROKE_PURE);
    }

    private static Icon createCenterIcon(final int size, final Color color) {
        return new Icon() {
            @Override
            public void paintIcon(Component c, Graphics g, int x, int y) {
                Graphics2D g2d = (Graphics2D) g.create();
                configureHighQualityRendering(g2d);
                g2d.translate(x, y);
                double s = size;
                float strokeW = (float)(s * 0.07);
                g2d.setColor(color);
                g2d.setStroke(new BasicStroke(strokeW, BasicStroke.CAP_ROUND, BasicStroke.JOIN_ROUND));
                double cx = s * 0.5, cy = s * 0.5, r = s * 0.32;
                g2d.draw(new java.awt.geom.Ellipse2D.Double(cx - r, cy - r, r * 2, r * 2));
                double lineLen = s * 0.18;
                g2d.draw(new java.awt.geom.Line2D.Double(cx, cy - r - lineLen * 0.5, cx, cy - r + lineLen));
                g2d.draw(new java.awt.geom.Line2D.Double(cx, cy + r - lineLen, cx, cy + r + lineLen * 0.5));
                g2d.draw(new java.awt.geom.Line2D.Double(cx - r - lineLen * 0.5, cy, cx - r + lineLen, cy));
                g2d.draw(new java.awt.geom.Line2D.Double(cx + r - lineLen, cy, cx + r + lineLen * 0.5, cy));
                double dotR = s * 0.06;
                g2d.fill(new java.awt.geom.Ellipse2D.Double(cx - dotR, cy - dotR, dotR * 2, dotR * 2));
                g2d.dispose();
            }
            @Override public int getIconWidth() { return size; }
            @Override public int getIconHeight() { return size; }
        };
    }

    private static Icon createRefreshIcon(final int size, final Color color) {
        return new Icon() {
            @Override
            public void paintIcon(Component c, Graphics g, int x, int y) {
                Graphics2D g2d = (Graphics2D) g.create();
                configureHighQualityRendering(g2d);
                g2d.translate(x, y);
                double s = size, cx = s / 2, cy = s / 2, radius = s * 0.32;
                float strokeW = (float)(s * 0.08);
                g2d.setColor(color);
                g2d.setStroke(new BasicStroke(strokeW, BasicStroke.CAP_ROUND, BasicStroke.JOIN_ROUND));
                g2d.draw(new java.awt.geom.Arc2D.Double(cx - radius, cy - radius, radius * 2, radius * 2, 45, 135, java.awt.geom.Arc2D.OPEN));
                g2d.draw(new java.awt.geom.Arc2D.Double(cx - radius, cy - radius, radius * 2, radius * 2, 225, 135, java.awt.geom.Arc2D.OPEN));
                double arrowSize = s * 0.15;
                double a1x = cx + radius * Math.cos(Math.toRadians(45));
                double a1y = cy - radius * Math.sin(Math.toRadians(45));
                GeneralPath arrow1 = new GeneralPath();
                arrow1.moveTo(a1x, a1y - arrowSize * 0.8);
                arrow1.lineTo(a1x + arrowSize * 0.7, a1y + arrowSize * 0.3);
                arrow1.lineTo(a1x - arrowSize * 0.5, a1y + arrowSize * 0.3);
                arrow1.closePath();
                g2d.fill(arrow1);
                double a2x = cx + radius * Math.cos(Math.toRadians(225));
                double a2y = cy - radius * Math.sin(Math.toRadians(225));
                GeneralPath arrow2 = new GeneralPath();
                arrow2.moveTo(a2x, a2y + arrowSize * 0.8);
                arrow2.lineTo(a2x - arrowSize * 0.7, a2y - arrowSize * 0.3);
                arrow2.lineTo(a2x + arrowSize * 0.5, a2y - arrowSize * 0.3);
                arrow2.closePath();
                g2d.fill(arrow2);
                g2d.dispose();
            }
            @Override public int getIconWidth() { return size; }
            @Override public int getIconHeight() { return size; }
        };
    }

    private static Icon createMapIcon(final int size, final Color color) {
        return new Icon() {
            @Override
            public void paintIcon(Component c, Graphics g, int x, int y) {
                Graphics2D g2d = (Graphics2D) g.create();
                configureHighQualityRendering(g2d);
                g2d.translate(x, y);
                double s = size, m = s * 0.12;
                float strokeW = (float)(s * 0.05);
                g2d.setColor(color);
                g2d.setStroke(new BasicStroke(strokeW, BasicStroke.CAP_ROUND, BasicStroke.JOIN_ROUND));
                g2d.draw(new java.awt.geom.RoundRectangle2D.Double(m, m, s - 2 * m, s - 2 * m, s * 0.1, s * 0.1));
                GeneralPath mountain = new GeneralPath();
                mountain.moveTo(s * 0.18, s * 0.75);
                mountain.lineTo(s * 0.35, s * 0.48);
                mountain.lineTo(s * 0.45, s * 0.58);
                mountain.lineTo(s * 0.62, s * 0.35);
                mountain.lineTo(s * 0.82, s * 0.75);
                mountain.closePath();
                g2d.fill(mountain);
                double sunR = s * 0.09;
                g2d.fill(new java.awt.geom.Ellipse2D.Double(s * 0.22, s * 0.22, sunR * 2, sunR * 2));
                g2d.dispose();
            }
            @Override public int getIconWidth() { return size; }
            @Override public int getIconHeight() { return size; }
        };
    }

    private static Icon createPlayIcon(final int size, final Color color) {
        return new Icon() {
            @Override
            public void paintIcon(Component c, Graphics g, int x, int y) {
                Graphics2D g2d = (Graphics2D) g.create();
                configureHighQualityRendering(g2d);
                g2d.translate(x, y);
                double s = size;
                GeneralPath triangle = new GeneralPath();
                triangle.moveTo(s * 0.25, s * 0.15);
                triangle.lineTo(s * 0.80, s * 0.50);
                triangle.lineTo(s * 0.25, s * 0.85);
                triangle.closePath();
                g2d.setColor(color);
                g2d.fill(triangle);
                g2d.setColor(color.darker().darker());
                g2d.setStroke(new BasicStroke((float)(s * 0.03), BasicStroke.CAP_ROUND, BasicStroke.JOIN_ROUND));
                g2d.draw(triangle);
                g2d.dispose();
            }
            @Override public int getIconWidth() { return size; }
            @Override public int getIconHeight() { return size; }
        };
    }

    private static Icon createStopIcon(final int size, final Color color) {
        return new Icon() {
            @Override
            public void paintIcon(Component c, Graphics g, int x, int y) {
                Graphics2D g2d = (Graphics2D) g.create();
                configureHighQualityRendering(g2d);
                g2d.translate(x, y);
                double s = size, margin = s * 0.22;
                java.awt.geom.RoundRectangle2D square = new java.awt.geom.RoundRectangle2D.Double(
                    margin, margin, s - 2 * margin, s - 2 * margin, s * 0.08, s * 0.08);
                g2d.setColor(color);
                g2d.fill(square);
                g2d.setColor(color.darker().darker());
                g2d.setStroke(new BasicStroke((float)(s * 0.03), BasicStroke.CAP_ROUND, BasicStroke.JOIN_ROUND));
                g2d.draw(square);
                g2d.dispose();
            }
            @Override public int getIconWidth() { return size; }
            @Override public int getIconHeight() { return size; }
        };
    }

    private Icon createExportIcon(int size, Color color) {
        return new Icon() {
            @Override
            public void paintIcon(Component c, Graphics g, int x, int y) {
                Graphics2D g2 = (Graphics2D) g.create();
                g2.setRenderingHint(RenderingHints.KEY_ANTIALIASING, RenderingHints.VALUE_ANTIALIAS_ON);
                g2.setColor(color);
                g2.setStroke(new BasicStroke(2f, BasicStroke.CAP_ROUND, BasicStroke.JOIN_ROUND));
                int cx = x + size / 2, cy = y + size / 2;
                g2.drawLine(cx, cy - 6, cx, cy + 4);
                g2.drawLine(cx - 4, cy, cx, cy + 4);
                g2.drawLine(cx + 4, cy, cx, cy + 4);
                g2.drawLine(cx - 6, cy + 7, cx + 6, cy + 7);
                g2.dispose();
            }
            @Override public int getIconWidth() { return size; }
            @Override public int getIconHeight() { return size; }
        };
    }
}
