package jline.io;

import jline.lang.*;
import jline.lang.constant.*;
import jline.lang.nodes.*;
import jline.lang.processes.*;
import jline.solvers.*;
// Matrix import will be added when needed
import jline.solvers.mva.SolverMVA;
import jline.solvers.jmt.SolverJMT;
import jline.solvers.mam.SolverMAM;
import jline.solvers.nc.SolverNC;
import jline.solvers.ctmc.SolverCTMC;
import jline.solvers.fluid.SolverFluid;
import jline.solvers.auto.SolverAUTO;

import java.io.*;
import java.util.ArrayList;
import java.util.Arrays;
import java.util.HashMap;
import java.util.List;
import java.util.Map;
import java.util.regex.Matcher;
import java.util.regex.Pattern;
import java.util.stream.Collectors;

/**
 * REPL (Read-Eval-Print Loop) for LINE solver
 * Provides an interactive environment to define models and execute solvers
 */
public class REPL {
    
    private static final String PROMPT = "LINE> ";
    private static final String VERSION = "1.0.0";
    
    private Network currentNetwork;
    private Map<String, Node> nodes;
    private Map<String, JobClass> jobClasses;
    private Map<String, Object> variables;
    private boolean running;
    private BufferedReader reader;
    private PrintWriter writer;
    private RoutingMatrix routingMatrix;
    
    public REPL() {
        this(new InputStreamReader(System.in), new OutputStreamWriter(System.out));
    }
    
    public REPL(Reader input, Writer output) {
        this.reader = new BufferedReader(input);
        this.writer = new PrintWriter(output, true);
        this.nodes = new HashMap<>();
        this.jobClasses = new HashMap<>();
        this.variables = new HashMap<>();
        this.running = true;
    }
    
    /**
     * Start the REPL
     */
    public void start() {
        printWelcome();
        
        while (running) {
            try {
                writer.print(PROMPT);
                writer.flush();
                
                String line = reader.readLine();
                if (line == null) {
                    break;
                }
                
                line = line.trim();
                if (line.isEmpty()) {
                    continue;
                }
                
                processCommand(line);
                
            } catch (IOException e) {
                writer.println("Error reading input: " + e.getMessage());
            } catch (Exception e) {
                writer.println("Error: " + e.getMessage());
                if (e.getCause() != null) {
                    writer.println("Caused by: " + e.getCause().getMessage());
                }
            }
        }
        
        writer.println("Goodbye!");
    }
    
    /**
     * Process a single command
     */
    private void processCommand(String command) {
        String[] tokens = tokenize(command);
        if (tokens.length == 0) {
            return;
        }
        
        String cmd = tokens[0].toLowerCase();
        
        switch (cmd) {
            case "help":
            case "?":
                showHelp(tokens.length > 1 ? tokens[1] : null);
                break;
                
            case "quit":
            case "exit":
                running = false;
                break;
                
            case "new":
                if (tokens.length > 1 && tokens[1].equalsIgnoreCase("network")) {
                    createNetwork(tokens);
                } else {
                    writer.println("Usage: new network [name]");
                }
                break;
                
            case "add":
                if (tokens.length > 1) {
                    addElement(Arrays.copyOfRange(tokens, 1, tokens.length));
                } else {
                    writer.println("Usage: add <queue|source|sink|delay|class> <name> [options]");
                }
                break;
                
            case "set":
                if (tokens.length > 1) {
                    setProperty(Arrays.copyOfRange(tokens, 1, tokens.length));
                } else {
                    writer.println("Usage: set <property> <value>");
                }
                break;
                
            case "link":
                if (tokens.length >= 3) {
                    linkNodes(Arrays.copyOfRange(tokens, 1, tokens.length));
                } else {
                    writer.println("Usage: link <from> <to> [probability] [class]");
                }
                break;
                
            case "solve":
                if (currentNetwork == null) {
                    writer.println("No network defined. Use 'new network' first.");
                } else {
                    solve(tokens.length > 1 ? tokens[1] : "auto");
                }
                break;
                
            case "show":
                if (tokens.length > 1) {
                    show(tokens[1]);
                } else {
                    show("network");
                }
                break;
                
            case "save":
                if (tokens.length > 1) {
                    saveModel(tokens[1]);
                } else {
                    writer.println("Usage: save <filename>");
                }
                break;
                
            case "load":
                if (tokens.length > 1) {
                    loadModel(tokens[1]);
                } else {
                    writer.println("Usage: load <filename>");
                }
                break;
                
            case "clear":
                clear();
                break;
                
            case "let":
                if (tokens.length >= 3) {
                    defineVariable(Arrays.copyOfRange(tokens, 1, tokens.length));
                } else {
                    writer.println("Usage: let <name> = <value>");
                }
                break;
                
            default:
                writer.println("Unknown command: " + cmd);
                writer.println("Type 'help' for available commands.");
        }
    }
    
    /**
     * Create a new network
     */
    private void createNetwork(String[] tokens) {
        String name = tokens.length > 2 ? tokens[2] : "Network";
        currentNetwork = new Network(name);
        nodes.clear();
        jobClasses.clear();
        writer.println("Created new network: " + name);
    }
    
    /**
     * Add elements to the network
     */
    private void addElement(String[] tokens) {
        if (currentNetwork == null) {
            writer.println("No network defined. Use 'new network' first.");
            return;
        }
        
        if (tokens.length < 2) {
            writer.println("Usage: add <type> <name> [options]");
            return;
        }
        
        String type = tokens[0].toLowerCase();
        String name = tokens[1];
        
        switch (type) {
            case "queue":
                addQueue(name, Arrays.copyOfRange(tokens, 2, tokens.length));
                break;
                
            case "source":
                addSource(name, Arrays.copyOfRange(tokens, 2, tokens.length));
                break;
                
            case "sink":
                addSink(name);
                break;
                
            case "delay":
                addDelay(name, Arrays.copyOfRange(tokens, 2, tokens.length));
                break;
                
            case "class":
                addClass(name, Arrays.copyOfRange(tokens, 2, tokens.length));
                break;
                
            default:
                writer.println("Unknown element type: " + type);
        }
    }
    
    /**
     * Add a queue to the network
     */
    private void addQueue(String name, String[] options) {
        int servers = 1;
        SchedStrategy sched = SchedStrategy.FCFS;
        
        for (int i = 0; i < options.length - 1; i += 2) {
            String opt = options[i].toLowerCase();
            String val = options[i + 1];
            
            switch (opt) {
                case "-servers":
                case "-s":
                    servers = Integer.parseInt(val);
                    break;
                    
                case "-sched":
                case "-strategy":
                    sched = parseSchedStrategy(val);
                    break;
            }
        }
        
        Queue queue = new Queue(currentNetwork, name, sched);
        queue.setNumberOfServers(servers);
        nodes.put(name, queue);
        writer.println("Added queue: " + name + " (servers=" + servers + ", strategy=" + sched + ")");
    }
    
    /**
     * Add a source to the network
     */
    private void addSource(String name, String[] options) {
        Source source = new Source(currentNetwork, name);
        nodes.put(name, source);
        writer.println("Added source: " + name);
    }
    
    /**
     * Add a sink to the network
     */
    private void addSink(String name) {
        Sink sink = new Sink(currentNetwork, name);
        nodes.put(name, sink);
        writer.println("Added sink: " + name);
    }
    
    /**
     * Add a delay station to the network
     */
    private void addDelay(String name, String[] options) {
        Delay delay = new Delay(currentNetwork, name);
        nodes.put(name, delay);
        writer.println("Added delay: " + name);
    }
    
    /**
     * Add a job class to the network
     */
    private void addClass(String name, String[] options) {
        String type = "closed";
        int population = 1;
        String refNode = null;
        int priority = 0;
        
        for (int i = 0; i < options.length - 1; i += 2) {
            String opt = options[i].toLowerCase();
            String val = options[i + 1];
            
            switch (opt) {
                case "-type":
                case "-t":
                    type = val.toLowerCase();
                    break;
                    
                case "-population":
                case "-pop":
                case "-n":
                    population = Integer.parseInt(val);
                    break;
                    
                case "-ref":
                case "-reference":
                    refNode = val;
                    break;
                    
                case "-priority":
                case "-p":
                    priority = Integer.parseInt(val);
                    break;
            }
        }
        
        JobClass jobClass;
        if (type.equals("open")) {
            if (refNode == null || !nodes.containsKey(refNode)) {
                writer.println("Open class requires a valid source as reference node");
                return;
            }
            Node ref = nodes.get(refNode);
            if (!(ref instanceof Source)) {
                writer.println("Reference node for open class must be a source");
                return;
            }
            jobClass = new OpenClass(currentNetwork, name, priority);
            // OpenClass reference is set at construction, not after
        } else {
            if (refNode != null && nodes.containsKey(refNode)) {
                Node ref = nodes.get(refNode);
                if (ref instanceof Station) {
                    jobClass = new ClosedClass(currentNetwork, name, (double) population, (Station) ref, priority);
                } else {
                    writer.println("Reference node for closed class must be a station (queue, delay, etc.)");
                    return;
                }
            } else {
                // Find first station to use as reference
                Station refStation = null;
                for (Node node : nodes.values()) {
                    if (node instanceof Station) {
                        refStation = (Station) node;
                        break;
                    }
                }
                jobClass = new ClosedClass(currentNetwork, name, (double) population, refStation, priority);
            }
        }
        
        jobClasses.put(name, jobClass);
        writer.println("Added " + type + " class: " + name + 
                      (type.equals("closed") ? " (population=" + population + ")" : ""));
    }
    
    /**
     * Set properties for nodes and classes
     */
    private void setProperty(String[] tokens) {
        if (tokens.length < 3) {
            writer.println("Usage: set <element>.<property> <value> [class]");
            return;
        }
        
        String[] parts = tokens[0].split("\\.");
        if (parts.length != 2) {
            writer.println("Property must be in format: element.property");
            return;
        }
        
        String elementName = parts[0];
        String property = parts[1].toLowerCase();
        String value = tokens[1];
        String className = tokens.length > 2 ? tokens[2] : null;
        
        if (nodes.containsKey(elementName)) {
            setNodeProperty(nodes.get(elementName), property, value, className);
        } else if (jobClasses.containsKey(elementName)) {
            setClassProperty(jobClasses.get(elementName), property, value);
        } else {
            writer.println("Unknown element: " + elementName);
        }
    }
    
    /**
     * Set node properties
     */
    private void setNodeProperty(Node node, String property, String value, String className) {
        JobClass jobClass = null;
        if (className != null) {
            if (!jobClasses.containsKey(className)) {
                writer.println("Unknown job class: " + className);
                return;
            }
            jobClass = jobClasses.get(className);
        }
        
        try {
            switch (property) {
                case "service":
                    if (node instanceof Queue && jobClass != null) {
                        Distribution dist = parseDistribution(value);
                        ((Queue) node).setService(jobClass, dist);
                        writer.println("Set service time for " + node.getName() + 
                                     " class " + className + ": " + value);
                    } else if (node instanceof Delay && jobClass != null) {
                        Distribution dist = parseDistribution(value);
                        ((Delay) node).setService(jobClass, dist);
                        writer.println("Set delay time for " + node.getName() + 
                                     " class " + className + ": " + value);
                    }
                    break;
                    
                case "arrival":
                    if (node instanceof Source && jobClass != null) {
                        Distribution dist = parseDistribution(value);
                        ((Source) node).setArrival(jobClass, dist);
                        writer.println("Set arrival rate for " + node.getName() + 
                                     " class " + className + ": " + value);
                    }
                    break;
                    
                case "servers":
                    if (node instanceof Queue) {
                        int servers = Integer.parseInt(value);
                        ((Queue) node).setNumberOfServers(servers);
                        writer.println("Set servers for " + node.getName() + ": " + servers);
                    }
                    break;
                    
                default:
                    writer.println("Unknown property: " + property);
            }
        } catch (Exception e) {
            writer.println("Error setting property: " + e.getMessage());
        }
    }
    
    /**
     * Set job class properties
     */
    private void setClassProperty(JobClass jobClass, String property, String value) {
        try {
            switch (property) {
                case "population":
                    if (jobClass instanceof ClosedClass) {
                        int pop = Integer.parseInt(value);
                        // Note: population is set at construction time for ClosedClass
                        writer.println("Note: Population for closed classes must be set at creation time");
                    }
                    break;
                    
                case "priority":
                    int priority = Integer.parseInt(value);
                    jobClass.setPriority(priority);
                    writer.println("Set priority for " + jobClass.getName() + ": " + priority);
                    break;
                    
                default:
                    writer.println("Unknown property: " + property);
            }
        } catch (Exception e) {
            writer.println("Error setting property: " + e.getMessage());
        }
    }
    
    /**
     * Link nodes with routing
     */
    private void linkNodes(String[] tokens) {
        if (tokens.length < 2) {
            writer.println("Usage: link <from> <to> [probability] [class]");
            return;
        }
        
        String fromName = tokens[0];
        String toName = tokens[1];
        double probability = tokens.length > 2 ? Double.parseDouble(tokens[2]) : 1.0;
        String className = tokens.length > 3 ? tokens[3] : null;
        
        if (!nodes.containsKey(fromName) || !nodes.containsKey(toName)) {
            writer.println("Both nodes must exist in the network");
            return;
        }
        
        Node from = nodes.get(fromName);
        Node to = nodes.get(toName);
        
        // Initialize routing matrix if needed
        if (this.routingMatrix == null) {
            this.routingMatrix = currentNetwork.initRoutingMatrix();
        }
        
        if (className != null) {
            if (!jobClasses.containsKey(className)) {
                writer.println("Unknown job class: " + className);
                return;
            }
            JobClass jobClass = jobClasses.get(className);
            
            // Set routing probability for specific class
            if (from instanceof Station) {
                ((Station) from).setProbRouting(jobClass, to, probability);
            }
            
            writer.println("Linked " + fromName + " -> " + toName + 
                         " (probability=" + probability + ", class=" + className + ")");
        } else {
            // Add link for all classes
            currentNetwork.addLink(from, to);
            writer.println("Linked " + fromName + " -> " + toName);
        }
    }
    
    /**
     * Solve the current network
     */
    private void solve(String solverName) {
        NetworkSolver solver;
        
        switch (solverName.toLowerCase()) {
            case "mva":
                solver = new SolverMVA(currentNetwork);
                break;
                
            case "jmt":
                solver = new SolverJMT(currentNetwork);
                break;
                
            case "mam":
                solver = new SolverMAM(currentNetwork);
                break;
                
            case "nc":
                solver = new SolverNC(currentNetwork);
                break;
                
            case "ctmc":
                solver = new SolverCTMC(currentNetwork);
                break;
                
            case "fluid":
                solver = new SolverFluid(currentNetwork);
                break;
                
            case "auto":
            default:
                solver = new SolverAUTO(currentNetwork);
                break;
        }
        
        writer.println("Running " + solver.getName() + " solver...");
        
        try {
            // Apply routing matrix if we've been using setProbRouting
            if (this.routingMatrix == null && currentNetwork.getConnectionMatrix().isEmpty()) {
                // Create simple serial routing for all nodes
                List<Node> nodeList = new ArrayList<>(nodes.values());
                for (int i = 0; i < nodeList.size() - 1; i++) {
                    currentNetwork.addLink(nodeList.get(i), nodeList.get(i + 1));
                }
            }
            
            // Get results (this automatically runs the analyzer)
            NetworkAvgTable avgTable = solver.getAvgTable();
            displayResults(solver, avgTable);
        } catch (Exception e) {
            writer.println("Solver error: " + e.getMessage());
            e.printStackTrace(writer);
        }
    }
    
    /**
     * Display solver results
     */
    private void displayResults(NetworkSolver solver, NetworkAvgTable avgTable) {
        writer.println("\n=== Performance Metrics ===");
        
        // Display the average table
        if (avgTable != null) {
            avgTable.print();
        }
        
        // Additional metrics can be displayed if needed
        writer.println("\nUse 'show detailed' for more performance metrics.");
    }
    
    /**
     * Show information about the current model
     */
    private void show(String what) {
        switch (what.toLowerCase()) {
            case "network":
                if (currentNetwork == null) {
                    writer.println("No network defined");
                } else {
                    writer.println("Network: " + currentNetwork.getName());
                    writer.println("Nodes: " + nodes.size());
                    writer.println("Classes: " + jobClasses.size());
                }
                break;
                
            case "nodes":
                if (nodes.isEmpty()) {
                    writer.println("No nodes defined");
                } else {
                    writer.println("Nodes:");
                    for (Map.Entry<String, Node> entry : nodes.entrySet()) {
                        Node node = entry.getValue();
                        writer.println("  " + entry.getKey() + " (" + 
                                     node.getClass().getSimpleName() + ")");
                    }
                }
                break;
                
            case "classes":
                if (jobClasses.isEmpty()) {
                    writer.println("No job classes defined");
                } else {
                    writer.println("Job Classes:");
                    for (Map.Entry<String, JobClass> entry : jobClasses.entrySet()) {
                        JobClass jc = entry.getValue();
                        String type = jc instanceof OpenClass ? "open" : "closed";
                        writer.println("  " + entry.getKey() + " (" + type + ")");
                    }
                }
                break;
                
            case "variables":
                if (variables.isEmpty()) {
                    writer.println("No variables defined");
                } else {
                    writer.println("Variables:");
                    for (Map.Entry<String, Object> entry : variables.entrySet()) {
                        writer.println("  " + entry.getKey() + " = " + entry.getValue());
                    }
                }
                break;
                
            default:
                writer.println("Unknown target: " + what);
                writer.println("Available: network, nodes, classes, variables");
        }
    }
    
    /**
     * Save the current model
     */
    private void saveModel(String filename) {
        writer.println("Save functionality not yet implemented");
    }
    
    /**
     * Load a model from file
     */
    private void loadModel(String filename) {
        writer.println("Load functionality not yet implemented");
    }
    
    /**
     * Clear the current model
     */
    private void clear() {
        currentNetwork = null;
        nodes.clear();
        jobClasses.clear();
        variables.clear();
        writer.println("Model cleared");
    }
    
    /**
     * Define a variable
     */
    private void defineVariable(String[] tokens) {
        if (tokens.length < 3 || !tokens[1].equals("=")) {
            writer.println("Usage: let <name> = <value>");
            return;
        }
        
        String name = tokens[0];
        String value = String.join(" ", Arrays.copyOfRange(tokens, 2, tokens.length));
        
        try {
            // Try to parse as number
            if (value.contains(".")) {
                variables.put(name, Double.parseDouble(value));
            } else {
                variables.put(name, Integer.parseInt(value));
            }
        } catch (NumberFormatException e) {
            // Store as string
            variables.put(name, value);
        }
        
        writer.println(name + " = " + variables.get(name));
    }
    
    /**
     * Parse a distribution string
     */
    private Distribution parseDistribution(String spec) {
        // Handle variable references
        if (spec.startsWith("$") && variables.containsKey(spec.substring(1))) {
            Object val = variables.get(spec.substring(1));
            if (val instanceof Number) {
                spec = val.toString();
            }
        }
        
        // Parse distribution
        if (spec.matches("\\d+(\\.\\d+)?")) {
            // Plain number - create exponential with rate
            double rate = Double.parseDouble(spec);
            return new Exp(1.0 / rate);
        } else if (spec.toLowerCase().startsWith("exp(")) {
            double param = extractParam(spec);
            return new Exp(param);
        } else if (spec.toLowerCase().startsWith("det(")) {
            double param = extractParam(spec);
            return new Det(param);
        } else if (spec.toLowerCase().startsWith("uniform(")) {
            String params = spec.substring(8, spec.length() - 1);
            String[] parts = params.split(",");
            if (parts.length == 2) {
                double min = Double.parseDouble(parts[0].trim());
                double max = Double.parseDouble(parts[1].trim());
                return new Uniform(min, max);
            }
        }
        
        throw new IllegalArgumentException("Unknown distribution: " + spec);
    }
    
    /**
     * Parse scheduling strategy
     */
    private SchedStrategy parseSchedStrategy(String name) {
        switch (name.toUpperCase()) {
            case "FCFS":
                return SchedStrategy.FCFS;
            case "PS":
                return SchedStrategy.PS;
            case "LCFS":
                return SchedStrategy.LCFS;
            case "SIRO":
                return SchedStrategy.SIRO;
            case "INF":
                return SchedStrategy.INF;
            case "HOL":
                return SchedStrategy.HOL;
            default:
                throw new IllegalArgumentException("Unknown scheduling strategy: " + name);
        }
    }
    
    /**
     * Extract parameter from distribution string
     */
    private double extractParam(String spec) {
        Pattern pattern = Pattern.compile("\\((.*?)\\)");
        Matcher matcher = pattern.matcher(spec);
        if (matcher.find()) {
            String param = matcher.group(1);
            if (param.startsWith("$") && variables.containsKey(param.substring(1))) {
                Object val = variables.get(param.substring(1));
                if (val instanceof Number) {
                    return ((Number) val).doubleValue();
                }
            }
            return Double.parseDouble(param);
        }
        throw new IllegalArgumentException("Invalid distribution specification: " + spec);
    }
    
    /**
     * Tokenize command line
     */
    private String[] tokenize(String line) {
        List<String> tokens = new ArrayList<>();
        Pattern pattern = Pattern.compile("\"([^\"]*)\"|'([^']*)'|(\\S+)");
        Matcher matcher = pattern.matcher(line);
        
        while (matcher.find()) {
            if (matcher.group(1) != null) {
                tokens.add(matcher.group(1));
            } else if (matcher.group(2) != null) {
                tokens.add(matcher.group(2));
            } else {
                tokens.add(matcher.group(3));
            }
        }
        
        return tokens.toArray(new String[0]);
    }
    
    /**
     * Format double for display
     */
    private String formatDouble(double value) {
        if (Double.isNaN(value)) {
            return "NaN";
        } else if (Double.isInfinite(value)) {
            return "Inf";
        } else if (value == 0.0) {
            return "0";
        } else if (Math.abs(value) < 0.0001 || Math.abs(value) > 10000) {
            return String.format("%.4e", value);
        } else {
            return String.format("%.4f", value);
        }
    }
    
    /**
     * Show help
     */
    private void showHelp(String topic) {
        if (topic == null) {
            writer.println("\n=== LINE REPL Commands ===");
            writer.println("");
            writer.println("Model Definition:");
            writer.println("  new network [name]              - Create a new network");
            writer.println("  add queue <name> [options]      - Add a queue node");
            writer.println("  add source <name>               - Add a source node");
            writer.println("  add sink <name>                 - Add a sink node");
            writer.println("  add delay <name>                - Add a delay station");
            writer.println("  add class <name> [options]      - Add a job class");
            writer.println("");
            writer.println("Configuration:");
            writer.println("  set <element>.<property> <value> [class] - Set element properties");
            writer.println("  link <from> <to> [prob] [class] - Link nodes with routing");
            writer.println("");
            writer.println("Solving:");
            writer.println("  solve [solver]                  - Solve the model (auto/mva/jmt/mam/nc/ctmc/fluid)");
            writer.println("");
            writer.println("Information:");
            writer.println("  show [network|nodes|classes|variables] - Display information");
            writer.println("  help [topic]                    - Show help");
            writer.println("");
            writer.println("Utilities:");
            writer.println("  let <name> = <value>            - Define a variable");
            writer.println("  clear                           - Clear the current model");
            writer.println("  save <filename>                 - Save model to file");
            writer.println("  load <filename>                 - Load model from file");
            writer.println("  quit/exit                       - Exit REPL");
            writer.println("");
            writer.println("Type 'help <command>' for detailed help on specific commands.");
            writer.println("");
            writer.println("Additional Help:");
            writer.println("  help syntax    - Complete language syntax reference");
            writer.println("  help quick     - Quick reference card");
            writer.println("  help examples  - Example models and use cases");
        } else {
            showDetailedHelp(topic);
        }
    }
    
    /**
     * Show detailed help for specific topics
     */
    private void showDetailedHelp(String topic) {
        switch (topic.toLowerCase()) {
            case "queue":
                writer.println("\n=== Adding Queues ===");
                writer.println("Usage: add queue <name> [options]");
                writer.println("");
                writer.println("Options:");
                writer.println("  -servers <n>    Number of servers (default: 1)");
                writer.println("  -sched <type>   Scheduling strategy (default: FCFS)");
                writer.println("");
                writer.println("Scheduling strategies: FCFS, PS, LCFS, SIRO, INF, HOL");
                writer.println("");
                writer.println("Example:");
                writer.println("  add queue CPU -servers 2 -sched PS");
                break;
                
            case "class":
                writer.println("\n=== Adding Job Classes ===");
                writer.println("Usage: add class <name> [options]");
                writer.println("");
                writer.println("Options:");
                writer.println("  -type <open|closed>     Class type (default: closed)");
                writer.println("  -population <n>         Population for closed class (default: 1)");
                writer.println("  -ref <node>             Reference node");
                writer.println("  -priority <n>           Priority level (default: 0)");
                writer.println("");
                writer.println("Example:");
                writer.println("  add class WebUsers -type closed -population 10 -ref Delay1");
                break;
                
            case "set":
                writer.println("\n=== Setting Properties ===");
                writer.println("Usage: set <element>.<property> <value> [class]");
                writer.println("");
                writer.println("Properties:");
                writer.println("  <queue>.service <dist> <class>   - Service time distribution");
                writer.println("  <source>.arrival <dist> <class>  - Arrival distribution");
                writer.println("  <delay>.service <dist> <class>   - Delay time distribution");
                writer.println("  <queue>.servers <n>              - Number of servers");
                writer.println("  <class>.population <n>           - Class population");
                writer.println("  <class>.priority <n>             - Class priority");
                writer.println("");
                writer.println("Distributions:");
                writer.println("  <number>         - Exponential with mean service time");
                writer.println("  exp(<rate>)      - Exponential distribution");
                writer.println("  det(<value>)     - Deterministic distribution");
                writer.println("  uniform(<min>,<max>) - Uniform distribution");
                writer.println("");
                writer.println("Example:");
                writer.println("  set CPU.service 0.1 WebUsers");
                writer.println("  set Source1.arrival exp(2.0) WebUsers");
                break;
                
            case "link":
                writer.println("\n=== Linking Nodes ===");
                writer.println("Usage: link <from> <to> [probability] [class]");
                writer.println("");
                writer.println("Creates routing between nodes.");
                writer.println("If probability is omitted, defaults to 1.0");
                writer.println("If class is omitted, applies to all classes");
                writer.println("");
                writer.println("Example:");
                writer.println("  link Source1 CPU");
                writer.println("  link CPU Disk 0.7 WebUsers");
                writer.println("  link CPU Sink 0.3 WebUsers");
                break;
                
            case "solve":
                writer.println("\n=== Solving Models ===");
                writer.println("Usage: solve [solver]");
                writer.println("");
                writer.println("Available solvers:");
                writer.println("  auto   - Automatic solver selection (default)");
                writer.println("  mva    - Mean Value Analysis");
                writer.println("  jmt    - Java Modelling Tools simulation");
                writer.println("  mam    - Matrix Analytic Methods");
                writer.println("  nc     - Normalizing Constant");
                writer.println("  ctmc   - Continuous Time Markov Chain");
                writer.println("  fluid  - Fluid approximation");
                writer.println("");
                writer.println("Example:");
                writer.println("  solve");
                writer.println("  solve mva");
                break;
                
            case "syntax":
                showSyntaxHelp();
                break;
                
            case "quick":
            case "quickref":
            case "reference":
                showQuickReference();
                break;
                
            case "examples":
                showExamples();
                break;
                
            default:
                writer.println("No detailed help available for: " + topic);
                writer.println("Try: help, help syntax, help quick, help examples");
        }
    }
    
    /**
     * Show language syntax help
     */
    private void showSyntaxHelp() {
        writer.println("\n=== LINE REPL Language Syntax Reference ===");
        writer.println("");
        writer.println("## MODEL DEFINITION SYNTAX");
        writer.println("");
        writer.println("### Network Creation");
        writer.println("  new network [<name>]");
        writer.println("    Creates a new queueing network model");
        writer.println("    Example: new network MyWebApp");
        writer.println("");
        writer.println("### Node Types");
        writer.println("  add queue <name> [-servers <n>] [-sched <strategy>]");
        writer.println("    -servers: number of servers (default: 1, use 'inf' for infinite)");
        writer.println("    -sched: FCFS, PS, LCFS, SIRO, INF, HOL (default: FCFS)");
        writer.println("");
        writer.println("  add source <name>");
        writer.println("    Entry point for open classes");
        writer.println("");
        writer.println("  add sink <name>");
        writer.println("    Exit point for open classes");
        writer.println("");
        writer.println("  add delay <name>");
        writer.println("    Infinite server station (think time)");
        writer.println("");
        writer.println("### Job Classes");
        writer.println("  add class <name> [-type <open|closed>] [-population <n>] [-ref <node>] [-priority <p>]");
        writer.println("    -type: open or closed (default: closed)");
        writer.println("    -population: number of jobs for closed class (default: 1)");
        writer.println("    -ref: reference node (required: source for open, any station for closed)");
        writer.println("    -priority: priority level (default: 0, higher = more priority)");
        writer.println("");
        writer.println("## CONFIGURATION SYNTAX");
        writer.println("");
        writer.println("### Service/Arrival Times");
        writer.println("  set <node>.<property> <distribution> <class>");
        writer.println("");
        writer.println("  Properties:");
        writer.println("    - queue.service: service time at queue");
        writer.println("    - delay.service: delay time at delay station");
        writer.println("    - source.arrival: interarrival time at source");
        writer.println("    - queue.servers: number of servers");
        writer.println("");
        writer.println("  Distributions:");
        writer.println("    <number>           - Exponential with mean <number>");
        writer.println("    exp(<rate>)        - Exponential with rate <rate>");
        writer.println("    det(<value>)       - Deterministic (constant) time");
        writer.println("    uniform(<a>,<b>)   - Uniform between a and b");
        writer.println("");
        writer.println("  Examples:");
        writer.println("    set CPU.service 0.01 WebUsers       # Exp with mean 0.01");
        writer.println("    set CPU.service exp(100) WebUsers   # Exp with rate 100");
        writer.println("    set Disk.service det(0.02) DBUsers  # Constant 0.02");
        writer.println("    set Think.service uniform(1,5) Users # Uniform [1,5]");
        writer.println("");
        writer.println("### Routing");
        writer.println("  link <from> <to> [<probability>] [<class>]");
        writer.println("    Connects nodes with optional routing probability");
        writer.println("    If probability omitted, defaults to 1.0");
        writer.println("    If class omitted, applies to all classes");
        writer.println("");
        writer.println("  Examples:");
        writer.println("    link Source1 CPU                    # Deterministic routing");
        writer.println("    link CPU Disk 0.3 WebUsers         # 30% to Disk");
        writer.println("    link CPU Network 0.7 WebUsers      # 70% to Network");
        writer.println("");
        writer.println("## VARIABLE SYNTAX");
        writer.println("");
        writer.println("  let <name> = <value>");
        writer.println("    Define a variable for reuse");
        writer.println("    Variables can be referenced with $name");
        writer.println("");
        writer.println("  Examples:");
        writer.println("    let thinkTime = 5");
        writer.println("    let cpuService = 0.01");
        writer.println("    set CPU.service $cpuService WebUsers");
        writer.println("    set Think.service $thinkTime Users");
        writer.println("");
        writer.println("## COMPLETE EXAMPLE");
        writer.println("");
        writer.println("  # Create a closed queueing network");
        writer.println("  new network WebApplication");
        writer.println("  ");
        writer.println("  # Add nodes");
        writer.println("  add queue CPU -servers 2 -sched PS");
        writer.println("  add queue Disk -servers 1 -sched FCFS");
        writer.println("  add delay Think");
        writer.println("  ");
        writer.println("  # Add job class");
        writer.println("  add class Users -type closed -population 50 -ref Think");
        writer.println("  ");
        writer.println("  # Configure service times");
        writer.println("  let cpuTime = 0.01");
        writer.println("  set CPU.service $cpuTime Users");
        writer.println("  set Disk.service 0.03 Users");
        writer.println("  set Think.service 5 Users");
        writer.println("  ");
        writer.println("  # Configure routing");
        writer.println("  link Think CPU");
        writer.println("  link CPU Disk 0.4 Users");
        writer.println("  link CPU Think 0.6 Users");
        writer.println("  link Disk CPU");
        writer.println("  ");
        writer.println("  # Solve");
        writer.println("  solve mva");
        writer.println("");
        writer.println("## CONVENTIONS");
        writer.println("");
        writer.println("  - Names are case-sensitive");
        writer.println("  - Quotes optional for single-word names");
        writer.println("  - Use quotes for names with spaces: \"Web Server\"");
        writer.println("  - Comments start with # (not implemented yet)");
        writer.println("  - Commands are case-insensitive");
        writer.println("  - Partial command matching supported (e.g., 'q' for 'quit')");
    }
    
    /**
     * Show quick reference card
     */
    private void showQuickReference() {
        writer.println("\n=== LINE REPL Quick Reference ===");
        writer.println("");
        writer.println("BASIC COMMANDS:");
        writer.println("  new network [name]          # Create network");
        writer.println("  add queue Q1 -servers 2     # Add queue (FCFS default)");
        writer.println("  add source S1               # Add source (open classes)");
        writer.println("  add sink K1                 # Add sink (open classes)");
        writer.println("  add delay D1                # Add delay (think time)");
        writer.println("  add class C1 -pop 10        # Add closed class");
        writer.println("  link A B [prob] [class]     # Connect nodes");
        writer.println("  solve [solver]              # Run solver");
        writer.println("");
        writer.println("SERVICE TIMES:");
        writer.println("  set Q1.service 0.1 C1       # Exponential mean 0.1");
        writer.println("  set Q1.service exp(10) C1   # Exponential rate 10");
        writer.println("  set Q1.service det(0.5) C1  # Deterministic 0.5");
        writer.println("");
        writer.println("SCHEDULING:");
        writer.println("  -sched FCFS    # First Come First Served");
        writer.println("  -sched PS      # Processor Sharing");
        writer.println("  -sched LCFS    # Last Come First Served");
        writer.println("  -sched INF     # Infinite Server");
        writer.println("");
        writer.println("SOLVERS:");
        writer.println("  solve          # Auto-select");
        writer.println("  solve mva      # Mean Value Analysis");
        writer.println("  solve jmt      # Simulation");
        writer.println("");
        writer.println("Type 'help syntax' for full syntax reference");
        writer.println("Type 'help examples' for complete examples");
    }
    
    /**
     * Show example models
     */
    private void showExamples() {
        writer.println("\n=== LINE REPL Examples ===");
        writer.println("");
        writer.println("## Example 1: Simple M/M/1 Queue (Open Network)");
        writer.println("  new network MM1");
        writer.println("  add source Source");
        writer.println("  add queue Server");
        writer.println("  add sink Sink");
        writer.println("  add class Customers -type open -ref Source");
        writer.println("  set Source.arrival 0.9 Customers    # λ = 0.9");
        writer.println("  set Server.service 1.0 Customers    # μ = 1.0");
        writer.println("  link Source Server");
        writer.println("  link Server Sink");
        writer.println("  solve");
        writer.println("");
        writer.println("## Example 2: Multi-Server Queue M/M/c");
        writer.println("  new network MMc");
        writer.println("  add source Arrivals");
        writer.println("  add queue CallCenter -servers 5");
        writer.println("  add sink Departures");
        writer.println("  add class Calls -type open -ref Arrivals");
        writer.println("  set Arrivals.arrival 4.5 Calls");
        writer.println("  set CallCenter.service 1.0 Calls");
        writer.println("  link Arrivals CallCenter");
        writer.println("  link CallCenter Departures");
        writer.println("  solve");
        writer.println("");
        writer.println("## Example 3: Closed Network (Computer System)");
        writer.println("  new network ComputerSystem");
        writer.println("  add queue CPU -sched PS");
        writer.println("  add queue Disk1");
        writer.println("  add queue Disk2");
        writer.println("  add delay Terminal");
        writer.println("  add class Users -type closed -pop 20 -ref Terminal");
        writer.println("  set CPU.service 0.005 Users");
        writer.println("  set Disk1.service 0.030 Users");
        writer.println("  set Disk2.service 0.027 Users");
        writer.println("  set Terminal.service 15.0 Users");
        writer.println("  link Terminal CPU");
        writer.println("  link CPU Disk1 0.5 Users");
        writer.println("  link CPU Disk2 0.3 Users");
        writer.println("  link CPU Terminal 0.2 Users");
        writer.println("  link Disk1 CPU");
        writer.println("  link Disk2 CPU");
        writer.println("  solve mva");
        writer.println("");
        writer.println("## Example 4: Multi-Class Network");
        writer.println("  new network WebServer");
        writer.println("  add queue Frontend -servers 2 -sched PS");
        writer.println("  add queue Backend -servers 4");
        writer.println("  add queue Database");
        writer.println("  add delay Browse");
        writer.println("  add delay Purchase");
        writer.println("  ");
        writer.println("  # Two types of users");
        writer.println("  add class Browsers -pop 100 -ref Browse");
        writer.println("  add class Buyers -pop 20 -ref Purchase -priority 1");
        writer.println("  ");
        writer.println("  # Service times for browsers");
        writer.println("  set Frontend.service 0.001 Browsers");
        writer.println("  set Backend.service 0.002 Browsers");
        writer.println("  set Database.service 0.001 Browsers");
        writer.println("  set Browse.service 5.0 Browsers");
        writer.println("  ");
        writer.println("  # Service times for buyers (slower)");
        writer.println("  set Frontend.service 0.002 Buyers");
        writer.println("  set Backend.service 0.005 Buyers");
        writer.println("  set Database.service 0.010 Buyers");
        writer.println("  set Purchase.service 10.0 Buyers");
        writer.println("  ");
        writer.println("  # Routing for browsers");
        writer.println("  link Browse Frontend");
        writer.println("  link Frontend Backend 0.7 Browsers");
        writer.println("  link Frontend Browse 0.3 Browsers");
        writer.println("  link Backend Database 0.5 Browsers");
        writer.println("  link Backend Browse 0.5 Browsers");
        writer.println("  link Database Browse");
        writer.println("  ");
        writer.println("  # Routing for buyers");
        writer.println("  link Purchase Frontend");
        writer.println("  link Frontend Backend 1.0 Buyers");
        writer.println("  link Backend Database 0.9 Buyers");
        writer.println("  link Backend Purchase 0.1 Buyers");
        writer.println("  link Database Backend 0.3 Buyers");
        writer.println("  link Database Purchase 0.7 Buyers");
        writer.println("  ");
        writer.println("  solve");
    }
    
    /**
     * Print welcome message
     */
    private void printWelcome() {
        writer.println("LINE REPL v" + VERSION);
        writer.println("Type 'help' for available commands, 'quit' to exit");
        writer.println("");
    }
    
    /**
     * Main method to start REPL
     */
    public static void main(String[] args) {
        REPL repl = new REPL();
        repl.start();
    }
}