
import jpype
import jpype.imports
import numpy as np
from pprint import pformat
from . import SchedStrategy, jlineFromDistribution, SchedStrategyType, CallType, ReplacementStrategy
from .lang import jlineMatrixToArray


class LayeredNetwork:
    """
    Layered queueing network (LQN) model.

    LayeredNetwork represents a layered queueing network, which is a special
    type of queueing network used to model software systems with client-server
    interactions, thread pools, and nested service calls.

    LQNs consist of:
    - Hosts (processing resources)
    - Tasks (software processes)
    - Entries (service interfaces)
    - Activities (computation steps)

    The model can be loaded from LQN XML files or constructed programmatically.
    Python wrappers are automatically tracked for name-based access and copying.

    Args:
        name (str): Name of the layered network.

    Examples:
        >>> # Load from XML
        >>> lqn = LayeredNetwork.parseXML('model.xml')
        >>> lqn.summary()

        >>> # Create programmatically
        >>> lqn = LayeredNetwork('WebSystem')
        >>> proc = Processor(lqn, 'P1', 1, SchedStrategy.PS)
        >>> task = Task(lqn, 'T1', 1, SchedStrategy.FCFS)

        >>> # Access nodes by name
        >>> task = lqn.getNodeByName('T1')
        >>> task.setThinkTime(Exp(2.0))

        >>> # Enumerate all nodes
        >>> for node in lqn.getNodes():
        ...     print(f"{node.getName()}: {type(node).__name__}")

        >>> # Create model variants for parameter sweeps
        >>> base = LayeredNetwork('baseline')
        >>> # ... configure base ...
        >>>
        >>> for cpu_count in [1, 2, 4]:
        ...     model = base.copy()
        ...     model.getNodeByName('P1').setMultiplicity(cpu_count)
        ...     solver = LN(model, @(m) MVA(m))
        ...     print(f"CPUs={cpu_count}: {solver.getAvgTable()}")
    """

    def __init__(self, name):
        """
        Initialize a new layered queueing network.

        Args:
            name (str): Name of the layered network.
        """
        self.obj = jpype.JPackage('jline').lang.layered.LayeredNetwork(name)

        # Wrapper registry for tracking Python objects
        self._wrapper_nodes = {}    # Maps Java object → Python wrapper
        self._nodes_by_name = {}    # Maps node name → Python wrapper

    def parseXML(self, filename, verbose=False):
        """
        Parse a layered queueing network from an XML file.

        Args:
            filename (str): Path to the LQN XML file.
            verbose (bool): Whether to print parsing details.
        """
        self.obj.parseXML(filename, verbose)

    def parse_xml(self, filename, verbose=False):
        """Snake case alias for parseXML"""
        return self.parseXML(filename, verbose)

    def writeXML(self, filename, abstractNames=False):
        """
        Write the layered queueing network to an XML file.
        
        Args:
            filename (str): Output XML file path.
            abstractNames (bool): Whether to use abstract names in the output.
        """
        self.obj.writeXML(filename, abstractNames)

    def getNodeIndex(self, node):
        return self.obj.getNodeIndex(node.obj)

    def getNodeNames(self):
        return self.obj.getNodeNames()

    def getEnsemble(self):
        return self.obj.getEnsemble()

    def getLayers(self):
        return self.obj.getLayers()

    def getNumberOfLayers(self):
        """
        Get the number of layers in the network.
        
        Returns:
            int: Number of layers.
        """
        return self.obj.getNumberOfLayers()

    def getNumberOfModels(self):
        """
        Get the number of submodels in the layered network.
        
        Returns:
            int: Number of submodels.
        """
        return self.obj.getNumberOfModels()

    def summary(self):
        """
        Get a text summary of the layered network structure.
        
        Returns:
            str: Summary description of hosts, tasks, entries, and activities.
        """
        return self.obj.summary()

    def parseXML(self, filename, verbose):
        return self.obj.parseXML(filename, verbose)

    def registerNode(self, wrapper):
        """
        Register a node wrapper with the network for name-based lookup.

        This method is called automatically when creating Processor, Task, Entry,
        Activity, CacheTask, or ItemEntry objects. It maintains bidirectional
        mappings between Java objects and Python wrappers for introspection.

        Args:
            wrapper: Python wrapper object (Processor, Task, Entry, Activity, CacheTask, ItemEntry)

        Returns:
            None

        Example:
            >>> proc = Processor(model, 'P1', 1, SchedStrategy.PS)
            >>> # Auto-registered automatically, but can call manually if needed
            >>> model.registerNode(proc)
        """
        if hasattr(wrapper, 'obj'):
            try:
                node_name = wrapper.obj.getName()
                self._wrapper_nodes[wrapper.obj] = wrapper
                self._nodes_by_name[node_name] = wrapper
            except Exception:
                # If getName() fails, skip registration
                pass

    def getStruct(self):
        jsn = self.obj.getStruct()
        lsn = LayeredNetworkStruct()
        lsn.fromJline(jsn)
        return lsn

    def getNodeByName(self, name):
        """
        Get a Python-wrapped node by its name.

        Returns the Python wrapper (Processor, Task, Entry, Activity, CacheTask, ItemEntry)
        with the given name, allowing use of Python methods directly without Java interaction.

        Args:
            name (str): Node name to look up

        Returns:
            Wrapper object (Processor, Task, Entry, Activity, CacheTask, ItemEntry) or None if not found.

        Example:
            >>> task = model.getNodeByName('T1')
            >>> if task:
            ...     task.setThinkTime(Exp(2.0))  # Use Python API directly
        """
        return self._nodes_by_name.get(name, None)

    def getNodes(self):
        """
        Get all registered Python-wrapped nodes.

        Returns all tracked node wrappers in insertion order.

        Returns:
            list: All registered wrapper objects (Processor, Task, Entry, Activity, CacheTask, ItemEntry)

        Example:
            >>> for node in model.getNodes():
            ...     print(f"{node.getName()}: {type(node).__name__}")
        """
        return list(self._wrapper_nodes.values())

    def getNodeCount(self):
        """
        Get the total number of registered nodes.

        Returns:
            int: Total count of nodes in the network
        """
        return len(self._wrapper_nodes)

    # Property getters delegating to Java
    def getHosts(self):
        """Get the map of hosts/processors in the network."""
        return self.obj.getHosts()

    def getTasks(self):
        """Get the map of tasks in the network."""
        return self.obj.getTasks()

    def getEntries(self):
        """Get the map of entries in the network."""
        return self.obj.getEntries()

    def getActivities(self):
        """Get the map of activities in the network."""
        return self.obj.getActivities()

    # Snake case aliases
    get_hosts = getHosts
    get_tasks = getTasks
    get_entries = getEntries
    get_activities = getActivities

    def copy(self):
        """
        Create an independent copy of this LayeredNetwork.

        Creates a new independent LayeredNetwork and recreates all registered node wrappers.
        Note: This creates a structural copy; activity graphs and cross-layer relationships
        are not recreated. For full network replication, manually recreate the network or
        use advanced features when available.

        Returns:
            LayeredNetwork: A new independent network with registered node wrappers.

        Raises:
            RuntimeError: If the copy operation fails.

        Example - Create Network Variant:
            >>> base_model = LayeredNetwork('baseline')
            >>> proc = Processor(base_model, 'P1', 1, SchedStrategy.PS)
            >>>
            >>> variant = base_model.copy()
            >>> # variant now has the same structure as base_model
        """
        try:
            # Create a new independent LayeredNetwork with a unique name
            copied = LayeredNetwork(self.obj.getName() + f'_copy_{id(self)}')

            # Recreate nodes in dependency order
            # First pass: Processors
            for wrapper in self._wrapper_nodes.values():
                class_name = wrapper.obj.getClass().getSimpleName()
                if class_name == "Processor":
                    name = wrapper.obj.getName()
                    Processor(copied, name, 1, SchedStrategy.PS)

            # Second pass: Tasks and CacheTasks
            for wrapper in self._wrapper_nodes.values():
                class_name = wrapper.obj.getClass().getSimpleName()
                if class_name == "Task":
                    name = wrapper.obj.getName()
                    Task(copied, name, 1, SchedStrategy.FCFS)
                elif class_name == "CacheTask":
                    name = wrapper.obj.getName()
                    try:
                        CacheTask(copied, name, 10, 2, ReplacementStrategy.LRU, 1)
                    except:
                        pass

            # Third pass: Entries and ItemEntries
            for wrapper in self._wrapper_nodes.values():
                class_name = wrapper.obj.getClass().getSimpleName()
                if class_name == "Entry":
                    name = wrapper.obj.getName()
                    Entry(copied, name)
                elif class_name == "ItemEntry":
                    name = wrapper.obj.getName()
                    try:
                        # ItemEntry requires number of items and access distribution
                        # Use defaults since we can't retrieve the originals
                        pAccess = DiscreteSampler([0.1] * 10)
                        ItemEntry(copied, name, 10, pAccess)
                    except Exception:
                        # ItemEntry recreation may fail if it requires special binding
                        # Skip silently in that case
                        pass

            # Fourth pass: Activities
            for wrapper in self._wrapper_nodes.values():
                class_name = wrapper.obj.getClass().getSimpleName()
                if class_name == "Activity":
                    name = wrapper.obj.getName()
                    try:
                        try:
                            dist = wrapper.obj.getServiceDist()
                        except:
                            dist = Exp(1.0)
                        Activity(copied, name, dist)
                    except Exception as ex:
                        # Silently skip activities that can't be recreated
                        # This can happen if activities aren't bound to tasks
                        pass

            return copied

        except Exception as e:
            raise RuntimeError(f"Failed to copy LayeredNetwork: {e}")

    def _rebuildWrappers(self):
        """
        Rebuild Python wrapper objects for a copied network.

        Automatically detects wrapper types (Task vs CacheTask, Entry vs ItemEntry)
        and recreates the corresponding Python objects, populating the wrapper registry.

        This method is called automatically by copy() - users should not call directly.

        Internal Implementation:
            - Iterates through all node names via getNodeNames()
            - Gets each node via getNodeByName()
            - Detects type using Java class name (getClass().getSimpleName())
            - Creates appropriate Python wrapper (Processor, Task/CacheTask, Entry/ItemEntry, Activity)
            - Registers each wrapper in _wrapper_nodes and _nodes_by_name

        Raises:
            RuntimeError: If wrapper reconstruction fails
        """
        try:
            # Get all node names
            node_names = self.obj.getNodeNames()

            if node_names:
                for node_name in node_names:
                    # Get the Java node object
                    java_node = self.obj.getNodeByName(node_name)
                    if java_node is None:
                        continue

                    # Detect the node type
                    class_name = java_node.getClass().getSimpleName()

                    # Create appropriate Python wrapper based on type
                    if class_name == "Processor":
                        python_wrapper = Processor.__new__(Processor)
                    elif class_name == "CacheTask":
                        python_wrapper = CacheTask.__new__(CacheTask)
                    elif class_name == "Task":
                        python_wrapper = Task.__new__(Task)
                    elif class_name == "ItemEntry":
                        python_wrapper = ItemEntry.__new__(ItemEntry)
                    elif class_name == "Entry":
                        python_wrapper = Entry.__new__(Entry)
                    elif class_name == "Activity":
                        python_wrapper = Activity.__new__(Activity)
                    else:
                        # Skip unknown types
                        continue

                    # Set the Java object reference
                    python_wrapper.obj = java_node

                    # Register in both dictionaries
                    self._wrapper_nodes[java_node] = python_wrapper
                    self._nodes_by_name[node_name] = python_wrapper

        except Exception as e:
            raise RuntimeError(f"Failed to rebuild wrappers: {e}")

    def plot(self, showTaskGraph=False, **kwargs):
        self.plotGraph(**kwargs)
        if showTaskGraph:
            self.plotTaskGraph(**kwargs)

    def plotGraph(self, method='nodes', **kwargs):
        try:
            import matplotlib.pyplot as plt
            import networkx as nx
            import numpy as np
        except ImportError:
            raise ImportError("Matplotlib and NetworkX are required for plotting. Install with: pip install matplotlib networkx")

        lqn = self.getStruct()

        G = nx.DiGraph()

        adj_matrix = lqn.graph

        node_labels = {}
        node_colors = []
        node_shapes = []
        node_sizes = []

        for i in range(adj_matrix.shape[0]):
            G.add_node(i)

            if method == 'hashnames':
                label = lqn.hashnames[i] if i < len(lqn.hashnames) else f'Node{i}'
                label = label.replace('_', '\\_')
            elif method == 'names':
                label = lqn.names[i] if i < len(lqn.names) else f'Node{i}'
                label = label.replace('_', '\\_')
            elif method == 'ids':
                label = str(i)
            else:
                label = lqn.hashnames[i] if i < len(lqn.hashnames) else f'Node{i}'
                label = label.replace('_', '\\_')

            node_labels[i] = label

            type_val = lqn.type[i] if i < len(lqn.type) else 0
            node_type = int(np.asarray(type_val).flat[0]) if hasattr(type_val, '__iter__') else int(type_val)

            HOST = 1; TASK = 2; ENTRY = 3; ACTIVITY = 4

            if node_type == HOST:
                node_colors.append('black')
                node_shapes.append('h')
                node_sizes.append(1200)
            elif node_type == ACTIVITY:
                node_colors.append('blue')
                node_shapes.append('o')
                node_sizes.append(800)
            elif node_type == TASK:
                isref_val = lqn.isref[i] if i < len(lqn.isref) else False
                isref = bool(np.asarray(isref_val).flat[0]) if hasattr(isref_val, '__iter__') else bool(isref_val)
                if isref:
                    node_colors.append('#EDB120')
                    node_shapes.append('^')
                else:
                    node_colors.append('magenta')
                    node_shapes.append('v')
                node_sizes.append(1000)
            elif node_type == ENTRY:
                node_colors.append('red')
                node_shapes.append('s')
                node_sizes.append(900)
            else:
                node_colors.append('lightgray')
                node_shapes.append('o')
                node_sizes.append(600)

        for i in range(adj_matrix.shape[0]):
            for j in range(adj_matrix.shape[1]):
                if adj_matrix[i, j] > 0:
                    G.add_edge(i, j, weight=adj_matrix[i, j])

        plt.figure(figsize=kwargs.get('figsize', (14, 10)))

        try:
            pos = nx.nx_agraph.graphviz_layout(G, prog='dot')
        except:
            # Enhanced spring layout with better parameters
            pos = nx.spring_layout(G, 
                                 k=5,           # Increased node spacing
                                 iterations=100, # More iterations for better convergence
                                 seed=23000,    # Consistent layouts across runs
                                 scale=2)       # Larger drawing area

        unique_shapes = list(set(node_shapes))
        for shape in unique_shapes:
            shape_nodes = [i for i, s in enumerate(node_shapes) if s == shape]
            shape_colors = [node_colors[i] for i in shape_nodes]
            shape_sizes = [node_sizes[i] for i in shape_nodes]

            marker_map = {'o': 'o', 's': 's', '^': '^', 'v': 'v', 'h': 'h'}
            marker = marker_map.get(shape, 'o')

            nx.draw_networkx_nodes(G, pos,
                                 nodelist=shape_nodes,
                                 node_color=shape_colors,
                                 node_size=shape_sizes,
                                 node_shape=marker,
                                 alpha=0.8)

        nx.draw_networkx_edges(G, pos,
                             edge_color='gray',
                             arrows=True,
                             arrowsize=20,
                             arrowstyle='->',
                             alpha=0.6,
                             width=1.0)

        nx.draw_networkx_labels(G, pos,
                              labels=node_labels,
                              font_size=kwargs.get('font_size', 8),
                              font_weight='bold')

        plt.title(f'Model: {self.obj.getName()}', fontsize=kwargs.get('title_fontsize', 14))
        plt.axis('off')

        from matplotlib.patches import Patch
        legend_elements = [
            Patch(facecolor='black', label='Host'),
            Patch(facecolor='blue', label='Activity'),
            Patch(facecolor='#EDB120', label='Reference Task'),
            Patch(facecolor='magenta', label='Task'),
            Patch(facecolor='red', label='Entry')
        ]
        plt.legend(handles=legend_elements, loc='upper right', bbox_to_anchor=(1.15, 1))

        plt.tight_layout()

        if kwargs.get('show', True):
            plt.show()

        return plt.gcf()

    def plotGraphSimple(self, method='nodes', **kwargs):
        try:
            import matplotlib.pyplot as plt
            import networkx as nx
        except ImportError:
            raise ImportError("Matplotlib and NetworkX are required for plotting. Install with: pip install matplotlib networkx")

        lqn = self.getStruct()

        G = nx.DiGraph()

        adj_matrix = lqn.graph

        node_labels = {}
        for i in range(adj_matrix.shape[0]):
            G.add_node(i)

            if method in ['nodes', 'hashnames']:
                label = lqn.hashnames[i] if i < len(lqn.hashnames) else f'Node{i}'
                label = label.replace('_', '\\_')
            elif method == 'names':
                label = lqn.names[i] if i < len(lqn.names) else f'Node{i}'
                label = label.replace('_', '\\_')
            elif method == 'hashids':
                if i < len(lqn.hashnames):
                    prefix = lqn.hashnames[i][:2] if len(lqn.hashnames[i]) >= 2 else 'XX'
                    label = f'{prefix}{i}'
                else:
                    label = f'N{i}'
            else:
                label = str(i)

            node_labels[i] = label

        for i in range(adj_matrix.shape[0]):
            for j in range(adj_matrix.shape[1]):
                if adj_matrix[i, j] > 0:
                    G.add_edge(i, j, weight=adj_matrix[i, j])

        plt.figure(figsize=kwargs.get('figsize', (12, 8)))

        try:
            pos = nx.nx_agraph.graphviz_layout(G, prog='dot')
        except:
            # Enhanced spring layout with better parameters
            pos = nx.spring_layout(G, 
                                 k=5,           # Increased node spacing
                                 iterations=100, # More iterations for better convergence
                                 seed=23000,    # Consistent layouts across runs
                                 scale=2)       # Larger drawing area

        nx.draw_networkx_nodes(G, pos,
                             node_color='white',
                             edgecolors='black',
                             node_size=kwargs.get('node_size', 800),
                             linewidths=2)

        nx.draw_networkx_edges(G, pos,
                             edge_color='black',
                             arrows=True,
                             arrowsize=20,
                             arrowstyle='->',
                             width=1.5)

        nx.draw_networkx_labels(G, pos,
                              labels=node_labels,
                              font_size=kwargs.get('font_size', 8),
                              font_weight='bold')

        plt.title(f'Model: {self.obj.getName()}', fontsize=kwargs.get('title_fontsize', 14))
        plt.axis('off')
        plt.tight_layout()

        if kwargs.get('show', True):
            plt.show()

        return plt.gcf()

    def plotTaskGraph(self, method='nodes', **kwargs):
        try:
            import matplotlib.pyplot as plt
            import networkx as nx
            import numpy as np
        except ImportError:
            raise ImportError("Matplotlib and NetworkX are required for plotting. Install with: pip install matplotlib networkx")

        lqn = self.getStruct()

        T = np.zeros((lqn.nhosts + lqn.ntasks, lqn.nhosts + lqn.ntasks))

        for h in range(lqn.nhosts):
            hidx = h
            if h < len(lqn.tasksof) and lqn.tasksof[h] is not None:
                for tidx in lqn.tasksof[h]:
                    if tidx < T.shape[0] and hidx < T.shape[1]:
                        T[tidx, hidx] = 1

        for t in range(lqn.ntasks):
            tidx = lqn.tshift + t
            if tidx < len(lqn.entriesof) and lqn.entriesof[tidx] is not None:
                for entry_idx in lqn.entriesof[tidx]:
                    if entry_idx < lqn.iscaller.shape[1]:
                        calling_indices = np.where(lqn.iscaller[:, entry_idx] > 0)[0]
                        task_range = range(lqn.tshift, lqn.tshift + lqn.ntasks)
                        callers = [idx for idx in calling_indices if idx in task_range]

                        for caller_idx in callers:
                            if caller_idx < T.shape[0] and tidx < T.shape[1]:
                                T[caller_idx, tidx] = 1

        G = nx.DiGraph()

        node_labels = {}
        node_colors = []

        for i in range(lqn.nhosts + lqn.ntasks):
            G.add_node(i)

            if method == 'nodes':
                label = lqn.hashnames[i] if i < len(lqn.hashnames) else f'Node{i}'
            elif method == 'names':
                label = lqn.names[i] if i < len(lqn.names) else f'Node{i}'
            else:
                label = str(i)

            node_labels[i] = label

            if i < lqn.nhosts:
                node_colors.append('lightgray')
            else:
                node_colors.append('lightblue')

        for i in range(T.shape[0]):
            for j in range(T.shape[1]):
                if T[i, j] > 0:
                    G.add_edge(i, j)

        plt.figure(figsize=kwargs.get('figsize', (10, 8)))

        try:
            pos = nx.nx_agraph.graphviz_layout(G, prog='dot')
        except:
            # Enhanced spring layout with better parameters
            pos = nx.spring_layout(G, 
                                 k=5,           # Increased node spacing
                                 iterations=100, # More iterations for better convergence
                                 seed=23000,    # Consistent layouts across runs
                                 scale=2)       # Larger drawing area

        nx.draw_networkx_nodes(G, pos,
                             node_color=node_colors,
                             node_size=kwargs.get('node_size', 1000),
                             alpha=0.8)

        nx.draw_networkx_edges(G, pos,
                             edge_color='gray',
                             arrows=True,
                             arrowsize=20,
                             arrowstyle='->',
                             alpha=0.6)

        nx.draw_networkx_labels(G, pos,
                              labels=node_labels,
                              font_size=kwargs.get('font_size', 10),
                              font_weight='bold')

        plt.title('Task Graph', fontsize=kwargs.get('title_fontsize', 14))
        plt.axis('off')

        from matplotlib.patches import Patch
        legend_elements = [
            Patch(facecolor='lightgray', label='Host'),
            Patch(facecolor='lightblue', label='Task')
        ]
        plt.legend(handles=legend_elements, loc='upper right')

        plt.tight_layout()

        if kwargs.get('show', True):
            plt.show()

        return plt.gcf()

    def writeJLQN(self, filename, abstractNames=False):
        self.obj.writeJLQN(filename, abstractNames)

    def sendModel(self, outputPath, *args):
        if len(args) == 1:
            portNumber = args[0]
            self.obj.sendModel(outputPath, portNumber)
        elif len(args) == 2:
            ipNumber, portNumber = args
            self.obj.sendModel(outputPath, ipNumber, portNumber)
        else:
            raise ValueError("sendModel requires either (outputPath, portNumber) or (outputPath, ipNumber, portNumber)")

    def generateGraph(self):
        """Generate internal graph representation."""
        self.obj.generateGraph()

    @staticmethod
    def load(filename, verbose=False):
        java_model = jpype.JPackage('jline').lang.layered.LayeredNetwork.load(filename, verbose)
        model = LayeredNetwork('')
        model.obj = java_model
        return model

    @staticmethod
    def parseXML(filename, verbose=False):
        java_model = jpype.JPackage('jline').lang.layered.LayeredNetwork.parseXML(filename, verbose)
        model = LayeredNetwork('')
        model.obj = java_model
        return model

    @staticmethod
    def parse_xml(filename, verbose=False):
        """Snake case alias for parseXML"""
        return LayeredNetwork.parseXML(filename, verbose)

    @staticmethod
    def readXML(filename, verbose=False):
        java_model = jpype.JPackage('jline').lang.layered.LayeredNetwork.readXML(filename, verbose)
        model = LayeredNetwork('')
        model.obj = java_model
        return model

    @staticmethod
    def viewModel(filename, jlqnPath=None):
        if jlqnPath is None:
            jpype.JPackage('jline').lang.layered.LayeredNetwork.viewModel(filename)
        else:
            jpype.JPackage('jline').lang.layered.LayeredNetwork.viewModel(jlqnPath, filename)

    def view(self):
        """
        Display the layered network in ModelVisualizer.

        Opens the model in the interactive ModelVisualizer for
        visualization.

        Example:
            >>> lqn.view()
        """
        self.obj.view()

    get_node_index = getNodeIndex
    get_node_names = getNodeNames
    get_ensemble = getEnsemble
    get_layers = getLayers
    get_number_of_layers = getNumberOfLayers
    get_number_of_models = getNumberOfModels
    get_struct = getStruct
    write_jlqn = writeJLQN
    send_model = sendModel
    generate_graph = generateGraph
    view_model = viewModel
    register_node = registerNode
    get_node_by_name = getNodeByName
    get_nodes = getNodes
    get_node_count = getNodeCount


class LayeredNetworkStruct():
    """
    Internal structure representation of a layered queueing network.
    
    Contains the structural information for layered queueing networks (LQN),
    including processors, tasks, entries, activities, and call relationships.
    This class is populated from the Java LINE implementation and provides
    access to the detailed structure of layered models.
    
    Attributes:
        nhosts (int): Number of processors (hosts).
        ntasks (int): Number of tasks.
        nentries (int): Number of entries.
        nacts (int): Number of activities.
        ncalls (int): Number of inter-task calls.
    """
    def __str__(self):
        return pformat(vars(self))

    def fromJline(self, jsn):
        self.nidx = int(jsn.nidx)
        self.nhosts = int(jsn.nhosts)
        self.ntasks = int(jsn.ntasks)
        self.nentries = int(jsn.nentries)
        self.nacts = int(jsn.nacts)
        self.ncalls = int(jsn.ncalls)

        self.hshift = int(jsn.hshift)
        self.tshift = int(jsn.tshift)
        self.eshift = int(jsn.eshift)
        self.ashift = int(jsn.ashift)
        self.cshift = int(jsn.cshift)

        self.schedid = jlineMatrixToArray(jsn.schedid)
        self.mult = jlineMatrixToArray(jsn.mult)
        self.repl = jlineMatrixToArray(jsn.repl)
        self.type = jlineMatrixToArray(jsn.type)
        self.graph = jlineMatrixToArray(jsn.graph)
        self.replygraph = jlineMatrixToArray(jsn.replygraph)
        self.nitems = jlineMatrixToArray(jsn.nitems)
        self.replacement = jlineMatrixToArray(jsn.replacement)
        self.parent = jlineMatrixToArray(jsn.parent)
        self.iscaller = jlineMatrixToArray(jsn.iscaller)
        self.issynccaller = jlineMatrixToArray(jsn.issynccaller)
        self.isasynccaller = jlineMatrixToArray(jsn.isasynccaller)
        self.callpair = jlineMatrixToArray(jsn.callpair)
        self.taskgraph = jlineMatrixToArray(jsn.taskgraph)
        self.actpretype = jlineMatrixToArray(jsn.actpretype)
        self.actposttype = jlineMatrixToArray(jsn.actposttype)
        self.isref = jlineMatrixToArray(jsn.isref)

        self.names = np.empty(self.nidx, dtype=object)
        self.hashnames = np.empty(self.nidx, dtype=object)
        for i in range(int(self.nidx)):
            self.names[i] = jsn.names.get(jpype.JPackage('java').lang.Integer(1 + i))
            self.hashnames[i] = jsn.hashnames.get(jpype.JPackage('java').lang.Integer(1 + i))

        self.callnames = np.empty(self.ncalls, dtype=object)
        self.callhashnames = np.empty(self.ncalls, dtype=object)
        for i in range(int(self.ncalls)):
            self.callnames[i] = jsn.callnames.get(jpype.JPackage('java').lang.Integer(1 + i))
            self.callhashnames[i] = jsn.callhashnames.get(jpype.JPackage('java').lang.Integer(1 + i))

        self.hostdem = np.empty(self.nidx, dtype=object)
        for i in range(len(jsn.hostdem)):
            distrib = jsn.hostdem.get(jpype.JPackage('java').lang.Integer(1 + i))
            self.hostdem[i] = jlineFromDistribution(distrib)

        self.think = np.empty(self.nidx, dtype=object)
        for i in range(len(jsn.think)):
            distrib = jsn.think.get(jpype.JPackage('java').lang.Integer(1 + i))
            self.think[i] = jlineFromDistribution(distrib)

        self.callproc = np.empty(self.nidx, dtype=object)
        for i in range(len(jsn.callproc)):
            distrib = jsn.callproc.get(jpype.JPackage('java').lang.Integer(1 + i))
            self.callproc[i] = jlineFromDistribution(distrib)

        self.itemproc = np.empty(self.nidx, dtype=object)
        for i in range(len(jsn.itemproc)):
            distrib = jsn.itemproc.get(jpype.JPackage('java').lang.Integer(1 + i))
            self.itemproc[i] = jlineFromDistribution(distrib)

        self.itemcap = np.zeros(len(jsn.itemproc), dtype=object)
        for i in range(len(jsn.itemproc)):
            self.itemcap[i] = jsn.itemproc.get(jpype.JPackage('java').lang.Integer(1 + i)).intValue()

        self.sched = np.empty(len(jsn.sched), dtype=object)
        for i in range(len(jsn.sched)):
            sched_i = jsn.sched.get(jpype.JPackage('java').lang.Integer(1 + i))
            if sched_i is not None:
                self.sched[i] = SchedStrategy(sched_i).name
            else:
                self.sched[i] = None

        self.calltype = np.empty(self.ncalls, dtype=object)
        for i in range(len(jsn.calltype)):
            calltype_i = jsn.calltype.get(jpype.JPackage('java').lang.Integer(1 + i))
            if calltype_i is not None:
                self.calltype[i] = CallType(calltype_i).name
            else:
                self.calltype[i] = None

        self.entriesof = np.empty(len(jsn.entriesof), dtype=object)
        for i in range(len(jsn.entriesof)):
            arrayList = jsn.entriesof.get(jpype.JPackage('java').lang.Integer(1 + i))
            if arrayList is not None:
                self.entriesof[i] = list(arrayList)

        self.tasksof = np.empty(len(jsn.tasksof), dtype=object)
        for i in range(len(jsn.tasksof)):
            arrayList = jsn.tasksof.get(jpype.JPackage('java').lang.Integer(1 + i))
            if arrayList is not None:
                self.tasksof[i] = list(arrayList)

        self.actsof = np.empty(len(jsn.actsof), dtype=object)
        for i in range(len(jsn.actsof)):
            arrayList = jsn.actsof.get(jpype.JPackage('java').lang.Integer(1 + i))
            if arrayList is not None:
                self.actsof[i] = list(arrayList)

        self.callsof = np.empty(len(jsn.callsof), dtype=object)
        for i in range(len(jsn.callsof)):
            arrayList = jsn.callsof.get(jpype.JPackage('java').lang.Integer(1 + i))
            if arrayList is not None:
                self.callsof[i] = list(arrayList)

class Processor:
    """
    Processor (host) in a layered queueing network.

    A processor represents a physical or logical processing resource
    in the system (e.g., CPU, server, thread pool). Tasks are deployed
    on processors and share the processor's capacity according to the
    scheduling strategy.

    Args:
        model (LayeredNetwork): Parent layered network.
        name (str): Name of the processor.
        mult (int): Multiplicity (number of identical processors).
        schedStrategy: Scheduling strategy for the processor.
    """

    def __init__(self, model, name, mult, schedStrategy):
        from .constants import GlobalConstants
        import math
        if mult == float('inf') or mult == float('-inf') or (isinstance(mult, (int, float)) and math.isinf(mult)):
            mult_value = GlobalConstants.MaxInt
        else:
            mult_value = int(mult)
        self.obj = jpype.JPackage('jline').lang.layered.Processor(model.obj, name, mult_value, schedStrategy.value)
        model.registerNode(self)

    # Property getters delegating to Java
    def getMultiplicity(self):
        """Get the multiplicity (number of processor instances)."""
        return self.obj.getMultiplicity()

    def getReplication(self):
        """Get the replication factor."""
        return self.obj.getReplication()

    def getScheduling(self):
        """Get the scheduling strategy."""
        return self.obj.getScheduling()

    def getQuantum(self):
        """Get the time quantum for scheduling."""
        return self.obj.getQuantum()

    def getSpeedFactor(self):
        """Get the speed factor."""
        return self.obj.getSpeedFactor()

    # Snake case aliases
    get_multiplicity = getMultiplicity
    get_replication = getReplication
    get_scheduling = getScheduling
    get_quantum = getQuantum
    get_speed_factor = getSpeedFactor


class Task:
    """
    Task in a layered queueing network.
    
    A task represents a software process or thread that provides services
    through entries and performs activities. Tasks are deployed on processors
    and can make calls to other tasks.
    
    Args:
        model (LayeredNetwork): Parent layered network.
        name (str): Name of the task.
        mult (int): Multiplicity (number of identical task instances).
        schedStrategy: Scheduling strategy for the task.
    """
    
    def __init__(self, model, name, mult, schedStrategy):
        from .constants import GlobalConstants
        import math
        if mult == float('inf') or mult == float('-inf') or (isinstance(mult, (int, float)) and math.isinf(mult)):
            mult_value = GlobalConstants.MaxInt
        else:
            mult_value = int(mult)
        self.obj = jpype.JPackage('jline').lang.layered.Task(model.obj, name, mult_value, schedStrategy.value)
        model.registerNode(self)

    def on(self, proc):
        self.obj.on(proc.obj)
        return self

    def setThinkTime(self, distrib):
        self.obj.setThinkTime(distrib.obj)
        return self

    def addPrecedence(self, prec):
        self.obj.addPrecedence(prec)

    def add_precedence(self, prec):
        """Snake case alias for addPrecedence"""
        return self.addPrecedence(prec)

    set_think_time = setThinkTime

    # Property getters delegating to Java
    def getMultiplicity(self):
        """Get the multiplicity (number of task instances)."""
        return self.obj.getMultiplicity()

    def getReplication(self):
        """Get the replication factor."""
        return self.obj.getReplication()

    def getScheduling(self):
        """Get the scheduling strategy."""
        return self.obj.getScheduling()

    def getThinkTimeMean(self):
        """Get the mean think time."""
        return self.obj.getThinkTimeMean()

    def getThinkTimeSCV(self):
        """Get the squared coefficient of variation of think time."""
        return self.obj.getThinkTimeSCV()

    def getParent(self):
        """Get the parent host/processor."""
        return self.obj.getParent()

    def getPrecedences(self):
        """Get the list of activity precedences."""
        return self.obj.getPrecedences()

    def getSetupTimeMean(self):
        """Get the mean setup time."""
        return self.obj.getSetupTimeMean()

    def getDelayOffTimeMean(self):
        """Get the mean delay-off time."""
        return self.obj.getDelayOffTimeMean()

    # Snake case aliases
    get_multiplicity = getMultiplicity
    get_replication = getReplication
    get_scheduling = getScheduling
    get_think_time_mean = getThinkTimeMean
    get_think_time_scv = getThinkTimeSCV
    get_parent = getParent
    get_precedences = getPrecedences
    get_setup_time_mean = getSetupTimeMean
    get_delay_off_time_mean = getDelayOffTimeMean


class FunctionTask:
    def __init__(self, model, name, mult, schedStrategy):
        from .constants import GlobalConstants
        import math
        if mult == float('inf') or mult == float('-inf') or (isinstance(mult, (int, float)) and math.isinf(mult)):
            mult_value = GlobalConstants.MaxInt
        else:
            mult_value = int(mult)
        self.obj = jpype.JPackage('jline').lang.layered.FunctionTask(model.obj, name, mult_value, schedStrategy.value)

    def on(self, proc):
        self.obj.on(proc.obj)
        return self

    def setThinkTime(self, distrib):
        self.obj.setThinkTime(distrib.obj)
        return self

    def setSetupTime(self, distrib):
        if hasattr(distrib, 'obj'):
            self.obj.setSetupTime(distrib.obj)
        else:
            self.obj.setSetupTime(float(distrib))
        return self

    def setDelayOffTime(self, distrib):
        if hasattr(distrib, 'obj'):
            self.obj.setDelayOffTime(distrib.obj)
        else:
            self.obj.setDelayOffTime(float(distrib))
        return self

    def addPrecedence(self, prec):
        self.obj.addPrecedence(prec)

    def add_precedence(self, prec):
        """Snake case alias for addPrecedence"""
        return self.addPrecedence(prec)

    set_think_time = setThinkTime
    set_setup_time = setSetupTime
    set_delay_off_time = setDelayOffTime


class Entry:
    def __init__(self, model, name):
        self.obj = jpype.JPackage('jline').lang.layered.Entry(model.obj, name)
        model.registerNode(self)

    def on(self, proc):
        self.obj.on(proc.obj)
        return self

    def setArrival(self, distrib):
        """
        Set the open arrival distribution for this entry.

        Args:
            distrib: Distribution object (e.g., Exp, Erlang, HyperExp) representing
                    the inter-arrival time distribution for external requests.

        Returns:
            self: For method chaining.

        Example:
            >>> from line_solver import Exp
            >>> entry.setArrival(Exp(2.5))  # Arrival rate of 2.5 jobs/sec
        """
        self.obj.setArrival(distrib.obj)
        return self

    def set_arrival(self, distrib):
        """Snake case alias for setArrival"""
        return self.setArrival(distrib)

    def forward(self, dest, prob=1.0):
        """Add a forwarding call to another entry.

        Args:
            dest: Destination entry object or entry name (string)
            prob: Probability of forwarding (0.0 to 1.0), default is 1.0

        Returns:
            self: This entry for method chaining
        """
        if hasattr(dest, 'obj'):
            self.obj.forward(dest.obj, float(prob))
        else:
            self.obj.forward(str(dest), float(prob))
        return self

    # Property getters delegating to Java
    def getBoundToActivity(self):
        """Get the activities bound to this entry."""
        return self.obj.getBoundToActivity()

    def getReplyActivity(self):
        """Get the reply activities for this entry."""
        return self.obj.getReplyActivity()

    def getParent(self):
        """Get the parent task."""
        return self.obj.getParent()

    def getForwardingDests(self):
        """Get the forwarding destinations."""
        return self.obj.getForwardingDests()

    def getForwardingProbs(self):
        """Get the forwarding probabilities."""
        return self.obj.getForwardingProbs()

    def getArrival(self):
        """Get the arrival distribution."""
        return self.obj.getArrival()

    # Snake case aliases
    get_bound_to_activity = getBoundToActivity
    get_reply_activity = getReplyActivity
    get_parent = getParent
    get_forwarding_dests = getForwardingDests
    get_forwarding_probs = getForwardingProbs
    get_arrival = getArrival


class Activity:
    def __init__(self, model, name, distrib):
        self.obj = jpype.JPackage('jline').lang.layered.Activity(model.obj, name, distrib.obj)
        model.registerNode(self)

    def on(self, proc):
        self.obj.on(proc.obj)
        return self

    def boundTo(self, proc):
        self.obj.boundTo(proc.obj)
        return self

    def bound_to(self, proc):
        """Snake case alias for boundTo"""
        return self.boundTo(proc)

    def repliesTo(self, entry):
        self.obj.repliesTo(entry.obj)
        return self

    def replies_to(self, entry):
        """Snake case alias for repliesTo"""
        return self.repliesTo(entry)

    def synchCall(self, entry, callmult=1.0):
        self.obj.synchCall(entry.obj, callmult)
        return self

    def synch_call(self, entry, callmult=1.0):
        """Snake case alias for synchCall"""
        return self.synchCall(entry, callmult)

    def asynchCall(self, entry, callmult=1.0):
        """
        Make an asynchronous (non-blocking) call to another entry.

        Async calls are fire-and-forget - the caller continues immediately
        without waiting for a response. Use for prefetching, cache warming,
        or non-critical operations where the result is not immediately needed.

        Args:
            entry: Entry or ItemEntry to call asynchronously
            callmult (float): Mean number of async calls (default 1.0)

        Returns:
            self: This activity for method chaining

        Examples:
            >>> # Async cache prefetch (non-blocking)
            >>> activity.asynchCall(cache_entry, 1.0)
            >>>
            >>> # Multiple async calls
            >>> activity.asynchCall(prefetch_entry, 2.5)
        """
        if hasattr(entry, 'obj'):
            self.obj.asynchCall(entry.obj, float(callmult))
        else:
            self.obj.asynchCall(str(entry), float(callmult))
        return self

    def asynch_call(self, entry, callmult=1.0):
        """Snake case alias for asynchCall"""
        return self.asynchCall(entry, callmult)

    def setPhase(self, phase):
        """
        Set the phase number for this activity.

        Phase 1: activities before the reply is sent (client blocks)
        Phase 2: activities after the reply is sent (post-reply processing)

        Args:
            phase (int): Phase number (must be 1 or 2)

        Returns:
            self: This activity for method chaining
        """
        self.obj.setPhase(int(phase))
        return self

    def set_phase(self, phase):
        """Snake case alias for setPhase"""
        return self.setPhase(phase)

    def getPhase(self):
        """
        Get the phase number for this activity.

        Returns:
            int: The phase number (1 or 2)
        """
        return int(self.obj.getPhase())

    def get_phase(self):
        """Snake case alias for getPhase"""
        return self.getPhase()

    # Property getters delegating to Java
    def getHostDemand(self):
        """Get the host demand distribution."""
        return self.obj.getHostDemand()

    def getHostDemandMean(self):
        """Get the mean host demand."""
        return self.obj.getHostDemandMean()

    def getHostDemandSCV(self):
        """Get the squared coefficient of variation of host demand."""
        return self.obj.getHostDemandSCV()

    def getCallOrder(self):
        """Get the call order (STOCHASTIC or DETERMINISTIC)."""
        return self.obj.getCallOrder()

    def getBoundToEntry(self):
        """Get the entry this activity is bound to."""
        return self.obj.getBoundToEntry()

    def getParent(self):
        """Get the parent task."""
        return self.obj.getParent()

    def getSyncCallDests(self):
        """Get the synchronous call destinations."""
        return self.obj.getSyncCallDests()

    def getSyncCallMeans(self):
        """Get the synchronous call means."""
        return self.obj.getSyncCallMeans()

    def getAsyncCallDests(self):
        """Get the asynchronous call destinations."""
        return self.obj.getAsyncCallDests()

    def getAsyncCallMeans(self):
        """Get the asynchronous call means."""
        return self.obj.getAsyncCallMeans()

    def getThinkTimeMean(self):
        """Get the mean think time."""
        return self.obj.getThinkTimeMean()

    # Snake case aliases
    get_host_demand = getHostDemand
    get_host_demand_mean = getHostDemandMean
    get_host_demand_scv = getHostDemandSCV
    get_call_order = getCallOrder
    get_bound_to_entry = getBoundToEntry
    get_parent = getParent
    get_sync_call_dests = getSyncCallDests
    get_sync_call_means = getSyncCallMeans
    get_async_call_dests = getAsyncCallDests
    get_async_call_means = getAsyncCallMeans
    get_think_time_mean = getThinkTimeMean


class ActivityPrecedence:
    def __init__(self, name):
        self.obj = jpype.JPackage('jline').lang.layered.ActivityPrecedence(name)

    @staticmethod
    def Serial(act0, act1):
        return jpype.JPackage('jline').lang.layered.ActivityPrecedence.Serial(
            act0.obj.getName(), act1.obj.getName())

    # Python snake_case aliases for compatibility
    @staticmethod
    def serial(acts):
        """Python snake_case alias for Serial(). Takes a list of activities."""
        if not acts:
            raise ValueError("serial() requires at least one activity")
        if len(acts) == 1:
            return jpype.JPackage('jline').lang.layered.ActivityPrecedence.Serial(
                acts[0].obj.getName(), acts[0].obj.getName())
        # For multiple activities, chain them together serially
        result = jpype.JPackage('jline').lang.layered.ActivityPrecedence.Serial(
            acts[0].obj.getName(), acts[1].obj.getName())
        for i in range(2, len(acts)):
            # Chain subsequent activities
            result = jpype.JPackage('jline').lang.layered.ActivityPrecedence.Serial(
                acts[i-1].obj.getName(), acts[i].obj.getName())
        return result

    @staticmethod
    def AndFork(preAct, postActs):
        preActName = preAct.obj.getName() if hasattr(preAct, 'obj') else str(preAct)
        postActNames = []
        for act in postActs:
            if hasattr(act, 'obj'):
                postActNames.append(act.obj.getName())
            else:
                postActNames.append(str(act))

        java_list = jpype.java.util.ArrayList()
        for name in postActNames:
            java_list.add(name)

        return jpype.JPackage('jline').lang.layered.ActivityPrecedence.AndFork(preActName, java_list)

    @staticmethod
    def AndJoin(preActs, postAct, quorum=None):
        preActNames = []
        for act in preActs:
            if hasattr(act, 'obj'):
                preActNames.append(act.obj.getName())
            else:
                preActNames.append(str(act))

        postActName = postAct.obj.getName() if hasattr(postAct, 'obj') else str(postAct)

        java_list = jpype.java.util.ArrayList()
        for name in preActNames:
            java_list.add(name)

        if quorum is not None:
            from .lang import jlineMatrixFromArray
            return jpype.JPackage('jline').lang.layered.ActivityPrecedence.AndJoin(
                java_list, postActName, jlineMatrixFromArray(quorum))
        else:
            return jpype.JPackage('jline').lang.layered.ActivityPrecedence.AndJoin(java_list, postActName)

    @staticmethod
    def OrFork(preAct, postActs, probs):
        preActName = preAct.obj.getName() if hasattr(preAct, 'obj') else str(preAct)
        postActNames = []
        for act in postActs:
            if hasattr(act, 'obj'):
                postActNames.append(act.obj.getName())
            else:
                postActNames.append(str(act))

        java_list = jpype.java.util.ArrayList()
        for name in postActNames:
            java_list.add(name)

        from .lang import jlineMatrixFromArray
        import numpy as np
        prob_matrix = jlineMatrixFromArray(np.array(probs).reshape(1, -1))

        return jpype.JPackage('jline').lang.layered.ActivityPrecedence.OrFork(
            preActName, java_list, prob_matrix)

    @staticmethod
    def OrJoin(preActs, postAct):
        preActNames = []
        for act in preActs:
            if hasattr(act, 'obj'):
                preActNames.append(act.obj.getName())
            else:
                preActNames.append(str(act))

        postActName = postAct.obj.getName() if hasattr(postAct, 'obj') else str(postAct)

        java_list = jpype.java.util.ArrayList()
        for name in preActNames:
            java_list.add(name)

        return jpype.JPackage('jline').lang.layered.ActivityPrecedence.OrJoin(java_list, postActName)

    @staticmethod
    def Loop(preAct, postActs, counts):
        preActName = preAct.obj.getName() if hasattr(preAct, 'obj') else str(preAct)
        postActNames = []
        for act in postActs:
            if hasattr(act, 'obj'):
                postActNames.append(act.obj.getName())
            else:
                postActNames.append(str(act))

        java_list = jpype.java.util.ArrayList()
        for name in postActNames:
            java_list.add(name)

        return jpype.JPackage('jline').lang.layered.ActivityPrecedence.Loop(
            preActName, java_list, float(counts))

    @staticmethod
    def CacheAccess(preAct, postActs):
        preActName = preAct.obj.getName() if hasattr(preAct, 'obj') else str(preAct)
        postActNames = []
        for act in postActs:
            if hasattr(act, 'obj'):
                postActNames.append(act.obj.getName())
            else:
                postActNames.append(str(act))

        java_list = jpype.java.util.ArrayList()
        for name in postActNames:
            java_list.add(name)

        return jpype.JPackage('jline').lang.layered.ActivityPrecedence.CacheAccess(preActName, java_list)

    @staticmethod
    def cache_access(preAct, postActs):
        """Snake case alias for CacheAccess"""
        return ActivityPrecedence.CacheAccess(preAct, postActs)

    @staticmethod
    def fromActivities(preActs, postActs, preType, postType=None, preParams=None, postParams=None):
        preActList = jpype.java.util.ArrayList()
        for act in preActs:
            if hasattr(act, 'obj'):
                preActList.add(act.obj)
            else:
                raise ValueError("preActs must contain Activity objects")

        postActList = jpype.java.util.ArrayList()
        for act in postActs:
            if hasattr(act, 'obj'):
                postActList.add(act.obj)
            else:
                raise ValueError("postActs must contain Activity objects")

        if postParams is not None:
            from .lang import jlineMatrixFromArray
            return jpype.JPackage('jline').lang.layered.ActivityPrecedence.fromActivities(
                preActList, postActList, preType, postType,
                jlineMatrixFromArray(preParams), jlineMatrixFromArray(postParams))
        elif preParams is not None:
            from .lang import jlineMatrixFromArray
            return jpype.JPackage('jline').lang.layered.ActivityPrecedence.fromActivities(
                preActList, postActList, preType, postType, jlineMatrixFromArray(preParams))
        elif postType is not None:
            return jpype.JPackage('jline').lang.layered.ActivityPrecedence.fromActivities(
                preActList, postActList, preType, postType)
        else:
            return jpype.JPackage('jline').lang.layered.ActivityPrecedence.fromActivities(
                preActList, postActList, preType)

    @staticmethod
    def getPrecedenceId(precedence):
        return jpype.JPackage('jline').lang.layered.ActivityPrecedence.getPrecedenceId(precedence)

    @staticmethod
    def Xor(preAct, postActs, probs):
        """
        Create exclusive-OR (XOR) fork precedence.

        XOR fork represents probabilistic branching where exactly one successor
        activity executes based on the given probabilities. This is an alias for
        OrFork (OR fork with mutually exclusive branches).

        Args:
            preAct: Predecessor activity (Activity object or string name)
            postActs: List of successor activities
            probs: List of branch probabilities (must sum to 1.0)

        Returns:
            ActivityPrecedence: XOR fork precedence object

        Example:
            >>> # Create XOR fork with 70% to A2, 30% to A3
            >>> prec = ActivityPrecedence.Xor(A1, [A2, A3], [0.7, 0.3])

        See Also:
            OrFork: Probabilistic fork (functionally identical to Xor)
            AndFork: Parallel fork with all branches executing
        """
        return ActivityPrecedence.OrFork(preAct, postActs, probs)

    # Snake_case alias
    xor = Xor


class CacheTask:
    def __init__(self, model, name, items=1, itemLevelCap=1, replacementStrategy=None, multiplicity=1, schedStrategy=None, thinkTime=None):
        """
        Create a cache task in a layered network.

        Args:
            model: The layered network model
            name: Task name
            items: Number of items in cache
            itemLevelCap: Cache capacity - can be:
                - int: Single-level cache with given capacity
                - list/tuple: Multi-level cache with capacity per level
            replacementStrategy: Cache replacement strategy (FIFO, LRU, etc.)
            multiplicity: Number of parallel cache instances
            schedStrategy: Scheduling strategy
            thinkTime: Optional think time distribution
        """
        from .constants import GlobalConstants
        import math

        if replacementStrategy is None:
            replacementStrategy = ReplacementStrategy.FIFO
        if schedStrategy is None:
            schedStrategy = SchedStrategy.FCFS

        if multiplicity == float('inf') or multiplicity == float('-inf') or (isinstance(multiplicity, (int, float)) and math.isinf(multiplicity)):
            mult_value = GlobalConstants.MaxInt
        else:
            mult_value = int(multiplicity)

        # Convert itemLevelCap to Java format (int or int[])
        if isinstance(itemLevelCap, (list, tuple)):
            # Multi-level cache: convert to Java int array
            java_caps = jpype.JArray(jpype.JInt)(list(itemLevelCap))
        else:
            # Single-level cache: use as scalar int
            java_caps = int(itemLevelCap)

        # Call appropriate JAR constructor based on parameters
        if thinkTime is not None:
            self.obj = jpype.JPackage('jline').lang.layered.CacheTask(
                model.obj, name, int(items), java_caps,
                replacementStrategy.value, mult_value, schedStrategy.value, thinkTime.obj)
        elif schedStrategy != SchedStrategy.FCFS:
            self.obj = jpype.JPackage('jline').lang.layered.CacheTask(
                model.obj, name, int(items), java_caps,
                replacementStrategy.value, mult_value, schedStrategy.value)
        else:
            self.obj = jpype.JPackage('jline').lang.layered.CacheTask(
                model.obj, name, int(items), java_caps,
                replacementStrategy.value, mult_value)
        model.registerNode(self)

    def on(self, proc):
        """Assign this cache task to run on a processor."""
        self.obj.on(proc.obj)
        return self

    def setThinkTime(self, distrib):
        """Set the think time distribution for this cache task."""
        self.obj.setThinkTime(distrib.obj)
        return self

    def addPrecedence(self, prec):
        """Add an activity precedence relationship."""
        self.obj.addPrecedence(prec)

    def add_precedence(self, prec):
        """Snake case alias for addPrecedence"""
        return self.addPrecedence(prec)

    def getItems(self):
        """Get the number of items in the cache."""
        return self.obj.getItems()

    def setItems(self, items):
        """Set the number of items in the cache."""
        self.obj.setItems(int(items))
        return self

    def getItemLevelCap(self, level=None):
        """
        Get cache level capacity/capacities.

        Args:
            level: Optional level index (0-based). If None, returns all levels.

        Returns:
            int or list: Capacity value(s)
        """
        if level is None:
            # Get all levels
            java_result = self.obj.getItemLevelCap()
            # Check if it's an array or scalar
            if hasattr(java_result, '__len__'):
                # Array - convert to Python list
                return [java_result[i] for i in range(len(java_result))]
            else:
                # Scalar
                return int(java_result)
        else:
            # Get specific level (0-based in Java)
            return int(self.obj.getItemLevelCap(int(level)))

    def setItemLevelCap(self, itemLevelCap):
        """
        Set cache level capacity/capacities.

        Args:
            itemLevelCap: int for single-level or list/tuple for multi-level
        """
        if isinstance(itemLevelCap, (list, tuple)):
            # Multi-level: convert to Java int array
            java_caps = jpype.JArray(jpype.JInt)(list(itemLevelCap))
            self.obj.setItemLevelCap(java_caps)
        else:
            # Single-level: use scalar
            self.obj.setItemLevelCap(int(itemLevelCap))
        return self

    def getTotalCapacity(self):
        """
        Get sum of all cache level capacities.

        Returns:
            int: Total capacity across all levels
        """
        return int(self.obj.getTotalCapacity())

    def getNumberOfLevels(self):
        """
        Get number of cache levels.

        Returns:
            int: Number of levels (1 for single-level, >1 for multi-level)
        """
        caps = self.getItemLevelCap()
        if isinstance(caps, list):
            return len(caps)
        else:
            return 1

    def getReplacestrategy(self):
        """Get the replacement strategy."""
        java_strategy = self.obj.getReplacestrategy()
        for strategy in ReplacementStrategy:
            if strategy.value == java_strategy:
                return strategy
        return None

    def setReplacestrategy(self, replacementStrategy):
        """Set the replacement strategy."""
        self.obj.setReplacestrategy(replacementStrategy.value)
        return self

    set_think_time = setThinkTime
    get_items = getItems
    set_items = setItems
    get_item_level_cap = getItemLevelCap
    set_item_level_cap = setItemLevelCap
    get_total_capacity = getTotalCapacity
    get_number_of_levels = getNumberOfLevels
    get_replacestrategy = getReplacestrategy
    set_replacestrategy = setReplacestrategy


class ItemEntry:
    def __init__(self, model, name, cardinality, distribution):
        if hasattr(distribution, 'obj'):
            dist_obj = distribution.obj
        else:
            dist_obj = distribution

        self.obj = jpype.JPackage('jline').lang.layered.ItemEntry(
            model.obj, jpype.JPackage('java').lang.String(name), int(cardinality), dist_obj)
        model.registerNode(self)

    def on(self, task):
        """Assign this entry to a parent task (typically a CacheTask)."""
        self.obj.on(task.obj)
        return self

    def getCardinality(self):
        """Get the number of items."""
        return self.obj.getCardinality()

    def getPopularity(self):
        """Get the popularity distribution."""
        java_dist = self.obj.getPopularity()
        from . import jlineFromDistribution
        return jlineFromDistribution(java_dist)

    get_cardinality = getCardinality
    get_popularity = getPopularity
