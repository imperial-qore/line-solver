
import os
import jpype
import jpype.imports
import numpy as np
import pandas as pd
from jpype import JArray

from line_solver import VerboseLevel, SolverType, jlineMatrixToArray, GlobalConstants, jlineMatrixFromArray, \
    jlineMapMatrixToArray, jlineMatrixCellToArray


# Track if JAR GlobalConstants have been initialized
_JAR_GLOBALS_INITIALIZED = False


def _initialize_jar_globals():
    """Initialize JAR global constants on first use."""
    global _JAR_GLOBALS_INITIALIZED

    if not _JAR_GLOBALS_INITIALIZED:
        try:
            jpype.JPackage('jline').lang.GlobalConstants.setLibraryAttributionShown(True)
            _JAR_GLOBALS_INITIALIZED = True
        except Exception:
            # Silently ignore if JAR is not available or GlobalConstants doesn't exist
            pass


def _get_java_network(model):
    """
    Get a Java Network object from a wrapper model.

    Args:
        model: Network model (wrapper)

    Returns:
        Java Network object
    """
    if hasattr(model, 'obj'):
        # Wrapper model - .obj is already Java
        return model.obj
    else:
        # Assume it's already a Java object
        return model


class SampleResult:
    """
    Container for sample-based simulation results.
    
    This class wraps results from simulation solvers that produce time-series
    samples or event traces. It provides access to timestamps, states, and
    event information from the simulation.
    
    Attributes:
        handle (str): Unique identifier for the result.
        t: Array of time points.
        state: State information at sampling points.
        event: List of events that occurred during simulation.
        isaggregate (bool): Whether results are aggregated across nodes.
        nodeIndex: Index of the node (for node-specific results).
        numSamples (int): Number of samples collected.
    """
    
    def __init__(self, java_result):
        self.handle = java_result.handle if hasattr(java_result, 'handle') else ""
        self.t = jlineMatrixToArray(java_result.t) if hasattr(java_result, 't') and java_result.t is not None else None
        self.state = java_result.state if hasattr(java_result, 'state') else None
        self.event = self._parse_events(java_result) if hasattr(java_result, 'event') else []
        self.isaggregate = java_result.isAggregate if hasattr(java_result, 'isAggregate') else False
        self.nodeIndex = java_result.nodeIndex if hasattr(java_result, 'nodeIndex') else None
        self.numSamples = java_result.numSamples if hasattr(java_result, 'numSamples') else 0

    def _parse_events(self, java_result):
        events = []
        if hasattr(java_result, 'event') and java_result.event is not None:
            if hasattr(java_result.event, '__iter__') and not hasattr(java_result.event, 'getNumRows'):
                for event_info in java_result.event:
                    event_obj = EventInfo()
                    event_obj.node = event_info.node if hasattr(event_info, 'node') else 0
                    event_obj.jobclass = event_info.jobclass if hasattr(event_info, 'jobclass') else 0
                    event_obj.t = event_info.t if hasattr(event_info, 't') else 0.0
                    event_obj.event = getattr(event_info, 'event', None)
                    events.append(event_obj)
            else:
                event_matrix = jlineMatrixToArray(java_result.event)
                if event_matrix is not None and len(event_matrix) > 0:
                    for row in event_matrix:
                        if len(row) >= 3:
                            event_obj = EventInfo()
                            event_obj.t = float(row[0])
                            event_obj.node = int(row[1])
                            event_obj.jobclass = 0

                            event_type_int = int(row[2])
                            if event_type_int == 1:
                                event_obj.event = "ARV"
                            elif event_type_int == 2:
                                event_obj.event = "DEP"
                            elif event_type_int == 3:
                                event_obj.event = "PHASE"
                            else:
                                event_obj.event = None

                            events.append(event_obj)
        return events


class EventInfo:
    """
    Information about a single event in a simulation trace.
    
    Attributes:
        node (int): Index of the node where the event occurred.
        jobclass (int): Index of the job class involved.
        t (float): Timestamp when the event occurred.
        event (str): Type of event ('ARV', 'DEP', 'PHASE', etc.).
    """
    
    def __init__(self):
        self.node = 0
        self.jobclass = 0
        self.t = 0.0
        self.event = None


class DistributionResult:
    """
    Container for cumulative distribution function (CDF) results.
    
    This class holds CDF data for response time and other performance metrics
    computed by analytical or simulation solvers. It provides methods to access
    distribution data for specific stations and job classes.
    
    Attributes:
        num_stations (int): Number of stations in the model.
        num_classes (int): Number of job classes.
        distribution_type (str): Type of distribution computed.
        is_transient (bool): Whether this is transient analysis.
        runtime (float): Solver runtime in seconds.
        time_points: Array of time points for transient analysis.
        cdf_data: Multi-dimensional array of CDF data.
    """
    
    def __init__(self, java_result):
        self.java_result = java_result
        self.num_stations = java_result.numStations if hasattr(java_result, 'numStations') else 0
        self.num_classes = java_result.numClasses if hasattr(java_result, 'numClasses') else 0
        self.distribution_type = java_result.distributionType if hasattr(java_result, 'distributionType') else ""
        self.is_transient = java_result.isTransient if hasattr(java_result, 'isTransient') else False
        self.runtime = java_result.runtime if hasattr(java_result, 'runtime') else 0.0
        self.time_points = jlineMatrixToArray(java_result.timePoints) if hasattr(java_result, 'timePoints') and java_result.timePoints is not None else None

        self.cdf_data = self._parse_cdf_data(java_result)

    def _parse_cdf_data(self, java_result):
        """Parse CDF data from Java result into Python format"""
        if not hasattr(java_result, 'cdfData') or java_result.cdfData is None:
            return []

        try:
            cdf_data = []
            for i in range(self.num_stations):
                station_cdfs = []
                for j in range(self.num_classes):
                    if java_result.hasCdf(i, j):
                        cdf_matrix = java_result.getCdf(i, j)
                        cdf_array = jlineMatrixToArray(cdf_matrix)
                        station_cdfs.append(cdf_array)
                    else:
                        station_cdfs.append(None)
                cdf_data.append(station_cdfs)
            return cdf_data
        except Exception as e:
            print(f"Error parsing CDF data: {e}")
            return []

    def get_cdf(self, station, job_class):
        """Get CDF for a specific station and job class"""
        if station < 0 or station >= self.num_stations or job_class < 0 or job_class >= self.num_classes:
            return None

        if self.cdf_data and station < len(self.cdf_data) and job_class < len(self.cdf_data[station]):
            return self.cdf_data[station][job_class]
        return None

    def has_cdf(self, station, job_class):
        """Check if CDF data is available for a specific station and job class"""
        if hasattr(self.java_result, 'hasCdf'):
            return self.java_result.hasCdf(station, job_class)
        return False

class ProbabilityResult:
    """
    Container for state probability results.
    
    This class holds steady-state or transient probability distributions
    over the system state space, typically computed by CTMC solvers.
    
    Attributes:
        probability: Matrix of state probabilities.
        log_normalizing_constant (float): Log of normalization constant.
        is_aggregated (bool): Whether probabilities are aggregated.
        node_index: Index of specific node (for node-specific results).
        state: State space specification.
    """
    
    def __init__(self, java_result):
        self.java_result = java_result
        self.probability = jlineMatrixToArray(java_result.probability) if hasattr(java_result, 'probability') and java_result.probability is not None else None
        self.log_normalizing_constant = java_result.logNormalizingConstant if hasattr(java_result, 'logNormalizingConstant') else 0.0
        self.is_aggregated = java_result.isAggregated if hasattr(java_result, 'isAggregated') else False
        self.node_index = java_result.nodeIndex if hasattr(java_result, 'nodeIndex') else None
        self.state = jlineMatrixToArray(java_result.state) if hasattr(java_result, 'state') and java_result.state is not None else None

    def get_probability(self):
        """Get the probability matrix as a numpy array"""
        return self.probability

    def get_log_normalizing_constant(self):
        """Get the logarithm of the normalizing constant"""
        return self.log_normalizing_constant

    def is_aggregated_result(self):
        """Check if this is an aggregated result"""
        return self.is_aggregated

    def get_node_index(self):
        """Get the node index (for node-specific results)"""
        return self.node_index

    def get_state(self):
        """Get the state specification"""
        return self.state

    def state(self):
        """Alias for get_state()"""
        return self.get_state()



class Solver:
    """
    Base class for all queueing network solvers.
    
    The Solver class provides the foundation for all analytical and simulation
    solvers in LINE. It handles common solver options like verbosity, iteration
    limits, and random seeds. Specific solvers like SolverMVA, SolverJMT, etc.
    inherit from this class.
    
    Attributes:
        solveropt: Solver options object containing algorithm parameters.
    """
    
    @staticmethod
    def defaultOptions():
        """
        Get default solver options.
        
        Returns:
            dict: Dictionary containing default solver options including:
                - keep (bool): Whether to keep intermediate results.
                - verbose (VerboseLevel): Output verbosity level.
                - cutoff (int): State space truncation parameter.
                - seed (int): Random number generator seed.
                - iter_max (int): Maximum number of iterations.
                - samples (int): Number of samples for simulation solvers.
                - method (str): Solution method ('default' for automatic selection).
        """
        options = {
            'keep': False,
            'verbose': VerboseLevel.STD,
            'cutoff': 10,
            'seed': 23000,
            'iter_max': 200,
            'samples': 10000,
            'method': 'default'
        }
        return options

    def __init__(self, options, *args, **kwargs):
        """
        Initialize a solver with the given options.
        
        Args:
            options: Solver options object.
            *args: Variable arguments for solver-specific parameters.
            **kwargs: Keyword arguments for solver options.
        """
        self.solveropt = options
        self._verbose_silent = False
        self._table_silent = False

        for key, value in kwargs.items():
            self._process_solver_option(key, value)

        if len(args) >= 1:
            ctr = 0
            while ctr < len(args):
                # Skip non-string args (e.g., numpy arrays passed by SolverENV)
                if not isinstance(args[ctr], str):
                    ctr += 1
                    continue
                if args[ctr] == 'cutoff':
                    self.solveropt.obj.cutoff(args[ctr + 1])
                    ctr += 2
                elif args[ctr] == 'method':
                    self.solveropt.obj.method(args[ctr + 1])
                    ctr += 2
                elif args[ctr] == 'exact':
                    self.solveropt.obj.method('exact')
                    ctr += 1
                elif args[ctr] == 'keep':
                    self.solveropt.obj.keep(args[ctr + 1])
                    ctr += 2
                elif args[ctr] == 'seed':
                    self.solveropt.obj.seed(args[ctr + 1])
                    ctr += 2
                elif args[ctr] == 'samples':
                    self.solveropt.obj.samples(args[ctr + 1])
                    ctr += 2
                elif args[ctr] == 'timespan':
                    self.solveropt.obj.timespan = JArray(jpype.JDouble)(args[ctr + 1])
                    ctr += 2
                elif args[ctr] == 'timestep':
                    self.solveropt.obj.timestep = jpype.JDouble(args[ctr + 1]) if args[ctr + 1] is not None else None
                    ctr += 2
                elif args[ctr] == 'verbose':
                    self._process_verbose_option(args[ctr + 1])
                    ctr += 2
                else:
                    ctr += 1

    def _process_solver_option(self, key, value):
        """Process a single solver option from keyword arguments"""
        if key == 'cutoff':
            if hasattr(value, '__iter__') and not isinstance(value, str):
                from line_solver import jlineMatrixFromArray
                self.solveropt.obj.cutoff(jlineMatrixFromArray(value))
            else:
                self.solveropt.obj.cutoff(value)
        elif key == 'method':
            self.solveropt.obj.method(str(value))
        elif key == 'keep':
            self.solveropt.obj.keep(bool(value))
        elif key == 'seed':
            self.solveropt.obj.seed(int(value))
        elif key == 'samples':
            self.solveropt.obj.samples(int(value))
        elif key == 'timespan':
            if hasattr(value, '__iter__'):
                self.solveropt.obj.timespan = JArray(jpype.JDouble)(value)
            else:
                self.solveropt.obj.timespan = value
        elif key == 'timestep':
            self.solveropt.obj.timestep = jpype.JDouble(value) if value is not None else None
        elif key == 'verbose':
            self._process_verbose_option(value)
        elif key == 'force':
            self.solveropt.obj.force(bool(value))
        elif key == 'cache':
            self.solveropt.obj.cache = bool(value)
        elif key == 'hide_immediate':
            self.solveropt.obj.hide_immediate = bool(value)
        elif key == 'iter_max':
            self.solveropt.obj.iter_max = int(value)
        elif key == 'iter_tol':
            self.solveropt.obj.iter_tol = float(value)
        elif key == 'tol':
            self.solveropt.obj.tol = float(value)
        elif key == 'lang':
            self.solveropt.obj.lang = str(value)
        elif key == 'remote':
            self.solveropt.obj.remote = bool(value)
        elif key == 'remote_endpoint':
            self.solveropt.obj.remote_endpoint = str(value)
        elif key == 'stiff':
            self.solveropt.obj.stiff = bool(value)
        elif key == 'init_sol':
            if hasattr(value, '__iter__'):
                from line_solver import jlineMatrixFromArray
                self.solveropt.obj.init_sol = jlineMatrixFromArray(value)
            else:
                self.solveropt.obj.init_sol = value
        else:
            if hasattr(self.solveropt.obj, key):
                try:
                    setattr(self.solveropt.obj, key, value)
                except (AttributeError, TypeError):
                    pass

    def _process_verbose_option(self, value):
        """Process verbose option handling both old and new formats"""
        if isinstance(value, bool):
            if value:
                self.solveropt.obj.verbose(jpype.JPackage('jline').VerboseLevel.STD)
            else:
                self.solveropt.obj.verbose(jpype.JPackage('jline').VerboseLevel.SILENT)
                self._verbose_silent = True
        else:
            if value is False:
                self.solveropt.obj.verbose(
                    jpype.JPackage('jline').VerboseLevel.SILENT)
                self._verbose_silent = True
            elif value is True:
                self.solveropt.obj.verbose(
                    jpype.JPackage('jline').VerboseLevel.STD)
                self._verbose_silent = False
            elif value == VerboseLevel.SILENT:
                self.solveropt.obj.verbose(
                    jpype.JPackage('jline').VerboseLevel.SILENT)
                self._verbose_silent = True
                self._table_silent = True
            elif value == VerboseLevel.STD:
                self.solveropt.obj.verbose(jpype.JPackage('jline').VerboseLevel.STD)
            elif value == VerboseLevel.DEBUG:
                self.solveropt.obj.verbose(jpype.JPackage('jline').VerboseLevel.DEBUG)

    def getName(self):
        """
        Get the name of this solver.
        
        Returns:
            str: The solver name (e.g., 'MVA', 'JMT', 'SSA').
        """
        return self.obj.getName()

    @staticmethod
    def defaultOptions():
        """Returns default solver options for the base Solver class."""
        java_options = jpype.JPackage('jline').solvers.mva.SolverMVA.defaultOptions()
        python_options = SolverOptions.__new__(SolverOptions)
        python_options.obj = java_options
        return python_options

    def supports(self, model):
        """
        Check if this solver supports the given model.
        
        Args:
            model: Network model to check for compatibility.
            
        Returns:
            bool: True if the solver can solve this model, False otherwise.
        """
        if hasattr(self, 'obj') and hasattr(self.obj, 'supports'):
            return self.obj.supports(model.obj if hasattr(model, 'obj') else model)
        return True

    get_name = getName
    default_options = defaultOptions

class EnsembleSolver(Solver):
    def __init__(self, options, *args, **kwargs):
        super().__init__(options, *args, **kwargs)
        pass

    def getNumberOfModels(self):
        """
        Get the number of models in the ensemble.

        Returns:
            int: Number of models in this ensemble solver
        """
        return self.obj.getNumberOfModels()

    def printEnsembleAvgTables(self):
        """
        Print average performance tables for all models in the ensemble.
        """
        self.obj.printEnsembleAvgTables()

    def printEnsembleAvgTs(self):
        """Short alias for printEnsembleAvgTables."""
        self.printEnsembleAvgTables()

    print_ensemble_avg_tables = printEnsembleAvgTables
    print_ensemble_avg_ts = printEnsembleAvgTs

    def numberOfModels(self):
        """Get the number of models in this ensemble solver.
        
        Returns:
            int: Number of models in the ensemble.
        """
        return self.getNumberOfModels()

    get_number_of_models = getNumberOfModels


class NetworkSolver(Solver):
    def __init__(self, options, *args, **kwargs):
        super().__init__(options, *args, **kwargs)
        pass

    def getAvgNodeTable(self):
        """
        Get average performance metrics per node.
        
        Returns:
            pandas.DataFrame: Performance metrics table with columns:
                - QLen: Average queue length per node
                - Util: Utilization per node
                - RespT: Response time per node
                - ResidT: Residence time per node
                - ArvR: Arrival rate per node
                - Tput: Throughput per node
        """
        table = self.obj.getAvgNodeTable()


        QLen = np.array(list(table.getQLen()))
        Util = np.array(list(table.getUtil()))
        RespT = np.array(list(table.getRespT()))
        ResidT = np.array(list(table.getResidT()))
        ArvR = np.array(list(table.getArvR()))
        Tput = np.array(list(table.getTput()))

        cols = ['QLen', 'Util', 'RespT', 'ResidT', 'ArvR', 'Tput']
        nodes = list(table.getNodeNames())
        nodenames = []
        for i in range(len(nodes)):
            nodenames.append(str(nodes[i]))
        jobclasses = list(table.getClassNames())

        classnames = []
        for i in range(len(jobclasses)):
            classnames.append(str(jobclasses[i]))
        AvgTable = pd.DataFrame(np.concatenate([[QLen, Util, RespT, ResidT, ArvR, Tput]]).T, columns=cols)
        tokeep = ~(AvgTable <= 0.0).all(axis=1)
        AvgTable.insert(0, "JobClass", classnames)
        AvgTable.insert(0, "Node", nodenames)
        AvgTable = AvgTable.loc[tokeep]
        if not self._table_silent:
            print(AvgTable)

        return AvgTable

    avg_node_table = getAvgNodeTable
    getNodeAvgT = getAvgNodeTable
    nodeAvgT = getAvgNodeTable
    aNT = getAvgNodeTable

    def getAvgChainTable(self):
        """
        Get average performance metrics per chain (customer class).
        
        Returns:
            pandas.DataFrame: Performance metrics table with columns:
                - QLen: Average queue length per chain
                - Util: Utilization per chain
                - RespT: Response time per chain
                - ResidT: Residence time per chain
                - ArvR: Arrival rate per chain
                - Tput: Throughput per chain
        """
        table = self.obj.getAvgChainTable()

        QLen = np.array(list(table.getQLen()))
        Util = np.array(list(table.getUtil()))
        RespT = np.array(list(table.getRespT()))
        ResidT = np.array(list(table.getResidT()))
        ArvR = np.array(list(table.getArvR()))
        Tput = np.array(list(table.getTput()))

        cols = ['QLen', 'Util', 'RespT', 'ResidT', 'ArvR', 'Tput']
        stations = list(table.getStationNames())
        statnames = []
        for i in range(len(stations)):
            statnames.append(str(stations[i]))
        jobchains = list(table.getChainNames())
        chainnames = []
        for i in range(len(jobchains)):
            chainnames.append(str(jobchains[i]))
        AvgChainTable = pd.DataFrame(np.concatenate([[QLen, Util, RespT, ResidT, ArvR, Tput]]).T, columns=cols)
        tokeep = ~(AvgChainTable <= 0.0).all(axis=1)
        AvgChainTable.insert(0, "Chain", chainnames)
        AvgChainTable.insert(0, "Station", statnames)
        AvgChainTable = AvgChainTable.loc[tokeep]
        if not self._table_silent:
            print(AvgChainTable)

        return AvgChainTable

    avg_chain_table = getAvgChainTable
    getChainAvgT = getAvgChainTable
    chainAvgT = getAvgChainTable

    def getAvgNodeChainTable(self):
        """
        Get average performance metrics per node and chain.

        Returns:
            pandas.DataFrame: Performance table with node and chain breakdown containing:
                - QLen: Average queue length per node-chain
                - Util: Utilization per node-chain
                - RespT: Response time per node-chain
                - ResidT: Residence time per node-chain  
                - ArvR: Arrival rate per node-chain
                - Tput: Throughput per node-chain
        """
        table = self.obj.getAvgNodeChainTable()

        QLen = np.array(list(table.getQLen()))
        Util = np.array(list(table.getUtil()))
        RespT = np.array(list(table.getRespT()))
        ResidT = np.array(list(table.getResidT()))
        ArvR = np.array(list(table.getArvR()))
        Tput = np.array(list(table.getTput()))

        cols = ['QLen', 'Util', 'RespT', 'ResidT', 'ArvR', 'Tput']
        nodes = list(table.getNodeNames())
        nodenames = []
        for i in range(len(nodes)):
            nodenames.append(str(nodes[i]))
        jobchains = list(table.getChainNames())
        chainnames = []
        for i in range(len(jobchains)):
            chainnames.append(str(jobchains[i]))
        AvgChainTable = pd.DataFrame(np.concatenate([[QLen, Util, RespT, ResidT, ArvR, Tput]]).T, columns=cols)
        tokeep = ~(AvgChainTable <= 0.0).all(axis=1)
        AvgChainTable.insert(0, "Chain", chainnames)
        AvgChainTable.insert(0, "Node", nodenames)
        AvgChainTable = AvgChainTable.loc[tokeep]
        if not self._table_silent:
            print(AvgChainTable)

        return AvgChainTable

    avg_node_chain_table = getAvgNodeChainTable
    getNodeChainAvgT = getAvgNodeChainTable
    nodeChainAvgT = getAvgNodeChainTable

    def getAvgTable(self):
        """
        Get comprehensive average performance metrics table.
        
        This method provides a detailed view of all performance metrics
        broken down by node and customer class combinations.
        
        Returns:
            pandas.DataFrame: Comprehensive performance metrics table with columns:
                - Node: Node/station name
                - JobClass: Customer class name
                - QLen: Average queue length
                - Util: Utilization
                - RespT: Response time
                - ResidT: Residence time
                - ArvR: Arrival rate
                - Tput: Throughput
        """
        table = self.obj.getAvgTable()


        QLen = np.array(list(table.getQLen()))
        Util = np.array(list(table.getUtil()))
        RespT = np.array(list(table.getRespT()))
        ResidT = np.array(list(table.getResidT()))
        ArvR = np.array(list(table.getArvR()))
        Tput = np.array(list(table.getTput()))

        cols = ['QLen', 'Util', 'RespT', 'ResidT', 'ArvR', 'Tput']

        stations = list(table.getStationNames())
        statnames = []
        for i in range(len(stations)):
            statnames.append(str(stations[i]))
        jobclasses = list(table.getClassNames())
        classnames = []
        for i in range(len(jobclasses)):
            classnames.append(str(jobclasses[i]))
        AvgTable = pd.DataFrame(np.concatenate([[QLen, Util, RespT, ResidT, ArvR, Tput]]).T, columns=cols)
        tokeep = ~(AvgTable <= 0.0).all(axis=1)
        AvgTable.insert(0, "JobClass", classnames)
        AvgTable.insert(0, "Station", statnames)
        AvgTable = AvgTable.loc[tokeep]
        if not self._table_silent:
            print(AvgTable)

        return AvgTable

    avg_table = getAvgTable
    getAvgT = getAvgTable
    avgT = getAvgTable
    aT = getAvgTable

    def getAvgSysTable(self):
        """
        Get system-wide performance metrics table.
        
        Returns a comprehensive table with system-level performance metrics
        including response times and throughputs for all job classes/chains.
        
        Returns:
            pandas.DataFrame: Table with columns ['Chain', 'JobClasses', 'SysRespT', 'SysTput']
                containing system response times and throughputs by chain and job class
        """
        table = self.obj.getAvgSysTable()

        SysRespT = np.array(list(table.getSysRespT()))
        SysTput = np.array(list(table.getSysTput()))

        cols = ['SysRespT', 'SysTput']
        jobchains = list(table.getChainNames())
        chains = []
        for i in range(len(jobchains)):
            chains.append(str(jobchains[i]))
        jobinchains = list(table.getInChainNames())
        inchains = []
        for i in range(len(jobinchains)):
            inchains.append(str(jobinchains[i]))
        AvgSysTable = pd.DataFrame(np.concatenate([[SysRespT, SysTput]]).T, columns=cols)
        tokeep = ~(AvgSysTable <= 0.0).all(axis=1)
        AvgSysTable.insert(0, "JobClasses", inchains)
        AvgSysTable.insert(0, "Chain", chains)
        AvgSysTable = AvgSysTable.loc[tokeep]
        if not self._table_silent:
            print(AvgSysTable)
        return AvgSysTable

    avg_sys_table = getAvgSysTable
    getSysAvgT = getAvgSysTable
    sysAvgT = getAvgSysTable

    def getDeadlineTable(self):
        """
        Get deadline-related metrics table (response time and tardiness).

        Returns:
            pandas.DataFrame: Deadline metrics table with columns:
                - Station: Node/station name
                - JobClass: Customer class name
                - RespT: Mean response time per station-class
                - Trdn: Mean tardiness per station-class
                - SysTrdn: Mean system tardiness per class
            Returns None if tardiness data is not available.
        """
        table = self.obj.getDeadlineTable()

        if table is None:
            return None

        RespT = np.array(list(table.getRespT()))
        Trdn = np.array(list(table.getTrdn()))
        SysTrdn = np.array(list(table.getSysTrdn()))

        if len(Trdn) == 0 and len(SysTrdn) == 0:
            return None

        cols = ['RespT', 'Trdn', 'SysTrdn']

        stations = list(table.getStationNames())
        statnames = []
        for i in range(len(stations)):
            statnames.append(str(stations[i]))
        jobclasses = list(table.getClassNames())
        classnames = []
        for i in range(len(jobclasses)):
            classnames.append(str(jobclasses[i]))

        DeadlineTable = pd.DataFrame(np.concatenate([[RespT, Trdn, SysTrdn]]).T, columns=cols)
        tokeep = ~(DeadlineTable <= 0.0).all(axis=1)
        DeadlineTable.insert(0, "JobClass", classnames)
        DeadlineTable.insert(0, "Station", statnames)
        DeadlineTable = DeadlineTable.loc[tokeep]
        if not self._table_silent:
            print(DeadlineTable)

        return DeadlineTable

    deadline_table = getDeadlineTable

    def getAvgTput(self):
        """
        Get average throughputs for all stations and classes.

        Returns:
            numpy.ndarray: Matrix of throughput values with shape (stations, classes)
        """
        Tput = jlineMatrixToArray(self.obj.getAvgTput())
        if not self._table_silent:
                print(Tput)
        return Tput

    avg_tput = getAvgTput

    def getAvgResidT(self):
        """
        Get average residence times for all stations and classes.

        Returns:
            numpy.ndarray: Matrix of residence time values with shape (stations, classes)
        """
        ResidT = jlineMatrixToArray(self.obj.getAvgResidT())
        if not self._table_silent:
            print(ResidT)
        return ResidT

    avg_residt = getAvgResidT

    def getAvgArvR(self):
        """
        Get average arrival rates for all stations and classes.

        Returns:
            numpy.ndarray: Matrix of arrival rate values with shape (stations, classes)
        """
        ArvR = jlineMatrixToArray(self.obj.getAvgArvR())
        if not self._table_silent:
            print(ArvR)
        return ArvR

    avg_arv_r = getAvgArvR

    def getAvgUtil(self):
        """
        Get average utilizations for all stations and classes.

        Returns:
            numpy.ndarray: Matrix of utilization values with shape (stations, classes)
        """
        Util = jlineMatrixToArray(self.obj.getAvgUtil())
        if not self._table_silent:
            print(Util)
        return Util

    avg_util = getAvgUtil

    def getAvgQLen(self):
        """
        Get average queue lengths for all stations and classes.

        Returns:
            numpy.ndarray: Matrix of queue length values with shape (stations, classes)
        """
        QLen = jlineMatrixToArray(self.obj.getAvgQLen())
        if not self._table_silent:
            print(QLen)
        return QLen

    avg_q_len = getAvgQLen

    def getAvgRespT(self):
        """
        Get average response times for all stations and classes.

        Returns:
            numpy.ndarray: Matrix of response time values with shape (stations, classes)
        """
        RespT = jlineMatrixToArray(self.obj.getAvgRespT())
        if not self._table_silent:
            print(RespT)
        return RespT

    avg_respt = getAvgRespT

    def getAvgWaitT(self):
        """
        Get average waiting times for all stations and classes.

        Returns:
            numpy.ndarray: Matrix of waiting time values with shape (stations, classes)
        """
        WaitT = jlineMatrixToArray(self.obj.getAvgWaitT())
        if not self._table_silent:
            print(WaitT)
        return WaitT

    avg_waitt = getAvgWaitT
    avg_wait_t = getAvgWaitT

    def getAvgSysTput(self):
        """
        Get average system throughput for all classes.

        Returns:
            numpy.ndarray: Vector of system throughput values for each class
        """
        SysTput = jlineMatrixToArray(self.obj.getAvgSysTput())
        if not self._table_silent:
            print(SysTput)
        return SysTput

    avg_sys_tput = getAvgSysTput

    def getAvgSysRespT(self):
        """
        Get average system response time for all classes.

        Returns:
            numpy.ndarray: Vector of system response time values for each class
        """
        SysRespT = jlineMatrixToArray(self.obj.getAvgSysRespT())
        if not self._table_silent:
            print(SysRespT)
        return SysRespT

    avg_sys_respt = getAvgSysRespT
    avg_sys_resp_t = getAvgSysRespT

    def getAvg(self):
        """
        Get comprehensive average performance metrics.
        
        Returns a matrix containing average performance metrics computed
        by the solver, including queue lengths, utilizations, and other
        performance measures depending on the solver type.
        
        Returns:
            numpy.ndarray: Matrix of average performance metrics
        """
        result = self.obj.getAvg()
        if hasattr(result, 'toArray2D'):
            avgRet = jlineMatrixToArray(result)
        else:
            avgRet = jlineMatrixToArray(result.QN)
        if not self._table_silent:
            print(avgRet)
        return avgRet

    avg = getAvg

    def getAvgArvRChain(self):
        """
        Get average arrival rates by routing chain.
        
        Returns arrival rates organized by routing chains, essential
        for multi-chain network analysis and chain-specific performance.
        
        Returns:
            numpy.ndarray: Matrix of average arrival rates by chain
        """
        ArvRChain = jlineMatrixToArray(self.obj.getAvgArvRChain())
        if not self._table_silent:
            print(ArvRChain)
        return ArvRChain

    avg_arv_r_chain = getAvgArvRChain

    def getAvgNodeArvRChain(self):
        """
        Get average node arrival rates per chain.

        Returns:
            numpy.ndarray: Node arrival rates by chain
        """
        NodeArvRChain = jlineMatrixToArray(self.obj.getAvgNodeArvRChain())
        if not self._table_silent:
            print(NodeArvRChain)
        return NodeArvRChain

    avg_node_arv_r_chain = getAvgNodeArvRChain

    def getAvgNodeQLenChain(self):
        """
        Get average node queue lengths per chain.

        Returns:
            numpy.ndarray: Node queue lengths by chain
        """
        NodeQLenChain = jlineMatrixToArray(self.obj.getAvgNodeQLenChain())
        if not self._table_silent:
            print(NodeQLenChain)
        return NodeQLenChain

    avg_node_q_len_chain = getAvgNodeQLenChain

    def getAvgNodeResidTChain(self):
        """
        Get average node residence times per chain.

        Returns:
            numpy.ndarray: Node residence times by chain
        """
        NodeResidTChain = jlineMatrixToArray(self.obj.getAvgNodeResidTChain())
        if not self._table_silent:
            print(NodeResidTChain)
        return NodeResidTChain

    avg_node_residt_chain = getAvgNodeResidTChain

    def getAvgNodeRespTChain(self):
        """
        Get average node response times per chain.

        Returns:
            numpy.ndarray: Node response times by chain
        """
        NodeRespTChain = jlineMatrixToArray(self.obj.getAvgNodeRespTChain())
        if not self._table_silent:
            print(NodeRespTChain)
        return NodeRespTChain

    avg_node_respt_chain = getAvgNodeRespTChain

    def getAvgNodeTputChain(self):
        """
        Get average node throughputs per chain.

        Returns:
            numpy.ndarray: Node throughputs by chain
        """
        NodeTputChain = jlineMatrixToArray(self.obj.getAvgNodeTputChain())
        if not self._table_silent:
            print(NodeTputChain)
        return NodeTputChain

    avg_node_tput_chain = getAvgNodeTputChain

    def getAvgNodeUtilChain(self):
        """
        Get average node utilizations per chain.

        Returns:
            numpy.ndarray: Node utilizations by chain
        """
        NodeUtilChain = jlineMatrixToArray(self.obj.getAvgNodeUtilChain())
        if not self._table_silent:
            print(NodeUtilChain)
        return NodeUtilChain

    avg_node_util_chain = getAvgNodeUtilChain

    def getAvgQLenChain(self):
        """
        Get average queue lengths per chain.

        Returns:
            numpy.ndarray: Queue lengths by chain
        """
        QLenChain = jlineMatrixToArray(self.obj.getAvgQLenChain())
        if not self._table_silent:
            print(QLenChain)
        return QLenChain

    avg_q_len_chain = getAvgQLenChain

    def getAvgResidTChain(self):
        """
        Get average residence times per chain.

        Returns:
            numpy.ndarray: Residence times by chain
        """
        ResidTChain = jlineMatrixToArray(self.obj.getAvgResidTChain())
        if not self._table_silent:
            print(ResidTChain)
        return ResidTChain

    avg_residt_chain = getAvgResidTChain

    def getAvgRespTChain(self):
        """
        Get average response times per chain.

        Returns:
            numpy.ndarray: Response times by chain
        """
        RespTChain = jlineMatrixToArray(self.obj.getAvgRespTChain())
        if not self._table_silent:
            print(RespTChain)
        return RespTChain

    avg_respt_chain = getAvgRespTChain

    def getAvgTputChain(self):
        """
        Get average throughputs per chain.

        Returns:
            numpy.ndarray: Throughputs by chain
        """
        TputChain = jlineMatrixToArray(self.obj.getAvgTputChain())
        if not self._table_silent:
            print(TputChain)
        return TputChain

    avg_tput_chain = getAvgTputChain

    def getAvgUtilChain(self):
        """
        Get average utilizations per chain.

        Returns:
            numpy.ndarray: Utilizations by chain
        """
        UtilChain = jlineMatrixToArray(self.obj.getAvgUtilChain())
        if not self._table_silent:
            print(UtilChain)
        return UtilChain

    avg_util_chain = getAvgUtilChain

    def getCdfRespT(self):
        """
        Get cumulative distribution function of response times.

        Returns:
            list: CDF data for response times per station and class
        """
        try:
            table = self.obj.cdfRespT()
            distribC = self.obj.fluidResult.distribC

            num_stations = self.model.get_number_of_stations()
            num_classes = self.model.get_number_of_classes()

            CdfRespT = []
            for i in range(num_stations):
                station_cdfs = []
                for c in range(num_classes):
                    if i < distribC.length and c < distribC[i].length:
                        F = jlineMatrixToArray(distribC[i][c])
                        station_cdfs.append(F)
                    else:
                        station_cdfs.append(None)
                CdfRespT.append(station_cdfs)

            return CdfRespT
        except:
            try:
                num_stations = self.model.get_number_of_stations()
                num_classes = self.model.get_number_of_classes()
                return [[None for _ in range(num_classes)] for _ in range(num_stations)]
            except:
                return [[]]

    cdf_respt = getCdfRespT
    getSjrnT = getCdfRespT
    sjrn_t = getSjrnT

    def getProbAggr(self, node, state=None):
        """
        Probability of a SPECIFIC per-class job distribution at a station.
        Returns P(n1 jobs of class 1, n2 jobs of class 2, ...) for given state.

        Compare with get_prob_marg: returns queue-length distribution for a
        single class, i.e., P(n jobs of class r) for n=0,1,...,N(r).

        Args:
            node: Node object or index
            state: Per-class job counts, e.g., [2,1] = 2 class-1, 1 class-2

        Returns:
            float: Scalar probability in [0,1]
        """
        if hasattr(node, 'obj'):
            if state is not None:
                java_result = self.obj.getProbAggr(node.obj, jlineMatrixFromArray(state))
            else:
                java_result = self.obj.getProbAggr(node.obj)
        else:
            if state is not None:
                java_result = self.obj.getProbAggr(int(node), jlineMatrixFromArray(state))
            else:
                java_result = self.obj.getProbAggr(int(node))

        if hasattr(java_result, 'getScalarProbability'):
            return java_result.getScalarProbability()
        elif hasattr(java_result, 'probability'):
            prob_matrix = java_result.probability
            if prob_matrix.getNumRows() == 1 and prob_matrix.getNumCols() == 1:
                return prob_matrix.get(0, 0)
            else:
                return jlineMatrixToArray(prob_matrix)
        else:
            return float(java_result)

    prob_aggr = getProbAggr

    def getProb(self, node, state=None):
        """
        Get state probability for a node.

        Args:
            node: Node object or index
            state: State specification (optional)

        Returns:
            float or numpy.ndarray: State probability
        """
        if hasattr(node, 'obj'):
            if state is not None:
                java_result = self.obj.prob(node.obj, jlineMatrixFromArray(state))
            else:
                java_result = self.obj.prob(node.obj)
        else:
            if state is not None:
                java_result = self.obj.prob(int(node), jlineMatrixFromArray(state))
            else:
                java_result = self.obj.prob(int(node))

        if hasattr(java_result, 'getScalarProbability'):
            return java_result.getScalarProbability()
        elif hasattr(java_result, 'probability'):
            prob_matrix = java_result.probability
            if prob_matrix.getNumRows() == 1 and prob_matrix.getNumCols() == 1:
                return prob_matrix.get(0, 0)
            else:
                return jlineMatrixToArray(prob_matrix)
        else:
            return float(java_result)

    prob = getProb

    def isSolved(self):
        """
        Check if the solver has results available.

        Returns:
            bool: True if solver has computed results, False otherwise
        """
        return hasattr(self, 'obj') and self.obj.hasResults()

    def getSolverType(self):
        """
        Get the type/name of this solver.

        Returns:
            str: Solver type name
        """
        if hasattr(self, 'obj'):
            return self.obj.getName()
        return self.__class__.__name__.replace('Solver', '')

    def reset(self):
        """
        Reset solver state and clear results.
        """
        if hasattr(self, 'obj') and hasattr(self.obj, 'reset'):
            self.obj.reset()

    def javaObj(self):
        """
        Get the underlying Java solver object.

        Returns:
            Java object or None: The Java solver object if available
        """
        return self.obj if hasattr(self, 'obj') else None

    def getUtil(self):
        """
        Get utilizations (alias for getAvgUtil).

        Returns:
            numpy.ndarray: Station utilizations
        """
        return self.getAvgUtil()

    def print(self):
        """
        Print solver information and results.
        """
        if hasattr(self, 'obj') and self.obj is not None:
            self.obj.print_()
        else:
            raise RuntimeError("No Java solver object available")

    def hasResults(self):
        """
        Check if solver has computed results.

        Returns:
            bool: True if results are available, False otherwise
        """
        if hasattr(self, 'obj') and self.obj is not None:
            try:
                return self.obj.hasResults()
            except:
                return False
        else:
            return False

    @classmethod
    def supportsModel(cls, model):
        """
        Check if this solver supports the given model type.

        Args:
            model: Model to check

        Returns:
            bool: True if model is supported, False otherwise
        """
        return True

    def avgNodeTable(self):
        """Get average performance metrics table organized by network nodes.
        
        Returns:
            DataFrame: Table with performance metrics (QLen, Util, RespT, etc.) for each node.
        """
        return self.avg_node_table()

    def avgChainTable(self):
        """Get average performance metrics table organized by routing chains.
        
        Returns:
            DataFrame: Table with performance metrics aggregated by routing chain.
        """
        return self.avg_chain_table()

    def avgNodeChainTable(self):
        """Get average performance metrics table organized by nodes and chains.
        
        Returns:
            DataFrame: Table with performance metrics for each node-chain combination.
        """
        return self.getAvgNodeChainTable()

    def avgTable(self):
        """Get the main average performance metrics table.
        
        Returns:
            DataFrame: Complete table with all average performance metrics.
        """
        return self.avg_table()

    def avgSysTable(self):
        """Get system-wide average performance metrics table.

        Returns:
            DataFrame: Table with aggregated system-level performance metrics.
        """
        return self.avg_sys_table()

    # Table -> T aliases
    def avgT(self):
        """Short alias for avgTable."""
        return self.avgTable()

    def avgSysT(self):
        """Short alias for avgSysTable."""
        return self.avgSysTable()

    def avgNodeT(self):
        """Short alias for avgNodeTable."""
        return self.avgNodeTable()

    def avgChainT(self):
        """Short alias for avgChainTable."""
        return self.avgChainTable()

    def avgNodeChainT(self):
        """Short alias for avgNodeChainTable."""
        return self.avgNodeChainTable()

    # Snake_case Table -> T aliases
    avg_t = avgT
    avg_sys_t = avgSysT
    avg_node_t = avgNodeT
    avg_chain_t = avgChainT
    avg_node_chain_t = avgNodeChainT

    # Short 4-letter aliases (consistent with MATLAB/JAR)
    def aT(self):
        """Short alias for getAvgTable (MATLAB-compatible)."""
        return self.avgTable()

    def aNT(self):
        """Short alias for getAvgNodeTable (MATLAB-compatible)."""
        return self.avgNodeTable()

    def aCT(self):
        """Short alias for getAvgChainTable (MATLAB-compatible)."""
        return self.avgChainTable()

    def aST(self):
        """Short alias for getAvgSysTable (MATLAB-compatible)."""
        return self.avgSysTable()

    def aNCT(self):
        """Short alias for getAvgNodeChainTable (MATLAB-compatible)."""
        return self.avgNodeChainTable()

    # Snake_case short aliases
    a_t = aT
    a_nt = aNT
    a_ct = aCT
    a_st = aST
    a_nct = aNCT

    def avgTput(self):
        """Get average throughput for all stations and job classes.
        
        Returns:
            ndarray: Matrix of average throughput values [stations x classes].
        """
        return self.getAvgTput()

    def avgResidT(self):
        """Get average residence time for all stations and job classes.
        
        Returns:
            ndarray: Matrix of average residence times [stations x classes].
        """
        return self.getAvgResidT()

    def avgArvR(self):
        """Get average arrival rates for all stations and job classes.
        
        Returns:
            ndarray: Matrix of average arrival rates [stations x classes].
        """
        return self.getAvgArvR()

    def avgUtil(self):
        """Get average utilization for all stations and job classes.
        
        Returns:
            ndarray: Matrix of average utilization values [stations x classes].
        """
        return self.getAvgUtil()

    def avgQLen(self):
        """Get average queue length for all stations and job classes.
        
        Returns:
            ndarray: Matrix of average queue lengths [stations x classes].
        """
        return self.getAvgQLen()

    def avgRespT(self):
        """Get average response time for all stations and job classes.
        
        Returns:
            ndarray: Matrix of average response times [stations x classes].
        """
        return self.getAvgRespT()

    def avgWaitT(self):
        """Get average waiting time for all stations and job classes.
        
        Returns:
            ndarray: Matrix of average waiting times [stations x classes].
        """
        return self.getAvgWaitT()

    def avgSysTput(self):
        """Get average system throughput for all job classes.
        
        Returns:
            ndarray: Vector of average system throughput values per class.
        """
        return self.getAvgSysTput()

    def avgSysRespT(self):
        """Get average system response time for all job classes.
        
        Returns:
            ndarray: Vector of average system response times per class.
        """
        return self.getAvgSysRespT()

    def avgArvRChain(self):
        """Get average arrival rates organized by routing chains.
        
        Returns:
            ndarray: Matrix of average arrival rates [stations x chains].
        """
        return self.getAvgArvRChain()

    def avgNodeArvRChain(self):
        """Get average node arrival rates organized by routing chains.
        
        Returns:
            ndarray: Matrix of average node arrival rates [nodes x chains].
        """
        return self.getAvgNodeArvRChain()

    def avgNodeQLenChain(self):
        """Get average node queue lengths organized by routing chains.
        
        Returns:
            ndarray: Matrix of average node queue lengths [nodes x chains].
        """
        return self.getAvgNodeQLenChain()

    def avgNodeResidTChain(self):
        """Get average node residence times organized by routing chains.
        
        Returns:
            ndarray: Matrix of average node residence times [nodes x chains].
        """
        return self.getAvgNodeResidTChain()

    def avgNodeRespTChain(self):
        """Get average node response times organized by routing chains.
        
        Returns:
            ndarray: Matrix of average node response times [nodes x chains].
        """
        return self.getAvgNodeRespTChain()

    def avgNodeTputChain(self):
        """Get average node throughput organized by routing chains.
        
        Returns:
            ndarray: Matrix of average node throughput [nodes x chains].
        """
        return self.getAvgNodeTputChain()

    def avgNodeUtilChain(self):
        """Get average node utilization organized by routing chains.
        
        Returns:
            ndarray: Matrix of average node utilization [nodes x chains].
        """
        return self.getAvgNodeUtilChain()

    def avgQLenChain(self):
        """Get average queue lengths organized by routing chains.
        
        Returns:
            ndarray: Matrix of average queue lengths [stations x chains].
        """
        return self.getAvgQLenChain()

    def avgResidTChain(self):
        """Get average residence times organized by routing chains.
        
        Returns:
            ndarray: Matrix of average residence times [stations x chains].
        """
        return self.getAvgResidTChain()

    def avgRespTChain(self):
        """Get average response times organized by routing chains.
        
        Returns:
            ndarray: Matrix of average response times [stations x chains].
        """
        return self.getAvgRespTChain()

    def avgTputChain(self):
        """Get average throughput organized by routing chains.
        
        Returns:
            ndarray: Matrix of average throughput [stations x chains].
        """
        return self.getAvgTputChain()

    def avgUtilChain(self):
        """Get average utilization organized by routing chains.
        
        Returns:
            ndarray: Matrix of average utilization [stations x chains].
        """
        return self.getAvgUtilChain()

    def util(self):
        """Get station utilizations.
        
        Returns:
            numpy.ndarray: Station utilization values.
        """
        return self.getUtil()

    def tput(self):
        """Get average throughputs for all stations and classes.
        
        Returns:
            numpy.ndarray: Matrix of throughput values with shape (stations, classes).
        """
        return self.getAvgTput()

    def respT(self):
        """Get average response times for all stations and classes.
        
        Returns:
            numpy.ndarray: Matrix of response time values with shape (stations, classes).
        """
        return self.getAvgRespT()

    def qLen(self):
        """Get average queue lengths for all stations and classes.
        
        Returns:
            numpy.ndarray: Matrix of queue length values with shape (stations, classes).
        """
        return self.getAvgQLen()

    def residT(self):
        """Get average residence times for all stations and classes.
        
        Returns:
            numpy.ndarray: Matrix of residence time values with shape (stations, classes).
        """
        return self.getAvgResidT()

    def waitT(self):
        """Get average waiting times for all stations and classes.
        
        Returns:
            numpy.ndarray: Matrix of waiting time values with shape (stations, classes).
        """
        return self.getAvgWaitT()

    def get_util(self):
        """Get station utilizations.
        
        Returns:
            numpy.ndarray: Station utilization values.
        """
        return self.getUtil()

    def get_tput(self):
        """Get average throughputs for all stations and classes.
        
        Returns:
            numpy.ndarray: Matrix of throughput values with shape (stations, classes).
        """
        return self.getAvgTput()

    def get_respt(self):
        """Get average response times for all stations and classes.
        
        Returns:
            numpy.ndarray: Matrix of response time values with shape (stations, classes).
        """
        return self.getAvgRespT()

    def get_q_len(self):
        """Get average queue lengths for all stations and classes.
        
        Returns:
            numpy.ndarray: Matrix of queue length values with shape (stations, classes).
        """
        return self.getAvgQLen()

    def get_residt(self):
        """Get average residence times for all stations and classes.
        
        Returns:
            numpy.ndarray: Matrix of residence time values with shape (stations, classes).
        """
        return self.getAvgResidT()

    def get_wait_t(self):
        """Get average waiting times for all stations and classes.
        
        Returns:
            numpy.ndarray: Matrix of waiting time values with shape (stations, classes).
        """
        return self.getAvgWaitT()

    def avg_node_table(self):
        """Get average performance metrics table aggregated by node.

        Returns:
            Performance metrics table organized by network nodes.
        """
        return self.getAvgNodeTable()

    def avg_chain_table(self):
        """Get average performance metrics table aggregated by chain.

        Returns:
            Performance metrics table organized by job chains.
        """
        return self.getAvgChainTable()

    def avg_node_chain_table(self):
        """Get average performance metrics table by node and chain.
        
        Returns:
            Performance metrics table organized by both nodes and chains.
        """
        return self.getAvgNodeChainTable()

    def avg_sys_table(self):
        """Get average system-level performance metrics table.

        Returns:
            Performance metrics table aggregated at the system level.
        """
        return self.getAvgSysTable()

    def avg_arv_r(self):
        """Get average arrival rates for all stations and classes.
        
        Returns:
            numpy.ndarray: Matrix of arrival rate values with shape (stations, classes).
        """
        return self.getAvgArvR()

    def avg_sys_tput(self):
        """Get average system throughput.
        
        Returns:
            float: System-level throughput value.
        """
        return self.getAvgSysTput()

    def avg_sys_respt(self):
        """Get average system response time.

        Returns:
            float: System-level response time value.
        """
        return self.getAvgSysRespT()

    def getAvgSys(self):
        """Get average system performance metrics.

        Returns:
            list: [response_time, throughput] - System-level performance values.
        """
        # Call the underlying Java method which computes system metrics
        self.obj.getAvgSys()
        # Return both response time and throughput
        return [self.getAvgSysRespT(), self.getAvgSysTput()]

    def avg_arv_r_chain(self):
        """Get average arrival rates per chain.
        
        Returns:
            numpy.ndarray: Arrival rates by chain.
        """
        return self.getAvgArvRChain()

    def avg_node_arv_r_chain(self):
        """Get average arrival rates per node and chain.
        
        Returns:
            numpy.ndarray: Matrix of arrival rates with shape (nodes, chains).
        """
        return self.getAvgNodeArvRChain()

    def avg_node_q_len_chain(self):
        """Get average queue lengths per node and chain.
        
        Returns:
            numpy.ndarray: Matrix of queue lengths with shape (nodes, chains).
        """
        return self.getAvgNodeQLenChain()

    def avg_node_residt_chain(self):
        """Get average residence times per node and chain.
        
        Returns:
            numpy.ndarray: Matrix of residence times with shape (nodes, chains).
        """
        return self.getAvgNodeResidTChain()

    def avg_node_respt_chain(self):
        """Get average response times per node and chain.
        
        Returns:
            numpy.ndarray: Matrix of response times with shape (nodes, chains).
        """
        return self.getAvgNodeRespTChain()

    def avg_node_tput_chain(self):
        """Get average throughputs per node and chain.
        
        Returns:
            numpy.ndarray: Matrix of throughputs with shape (nodes, chains).
        """
        return self.getAvgNodeTputChain()

    def avg_node_util_chain(self):
        """Get average utilizations per node and chain.
        
        Returns:
            numpy.ndarray: Matrix of utilizations with shape (nodes, chains).
        """
        return self.getAvgNodeUtilChain()

    def avg_q_len_chain(self):
        """Get average queue lengths per chain.
        
        Returns:
            numpy.ndarray: Queue lengths by chain.
        """
        return self.getAvgQLenChain()

    def avg_residt_chain(self):
        """Get average residence times per chain.
        
        Returns:
            numpy.ndarray: Residence times by chain.
        """
        return self.getAvgResidTChain()

    def avg_respt(self):
        """Get average response times per chain.

        Returns:
            numpy.ndarray: Response times by chain.
        """
        return self.getAvgRespT()

    def avg_respt_chain(self):
        """Get average response times per chain.
        
        Returns:
            numpy.ndarray: Response times by chain.
        """
        return self.getAvgRespTChain()

    def avg_tput_chain(self):
        """Get average throughputs per chain.
        
        Returns:
            numpy.ndarray: Throughputs by chain.
        """
        return self.getAvgTputChain()

    def avg_util_chain(self):
        """Get average utilizations per chain.
        
        Returns:
            numpy.ndarray: Utilizations by chain.
        """
        return self.getAvgUtilChain()

    def tran_prob(self, node):
        """Get transient state probabilities.
        
        Returns:
            numpy.ndarray: Transient probability values.
        """
        return self.getTranProb(node)

    def tran_prob_aggr(self, node):
        """Get aggregated transient state probabilities.
        
        Returns:
            numpy.ndarray: Aggregated transient probability values.
        """
        return self.getTranProbAggr(node)

    def tran_prob_sys(self):
        """Get system-level transient state probabilities.
        
        Returns:
            numpy.ndarray: System transient probability values.
        """
        return self.getTranProbSys()

    def tran_prob_sys_aggr(self):
        """Get aggregated system-level transient state probabilities.

        Returns:
            numpy.ndarray: Aggregated system transient probability values.
        """
        return self.getTranProbSysAggr()

    def tran_avg(self):
        """Get transient average performance metrics.
        
        Returns:
            Transient performance metrics.
        """
        return self.getTranAvg()

    def tran_cdf_respt(self, R=None):
        """Get transient cumulative distribution function of response times.
        
        Returns:
            numpy.ndarray: CDF values for response times.
        """
        return self.getTranCdfRespT(R)

    def tran_cdf_passt(self, R=None):
        """Get transient cumulative distribution function of passage times.
        
        Returns:
            numpy.ndarray: CDF values for passage times.
        """
        return self.getTranCdfPassT(R)

    def distrib_respt(self):
        """Get response time distributions.
        
        Returns:
            Response time distribution data.
        """
        return self.getDistribRespT()

    def distrib_respt_chain(self):
        """Get response time distributions per chain.
        
        Returns:
            Response time distribution data by chain.
        """
        return self.getDistribRespTChain()

    def distrib_respt_node(self):
        """Get response time distributions per node.
        
        Returns:
            Response time distribution data by node.
        """
        return self.getDistribRespTNode()

    def distrib_respt_node_chain(self):
        """Get response time distributions per node and chain.
        
        Returns:
            Response time distribution data by node and chain.
        """
        return self.getDistribRespTNodeChain()

    def tranProb(self, node):
        """Get transient state probabilities for a specific node.
        
        Args:
            node: The node to get probabilities for.
            
        Returns:
            ProbabilityResult: Transient state probability distribution.
        """
        return self.getTranProb(node)

    def tranProbAggr(self, node):
        """Get aggregated transient state probabilities for a specific node.
        
        Args:
            node: The node to get aggregated probabilities for.
            
        Returns:
            ProbabilityResult: Aggregated transient state probability distribution.
        """
        return self.getTranProbAggr(node)

    def tranProbSys(self):
        """Get system-wide transient state probabilities.
        
        Returns:
            ProbabilityResult: System-level transient state probability distribution.
        """
        return self.getTranProbSys()

    def tranProbSysAggr(self):
        """Get aggregated system-wide transient state probabilities.
        
        Returns:
            ProbabilityResult: Aggregated system-level transient probability distribution.
        """
        return self.tran_prob_sys_aggr()

    def tranAvg(self):
        """Get transient average performance metrics.
        
        Returns:
            dict: Dictionary containing transient performance metrics over time.
        """
        return self.getTranAvg()

    def tranCdfRespT(self, R=None):
        """Get transient cumulative distribution function of response times.
        
        Args:
            R: Optional response time handles.
            
        Returns:
            DistributionResult: Transient response time CDF.
        """
        return self.getTranCdfRespT(R)

    def tranCdfPassT(self, R=None):
        """Get transient cumulative distribution function of passage times.
        
        Args:
            R: Optional passage time handles.
            
        Returns:
            DistributionResult: Transient passage time CDF.
        """
        return self.getTranCdfPassT(R)

    def prob(self, node):
        """Get transient state probabilities for a specific node (short alias).
        
        Args:
            node: The node to get probabilities for.
            
        Returns:
            ProbabilityResult: Transient state probability distribution.
        """
        return self.getTranProb(node)

    def probAggr(self, node):
        """Get aggregated transient state probabilities for a specific node (short alias).
        
        Args:
            node: The node to get aggregated probabilities for.
            
        Returns:
            ProbabilityResult: Aggregated transient state probability distribution.
        """
        return self.getTranProbAggr(node)

    def probSys(self):
        """Get system-wide transient state probabilities (short alias).
        
        Returns:
            ProbabilityResult: System-level transient state probability distribution.
        """
        return self.getTranProbSys()

    def probSysAggr(self):
        """Get aggregated system-wide transient state probabilities (short alias).
        
        Returns:
            ProbabilityResult: Aggregated system-level transient probability distribution.
        """
        return self.tran_prob_sys_aggr()

    def name(self):
        """Get the name of this solver.
        
        Returns:
            str: Name of the solver.
        """
        return self.get_name()

    def numberOfModels(self):
        """Get the number of models in this solver (duplicate alias for compatibility)."""
        return self.getNumberOfModels()

    get_avg_node_table = getAvgNodeTable
    get_avg_chain_table = getAvgChainTable
    get_avg_node_chain_table = getAvgNodeChainTable
    get_avg_table = getAvgTable
    get_avg_sys_table = getAvgSysTable
    get_avg_tput = getAvgTput
    get_avg_residt = getAvgResidT
    get_avg_arv_r = getAvgArvR
    get_avg_util = getAvgUtil
    get_avg_q_len = getAvgQLen
    get_avg_respt = getAvgRespT
    get_avg_wait_t = getAvgWaitT
    get_avg_sys_tput = getAvgSysTput
    get_avg_sys_respt = getAvgSysRespT
    get_avg_sys = getAvgSys
    avg_sys = getAvgSys
    get_avg = getAvg
    get_avg_arv_r_chain = getAvgArvRChain
    get_avg_node_arv_r_chain = getAvgNodeArvRChain
    get_avg_node_q_len_chain = getAvgNodeQLenChain
    get_avg_node_residt_chain = getAvgNodeResidTChain
    get_avg_node_respt_chain = getAvgNodeRespTChain
    get_avg_node_tput_chain = getAvgNodeTputChain
    get_avg_node_util_chain = getAvgNodeUtilChain
    get_avg_q_len_chain = getAvgQLenChain
    get_avg_residt_chain = getAvgResidTChain
    get_avg_respt_chain = getAvgRespTChain
    get_avg_tput_chain = getAvgTputChain
    get_avg_util_chain = getAvgUtilChain
    get_cdf_respt = getCdfRespT
    get_prob_aggr = getProbAggr
    get_prob = getProb
    get_solver_type = getSolverType
    get_util = getUtil

class SolverCTMC(NetworkSolver):
    """
    Continuous-Time Markov Chain (CTMC) solver.

    SolverCTMC analyzes queueing networks by modeling them as continuous-time
    Markov chains and computing steady-state and transient solutions through
    matrix-based methods. This solver provides exact results for networks with
    exponential and phase-type distributions.

    The solver supports:
    - Open, closed, and mixed networks
    - Exponential and phase-type service distributions
    - Markovian arrival processes (MAP, MMPP)
    - State-dependent routing
    - Load-dependent service rates
    - Small to medium state spaces (typically < 10^6 states)

    Args:
        model: Network model to solve
        lang: Language for solver execution ('python' for native, 'java' for wrapper).
              Default is 'python' which uses native Python implementation.

    The CTMC solver constructs the infinitesimal generator matrix and solves
    for the steady-state probability vector. It also supports state probability
    queries and passage time distributions.

    Note:
        State space grows exponentially with network size. For large networks,
        consider using approximate solvers like SolverFluid or SolverMVA.
        Best suited for models with small to medium state spaces where exact
        results are required.
    """

    def __init__(self, *args, **kwargs):
        self._model = args[0]

        _initialize_jar_globals()
        if len(args) > 1 and hasattr(args[1], 'obj'):
            options = args[1]
            super().__init__(options, *args[2:], **kwargs)
        else:
            options = SolverOptions(jpype.JPackage('jline').lang.constant.SolverType.CTMC)
            super().__init__(options, *args[1:], **kwargs)
        model = args[0]
        java_network = _get_java_network(model)
        self.obj = jpype.JPackage('jline').solvers.ctmc.SolverCTMC(java_network, self.solveropt.obj)

    def runAnalyzer(self):
        """Run the CTMC analysis."""
        self.obj.runAnalyzer()
        return self

    def getAvgTable(self):
        """Get comprehensive average performance metrics table."""
        return super().getAvgTable()

    avg_table = getAvgTable
    getAvgT = getAvgTable
    avgT = getAvgTable
    aT = getAvgTable

    def getStateSpace(self):
        """
        Get state space representation for CTMC analysis.
        
        Returns the complete state space representation of the continuous-time
        Markov chain model, including global and local state spaces.
        
        Returns:
            tuple: (state_space, local_state_space) where state_space is the
                global state matrix and local_state_space contains station-specific states
        """
        StateSpace = self.obj.getStateSpace(self.solveropt.obj)
        return jlineMatrixToArray(StateSpace.stateSpace), jlineMatrixCellToArray(StateSpace.localStateSpace)

    def getGenerator(self):
        """
        Get infinitesimal generator matrix for CTMC.

        Returns the infinitesimal generator matrix that defines the
        transition rates between states in the continuous-time Markov chain.

        Returns:
            tuple: (generator_matrix, event_filters) where generator_matrix is the
                infinitesimal generator and event_filters contains event-specific filters
        """
        generatorResult = self.obj.getGenerator()
        return jlineMatrixToArray(generatorResult.infGen), jlineMapMatrixToArray(generatorResult.eventFilt.toMap())

    generator = getGenerator
    get_generator = getGenerator

    @staticmethod
    def printInfGen(infGen, stateSpace):
        jpype.JPackage('jline').solvers.ctmc.SolverCTMC.printInfGen(jlineMatrixFromArray(infGen), jlineMatrixFromArray(stateSpace))

    print_inf_gen = printInfGen

    def getTranProb(self, node):
        try:
            java_result = self.obj.getTranProb(node.obj if hasattr(node, 'obj') else node)
            return ProbabilityResult(java_result)
        except Exception as e:
            if not self._verbose_silent:
                print(f"CTMC getTranProb failed: {e}")
            return None

    tran_prob = getTranProb

    def getTranProbAggr(self, node):
        try:
            java_result = self.obj.getTranProbAggr(node.obj if hasattr(node, 'obj') else node)
            return ProbabilityResult(java_result)
        except Exception as e:
            if not self._verbose_silent:
                print(f"CTMC getTranProbAggr failed: {e}")
            return None

    tran_prob_aggr = getTranProbAggr

    def getTranProbSys(self):
        try:
            java_result = self.obj.getTranProbSys()
            return ProbabilityResult(java_result)
        except Exception as e:
            if not self._verbose_silent:
                print(f"CTMC getTranProbSys failed: {e}")
            return None

    tran_prob_sys = getTranProbSys

    def getTranProbSysAggr(self):
        try:
            java_result = self.obj.tranProbSysAggr()
            return ProbabilityResult(java_result)
        except Exception as e:
            if not self._verbose_silent:
                print(f"CTMC getTranProbSysAggr failed: {e}")
            return None

    tran_prob_sys_aggr = getTranProbSysAggr

    def getProbSysAggr(self):
        """
        Get aggregated system-wide state probabilities.
        
        Returns probability distribution over aggregated system states
        for the CTMC model.
        
        Returns:
            ProbabilityResult: Aggregated system state probabilities
        """
        try:
            java_result = self.obj.probSysAggr()
            return ProbabilityResult(java_result)
        except Exception as e:
            if not self._verbose_silent:
                print(f"CTMC getProbSysAggr failed: {e}")
            return None

    prob_sys_aggr = getProbSysAggr

    def sample(self, node, numSamples):
        """
        Sample performance metrics at a specific node.
        
        Generates random samples of performance metrics (queue lengths,
        response times, etc.) at the specified node using CTMC simulation.
        
        Args:
            node: Network node to sample
            numSamples: Number of samples to generate
            
        Returns:
            SampleResult: Object containing sampled performance metrics
        """
        java_result = self.obj.sample(node.obj, numSamples)
        return SampleResult(java_result) if java_result is not None else None

    def sampleAggr(self, node, numSamples):
        """
        Sample aggregated performance metrics at a specific node.
        
        Generates random samples of performance metrics aggregated across
        all job classes at the specified node using CTMC simulation.
        
        Args:
            node: Network node to sample
            numSamples: Number of samples to generate
            
        Returns:
            SampleResult: Object containing aggregated sampled metrics
        """
        java_result = self.obj.sampleAggr(node.obj, numSamples)
        return SampleResult(java_result) if java_result is not None else None

    def sampleSys(self, numEvents):
        """
        Sample system-wide performance metrics.
        
        Generates random samples of performance metrics across the entire
        system for a specified number of events using CTMC simulation.
        
        Args:
            numEvents: Number of events to simulate
            
        Returns:
            SampleResult: Object containing system-wide sampled metrics
        """
        java_result = self.obj.sampleSys(numEvents)
        return SampleResult(java_result) if java_result is not None else None

    def sampleSysAggr(self, numEvents):
        """
        Sample aggregated system-wide performance metrics.
        
        Generates random samples of performance metrics aggregated across
        all nodes and job classes in the system for a specified number of events.
        
        Args:
            numEvents: Number of events to simulate
            
        Returns:
            SampleResult: Object containing aggregated system-wide samples
        """
        java_result = self.obj.sampleSysAggr(numEvents)
        return SampleResult(java_result) if java_result is not None else None

    def getAvgReward(self, reward_name=None):
        """
        Get steady-state expected reward values.

        Computes the steady-state expected reward for reward functions
        previously defined using model.setReward().

        Args:
            reward_name (str, optional): Name of specific reward to get.
                                        If None, returns all rewards as a dict.

        Returns:
            dict or float: If reward_name is None, returns a dictionary mapping
                          reward names to their expected values.
                          If reward_name is specified, returns the float value
                          for that specific reward.

        Raises:
            Exception: If reward computation fails or reward name not found.

        Example:
            # Get all rewards
            rewards = solver.getAvgReward()
            print(rewards)  # {'QueueLength': 1.5, 'Utilization': 0.75}

            # Get specific reward
            qlen = solver.getAvgReward('QueueLength')  # Returns 1.5
        """
        try:
            if reward_name is not None:
                return float(self.obj.getAvgReward(reward_name))
            else:
                java_map = self.obj.getAvgReward()
                result = {}
                for entry in java_map.entrySet():
                    result[str(entry.getKey())] = float(entry.getValue())
                return result
        except Exception as e:
            if not self._verbose_silent:
                print(f"CTMC getAvgReward failed: {e}")
            raise

    get_avg_reward = getAvgReward

    @staticmethod
    def defaultOptions():
        """Returns default solver options for the CTMC solver."""
        java_options = jpype.JPackage('jline').solvers.ctmc.SolverCTMC.defaultOptions()
        python_options = SolverOptions.__new__(SolverOptions)
        python_options.obj = java_options
        return python_options

    default_options = defaultOptions

    @staticmethod
    def supportsModel(model):
        java_solver_class = jpype.JPackage('jline').solvers.ctmc.SolverCTMC
        return java_solver_class.supports(model.obj if hasattr(model, 'obj') else model)

    get_state_space = getStateSpace
    state_space = getStateSpace

    def getCdfSysRespT(self):
        """
        Get cumulative distribution function of system response times.

        Returns:
            numpy.ndarray: System response time CDF data
        """
        try:
            java_result = self.obj.getCdfSysRespT()
            if java_result is not None:
                return jlineMatrixToArray(java_result)
            else:
                return None
        except Exception as e:
            if not self._verbose_silent:
                print(f"CTMC getCdfSysRespT failed: {e}")
            return None

    cdf_sys_respt = getCdfSysRespT
    get_cdf_sys_respt = getCdfSysRespT

    get_generator = getGenerator
    get_tran_prob = getTranProb
    get_tran_prob_aggr = getTranProbAggr
    get_tran_prob_sys = getTranProbSys
    get_tran_prob_sys_aggr = getTranProbSysAggr
    get_prob_sys_aggr = getProbSysAggr
    sample_sys_aggr = sampleSysAggr
    sample_sys = sampleSys

class SolverENV(EnsembleSolver):
    """
    ENV - Ensemble environment solver for models with regime switching.

    SolverENV analyzes queueing networks operating in randomly varying environments.
    The environment alternates between different regimes (modes), each with its own
    service and arrival characteristics. The solver combines multiple sub-solvers,
    one for each regime, and aggregates results using regime probabilities.

    The solver supports:
    - Networks with random environment (regime switching)
    - Multiple operational modes with different parameters
    - Markovian regime transitions
    - Arbitrary sub-solver for each regime
    - Ensemble averaging across regimes

    Args:
        model: Network model to solve
        lang: Language for solver execution ('python' for native, 'java' for wrapper).
              Default is 'python' which uses native Python implementation.

    This is particularly useful for modeling:
    - Systems with time-varying workloads
    - Networks with failure/repair modes
    - Performance under different operational conditions

    Note:
        Requires specifying a solver for each environment regime. The environment
        transitions are modeled as a Markov chain.
    """
    def __init__(self, *args, **kwargs):
        self._model = args[0]

        _initialize_jar_globals()
        options = SolverOptions(jpype.JPackage('jline').lang.constant.SolverType.ENV)
        super().__init__(options, *args[1:], **kwargs)
        model = args[0]
        solvers = jpype.JPackage('jline').solvers.NetworkSolver[len(args[1])]
        for i in range(len(solvers)):
            solvers[i] = args[1][i].obj
        java_network = _get_java_network(model)
        self.obj = jpype.JPackage('jline').solvers.env.SolverENV(java_network, solvers, self.solveropt.obj)

    def getEnsembleAvg(self):
        """
        Get ensemble average of performance metrics.
        
        Computes the average performance metrics across multiple
        solver instances in the ensemble.
        
        Returns:
            Ensemble average results
        """
        return self.obj.ensembleAvg()

    def printAvgTable(self):
        self.obj.printAvgTable()

    def printAvgT(self):
        """Short alias for printAvgTable."""
        self.printAvgTable()

    print_avg_table = printAvgTable
    print_avg_t = printAvgT

    def runAnalyzer(self):
        """
        Run the performance analysis.
        
        Executes the ensemble solver to analyze the network model
        and compute performance metrics across all solver instances.
        """
        self.obj.runAnalyzer()

    def getName(self):
        """
        Get the name identifier of this solver.

        Returns the string identifier for this solver instance,
        useful for solver identification and debugging.

        Returns:
            str: Solver name/identifier
        """
        return self.obj.getName()

    def getSamplePathTable(self, sample_path):
        """
        Compute transient performance metrics for a sample path through environment states.

        The method runs transient analysis for each segment of the sample path and
        extracts initial and final metric values (QLen, Util, Tput).

        Args:
            sample_path: List of tuples (stage, duration) where:
                - stage: int (0-based index) or str (stage name)
                - duration: float (time spent in stage, must be positive)

        Returns:
            pd.DataFrame: Table with columns:
                - Segment: segment number (1-based)
                - Stage: stage name
                - Duration: time spent in this segment
                - Station: station name
                - JobClass: job class name
                - InitQLen, InitUtil, InitTput: metrics at segment start
                - FinalQLen, FinalUtil, FinalTput: metrics at segment end

        Example:
            >>> sample_path = [('Fast', 5.0), ('Slow', 10.0), ('Fast', 3.0)]
            >>> table = solver.getSamplePathTable(sample_path)
        """
        import pandas as pd
        import numpy as np

        if not sample_path:
            raise ValueError("Sample path cannot be empty.")

        # Build Java sample path list
        ArrayList = jpype.JPackage('java').util.ArrayList
        java_path = ArrayList()
        for stage, duration in sample_path:
            entry = jpype.JArray(jpype.JObject)(2)
            if isinstance(stage, str):
                entry[0] = jpype.JString(stage)
            else:
                entry[0] = jpype.JInt(int(stage))
            entry[1] = jpype.JDouble(float(duration))
            java_path.add(entry)

        # Call Java method
        java_result = self.obj.getSamplePathTable(java_path)

        # Get model structure for names
        sn = self._model.getStruct(force=True)
        M = sn.nstations
        K = sn.nclasses
        stationToNode = np.asarray(sn.stationToNode).flatten()
        station_names = [str(sn.nodenames[int(stationToNode[i])]) for i in range(M)]
        class_names = list(sn.classnames)

        # Build DataFrame from Java result
        rows = []
        for seg in java_result.segments:
            seg_idx = seg.segmentIndex + 1  # 1-based for output
            stage_name = str(seg.stageName)
            duration = seg.duration

            for i in range(M):
                for r in range(K):
                    init_q = seg.initialQ.get(i, r) if seg.initialQ else 0.0
                    init_u = seg.initialU.get(i, r) if seg.initialU else 0.0
                    init_t = seg.initialT.get(i, r) if seg.initialT else 0.0
                    final_q = seg.finalQ.get(i, r) if seg.finalQ else 0.0
                    final_u = seg.finalU.get(i, r) if seg.finalU else 0.0
                    final_t = seg.finalT.get(i, r) if seg.finalT else 0.0

                    # Only include rows with non-zero metrics
                    if any([init_q, init_u, init_t, final_q, final_u, final_t]):
                        rows.append({
                            'Segment': seg_idx,
                            'Stage': stage_name,
                            'Duration': duration,
                            'Station': station_names[i],
                            'JobClass': class_names[r],
                            'InitQLen': init_q,
                            'InitUtil': init_u,
                            'InitTput': init_t,
                            'FinalQLen': final_q,
                            'FinalUtil': final_u,
                            'FinalTput': final_t,
                        })

        df = pd.DataFrame(rows)

        if not self._table_silent if hasattr(self, '_table_silent') else True:
            if not df.empty:
                print(df.to_string(index=False))

        return df

    get_sample_path_table = getSamplePathTable

    @staticmethod
    def defaultOptions():
        """Returns default solver options for the Env solver."""
        java_options = jpype.JPackage('jline').solvers.env.SolverENV.defaultOptions()
        python_options = SolverOptions.__new__(SolverOptions)
        python_options.obj = java_options
        return python_options

    default_options = defaultOptions
    get_ensemble_avg = getEnsembleAvg
    get_name = getName

class ENV(SolverENV):
    """
    Alias for SolverENV.
    """
    def __init__(self, *args, **kwargs):
        super().__init__(*args, **kwargs)



class SolverFluid(NetworkSolver):
    """
    FLD - Fluid/Mean-Field Approximation solver.

    SolverFluid uses fluid (mean-field) approximations to analyze queueing
    networks with large population sizes. The fluid approximation treats
    discrete job populations as continuous fluid flows, providing good
    approximations for heavily loaded systems.

    The solver supports:
    - Large population closed networks
    - Mixed networks with high loads
    - Multi-class networks
    - Load-dependent service rates
    - Fast computation for large-scale systems

    Methods:
        Set via solver.solveropt.method(name):

        - 'default': Standard fluid approximation using ODE integration
          * Suitable for general closed and open networks
          * Recommended for most use cases
          * Works well for systems with moderate to high loads

        - 'statedep': State-dependent fluid analysis
          * For systems with state-dependent arrival/service rates
          * Handles load-dependent scheduling policies

        - 'closing': Closing method for open networks
          * Analyzes open queueing networks with exponential arrivals
          * Approximates behavior by introducing feedback loops

        - 'diffusion': Diffusion approximation using Euler-Maruyama SDE
          * Closed networks only
          * Supports PS, FCFS, INF, SIRO scheduling disciplines
          * Uses stochastic differential equations for transient analysis
          * Particularly accurate for heavily loaded systems
          * Provides insight into system variability and stability

          Advantages:
          - More accurate than ODE for small to moderate populations
          - Captures stochastic effects (variance, tail behavior)
          - Useful for analyzing transient dynamics

          Parameters:
          - Can set parameters via solver.solveropt attributes if available
          - Typically requires population N >= 5 for reasonable accuracy

          When to use:
          - System population is moderate to large (N >= 10 recommended)
          - Need to understand stochastic behavior, not just means
          - Studying system stability and response to perturbations
          - Heavily loaded systems (utilization > 0.7)

    Args:
        model: Network model to solve
        lang: Language for solver execution ('python' for native, 'java' for wrapper).
              Default is 'python' which uses native Python implementation.

    The fluid approximation becomes more accurate as population sizes increase
    and is particularly effective for systems operating near saturation.

    Examples:
        >>> # Basic fluid analysis with default method
        >>> solver = SolverFluid(model)
        >>> solver.runAnalyzer()
        >>> results = solver.getAvgTable()

        >>> # Using diffusion approximation for closed network
        >>> solver = SolverFluid(model)
        >>> solver.solveropt.method('diffusion')  # Use diffusion approximation
        >>> results = solver.getAvgTable()

        >>> # Analyze using state-dependent method
        >>> solver = SolverFluid(model)
        >>> solver.solveropt.method('statedep')
        >>> results = solver.getAvgTable()

    Note:
        Most accurate for heavily loaded systems with large populations.
        For lightly loaded systems, other solvers (e.g., MVA, JMT) may be more appropriate.

        The diffusion method is recommended for:
        - Closed networks with N >= 10 jobs
        - Systems with utilization >= 0.7
        - When stochastic effects are important
    """

    def __init__(self, *args, **kwargs):
        _initialize_jar_globals()
        if len(args) > 1 and hasattr(args[1], 'obj'):
            options = args[1]
            super().__init__(options, *args[2:], **kwargs)
        else:
            options = SolverOptions(jpype.JPackage('jline').lang.constant.SolverType.FLUID)
            super().__init__(options, *args[1:], **kwargs)
        self.model = args[0]
        self._model = self.model  # For consistency with other solvers
        java_network = _get_java_network(self.model)
        self.obj = jpype.JPackage('jline').solvers.fluid.SolverFluid(java_network, self.solveropt.obj)

    def runAnalyzer(self):
        """Run the Fluid analysis."""
        self.obj.runAnalyzer()
        return self

    def getAvgTable(self):
        """Get comprehensive average performance metrics table."""
        return super().getAvgTable()

    avg_table = getAvgTable
    getAvgT = getAvgTable
    avgT = getAvgTable
    aT = getAvgTable

    def getTranAvg(self):

        result = self.obj.result

        M = result.QNt.length
        K = result.QNt[0].length

        def extract(metrics):
            return [[jlineMatrixToArray(metrics[i][k]) for k in range(K)] for i in range(M)]

        return {
            'QNt': extract(result.QNt),
            'UNt': extract(result.UNt),
            'TNt': extract(result.TNt),
        }

    tran_avg = getTranAvg

    def getCdfRespT(self, R=None):
        try:
            java_result = self.obj.cdfRespT()

            if hasattr(self.obj, 'result') and self.obj.result is not None:
                fluid_result = self.obj.result
                if hasattr(fluid_result, 'distribC') and fluid_result.distribC is not None:
                    distribC = fluid_result.distribC

                    num_stations = self.model.get_number_of_stations()
                    num_classes = self.model.get_number_of_classes()

                    CdfRespT = []
                    for i in range(num_stations):
                        station_cdfs = []
                        for c in range(num_classes):
                            if (i < distribC.length and
                                distribC[i] is not None and
                                c < distribC[i].length and
                                distribC[i][c] is not None):
                                cdf_array = jlineMatrixToArray(distribC[i][c])
                                station_cdfs.append(cdf_array)
                            else:
                                station_cdfs.append(None)
                        CdfRespT.append(station_cdfs)

                    return CdfRespT

            num_stations = self.model.get_number_of_stations()
            num_classes = self.model.get_number_of_classes()
            return [[None for _ in range(num_classes)] for _ in range(num_stations)]

        except Exception as e:
            if not self._verbose_silent:
                print(f"Fluid getCdfRespT failed: {e}")
            try:
                num_stations = self.model.get_number_of_stations()
                num_classes = self.model.get_number_of_classes()
                return [[None for _ in range(num_classes)] for _ in range(num_stations)]
            except:
                return [[]]

    cdf_respt = getCdfRespT
    cdf_resp_t = getCdfRespT
    getSjrnT = getCdfRespT
    sjrn_t = getSjrnT

    @staticmethod
    def defaultOptions():
        """Returns default solver options for the Fluid solver."""
        java_options = jpype.JPackage('jline').solvers.fluid.SolverFluid.defaultOptions()
        python_options = SolverOptions.__new__(SolverOptions)
        python_options.obj = java_options
        return python_options

    default_options = defaultOptions

    @staticmethod
    def supportsModel(model):
        java_solver_class = jpype.JPackage('jline').solvers.fluid.SolverFluid
        return java_solver_class.supports(model.obj if hasattr(model, 'obj') else model)

    def getCdfPassT(self, R=None):
        """
        Get cumulative distribution function of passage times.

        Args:
            R: Optional response time handle

        Returns:
            list: CDF data for passage times per station and class
        """
        try:
            if R is None:
                java_result = self.obj.getCdfPassT()
            else:
                java_result = self.obj.getCdfPassT(R.obj if hasattr(R, 'obj') else R)

            if java_result is None:
                num_stations = self.model.get_number_of_stations()
                num_classes = self.model.get_number_of_classes()
                return [[None for _ in range(num_classes)] for _ in range(num_stations)]

            num_stations = java_result.numStations
            num_classes = java_result.numClasses

            CdfPassT = []
            for i in range(num_stations):
                station_cdfs = []
                for c in range(num_classes):
                    if java_result.hasCdf(i, c):
                        cdf_matrix = java_result.getCdf(i, c)
                        if cdf_matrix is not None:
                            cdf_array = jlineMatrixToArray(cdf_matrix)
                            station_cdfs.append(cdf_array)
                        else:
                            station_cdfs.append(None)
                    else:
                        station_cdfs.append(None)
                CdfPassT.append(station_cdfs)

            return CdfPassT
        except Exception as e:
            if not self._verbose_silent:
                print(f"Fluid getCdfPassT failed: {e}")
            try:
                num_stations = self.model.get_number_of_stations()
                num_classes = self.model.get_number_of_classes()
                return [[None for _ in range(num_classes)] for _ in range(num_stations)]
            except:
                return [[]]

    cdf_passt = getCdfPassT

    get_tran_avg = getTranAvg
    get_cdf_respt = getCdfRespT
    get_cdf_passt = getCdfPassT


class SolverJMT(NetworkSolver):
    """
    Java Modelling Tools (JMT) simulation solver.

    SolverJMT provides discrete-event simulation capabilities through integration
    with Java Modelling Tools. This solver can handle arbitrary network topologies
    and non-product-form features that cannot be solved analytically.

    The solver supports:
    - Open, closed, and mixed networks
    - Arbitrary service time distributions
    - Complex routing strategies
    - Fork-join structures
    - Load-dependent and state-dependent behavior
    - Non-product-form scheduling disciplines
    - Confidence intervals for statistical results

    Results include mean performance measures with confidence intervals.
    The simulation length is controlled by the 'samples' parameter.

    Note:
        Requires JMT.jar in the classpath. Simulation results are statistical
        approximations with confidence intervals.

    Args:
        model: Network model to solve
    """

    def __init__(self, *args, **kwargs):
        _initialize_jar_globals()
        if len(args) > 1 and hasattr(args[1], 'obj'):
            options = args[1]
            super().__init__(options, *args[2:], **kwargs)
        else:
            options = SolverOptions(jpype.JPackage('jline').lang.constant.SolverType.JMT)
            super().__init__(options, *args[1:], **kwargs)
        self.model = args[0]

        package_dir = os.path.dirname(os.path.abspath(__file__))
        python_dir = os.path.dirname(package_dir)
        root_dir = os.path.dirname(python_dir)
        common_dir = os.path.join(root_dir, 'common')
        jmt_jar_path = os.path.join(common_dir, 'JMT.jar')

        if not os.path.isfile(jmt_jar_path):
            print("JMT.jar not found in", common_dir)
            print("Attempting to download JMT.jar...")
            try:
                from urllib.request import urlretrieve
                jmt_url = 'https://line-solver.sourceforge.net/latest/JMT.jar'
                urlretrieve(jmt_url, jmt_jar_path)
                print("Successfully downloaded JMT.jar to", common_dir)
            except Exception as e:
                raise RuntimeError(f"Failed to download JMT.jar: {e}\n"
                                 f"Please manually download https://line-solver.sourceforge.net/latest/JMT.jar "
                                 f"and place it in {common_dir}")

        self.jmtPath = jpype.JPackage('java').lang.String(jmt_jar_path)
        java_network = _get_java_network(self.model)
        self.obj = jpype.JPackage('jline').solvers.jmt.SolverJMT(java_network, self.solveropt.obj, self.jmtPath)

    def runAnalyzer(self):
        """Run the JMT analyzer."""
        self.obj.runAnalyzer()
        return self

    def getAvg(self):
        """Get average performance metrics."""
        return super().getAvg()

    def jsimwView(self):
        self.obj.jsimwView(self.jmtPath)

    def jsimgView(self):
        self.obj.jsimgView(self.jmtPath)

    jsimg_view = jsimgView
    jsimw_view = jsimwView

    def sampleAggr(self, node, numEvents=None, markActivePassive=False):
        """Sample from a specific node using aggregated states"""
        try:
            if numEvents is None:
                java_result = self.obj.sampleAggr(node.obj)
            else:
                java_result = self.obj.sampleAggr(node.obj, numEvents, markActivePassive)
            return java_result
        except Exception as e:
            if not self._verbose_silent:
                print(f"JMT aggregated sampling failed: {e}")
            return None

    def sampleSysAggr(self, numEvents=None, markActivePassive=False):
        """Sample system-wide using aggregated states"""
        try:
            if numEvents is None:
                java_result = self.obj.sampleSysAggr()
            else:
                java_result = self.obj.sampleSysAggr(numEvents, markActivePassive)
            return java_result
        except Exception as e:
            if not self._verbose_silent:
                print(f"JMT system aggregated sampling failed: {e}")
            return None

    def getProbSysAggr(self):
        try:
            java_result = self.obj.probSysAggr()
            return ProbabilityResult(java_result)
        except Exception as e:
            if not self._verbose_silent:
                print(f"JMT getProbSysAggr failed: {e}")
            return None

    def getCdfRespT(self, R=None):
        try:
            if R is None:
                java_result = self.obj.cdfRespT()
            else:
                java_result = self.obj.getCdfRespT(R.obj if hasattr(R, 'obj') else R)

            if java_result is None:
                num_stations = self.model.get_number_of_stations()
                num_classes = self.model.get_number_of_classes()
                return [[None for _ in range(num_classes)] for _ in range(num_stations)]

            num_stations = java_result.numStations
            num_classes = java_result.numClasses

            CdfRespT = []
            for i in range(num_stations):
                station_cdfs = []
                for c in range(num_classes):
                    if java_result.hasCdf(i, c):
                        cdf_matrix = java_result.getCdf(i, c)
                        if cdf_matrix is not None:
                            cdf_array = jlineMatrixToArray(cdf_matrix)
                            station_cdfs.append(cdf_array)
                        else:
                            station_cdfs.append(None)
                    else:
                        station_cdfs.append(None)
                CdfRespT.append(station_cdfs)

            return CdfRespT
        except Exception as e:
            if not self._verbose_silent:
                print(f"JMT getCdfRespT failed: {e}")
            try:
                num_stations = self.model.get_number_of_stations()
                num_classes = self.model.get_number_of_classes()
                return [[None for _ in range(num_classes)] for _ in range(num_stations)]
            except:
                return [[]]

    getSjrnT = getCdfRespT
    sjrn_t = getSjrnT

    def getTranCdfRespT(self, R=None):
        try:
            if R is None:
                java_result = self.obj.getTranCdfRespT()
            else:
                java_result = self.obj.getTranCdfRespT(R.obj if hasattr(R, 'obj') else R)

            if java_result is None:
                num_stations = self.model.get_number_of_stations()
                num_classes = self.model.get_number_of_classes()
                return [[None for _ in range(num_classes)] for _ in range(num_stations)]

            num_stations = java_result.numStations
            num_classes = java_result.numClasses

            CdfRespT = []
            for i in range(num_stations):
                station_cdfs = []
                for c in range(num_classes):
                    if java_result.hasCdf(i, c):
                        cdf_matrix = java_result.getCdf(i, c)
                        if cdf_matrix is not None:
                            cdf_array = jlineMatrixToArray(cdf_matrix)
                            station_cdfs.append(cdf_array)
                        else:
                            station_cdfs.append(None)
                    else:
                        station_cdfs.append(None)
                CdfRespT.append(station_cdfs)

            return CdfRespT
        except Exception as e:
            if not self._verbose_silent:
                print(f"JMT getTranCdfRespT failed: {e}")
            try:
                num_stations = self.model.get_number_of_stations()
                num_classes = self.model.get_number_of_classes()
                return [[None for _ in range(num_classes)] for _ in range(num_stations)]
            except:
                return [[]]

    tran_cdf_respt = getTranCdfRespT

    def getTranCdfPassT(self, R=None):
        try:
            if R is None:
                java_result = self.obj.getTranCdfPassT()
            else:
                java_result = self.obj.getTranCdfPassT(R.obj if hasattr(R, 'obj') else R)
            return DistributionResult(java_result)
        except Exception as e:
            if not self._verbose_silent:
                print(f"JMT getTranCdfPassT failed: {e}")
            return None

    tran_cdf_passt = getTranCdfPassT

    def getTranAvg(self, Qt=None, Ut=None, Tt=None):
        """Get transient average station metrics.
        
        Args:
            Qt: Optional queue length handles 
            Ut: Optional utilization handles
            Tt: Optional throughput handles
            
        Returns:
            Tuple of (QNclass_t, UNclass_t, TNclass_t) containing transient metrics
            Each element is a list of lists [station][class] containing time series data
        """
        try:
            # Check that timespan is finite for transient analysis
            options = self.obj.getOptions()
            timespan_end = options.timespan[1]
            if timespan_end == float('inf'):
                raise RuntimeError("Transient analysis requires finite timespan. "
                                 "Use SolverJMT(model, timespan=[0, T]) to specify timespan.")
            
            # Call Java getTranAvg method to populate transient results
            self.obj.getTranAvg()
            
            # Check if transient results are available
            if not hasattr(self.obj, 'result') or self.obj.result is None:
                print("Warning: Transient results not available. Check solver execution.")
                return None, None, None
            
            result = self.obj.result
            if not hasattr(result, 'QNt') or result.QNt is None:
                print("Warning: Transient data not available. Check solver execution.")
                return None, None, None
            
            # Get network dimensions
            M = self.model.get_number_of_stations()
            K = self.model.get_number_of_classes()
            
            # Initialize output lists
            QNclass_t = [[None for _ in range(K)] for _ in range(M)]
            UNclass_t = [[None for _ in range(K)] for _ in range(M)]
            TNclass_t = [[None for _ in range(K)] for _ in range(M)]
            
            # Extract transient results from Java solver
            for k in range(K):
                for ist in range(M):
                    # Queue length transients
                    try:
                        java_result = self.obj.getTranQLen()
                        if java_result is not None and len(java_result) > ist and len(java_result[ist]) > k:
                            if java_result[ist][k] is not None:
                                ret = jlineMatrixToArray(java_result[ist][k])
                                if ret is not None and len(ret) > 0 and len(ret[0]) >= 2:
                                    QNclass_t[ist][k] = {
                                        'handle': (self.model.getStations()[ist], self.model.getClasses()[k]),
                                        't': [row[1] for row in ret],  # time points
                                        'metric': [row[0] for row in ret],  # queue length values
                                        'isaggregate': True
                                    }
                    except:
                        pass
                    
                    # Utilization transients
                    try:
                        java_result = self.obj.getTranUtil()
                        if java_result is not None and len(java_result) > ist and len(java_result[ist]) > k:
                            if java_result[ist][k] is not None:
                                ret = jlineMatrixToArray(java_result[ist][k])
                                if ret is not None and len(ret) > 0 and len(ret[0]) >= 2:
                                    UNclass_t[ist][k] = {
                                        'handle': (self.model.getStations()[ist], self.model.getClasses()[k]),
                                        't': [row[1] for row in ret],  # time points  
                                        'metric': [row[0] for row in ret],  # utilization values
                                        'isaggregate': True
                                    }
                    except:
                        pass
                    
                    # Throughput transients
                    try:
                        java_result = self.obj.getTranTput()
                        if java_result is not None and len(java_result) > ist and len(java_result[ist]) > k:
                            if java_result[ist][k] is not None:
                                ret = jlineMatrixToArray(java_result[ist][k])
                                if ret is not None and len(ret) > 0 and len(ret[0]) >= 2:
                                    TNclass_t[ist][k] = {
                                        'handle': (self.model.getStations()[ist], self.model.getClasses()[k]),
                                        't': [row[1] for row in ret],  # time points
                                        'metric': [row[0] for row in ret],  # throughput values
                                        'isaggregate': True
                                    }
                    except:
                        pass
            
            return QNclass_t, UNclass_t, TNclass_t
            
        except Exception as e:
            raise RuntimeError(f"Failed to compute transient metrics: {e}")

    @staticmethod
    def defaultOptions():
        """Returns default solver options for the JMT solver."""
        java_options = jpype.JPackage('jline').solvers.jmt.SolverJMT.defaultOptions()
        python_options = SolverOptions.__new__(SolverOptions)
        python_options.obj = java_options
        return python_options

    default_options = defaultOptions

    @staticmethod
    def supportsModel(model):
        java_solver_class = jpype.JPackage('jline').solvers.jmt.SolverJMT
        return java_solver_class.supports(model.obj if hasattr(model, 'obj') else model)

    def QN2JSIMG(self, sn=None, outputFileName=None):
        """
        Writes queueing network model to JMT JSIMG format.
        Delegates to the standalone QN2JSIMG function in the io package.

        Args:
            sn: Network structure (optional, ignored - uses model from solver)
            outputFileName: Output file name (optional)

        Returns:
            Path to the JSIM file
        """
        from .io import QN2JSIMG as io_qn2jsimg
        return io_qn2jsimg(self.model, outputFileName)

    get_prob_sys_aggr = getProbSysAggr
    prob_sys_aggr = getProbSysAggr
    get_cdf_respt = getCdfRespT
    cdf_respt = getCdfRespT
    cdf_resp_t = getCdfRespT
    get_tran_cdf_respt = getTranCdfRespT
    get_tran_cdf_passt = getTranCdfPassT
    sample_sys_aggr = sampleSysAggr
    #sample_sys = sampleSys


class SolverMAM(NetworkSolver):
    """
    Matrix-Analytic Methods (MAM) solver.

    SolverMAM implements matrix-analytic methods for analyzing queueing systems
    with structured Markov chains, including quasi-birth-death (QBD) processes,
    Markovian arrival processes (MAP), and rational arrival processes (RAP).

    The solver supports:
    - MAP/MAP/1 and related queues
    - QBD-type Markov chains
    - Phase-type service distributions
    - Batch arrivals and services
    - Multi-server queues with structured arrivals
    - Response time distributions (via Laplace transforms)

    Args:
        model: Network model to solve
        lang: Language for solver execution ('python' for native, 'java' for wrapper).
              Default is 'python' which uses native Python implementation.

    Matrix-analytic methods are particularly effective for:
    - Queues with correlated arrivals (MAP, MMPP)
    - Systems with complex service mechanisms
    - Models requiring transient or passage time analysis

    The solver computes exact solutions through iterative algorithms for
    rate matrices (R, G matrices) and probability vectors.

    Note:
        Requires phase-type or MAP distributions for accurate analysis.
        More computationally intensive than MVA but handles broader model classes.
    """
    def __init__(self, *args, **kwargs):
        self._model = args[0]

        _initialize_jar_globals()
        if len(args) > 1 and hasattr(args[1], 'obj'):
            options = args[1]
            super().__init__(options, *args[2:], **kwargs)
        else:
            options = SolverOptions(jpype.JPackage('jline').lang.constant.SolverType.MAM)
            super().__init__(options, *args[1:], **kwargs)
        model = args[0]
        java_network = _get_java_network(model)
        self.obj = jpype.JPackage('jline').solvers.mam.SolverMAM(java_network, self.solveropt.obj)

    def runAnalyzer(self):
        """Run the MAM analysis."""
        self.obj.runAnalyzer()
        return self

    def getAvgTable(self):
        """Get comprehensive average performance metrics table."""
        return super().getAvgTable()

    avg_table = getAvgTable
    getAvgT = getAvgTable
    avgT = getAvgTable
    aT = getAvgTable

    @staticmethod
    def defaultOptions():
        """Returns default solver options for the MAM solver."""
        java_options = jpype.JPackage('jline').solvers.mam.SolverMAM.defaultOptions()
        python_options = SolverOptions.__new__(SolverOptions)
        python_options.obj = java_options
        return python_options

    default_options = defaultOptions

    def getCdfPassT(self, R=None):
        """
        Get cumulative distribution function of passage times.

        Args:
            R: Optional response time handle

        Returns:
            list: CDF data for passage times per station and class
        """
        try:
            if R is None:
                java_result = self.obj.getCdfPassT()
            else:
                java_result = self.obj.getCdfPassT(R.obj if hasattr(R, 'obj') else R)

            if java_result is None:
                num_stations = self.model.get_number_of_stations()
                num_classes = self.model.get_number_of_classes()
                return [[None for _ in range(num_classes)] for _ in range(num_stations)]

            num_stations = java_result.numStations
            num_classes = java_result.numClasses

            CdfPassT = []
            for i in range(num_stations):
                station_cdfs = []
                for c in range(num_classes):
                    if java_result.hasCdf(i, c):
                        cdf_matrix = java_result.getCdf(i, c)
                        if cdf_matrix is not None:
                            cdf_array = jlineMatrixToArray(cdf_matrix)
                            station_cdfs.append(cdf_array)
                        else:
                            station_cdfs.append(None)
                    else:
                        station_cdfs.append(None)
                CdfPassT.append(station_cdfs)

            return CdfPassT
        except Exception as e:
            if not self._verbose_silent:
                print(f"MAM getCdfPassT failed: {e}")
            try:
                num_stations = self.model.get_number_of_stations()
                num_classes = self.model.get_number_of_classes()
                return [[None for _ in range(num_classes)] for _ in range(num_stations)]
            except:
                return [[]]

    cdf_passt = getCdfPassT
    get_cdf_passt = getCdfPassT

    def getCdfRespT(self, R=None):
        """
        Get cumulative distribution function of response times.

        Args:
            R: Optional response time handle

        Returns:
            list: CDF data for response times per station and class
        """
        try:
            if R is None:
                java_result = self.obj.cdfRespT()
            else:
                java_result = self.obj.getCdfRespT(R.obj if hasattr(R, 'obj') else R)

            if java_result is None:
                num_stations = self.model.get_number_of_stations()
                num_classes = self.model.get_number_of_classes()
                return [[None for _ in range(num_classes)] for _ in range(num_stations)]

            num_stations = java_result.numStations
            num_classes = java_result.numClasses

            CdfRespT = []
            for i in range(num_stations):
                station_cdfs = []
                for c in range(num_classes):
                    if java_result.hasCdf(i, c):
                        cdf_matrix = java_result.getCdf(i, c)
                        if cdf_matrix is not None:
                            cdf_array = jlineMatrixToArray(cdf_matrix)
                            station_cdfs.append(cdf_array)
                        else:
                            station_cdfs.append(None)
                    else:
                        station_cdfs.append(None)
                CdfRespT.append(station_cdfs)

            return CdfRespT
        except Exception as e:
            if not self._verbose_silent:
                print(f"MAM getCdfRespT failed: {e}")
            try:
                num_stations = self.model.get_number_of_stations()
                num_classes = self.model.get_number_of_classes()
                return [[None for _ in range(num_classes)] for _ in range(num_stations)]
            except:
                return [[]]

    cdf_respt = getCdfRespT
    get_cdf_respt = getCdfRespT
    getSjrnT = getCdfRespT
    sjrn_t = getSjrnT

    def getPerctRespT(self, percentiles=None):
        """
        Get response time percentiles from Fork-Join analysis.

        This method extracts response time percentiles for Fork-Join queueing models.
        It uses the FJ_codes algorithm to compute accurate tail percentiles for systems
        with parallel service queues. The method automatically detects valid Fork-Join
        topologies and applies the appropriate analysis.

        Only works for Fork-Join topologies with the structure:
        Source → Fork → K parallel Queues → Join → Sink

        Args:
            percentiles: List or array of percentile levels (0-100 scale).
                        Default: [50, 90, 95, 99]
                        Examples: [50, 90, 95, 99] for median, 90th, 95th, 99th
                                 [75, 95, 99.9] for 75th, 95th, 99.9th

        Returns:
            list: List of dictionaries, one per job class, each containing:
                - 'class': int - Job class index (0-based)
                - 'percentiles': numpy array - Requested percentile levels (0-100)
                - 'values': numpy array - Percentile values (response times)
                - 'K': int - Number of parallel queues in the Fork-Join system
                - 'method': str - Analysis method used ("FJ_NARE" or "FJ_Sylves")

        Raises:
            RuntimeError: If the model is not a valid Fork-Join topology,
                         or if the solver has not been run yet, or if
                         percentile results are unavailable.

        Example:
            >>> model = create_fork_join_model(K=10, lam=0.5, mu=1.0)
            >>> solver = SolverMAM(model)
            >>> solver.run()
            >>> results = solver.getPerctRespT([50, 90, 95, 99])
            >>> for r in results:
            ...     print(f"Class {r['class']}: median={r['values'][0]:.3f}, "
            ...           f"95th={r['values'][2]:.3f}")
            Class 0: median=1.234, 95th=3.567

        Note:
            The Fork-Join topology must satisfy:
            - Open classes only
            - Homogeneous service distributions across parallel queues
            - Supported distributions: Exp, HyperExp(2), Erlang(2), MAP(2)
            - FCFS or PS scheduling at queues

        References:
            Z. Qiu, J.F. Pérez, and P. Harrison, "Beyond the Mean in Fork-Join
            Queues: Efficient Approximation for Response-Time Tails",
            IFIP Performance 2015.
        """
        if percentiles is None:
            percentiles = [50.0, 90.0, 95.0, 99.0]

        try:
            # Convert Python list/array to Java double array
            java_percentiles = jpype.JArray(jpype.JDouble)(len(percentiles))
            for i, p in enumerate(percentiles):
                java_percentiles[i] = float(p)

            # Call Java method
            java_results = self.obj.getPerctRespT(java_percentiles)

            if java_results is None or len(java_results) == 0:
                raise RuntimeError(
                    "No percentile results available. "
                    "Ensure the model is a valid Fork-Join topology "
                    "and the solver has been run successfully."
                )

            # Convert Java List<FJPercentileResult> to Python list of dicts
            results = []
            for java_result in java_results:
                # Convert Java arrays to numpy arrays
                perc_array = []
                for i in range(len(java_result.getPercentiles())):
                    perc_array.append(java_result.getPercentiles()[i])

                vals_array = []
                for i in range(len(java_result.getRTp())):
                    vals_array.append(java_result.getRTp()[i])

                results.append({
                    'class': java_result.getJobClass(),
                    'percentiles': jpype.numpy.array(perc_array),
                    'values': jpype.numpy.array(vals_array),
                    'K': java_result.getK(),
                    'method': str(java_result.getMethod())
                })

            return results

        except jpype.JException as e:
            raise RuntimeError(
                f"Failed to get percentiles: {str(e)}. "
                "Ensure the model is a valid Fork-Join topology "
                "(Source → Fork → K Queues → Join → Sink) "
                "and the solver has been run."
            ) from e
        except Exception as e:
            raise RuntimeError(
                f"Unexpected error getting percentiles: {str(e)}"
            ) from e

    get_perct_resp_t = getPerctRespT

    def me_open(self, tol=1e-6, max_iter=1000, verbose=False):
        """
        Maximum Entropy algorithm for Open Queueing Networks.

        Applies the ME algorithm from Kouvatsos (1994) to the model.
        Only supports open queueing networks (no closed classes).

        Args:
            tol: Convergence tolerance (default: 1e-6)
            max_iter: Maximum number of iterations (default: 1000)
            verbose: Print iteration information (default: False)

        Returns:
            dict: Dictionary containing:
                - 'L': Mean queue lengths (numpy array [M x R])
                - 'W': Mean waiting times (numpy array [M x R])
                - 'Ca': Arrival SCV at each queue (numpy array [M x R])
                - 'Cd': Departure SCV at each queue (numpy array [M x R])
                - 'lambda': Total arrival rates (numpy array [M x R])
                - 'rho': Utilizations (numpy array [M x R])

        Reference:
            Kouvatsos (1994) "Entropy Maximisation and Queueing Network Models"
        """
        try:
            # Create ME options
            MeOqnOptions = jpype.JPackage('jline').api.mam.MeOqnOptions
            me_options = MeOqnOptions(tol, max_iter, verbose)

            # Call Java meOpen method
            java_result = self.obj.meOpen(me_options)

            if java_result is None:
                return None

            # Convert Java matrices to Python arrays
            result = {
                'L': jlineMatrixToArray(java_result.QN),
                'W': jlineMatrixToArray(java_result.RN),
                'lambda': jlineMatrixToArray(java_result.TN),
                'rho': jlineMatrixToArray(java_result.UN)
            }

            return result

        except Exception as e:
            if not self._verbose_silent:
                print(f"MAM me_open failed: {e}")
            return None

    def getProbMarg(self, node, jobclass, state_m=None):
        """
        Probability distribution for queue length of a SINGLE class at a station.
        Returns P(n jobs of class r) for n=0,1,...,N(r).

        Compare with get_prob_aggr: returns probability of a specific per-class
        distribution, e.g., P(2 class-1, 1 class-2) as a scalar.

        Args:
            node: Station/node index (0-based) or node object
            jobclass: Job class index (0-based)
            state_m: Optional state levels to query (None for all states).
                     If provided, returns probabilities only for the specified states.

        Returns:
            numpy.ndarray: Probability distribution where prob[n] = P(n jobs of this class).
                          For closed models, length is population+1.
                          For open models with cutoff, length is cutoff.

        Raises:
            RuntimeError: If the solver fails to compute probabilities
            ValueError: If node or jobclass indices are invalid

        Example:
            >>> solver = SolverMAM(model)
            >>> solver.run()
            >>> prob = solver.getProbMarg(0, 0)  # Station 0, Class 0
            >>> print(f"P(0 jobs) = {prob[0]:.4f}")
            >>> print(f"P(1 job) = {prob[1]:.4f}")

        Note:
            For networks with multiple queues, this method is only supported
            for single-queue models. Use SolverCTMC for multi-queue networks.
        """
        # Convert node to index if needed
        if hasattr(node, 'obj'):
            node_idx = self.model.getNodes().indexOf(node.obj)
        else:
            node_idx = int(node)

        # Convert state_m to Java Matrix if provided
        java_state = None
        if state_m is not None:
            java_state = jlineMatrixFromArray(np.array(state_m).reshape(1, -1))

        try:
            result = self.obj.getProbMarg(node_idx, int(jobclass), java_state)

            # Extract probability distribution from ProbabilityResult
            if result.hasDistribution():
                prob_vector = result.getProbabilityVector()
                return jlineMatrixToArray(prob_vector).flatten()
            else:
                # Scalar result - return as single-element array
                return np.array([result.probability])

        except jpype.JException as e:
            error_msg = str(e)
            if "UnsupportedOperationException" in error_msg:
                raise RuntimeError(
                    f"getProbMarg is not supported for this model configuration: {error_msg}. "
                    "For networks with multiple queues, use SolverCTMC instead."
                ) from e
            elif "IllegalArgumentException" in error_msg:
                raise ValueError(f"Invalid node or jobclass index: {error_msg}") from e
            else:
                raise RuntimeError(f"Failed to compute marginal probabilities: {error_msg}") from e
        except Exception as e:
            raise RuntimeError(f"Unexpected error in getProbMarg: {e}") from e

    get_prob_marg = getProbMarg

    @staticmethod
    def supportsModel(model):
        java_solver_class = jpype.JPackage('jline').solvers.mam.SolverMAM
        return java_solver_class.supports(model.obj if hasattr(model, 'obj') else model)

class SolverMVA(NetworkSolver):
    """
    Mean Value Analysis (MVA) solver.

    SolverMVA implements the Mean Value Analysis algorithm for solving
    product-form queueing networks. MVA provides exact results for networks
    that satisfy product-form conditions and is generally efficient for
    moderate-sized closed networks.

    The solver supports:
    - Closed queueing networks
    - Mixed networks (both open and closed classes)
    - Load-dependent service rates
    - Multiple job classes
    - Priority scheduling at some nodes

    The algorithm works by iteratively computing mean performance measures
    using the arrival theorem and Little's law.

    Args:
        model: Network model to solve
        lang: Language for solver execution ('python' for native, 'java' for wrapper).
              Default is 'python' which uses native Python implementation.

    Note:
        Limited to product-form networks. Non-product-form features require
        other solvers like SolverSSA or SolverJMT.
    """

    def __init__(self, *args, **kwargs):
        self._model = args[0]

        _initialize_jar_globals()
        if len(args) > 1 and hasattr(args[1], 'obj'):
            options = args[1]
            super().__init__(options, *args[2:], **kwargs)
        else:
            try:
                options = SolverMVA.defaultOptions()
            except:
                import line_solver
                options = line_solver.lineDefaults()
            super().__init__(options, *args[1:], **kwargs)
        model = args[0]

        # Get Java network - handles both wrapper and pure Python models
        java_network = _get_java_network(model)

        try:
            self.obj = jpype.JPackage('jline').solvers.mva.SolverMVA(java_network, self.solveropt.obj)
        except jpype.JException as e:
            if "Outside of matrix bounds" in str(e):
                import line_solver
                working_options = line_solver.lineDefaults()
                self.solveropt = working_options
                self.obj = jpype.JPackage('jline').solvers.mva.SolverMVA(java_network, self.solveropt.obj)
            else:
                raise

    def runAnalyzer(self):
        """Run the MVA analysis."""
        self.obj.runAnalyzer()
        return self

    def getAvgTable(self):
        """Get comprehensive average performance metrics table."""
        return super().getAvgTable()

    avg_table = getAvgTable
    getAvgT = getAvgTable
    avgT = getAvgTable
    aT = getAvgTable

    @staticmethod
    def defaultOptions():
        """Returns default solver options for the MVA solver."""
        java_options = jpype.JPackage('jline').solvers.mva.SolverMVA.defaultOptions()
        python_options = SolverOptions.__new__(SolverOptions)
        python_options.obj = java_options
        return python_options

    default_options = defaultOptions

    @staticmethod
    def supportsModel(model):
        java_solver_class = jpype.JPackage('jline').solvers.mva.SolverMVA
        return java_solver_class.supports(model.obj if hasattr(model, 'obj') else model)

class SolverQNS(NetworkSolver):
    """
    Queueing Network Solver (QNS) - wrapper for LQNS's qns utility.

    SolverQNS provides access to the qns solver from the Layered Queueing Network
    Solver (LQNS) toolkit. It implements various approximation methods for
    queueing networks, including Conway's method and variations.

    The solver supports:
    - Closed and mixed queueing networks
    - Multiple job classes
    - Various approximation methods (Conway, Rolia, Zhou, etc.)
    - Product-form and non-product-form features

    Available approximation methods:
    - conway: Conway's approximation
    - rolia: Rolia's method
    - zhou: Zhou's approximation
    - suri: Suri's approximation
    - reiser: Reiser's method
    - schmidt: Schmidt's method

    Args:
        model: Network model to analyze
        method: Approximation method to use

    Note:
        Requires the 'qnsolver' command-line utility to be installed and available
        in the system PATH. Check availability with SolverQNS.isAvailable().
    """
    def __init__(self, *args, method='default', **kwargs):
        self._model = args[0] if args else kwargs.get('model')
        self._method = method
        self._table_silent = kwargs.get('table_silent', False)

        _initialize_jar_globals()
        options = SolverOptions(jpype.JPackage('jline').lang.constant.SolverType.QNS)
        super().__init__(options, *args[1:], **kwargs)
        model = args[0]
        java_network = _get_java_network(model)
        self.obj = jpype.JPackage('jline').solvers.qns.SolverQNS(java_network, self.solveropt.obj)

    def runAnalyzer(self):
        """Run the QNS analysis."""
        return super().runAnalyzer()

    def getName(self):
        """Get the name of this solver."""
        return self.obj.getName()

    def getAvgTable(self):
        """Get comprehensive average performance metrics table."""
        return super().getAvgTable()

    def listValidMethods(self):
        """List valid methods for the QNS solver"""
        return ['default', 'conway', 'rolia', 'zhou', 'suri', 'reiser', 'schmidt']

    @staticmethod
    def isAvailable():
        """Check if the QNS solver is available"""
        try:
            return jpype.JPackage('jline').solvers.qns.SolverQNS.isAvailable()
        except:
            return False

    @staticmethod
    def defaultOptions():
        """Returns default solver options for the QNS solver."""
        java_options = jpype.JPackage('jline').solvers.qns.SolverQNS.defaultOptions()
        python_options = SolverOptions.__new__(SolverOptions)
        python_options.obj = java_options
        return python_options

    default_options = defaultOptions
    run_analyzer = runAnalyzer
    get_avg_table = getAvgTable
    get_name = getName

    @staticmethod
    def supportsModel(model):
        java_solver_class = jpype.JPackage('jline').solvers.qns.SolverQNS
        return java_solver_class.supports(model.obj if hasattr(model, 'obj') else model)

class SolverLQNS(Solver):
    """
    Layered Queueing Network Solver (LQNS) - wrapper for external LQNS tool.

    SolverLQNS provides integration with the external Layered Queueing Network
    Solver, a specialized tool for analyzing layered queueing networks (LQNs)
    commonly used in software performance engineering.

    The solver supports:
    - Layered queueing networks (client-server systems)
    - Processor-task-entry hierarchies
    - Synchronous and asynchronous calls
    - Multi-tier architectures
    - Activity graphs
    - Precedence constraints

    Layered queueing networks are particularly useful for modeling:
    - Distributed systems and microservices
    - Client-server applications
    - Multi-tier web applications
    - Software with nested service calls

    Args:
        model: LayeredNetwork model to analyze
        method: Analysis method (default, lqns, srvn, exactmva, sim, lqsim)

    Note:
        Requires the 'lqns' and 'lqsim' command-line utilities to be installed
        and available in the system PATH. Check availability with
        SolverLQNS.isAvailable() before creating solver instances.
    """
    def __init__(self, *args, method='default', **kwargs):
        self._model = args[0] if args else kwargs.get('model')
        self._method = method
        self._table_silent = kwargs.get('table_silent', False)

        _initialize_jar_globals()
        options = SolverOptions(jpype.JPackage('jline').lang.constant.SolverType.LQNS)
        super().__init__(options, *args[1:], **kwargs)

        try:
            java_network = _get_java_network(self._model)
            self.obj = jpype.JPackage('jline').solvers.lqns.SolverLQNS(java_network, self.solveropt.obj)
        except Exception as e:
            error_msg = str(e)
            if "lqns" in error_msg.lower() and "lqsim" in error_msg.lower():
                raise RuntimeError(
                    "SolverLQNS requires the 'lqns' and 'lqsim' commands to be available in your system PATH.\n"
                    "You can install them from: http://www.sce.carleton.ca/rads/lqns/\n\n"
                    "Alternatively, use remote execution via Docker:\n"
                    "  1. Pull and run: docker run -d -p 8080:8080 imperialqore/lqns-rest:latest\n"
                    "  2. Configure remote execution:\n"
                    "     options = SolverOptions()\n"
                    "     options.config.remote = True\n"
                    "     options.config.remote_url = 'http://localhost:8080'\n"
                    "     solver = SolverLQNS(model, options)\n\n"
                    "To skip this solver in advance, check if LQNS is available before creating the solver:\n"
                    "if SolverLQNS.isAvailable():\n"
                    "    solver = SolverLQNS(model)"
                ) from e
            else:
                raise e

    def runAnalyzer(self):
        """Run the LQNS analysis."""
        return super().runAnalyzer()

    @staticmethod
    def isAvailable():
        """Check if the LQNS solver is available."""
        try:
            return jpype.JPackage('jline').solvers.lqns.SolverLQNS.isAvailable()
        except Exception:
            return False

    def getAvgTable(self):
        """Get performance metrics table."""
        table = self.obj.getAvgTable()

        QLen = np.array(list(table.getQLen()))
        Util = np.array(list(table.getUtil()))
        RespT = np.array(list(table.getRespT()))
        ResidT = np.array(list(table.getResidT()))
        ArvR = np.array(list(table.getArvR()))
        Tput = np.array(list(table.getTput()))

        cols = ['QLen', 'Util', 'RespT', 'ResidT', 'ArvR', 'Tput']
        nodenames = list(table.getNodeNames())
        mynodenames = []
        for i in range(len(nodenames)):
            mynodenames.append(str(nodenames[i]))
        nodetypes = list(table.getNodeTypes())
        mynodetypes = []
        for i in range(len(nodetypes)):
            mynodetypes.append(str(nodetypes[i]))
        AvgTable = pd.DataFrame(np.concatenate([[QLen, Util, RespT, ResidT, ArvR, Tput]]).T, columns=cols)
        tokeep = ~(AvgTable <= 0.0).all(axis=1)
        AvgTable.insert(0, "NodeType", mynodetypes)
        AvgTable.insert(0, "Node", mynodenames)
        AvgTable = AvgTable.loc[tokeep]
        if not self._table_silent:
            print(AvgTable)

        return AvgTable

    avg_table = getAvgTable

    def avgT(self):
        """Short alias for avgTable/getAvgTable."""
        return self.getAvgTable()

    aT = avgT
    avg_t = avgT
    run_analyzer = runAnalyzer
    get_avg_table = getAvgTable

    @staticmethod
    def defaultOptions():
        """Returns default solver options for the LQNS solver."""
        java_options = jpype.JPackage('jline').solvers.lqns.SolverLQNS.defaultOptions()
        python_options = SolverOptions.__new__(SolverOptions)
        python_options.obj = java_options
        return python_options

    default_options = defaultOptions

    @staticmethod
    def supportsModel(model):
        java_solver_class = jpype.JPackage('jline').solvers.lqns.SolverLQNS
        return java_solver_class.supports(model.obj if hasattr(model, 'obj') else model)

class SolverLN(EnsembleSolver):
    """
    Layered Network (LN) solver for layered queueing network models.

    SolverLN is LINE's built-in solver for layered queueing networks (LQNs).
    It analyzes hierarchical client-server systems where processors host tasks
    that provide service through entries, and tasks can make calls to other
    tasks creating a layered service structure.

    The solver supports:
    - Multi-layered architectures
    - Synchronous and asynchronous task calls
    - Activity precedence graphs
    - Think time delays
    - Multiple job classes
    - Processor sharing and priority scheduling

    The solver uses iterative fixed-point methods to handle the interdependencies
    between layers, computing response times and throughputs for each entry and
    activity in the system.

    Args:
        model: LayeredNetwork model to solve
        lang: Language for solver execution ('python' for native, 'java' for wrapper).
              Default is 'python' which uses native Python implementation.

    Typical applications:
    - Software performance modeling
    - Distributed system analysis
    - Service-oriented architectures
    - Multi-tier web applications

    Note:
        This is LINE's native LQN solver. For compatibility with external LQN
        models, use SolverLQNS which wraps the external LQNS tool.
    """
    def __init__(self, *args, **kwargs):
        self._model = args[0]

        _initialize_jar_globals()
        options = SolverOptions(jpype.JPackage('jline').lang.constant.SolverType.LN)
        super().__init__(options, *args[1:], **kwargs)
        model = args[0]
        java_network = _get_java_network(model)
        self.obj = jpype.JPackage('jline').solvers.ln.SolverLN(java_network, self.solveropt.obj)

    def runAnalyzer(self):
        """Run the LN analysis."""
        self.obj.runAnalyzer()
        return self

    def getAvgTable(self):
        """Get comprehensive average performance metrics table."""
        table = self.obj.getAvgTable()

        # Check if solver returned null (e.g., due to unsupported features like phase-2)
        if table is None:
            return pd.DataFrame()

        QLen = np.array(list(table.getQLen()))
        Util = np.array(list(table.getUtil()))
        RespT = np.array(list(table.getRespT()))
        ResidT = np.array(list(table.getResidT()))
        ArvR = np.array(list(table.getArvR()))
        Tput = np.array(list(table.getTput()))

        cols = ['QLen', 'Util', 'RespT', 'ResidT', 'ArvR', 'Tput']
        nodenames = list(table.getNodeNames())
        mynodenames = []
        for i in range(len(nodenames)):
            mynodenames.append(str(nodenames[i]))
        nodetypes = list(table.getNodeTypes())
        mynodetypes = []
        for i in range(len(nodetypes)):
            mynodetypes.append(str(nodetypes[i]))
        AvgTable = pd.DataFrame(np.concatenate([[QLen, Util, RespT, ResidT, ArvR, Tput]]).T, columns=cols)
        tokeep = ~(AvgTable <= 0.0).all(axis=1)
        AvgTable.insert(0, "NodeType", mynodetypes)
        AvgTable.insert(0, "Node", mynodenames)
        AvgTable = AvgTable.loc[tokeep]
        if not self._table_silent:
            print(AvgTable)

        return AvgTable

    avg_table = getAvgTable

    def avgT(self):
        """Short alias for avgTable/getAvgTable."""
        return self.getAvgTable()

    aT = avgT
    avg_t = avgT

    @staticmethod
    def defaultOptions():
        """Returns default solver options for the LN solver."""
        java_options = jpype.JPackage('jline').solvers.ln.SolverLN.defaultOptions()
        python_options = SolverOptions.__new__(SolverOptions)
        python_options.obj = java_options
        return python_options

    default_options = defaultOptions

    @staticmethod
    def supportsModel(model):
        java_solver_class = jpype.JPackage('jline').solvers.ln.SolverLN
        return java_solver_class.supports(model.obj if hasattr(model, 'obj') else model)

    get_avg_table = getAvgTable

    def get_state(self):
        """
        Export current solver state for continuation.

        Returns:
            LNState: Object containing the current solution state, which can be
                used to continue iteration with a different solver via set_state().

        Example:
            >>> solver1 = SolverLN(model)
            >>> solver1.getAvgTable()
            >>> state = solver1.get_state()
            >>>
            >>> solver2 = SolverLN(model)
            >>> solver2.set_state(state)
            >>> solver2.getAvgTable()  # Continues from solver1's solution
        """
        return self.obj.getState()

    getState = get_state

    def set_state(self, state):
        """
        Import solution state for continuation.

        Initializes the solver with a previously exported state, allowing
        iteration to continue from where a previous solver left off.

        Args:
            state: LNState object from get_state()

        Example:
            >>> solver1 = SolverLN(model)
            >>> solver1.getAvgTable()
            >>> state = solver1.get_state()
            >>>
            >>> solver2 = SolverLN(model)
            >>> solver2.set_state(state)
            >>> solver2.getAvgTable()  # Continues from solver1's solution
        """
        self.obj.setState(state)

    setState = set_state

    def update_solver(self, solver_factory):
        """
        Change the solver for all layers.

        Replaces all layer solvers with new solvers created by the given factory.
        This allows switching between different solving methods (e.g., from MVA
        to JMT) while preserving the current solution state.

        Args:
            solver_factory: Factory function to create new layer solvers.
                Should be a Java SolverFactory instance.

        Example:
            >>> solver = SolverLN(model)
            >>> solver.getAvgTable()  # Fast MVA solution
            >>> # Switch to simulation for refinement
            >>> solver.update_solver(jline.solvers.ln.SolverLN.createSolverFactory(
            ...     jline.lang.constant.SolverType.JMT))
            >>> solver.getAvgTable()  # Continue with JMT
        """
        self.obj.updateSolver(solver_factory)

    updateSolver = update_solver

class SolverNC(NetworkSolver):
    """
    Normalizing Constant (NC) solver for closed queueing networks.

    SolverNC implements the convolution algorithm for closed product-form
    queueing networks. It computes the normalizing constant and uses it to
    derive exact performance measures including queue lengths, throughputs,
    and utilizations.

    The solver supports:
    - Closed queueing networks with fixed populations
    - Multiple job classes
    - Load-dependent service rates
    - Product-form scheduling disciplines (FCFS, PS, INF, LCFS-PR)
    - State probability computations

    The convolution algorithm works by recursively computing normalizing
    constants for increasing population levels, making it particularly
    efficient for networks with many stations but moderate populations.

    Args:
        model: Network model to solve
        lang: Language for solver execution ('python' for native, 'java' for wrapper).
              Default is 'python' which uses native Python implementation.

    Advantages:
    - More numerically stable than MVA for large populations
    - Provides access to state probabilities
    - Exact results for product-form networks

    Note:
        Limited to closed networks satisfying product-form conditions (BCMP theorem).
        For open or mixed networks, use SolverMVA or other appropriate solvers.
    """
    def __init__(self, *args, **kwargs):
        self._model = args[0]

        _initialize_jar_globals()
        if len(args) > 1 and hasattr(args[1], 'obj'):
            options = args[1]
            super().__init__(options, *args[2:], **kwargs)
        else:
            options = SolverOptions(jpype.JPackage('jline').lang.constant.SolverType.NC)
            super().__init__(options, *args[1:], **kwargs)
        model = args[0]
        java_network = _get_java_network(model)
        self.obj = jpype.JPackage('jline').solvers.nc.SolverNC(java_network, self.solveropt.obj)

    def runAnalyzer(self):
        """Run the NC analysis."""
        self.obj.runAnalyzer()
        return self

    def getAvgTable(self):
        """Get comprehensive average performance metrics table."""
        return super().getAvgTable()

    avg_table = getAvgTable
    getAvgT = getAvgTable
    avgT = getAvgTable
    aT = getAvgTable

    def getProbSysAggr(self):
        try:
            java_result = self.obj.probSysAggr()
            return ProbabilityResult(java_result)
        except Exception as e:
            if not self._verbose_silent:
                print(f"NC getProbSysAggr failed: {e}")
            return None

    def getProbMarg(self, node):
        """
        Probability distribution for TOTAL queue length at a station.
        Returns P(n total jobs) for n=0,1,...,N, summing over all class combinations.

        Compare with get_prob_aggr: returns probability of a specific per-class
        distribution, e.g., P(2 class-1, 1 class-2) as a scalar.

        Args:
            node: Node object or index

        Returns:
            tuple: (prob, log_prob) where prob[n] = P(n total jobs at node)
        """
        try:
            if hasattr(node, 'obj'):
                java_node = node.obj
            else:
                # Convert index to node object
                java_node = self.obj.getStruct().nodes[int(node) - 1]

            # Call JAR method - returns Pair<Matrix, Matrix> of (prob, logprob)
            java_result = self.obj.getProbMarg(java_node)

            # Extract the two matrices
            prob = jlineMatrixToArray(java_result.getKey())
            log_prob = jlineMatrixToArray(java_result.getValue())

            return (prob, log_prob)
        except Exception as e:
            if not self._verbose_silent:
                print(f"NC getProbMarg failed: {e}")
            return (None, None)

    get_prob_marg = getProbMarg

    def getProbSys(self):
        """
        Get system-wide joint probability.

        Returns:
            ProbabilityResult: System-wide joint probability value with metadata
        """
        try:
            java_result = self.obj.getProbSys()
            return ProbabilityResult(java_result)
        except Exception as e:
            if not self._verbose_silent:
                print(f"NC getProbSys failed: {e}")
            return None

    get_prob_sys = getProbSys

    @staticmethod
    def defaultOptions():
        """Returns default solver options for the NC solver."""
        java_options = jpype.JPackage('jline').solvers.nc.SolverNC.defaultOptions()
        python_options = SolverOptions.__new__(SolverOptions)
        python_options.obj = java_options
        return python_options

    default_options = defaultOptions

    @staticmethod
    def supportsModel(model):
        java_solver_class = jpype.JPackage('jline').solvers.nc.SolverNC
        return java_solver_class.supports(model.obj if hasattr(model, 'obj') else model)

    get_prob_sys_aggr = getProbSysAggr
    prob_sys_aggr = getProbSysAggr


class SolverSSA(NetworkSolver):
    """
    State-Space Analysis (SSA) solver using stochastic simulation.

    SolverSSA implements discrete-event simulation for queueing networks using
    a state-space approach. It can handle arbitrary network topologies and
    distributions, providing statistical estimates of performance metrics
    along with sample paths and event traces.

    The solver supports:
    - Open, closed, and mixed networks
    - Arbitrary service time distributions
    - General arrival processes
    - Fork-join synchronization
    - State-dependent routing and service
    - Any scheduling discipline
    - Event logging and trace analysis

    Features:
    - Sample path generation
    - State probability estimation
    - Transient and steady-state analysis
    - Configurable sampling strategies
    - Event-by-event traces

    Args:
        model: Network model to solve
        lang: Language for solver execution ('python' for native, 'java' for wrapper).
              Default is 'python' which uses native Python implementation.

    The SSA solver uses the Stochastic Simulation Algorithm (Gillespie's algorithm)
    to generate exact sample paths of the stochastic process, from which
    performance metrics are estimated.

    Note:
        Results are statistical estimates with sampling variability. For exact
        analytical results on product-form networks, consider SolverMVA or
        SolverNC. Use larger sample sizes for higher accuracy.
    """
    def __init__(self, *args, **kwargs):
        self._samples = kwargs.pop('samples', 10000)
        self._seed = kwargs.pop('seed', 0)
        self._model = args[0]

        _initialize_jar_globals()
        if len(args) > 1 and hasattr(args[1], 'obj'):
            options = args[1]
            super().__init__(options, *args[2:], **kwargs)
        else:
            options = SolverOptions(jpype.JPackage('jline').lang.constant.SolverType.SSA)
            super().__init__(options, *args[1:], **kwargs)

        # Set seed and samples on Java solver options (they were popped from kwargs)
        self.solveropt.obj.seed(int(self._seed))
        self.solveropt.obj.samples(int(self._samples))

        model = args[0]
        java_network = _get_java_network(model)
        self.obj = jpype.JPackage('jline').solvers.ssa.SolverSSA(java_network, self.solveropt.obj)

    def runAnalyzer(self):
        """Run the SSA analysis."""
        self.obj.runAnalyzer()
        return self

    def getAvgTable(self):
        """Get comprehensive average performance metrics table."""
        return super().getAvgTable()

    avg_table = getAvgTable
    getAvgT = getAvgTable
    avgT = getAvgTable
    aT = getAvgTable

    def sample(self, node, numSamples=10000, markActivePassive=False):
        """Sample from a specific node"""
        try:
            if numSamples is None:
                java_result = self.obj.sample(node.obj)
            else:
                java_result = self.obj.sample(node.obj, numSamples, markActivePassive)
            return SampleResult(java_result) if java_result is not None else None
        except Exception as e:
            if not self._verbose_silent:
                print(f"SSA sampling failed: {e}")
            return None

    def sampleAggr(self, node, numSamples=10000, markActivePassive=False):
        """Sample from a specific node using aggregated states"""
        try:
            if numSamples is None:
                java_result = self.obj.sampleAggr(node.obj)
            else:
                java_result = self.obj.sampleAggr(node.obj, numSamples, markActivePassive)
            return SampleResult(java_result) if java_result is not None else None
        except Exception as e:
            if not self._verbose_silent:
                print(f"SSA aggregated sampling failed: {e}")
            return None

    def sampleSys(self, numSamples=10000):
        """Sample system-wide"""
        try:
            if numSamples is None:
                java_result = self.obj.sampleSys()
            else:
                java_result = self.obj.sampleSys(numSamples)
            return SampleResult(java_result) if java_result is not None else None
        except Exception as e:
            if not self._verbose_silent:
                print(f"SSA system sampling failed: {e}")
            return None

    def sampleSysAggr(self, numSamples=10000):
        """Sample system-wide using aggregated states"""
        try:
            if numSamples is None:
                java_result = self.obj.sampleSysAggr()
            else:
                java_result = self.obj.sampleSysAggr(numSamples)
            return SampleResult(java_result) if java_result is not None else None
        except Exception as e:
            if not self._verbose_silent:
                print(f"SSA system aggregated sampling failed: {e}")
            return None

    @staticmethod
    def defaultOptions():
        """Returns default solver options for the SSA solver."""
        java_options = jpype.JPackage('jline').solvers.ssa.SolverSSA.defaultOptions()
        python_options = SolverOptions.__new__(SolverOptions)
        python_options.obj = java_options
        return python_options

    default_options = defaultOptions

    @staticmethod
    def supportsModel(model):
        java_solver_class = jpype.JPackage('jline').solvers.ssa.SolverSSA
        return java_solver_class.supports(model.obj if hasattr(model, 'obj') else model)

    def getTranProb(self, node):
        """
        Get transient state probabilities for a specific node.

        Args:
            node: The node to get probabilities for.

        Returns:
            ProbabilityResult: Transient state probability distribution.
        """
        try:
            java_result = self.obj.getTranProb(node.obj if hasattr(node, 'obj') else node)
            return ProbabilityResult(java_result)
        except Exception as e:
            if not self._verbose_silent:
                print(f"SSA getTranProb failed: {e}")
            return None

    def getTranProbAggr(self, node):
        """
        Get aggregated transient state probabilities for a node.

        Args:
            node: The node to get probabilities for.

        Returns:
            ProbabilityResult: Aggregated transient state probabilities.
        """
        try:
            java_result = self.obj.getTranProbAggr(node.obj if hasattr(node, 'obj') else node)
            return ProbabilityResult(java_result)
        except Exception as e:
            if not self._verbose_silent:
                print(f"SSA getTranProbAggr failed: {e}")
            return None

    def getTranProbSys(self):
        """
        Get system-wide transient state probabilities.

        Returns:
            ProbabilityResult: System-wide transient probabilities.
        """
        try:
            java_result = self.obj.getTranProbSys()
            return ProbabilityResult(java_result)
        except Exception as e:
            if not self._verbose_silent:
                print(f"SSA getTranProbSys failed: {e}")
            return None

    def getTranProbSysAggr(self):
        """
        Get aggregated system-wide transient state probabilities.

        Returns:
            ProbabilityResult: Aggregated system-wide transient probabilities.
        """
        try:
            java_result = self.obj.getTranProbSysAggr()
            return ProbabilityResult(java_result)
        except Exception as e:
            if not self._verbose_silent:
                print(f"SSA getTranProbSysAggr failed: {e}")
            return None

    sample_sys_aggr = sampleSysAggr
    sample_sys = sampleSys
    tran_prob = getTranProb
    get_tran_prob = getTranProb
    tran_prob_aggr = getTranProbAggr
    get_tran_prob_aggr = getTranProbAggr
    tran_prob_sys = getTranProbSys
    get_tran_prob_sys = getTranProbSys
    tran_prob_sys_aggr = getTranProbSysAggr
    get_tran_prob_sys_aggr = getTranProbSysAggr


class SolverDES(NetworkSolver):
    """
    Discrete Event Simulation (DES) solver using SimPy.

    SolverDES implements discrete-event simulation for queueing networks
    using SimPy (native Python) or SSJ (Java backend). It supports open
    and closed networks with various service distributions, scheduling
    strategies, and advanced node types.

    The solver supports:
    - Node types: Queue, Delay, Fork, Join, Router, ClassSwitch, Logger
    - Distributions: Exp, Erlang, HyperExp, PH, APH, Coxian, Replayer, Immediate, Disabled
    - Scheduling: FCFS, PS, DPS, GPS, LCFS variants, SIRO, SJF, LJF, SEPT, LEPT, HOL (priority)
    - Routing: Probabilistic, Random, Round-Robin, Weighted RR, K-Choices
    - Load-dependent service rates
    - Fork-join parallelism
    - Class switching and routing
    - Finite buffer capacity and capacity regions
    - Open and closed classes
    - Multiclass workloads

    Features:
    - High-performance discrete-event simulation
    - Statistical confidence estimation via MSER-5 and OBM
    - Configurable sample size for accuracy control
    - Reproducible results via seed control
    - Support for complex queueing networks with advanced features

    Args:
        model: Network model to solve
        lang: Language for solver execution ('python' for native SimPy, 'java' for SSJ).
              Default is 'python' which uses native Python implementation.

    Note:
        Results are statistical estimates with sampling variability. For exact
        analytical results on product-form networks, consider SolverMVA.
        The DES solver supports a wide range of network types and features.
    """
    def __init__(self, *args, **kwargs):
        self._model = args[0]
        self._verbose_silent = kwargs.get('verbose', VerboseLevel.SILENT) == VerboseLevel.SILENT

        _initialize_jar_globals()
        if len(args) > 1 and hasattr(args[1], 'obj'):
            options = args[1]
            super().__init__(options, *args[2:], **kwargs)
        else:
            try:
                options = SolverDES.defaultOptions()
            except:
                options = SolverOptions(jpype.JPackage('jline').lang.constant.SolverType.DES)
            super().__init__(options, *args[1:], **kwargs)
        model = args[0]
        java_network = _get_java_network(model)
        self.obj = jpype.JPackage('jline').solvers.des.SolverDES(java_network, self.solveropt.obj)

    def runAnalyzer(self):
        """Run the DES analysis."""
        # Java DES solver auto-runs simulation when getAvgTable() is called
        return self

    run_analyzer = runAnalyzer

    def getAvgTable(self):
        """Get comprehensive average performance metrics table."""
        return super().getAvgTable()

    get_avg_table = getAvgTable
    avg_table = getAvgTable

    @staticmethod
    def defaultOptions():
        """Returns default solver options for the DES solver."""
        java_options = jpype.JPackage('jline').solvers.des.SolverDES.defaultOptions()
        python_options = SolverOptions.__new__(SolverOptions)
        python_options.obj = java_options
        return python_options

    default_options = defaultOptions

    @staticmethod
    def supportsModel(model):
        """Check if the model is supported by the DES solver."""
        java_solver_class = jpype.JPackage('jline').solvers.des.SolverDES
        return java_solver_class.supports(model.obj if hasattr(model, 'obj') else model)

    supports_model = supportsModel


class SolverAuto(NetworkSolver):
    """
    Automatic solver selection for queueing network models.

    SolverAuto analyzes the network model and automatically selects the most
    appropriate solver based on model characteristics, feature support, and
    performance considerations. This provides a convenient high-level interface
    when the optimal solver choice is not immediately obvious.

    Selection criteria include:
    - Network type (open, closed, mixed)
    - Distribution types (exponential, phase-type, general)
    - Scheduling disciplines
    - Network size and complexity
    - Desired accuracy vs computation time tradeoffs

    Available selection methods:
    - default: Rule-based selection using model features
    - heur: Heuristic-based selection
    - ai: AI-assisted selection (if available)
    - nn: Neural network-based selection (if available)

    The solver maintains a list of candidate solvers and selects the best one
    according to the specified selection method. You can query which solver
    was selected and see alternative candidates.

    Args:
        model: Network model to analyze
        method: Solver selection method ('default', 'heur', 'ai')
        **kwargs: Additional options passed to the selected solver

    Usage:
        solver = SolverAuto(model)
        results = solver.avg_table()
        print(f"Selected solver: {solver.get_selected_solver_name()}")
        print(f"Candidates: {solver.get_candidate_solver_names()}")

    Note:
        Can be forced to use a specific solver via set_forced_solver().
        The LINE class is an alias for SolverAuto.
    """
    def __init__(self, *args, method='default', **kwargs):
        self._model = args[0] if args else kwargs.get('model')
        self._method = method
        self._kwargs = kwargs

        _initialize_jar_globals()
        options = SolverOptions(jpype.JPackage('jline').lang.constant.SolverType.MVA)
        super().__init__(options, *args[1:], **kwargs)
        model = args[0]

        # Set selection method on options before passing to Java
        self.solveropt.method(method)

        SolverAutoClass = jpype.JClass('jline.solvers.auto.SolverAUTO')
        java_network = _get_java_network(model)
        # Pass options object to Java SolverAuto constructor
        self.obj = SolverAutoClass(java_network, self.solveropt.obj)

    def runAnalyzer(self):
        """Run the analysis with the automatically selected solver."""
        # Java implementation auto-runs when getting results
        return self

    def getName(self):
        """Get the name of this solver."""
        return self.obj.getName()

    def getAvgTable(self):
        """Get comprehensive average performance metrics table."""
        return super().getAvgTable()

    def getSelectedSolverName(self):
        """Get the name of the solver that was automatically selected"""
        return str(self.obj.getSelectedSolverName())

    def getCandidateSolverNames(self):
        """Get list of candidate solver names that could be used for this model"""
        names = self.obj.getCandidateSolverNames()
        return [str(name) for name in names]

    def setSelectionMethod(self, method):
        """Set the solver selection method: 'default', 'heur', 'ai', or 'nn'"""
        self.solveropt.obj.method = method

    def setForcedSolver(self, solver_name):
        """Force a specific solver: 'mva', 'nc', 'mam', 'fluid', 'jmt', 'ssa', 'ctmc'"""
        auto_opts = AutoOptions()
        auto_opts.setForcedSolver(solver_name)

    @staticmethod
    def supportsModel(model):
        SolverAutoClass = jpype.JClass('jline.solvers.auto.SolverAUTO')
        java_network = _get_java_network(model)
        return SolverAutoClass.supports(java_network)

    # ========== Basic Metric Methods ==========

    def getAvg(self):
        """Get all average metrics."""
        return super().getAvg()

    def getAvgQLen(self):
        """Get average queue lengths at steady-state."""
        return super().getAvgQLen()

    def getAvgUtil(self):
        """Get average utilizations at steady-state."""
        return super().getAvgUtil()

    def getAvgRespT(self):
        """Get average response times at steady-state."""
        return super().getAvgRespT()

    def getAvgResidT(self):
        """Get average residence times at steady-state."""
        return super().getAvgResidT()

    def getAvgTput(self):
        """Get average throughputs at steady-state."""
        return super().getAvgTput()

    def getAvgArvR(self):
        """Get average arrival rates at steady-state."""
        return super().getAvgArvR()

    def getAvgWaitT(self):
        """Get average waiting times."""
        return super().getAvgWaitT()

    # ========== System Metric Methods ==========

    def getAvgSysRespT(self):
        """Get system average response time."""
        return super().getAvgSysRespT()

    def getAvgSysTput(self):
        """Get system average throughput."""
        return super().getAvgSysTput()

    # ========== Chain Methods ==========

    def getAvgArvRChain(self):
        """Get average arrival rates by chain."""
        return super().getAvgArvRChain()

    def getAvgQLenChain(self):
        """Get average queue lengths by chain."""
        return super().getAvgQLenChain()

    def getAvgResidTChain(self):
        """Get average residence times by chain."""
        return super().getAvgResidTChain()

    def getAvgRespTChain(self):
        """Get average response times by chain."""
        return super().getAvgRespTChain()

    def getAvgTputChain(self):
        """Get average throughputs by chain."""
        return super().getAvgTputChain()

    def getAvgUtilChain(self):
        """Get average utilizations by chain."""
        return super().getAvgUtilChain()

    # ========== Node Chain Methods ==========

    def getAvgNodeArvRChain(self):
        """Get average node arrival rates by chain."""
        return super().getAvgNodeArvRChain()

    def getAvgNodeQLenChain(self):
        """Get average node queue lengths by chain."""
        return super().getAvgNodeQLenChain()

    def getAvgNodeResidTChain(self):
        """Get average node residence times by chain."""
        return super().getAvgNodeResidTChain()

    def getAvgNodeRespTChain(self):
        """Get average node response times by chain."""
        return super().getAvgNodeRespTChain()

    def getAvgNodeTputChain(self):
        """Get average node throughputs by chain."""
        return super().getAvgNodeTputChain()

    def getAvgNodeUtilChain(self):
        """Get average node utilizations by chain."""
        return super().getAvgNodeUtilChain()

    # ========== Table Methods ==========

    def getAvgNodeTable(self):
        """Get average metrics table by node."""
        return super().getAvgNodeTable()

    def getAvgChainTable(self):
        """Get average metrics table by chain."""
        return super().getAvgChainTable()

    def getAvgSysTable(self):
        """Get system average metrics table."""
        return super().getAvgSysTable()

    def getAvgNodeChainTable(self):
        """Get average metrics table by node and chain."""
        return super().getAvgNodeChainTable()

    # ========== Distribution Methods ==========

    def getCdfRespT(self):
        """Get CDF of response times at steady-state."""
        return super().getCdfRespT()

    # ========== Utility Methods ==========

    def hasResults(self):
        """Check if solver has results available."""
        return super().hasResults()

    def reset(self):
        """Reset solver state."""
        super().reset()

    run_analyzer = runAnalyzer
    get_name = getName
    get_selected_solver_name = getSelectedSolverName
    get_candidate_solver_names = getCandidateSolverNames
    set_selection_method = setSelectionMethod
    set_forced_solver = setForcedSolver

    # Snake_case aliases for basic metrics
    get_avg = getAvg
    get_avg_qlen = getAvgQLen
    get_avg_util = getAvgUtil
    get_avg_resp_t = getAvgRespT
    get_avg_resid_t = getAvgResidT
    get_avg_tput = getAvgTput
    get_avg_arv_r = getAvgArvR
    get_avg_wait_t = getAvgWaitT

    # Snake_case aliases for system metrics
    get_avg_sys_resp_t = getAvgSysRespT
    get_avg_sys_tput = getAvgSysTput

    # Snake_case aliases for chain metrics
    get_avg_arv_r_chain = getAvgArvRChain
    get_avg_qlen_chain = getAvgQLenChain
    get_avg_resid_t_chain = getAvgResidTChain
    get_avg_resp_t_chain = getAvgRespTChain
    get_avg_tput_chain = getAvgTputChain
    get_avg_util_chain = getAvgUtilChain

    # Snake_case aliases for node chain metrics
    get_avg_node_arv_r_chain = getAvgNodeArvRChain
    get_avg_node_qlen_chain = getAvgNodeQLenChain
    get_avg_node_resid_t_chain = getAvgNodeResidTChain
    get_avg_node_resp_t_chain = getAvgNodeRespTChain
    get_avg_node_tput_chain = getAvgNodeTputChain
    get_avg_node_util_chain = getAvgNodeUtilChain

    # Snake_case aliases for table methods
    get_avg_table = getAvgTable
    get_avg_node_table = getAvgNodeTable
    get_avg_chain_table = getAvgChainTable
    get_avg_sys_table = getAvgSysTable
    get_avg_node_chain_table = getAvgNodeChainTable

    # Snake_case aliases for distribution methods
    get_cdf_resp_t = getCdfRespT

    # Snake_case aliases for utility methods
    has_results = hasResults


class LINE(SolverAuto):
    """Alias for SolverAuto - automatic solver selection.

    Also provides static methods for loading models from files.
    """
    def __init__(self, *args, **kwargs):
        super().__init__(*args, **kwargs)

    @staticmethod
    def load(filename, verbose=False):
        """Load a Network model from a file.

        Detects the file format and loads the appropriate model type.
        Currently supports:
        - .jsim, .jsimg, .jsimw: JSIM/JMT format (queueing networks)
        - .xml, .lqn, .lqnx: Layered Queueing Network (LQN) models
        - .mat: MATLAB format
        - .pkl: Python pickle format

        Args:
            filename (str): The path to the file to load
            verbose (bool): Whether to print verbose loading information

        Returns:
            Network or LayeredNetwork: The loaded model

        Raises:
            ValueError: If the file format is unsupported
            RuntimeError: If loading fails

        Examples:
            >>> model = LINE.load('mymodel.lqn')
            >>> model = LINE.load('network.xml', verbose=True)
            >>> model = LINE.load('network.jsimg')
        """
        import os
        from .layered import LayeredNetwork
        from .lang import M2M

        # Get file extension
        _, ext = os.path.splitext(filename.lower())

        if ext in ['.xml', '.lqn', '.lqnx']:
            # Load as LayeredNetwork (LQN model)
            try:
                lqn_model = LayeredNetwork.load(filename, verbose)
                if verbose:
                    print(f"Loaded LQN model from: {filename}")
                # Return the underlying Network model
                return lqn_model.get_model()
            except Exception as e:
                raise RuntimeError(f"Failed to load LQN model from: {filename}") from e

        elif ext in ['.jsim', '.jsimg', '.jsimw']:
            # JSIM/JMT format - convert to LINE Network using M2M
            try:
                m2m = M2M()
                jsim_model = m2m.JSIM2LINE(filename)
                if verbose:
                    print(f"Loaded JSIM model from: {filename}")
                return jsim_model
            except Exception as e:
                raise RuntimeError(f"Failed to load JSIM model from: {filename}") from e

        elif ext == '.jmva':
            # JMVA format - would need JMVA2LINE implementation
            raise NotImplementedError(
                "JMVA format loading not yet implemented in Python. "
                "Use MATLAB LINE.load() for JMVA files."
            )

        elif ext == '.pkl':
            # Python pickle format for saved models
            import pickle
            try:
                with open(filename, 'rb') as f:
                    model = pickle.load(f)
                if verbose:
                    print(f"Loaded model from pickle file: {filename}")
                return model
            except Exception as e:
                raise RuntimeError(f"Failed to load pickled model from: {filename}") from e

        elif ext == '.mat':
            # MATLAB .mat format - convert to LINE Network using M2M
            try:
                m2m = M2M()
                mat_model = m2m.MAT2LINE(filename)
                if verbose:
                    print(f"Loaded MATLAB model from: {filename}")
                return mat_model
            except Exception as e:
                raise RuntimeError(f"Failed to load MATLAB model from: {filename}") from e

        else:
            raise ValueError(
                f"Unsupported file format: {ext}. "
                f"Supported formats: .jsim, .jsimg, .jsimw (JSIM models), "
                f".xml, .lqn, .lqnx (LQN models), .mat (MATLAB), .pkl (pickled models)"
            )


# Alias for SolverAuto - matches MATLAB naming convention
AUTO = SolverAuto


class SolverConfig():
    """Python wrapper for SolverOptions.Config class"""
    def __init__(self, java_config):
        self.obj = java_config

    @property
    def multiserver(self):
        return self.obj.multiserver

    @multiserver.setter
    def multiserver(self, value):
        self.obj.multiserver = value

    @property
    def highvar(self):
        return self.obj.highvar

    @highvar.setter
    def highvar(self, value):
        self.obj.highvar = value

    @property
    def np_priority(self):
        return self.obj.npPriority

    @np_priority.setter
    def np_priority(self, value):
        self.obj.npPriority = value

    @property
    def fork_join(self):
        return self.obj.fork_join

    @fork_join.setter
    def fork_join(self, value):
        self.obj.fork_join = value

    @property
    def eventcache(self):
        return self.obj.eventcache

    @eventcache.setter
    def eventcache(self, value):
        self.obj.eventcache = value

    @property
    def nonmkv(self):
        return self.obj.nonmkv

    @nonmkv.setter
    def nonmkv(self, value):
        self.obj.nonmkv = value

    @property
    def nonmkvorder(self):
        return self.obj.nonmkvorder

    @nonmkvorder.setter
    def nonmkvorder(self, value):
        self.obj.nonmkvorder = value

class SolverOptions():
    def __init__(self, solvertype):
        self.obj = jpype.JPackage('jline').solvers.SolverOptions(solvertype)
        self._config = None

    def method(self, value):
        self.obj.method(value)

    def samples(self, value):
        self.obj.samples(value)

    def seed(self, value):
        self.obj.seed(value)

    def verbose(self, value):
        if hasattr(value, 'value'):
            self.obj.verbose(value.value)
        else:
            self.obj.verbose(value)

    @property
    def config(self):
        """Access to advanced configuration options"""
        if not hasattr(self, '_config') or self._config is None:
            self._config = SolverConfig(self.obj.config)
        return self._config

    @property
    def iter_max(self):
        return self.obj.iter_max

    @iter_max.setter
    def iter_max(self, value):
        self.obj.iter_max = value

    @property
    def iter_tol(self):
        return self.obj.iter_tol

    @iter_tol.setter
    def iter_tol(self, value):
        self.obj.iter_tol = value

    @property
    def tol(self):
        return self.obj.tol

    @tol.setter
    def tol(self, value):
        self.obj.tol = value

    def __setitem__(self, key, value):
        """Support dictionary-style assignment for solver options."""
        if key == 'keep':
            self.obj.keep(bool(value))
        elif key == 'verbose':
            if hasattr(value, 'value'):
                self.obj.verbose(value.value)
            else:
                self.obj.verbose(value)
        elif key == 'cutoff':
            if hasattr(value, '__iter__') and not isinstance(value, str):
                from line_solver import jlineMatrixFromArray
                self.obj.cutoff(jlineMatrixFromArray(value))
            else:
                self.obj.cutoff(value)
        elif key == 'seed':
            self.obj.seed(int(value))
        elif key == 'samples':
            self.obj.samples(int(value))
        elif key == 'method':
            self.obj.method(str(value))
        elif key == 'force':
            self.obj.force(bool(value))
        elif key in ['iter_max', 'iter_tol', 'tol']:
            setattr(self.obj, key, value)
        elif key == 'timespan':
            if hasattr(value, '__iter__'):
                from jpype import JArray, JDouble
                self.obj.timespan = JArray(JDouble)(value)
            else:
                self.obj.timespan = value
        elif key == 'timestep':
            from jpype import JDouble
            self.obj.timestep = JDouble(value) if value is not None else None
        elif hasattr(self.obj, key):
            try:
                setattr(self.obj, key, value)
            except AttributeError:
                raise KeyError(f"Solver option '{key}' is not settable")
        else:
            raise KeyError(f"Unknown solver option: '{key}'")

    def __getitem__(self, key):
        """Support dictionary-style access for solver options."""
        if hasattr(self.obj, key):
            value = getattr(self.obj, key)

            if hasattr(value, '__call__'):
                raise KeyError(f"Cannot read solver option '{key}' - use property-style access instead (options.{key})")

            if hasattr(value, 'value'):
                return value.value
            elif str(type(value)).startswith('<java'):
                try:
                    if hasattr(value, 'toString'):
                        return value.toString()
                    else:
                        return str(value)
                except:
                    return value
            else:
                return value
        else:
            raise KeyError(f"Unknown solver option: '{key}'")

class CTMCOptions():
    """
    Configuration options for CTMC (Continuous Time Markov Chain) solver.
    
    Provides access to solver-specific options for CTMC analysis including
    numerical tolerances, iteration limits, and analysis methods.
    """
    def __init__(self):
        self.obj = jpype.JPackage('jline').solvers.CTMCOptions()

class EnvOptions():
    """
    Configuration options for Environment solver.
    
    Provides access to solver-specific options for ensemble analysis
    and environmental modeling scenarios.
    """
    def __init__(self):
        self.obj = jpype.JPackage('jline').solvers.EnvOptions()

class FluidOptions():
    """
    Configuration options for Fluid solver.
    
    Provides access to solver-specific options for fluid approximation
    analysis of queueing networks.
    """
    def __init__(self):
        self.obj = jpype.JPackage('jline').solvers.FluidOptions()

class JMTOptions():
    """
    Configuration options for JMT (Java Modelling Tools) solver.
    
    Provides access to solver-specific options for discrete-event
    simulation using the JMT simulation engine.
    """
    def __init__(self):
        self.obj = jpype.JPackage('jline').solvers.JMTOptions()

class LNOptions():
    """
    Configuration options for LN (Layered Network) solver.
    
    Provides access to solver-specific options for layered queueing
    network analysis.
    """
    def __init__(self):
        self.obj = jpype.JPackage('jline').solvers.LNOptions()

class LQNSOptions():
    """
    Configuration options for LQNS (Layered Queueing Network Solver).
    
    Provides access to solver-specific options for layered queueing
    network analysis using the LQNS solver engine.
    """
    def __init__(self):
        self.obj = jpype.JPackage('jline').solvers.LQNSOptions()

class MAMOptions():
    """
    Configuration options for MAM (Matrix Analytic Methods) solver.
    
    Provides access to solver-specific options for matrix-analytic
    methods including MAP/PH analysis and QBD processes.
    """
    def __init__(self):
        self.obj = jpype.JPackage('jline').solvers.MAMOptions()

class MVAOptions():
    """
    Configuration options for MVA (Mean Value Analysis) solver.
    
    Provides access to solver-specific options for mean value analysis
    of product-form queueing networks.
    """
    def __init__(self):
        self.obj = jpype.JPackage('jline').solvers.MVAOptions()

class NCOptions():
    """
    Configuration options for NC (Normalizing Constant) solver.
    
    Provides access to solver-specific options for normalizing constant
    algorithms in product-form queueing networks.
    """
    def __init__(self):
        self.obj = jpype.JPackage('jline').solvers.NCOptions()

class QNSOptions():
    """
    Configuration options for QNS (Queueing Network Solver).
    
    Provides access to solver-specific options for general
    queueing network analysis.
    """
    def __init__(self):
        self.obj = jpype.JPackage('jline').solvers.qns.SolverQNS.defaultOptions()

class SSAOptions():
    """
    Configuration options for SSA (Stochastic State-space Analysis) solver.
    
    Provides access to solver-specific options for state-space analysis
    and discrete-event simulation methods.
    """
    def __init__(self):
        self.obj = jpype.JPackage('jline').solvers.SSAOptions()

class AutoOptions():
    def __init__(self):
        self.obj = jpype.JPackage('jline').solvers.auto.AutoOptions()

    def setSelectionMethod(self, method):
        """Set the solver selection method: 'default', 'heur', 'ai', or 'nn'"""
        self.obj.selectionMethod = method

    def setForcedSolver(self, solver_name):
        """Force a specific solver: 'mva', 'nc', 'mam', 'fluid', 'jmt', 'ssa', 'ctmc'"""
        self.obj.forceSolver = solver_name

    set_selection_method = setSelectionMethod
    set_forced_solver = setForcedSolver


# Solver Aliases - Short names for convenience
MVA = SolverMVA
"""Alias for SolverMVA (Mean Value Analysis solver)."""

MAM = SolverMAM
"""Alias for SolverMAM (Matrix Analytic Methods solver)."""

NC = SolverNC
"""Alias for SolverNC (Normalizing Constant solver)."""

CTMC = SolverCTMC
"""Alias for SolverCTMC (Continuous Time Markov Chain solver)."""

SSA = SolverSSA
"""Alias for SolverSSA (Stochastic State-space Analysis solver)."""

JMT = SolverJMT
"""Alias for SolverJMT (Java Modelling Tools solver)."""

DES = SolverDES
"""Alias for SolverDES (Discrete Event Simulation solver)."""

Fluid = SolverFluid
"""Alias for SolverFluid (Fluid/Mean-Field Approximation solver)."""

FLD = SolverFluid
"""Alias for SolverFluid (Fluid/Mean-Field Approximation solver)."""

LN = SolverLN
"""Alias for SolverLN (Layered Network solver)."""

LQNS = SolverLQNS
"""Alias for SolverLQNS (Layered Queueing Network Solver)."""

QNS = SolverQNS
"""Alias for SolverQNS (Queueing Network Solver)."""

Auto = SolverAuto
"""Alias for SolverAuto (Automatic solver selection)."""

AUTO = SolverAuto
"""Alias for SolverAuto (Automatic solver selection)."""

SolverAUTO = SolverAuto
"""Alias for SolverAuto (Automatic solver selection)."""

ENV = SolverENV
"""Alias for SolverENV (Ensemble environment solver)."""

Posterior = None  # Will be assigned after SolverPosterior class definition
"""Alias for SolverPosterior (Bayesian posterior analysis solver)."""


class SolverPosterior(EnsembleSolver):
    """
    Posterior solver for Bayesian-style parameter uncertainty analysis.

    SolverPosterior handles Prior distributions by expanding the model into a family
    of networks, one for each alternative in the Prior. Results are aggregated using
    prior-weighted expectations.

    The solver:
    - Detects Prior distributions in service/arrival processes
    - Creates a copy of the model for each alternative
    - Runs the wrapped solver on each model
    - Aggregates results using prior probabilities as weights

    This is useful for:
    - Sensitivity analysis with parameter uncertainty
    - Bayesian inference on model parameters
    - Computing expected performance under parameter distributions

    Args:
        model: Network model containing Prior distributions
        solver_factory: Function that creates a solver for a given model
        options: Optional solver options

    Examples:
        >>> from line_solver import Network, Queue, Source, Sink, OpenClass
        >>> from line_solver.distributions import Exp, Prior
        >>> from line_solver import SolverPosterior, SolverMVA
        >>>
        >>> model = Network('UncertainMM1')
        >>> source = Source(model, 'Source')
        >>> queue = Queue(model, 'Queue', SchedStrategy.FCFS)
        >>> sink = Sink(model, 'Sink')
        >>> job_class = OpenClass(model, 'Jobs')
        >>> source.setArrival(job_class, Exp(0.5))
        >>> prior = Prior([Exp(1.0), Exp(2.0)], [0.6, 0.4])
        >>> queue.setService(job_class, prior)
        >>> model.link(model.serialRouting(source, queue, sink))
        >>>
        >>> solver = SolverPosterior(model, lambda m: SolverMVA(m))
        >>> avg_table = solver.getAvgTable()  # Prior-weighted expectations
    """

    def __init__(self, model, solver_factory, options=None):
        """
        Create a Posterior solver.

        Args:
            model: Network model containing Prior distributions
            solver_factory: Callable that takes a Network and returns a NetworkSolver
            options: Optional SolverOptions
        """
        # Create default options if not provided
        if options is None:
            options = SolverOptions(jpype.JPackage('jline').lang.constant.SolverType.POSTERIOR)
        super().__init__(options)

        # Store references
        self.model = model
        self.solver_factory = solver_factory

        # Create Java solver factory as a lambda
        # We need to wrap the Python solver factory for Java
        @jpype.JImplements('jline.solvers.posterior.SolverPosterior$SolverFactory')
        class PythonSolverFactory:
            def __init__(self, py_factory):
                self._py_factory = py_factory

            @jpype.JOverride
            def create(self, java_model):
                # Wrap Java model in Python Network
                from line_solver import Network
                py_model = Network.__new__(Network)
                py_model.obj = java_model
                # Create solver using Python factory
                py_solver = self._py_factory(py_model)
                # Return the underlying Java solver
                return py_solver.obj

        java_factory = PythonSolverFactory(solver_factory)

        # Create Java Posterior solver
        self.obj = jpype.JPackage('jline').solvers.posterior.SolverPosterior(
            model.obj, java_factory
        )

    def runAnalyzer(self):
        """
        Run the posterior analysis.

        Executes the solver for each alternative in the Prior distribution
        and aggregates results.
        """
        self.obj.runAnalyzer()

    def getAvgTable(self):
        """
        Get prior-weighted average performance metrics.

        Returns:
            pandas.DataFrame: Performance metrics table with prior-weighted
                expectations for QLen, Util, RespT, ResidT, ArvR, Tput.
        """
        table = self.obj.getAvgTable()

        QLen = np.array(list(table.getQLen()))
        Util = np.array(list(table.getUtil()))
        RespT = np.array(list(table.getRespT()))
        ResidT = np.array(list(table.getResidT()))
        ArvR = np.array(list(table.getArvR()))
        Tput = np.array(list(table.getTput()))

        cols = ['QLen', 'Util', 'RespT', 'ResidT', 'ArvR', 'Tput']
        stations = list(table.getStationNames())
        stationnames = []
        for i in range(len(stations)):
            stationnames.append(str(stations[i]))
        jobclasses = list(table.getClassNames())
        classnames = []
        for i in range(len(jobclasses)):
            classnames.append(str(jobclasses[i]))

        AvgTable = pd.DataFrame(
            np.concatenate([[QLen, Util, RespT, ResidT, ArvR, Tput]]).T,
            columns=cols
        )
        tokeep = ~(AvgTable <= 0.0).all(axis=1)
        AvgTable.insert(0, "JobClass", classnames)
        AvgTable.insert(0, "Station", stationnames)
        AvgTable = AvgTable.loc[tokeep]

        if not self._table_silent:
            print(AvgTable)

        return AvgTable

    def getPosteriorTable(self):
        """
        Get per-alternative results with probabilities.

        Returns:
            pandas.DataFrame: Table with columns:
                - Alternative: Index of the alternative (0-based)
                - Probability: Prior probability of this alternative
                - Station: Station name
                - JobClass: Job class name
                - Q, U, R, T, A: Performance metrics
        """
        table = self.obj.getPosteriorTable()
        rows = table.rows

        data = []
        for row in rows:
            data.append({
                'Alternative': row.alternativeIdx,
                'Probability': row.probability,
                'Station': row.station,
                'JobClass': row.jobClass,
                'Q': row.Q,
                'U': row.U,
                'R': row.R,
                'T': row.T,
                'A': row.A
            })

        df = pd.DataFrame(data)
        if not self._table_silent:
            print(df)
        return df

    def getPosteriorDist(self, metric, station, job_class):
        """
        Get the posterior distribution for a specific metric.

        Args:
            metric (str): Metric name ('Q', 'U', 'R', 'T', 'A')
            station: Station object
            job_class: JobClass object

        Returns:
            dict: Dictionary with 'values', 'probabilities', 'cdf' arrays
        """
        java_cdf = self.obj.getPosteriorDist(metric, station.obj, job_class.obj)

        return {
            'values': np.array([float(v) for v in java_cdf.values]),
            'probabilities': np.array([float(p) for p in java_cdf.probabilities]),
            'cdf': np.array([float(c) for c in java_cdf.cdf]),
            'mean': float(java_cdf.getMean())
        }

    def hasPriorDistribution(self):
        """
        Check if the model has a Prior distribution.

        Returns:
            bool: True if the model contains a Prior distribution
        """
        return bool(self.obj.hasPriorDistribution())

    def getNumAlternatives(self):
        """
        Get the number of alternatives in the Prior.

        Returns:
            int: Number of alternatives
        """
        return int(self.obj.getNumAlternatives())

    @staticmethod
    def defaultOptions():
        """Returns default solver options for the Posterior solver."""
        java_options = jpype.JPackage('jline').solvers.posterior.SolverPosterior.defaultOptions()
        python_options = SolverOptions.__new__(SolverOptions)
        python_options.obj = java_options
        return python_options

    # Snake_case aliases
    run_analyzer = runAnalyzer
    get_avg_table = getAvgTable
    get_posterior_table = getPosteriorTable
    get_posterior_dist = getPosteriorDist
    has_prior_distribution = hasPriorDistribution
    get_num_alternatives = getNumAlternatives
    default_options = defaultOptions


# Update the Posterior alias
Posterior = SolverPosterior
