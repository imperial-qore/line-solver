import os
import jpype
import jpype.imports
import numpy as np
import pandas as pd
from jpype import JArray

from . import jlineMatrixToArray, is_interactive
from .constants import SolverType, VerboseLevel, GlobalConstants


class Solver:
    def __init__(self, options, args):
        self.solveropt = options
        if len(args) > 1:
            ctr = 0
            for ctr in range(len(args) - 1):
                match args[ctr]:
                    case 'cutoff':
                        self.solveropt.obj.cutoff = args[ctr + 1]
                    case 'method':
                        self.solveropt.obj.method = args[ctr + 1]
                    case 'keep':
                        self.solveropt.obj.keep = args[ctr + 1]
                    case 'seed':
                        self.solveropt.obj.seed = args[ctr + 1]
                    case 'samples':
                        self.solveropt.obj.samples = args[ctr + 1]
                    case 'verbose':
                        if isinstance(args[ctr + 1], bool):
                            if args[ctr + 1]:
                                self.solveropt.obj.verbose = VerboseLevel.STD.value
                            else:
                                self.solveropt.obj.verbose = VerboseLevel.SILENT.value
                        else:
                            self.solveropt.obj.verbose = args[ctr + 1].value
                ctr += 2

    def getName(self):
        return self.obj.getName()

    def getAvgTable(self):
        table = self.obj.getAvgTable()
        # convert to NumPy

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
        AvgTable = AvgTable.loc[tokeep]  # eliminate zero rows
        if not (
                GlobalConstants.getVerbose() == VerboseLevel.SILENT) and not self.solveropt.obj.verbose == VerboseLevel.SILENT.value:  # and not is_interactive():
            print(AvgTable)

        return AvgTable

    def getAvgRespTTable(self):
        table = self.obj.getAvgTable()
        # convert to NumPy
        RespT = np.array(list(table.getRespT()))

        cols = ['RespT']
        stations = list(table.getStationNames())
        statnames = []
        for i in range(len(stations)):
            statnames.append(str(stations[i]))
        jobclasses = list(table.getClassNames())
        classnames = []
        for i in range(len(jobclasses)):
            classnames.append(str(jobclasses[i]))
        AvgTable = pd.DataFrame(np.concatenate([[RespT]]).T, columns=cols)
        tokeep = ~(AvgTable <= 0.0).all(axis=1)
        AvgTable.insert(0, "JobClass", classnames)
        AvgTable.insert(0, "Station", statnames)
        AvgTable = AvgTable.loc[tokeep]  # eliminate zero rows
        if not (
                GlobalConstants.getVerbose() == VerboseLevel.SILENT) and not self.solveropt.obj.verbose == VerboseLevel.SILENT.value:  # and not is_interactive():
            print(AvgTable)

        return AvgTable

    def getAvgResidTTable(self):
        table = self.obj.getAvgTable()

        # convert to NumPy
        ResidT = np.array(list(table.getResidT()))

        cols = ['ResidT']
        stations = list(table.getStationNames())
        statnames = []
        for i in range(len(stations)):
            statnames.append(str(stations[i]))
        jobclasses = list(table.getClassNames())
        classnames = []
        for i in range(len(jobclasses)):
            classnames.append(str(jobclasses[i]))
        AvgTable = pd.DataFrame(np.concatenate([[ResidT]]).T, columns=cols)
        tokeep = ~(AvgTable <= 0.0).all(axis=1)
        AvgTable.insert(0, "JobClass", classnames)
        AvgTable.insert(0, "Station", statnames)
        AvgTable = AvgTable.loc[tokeep]  # eliminate zero rows
        if not (GlobalConstants.getVerbose() == VerboseLevel.SILENT) and not self.solveropt.obj.verbose == VerboseLevel.SILENT.value:  # and not is_interactive():
            print(AvgTable)

        return AvgTable

    def getAvgUtilTable(self):
        table = self.obj.getAvgTable()

        # convert to NumPy
        Util = np.array(list(table.getUtil()))

        cols = ['Util']
        stations = list(table.getStationNames())
        statnames = []
        for i in range(len(stations)):
            statnames.append(str(stations[i]))
        jobclasses = list(table.getClassNames())
        classnames = []
        for i in range(len(jobclasses)):
            classnames.append(str(jobclasses[i]))
        AvgTable = pd.DataFrame(np.concatenate([[Util]]).T, columns=cols)
        tokeep = ~(AvgTable <= 0.0).all(axis=1)
        AvgTable.insert(0, "JobClass", classnames)
        AvgTable.insert(0, "Station", statnames)
        AvgTable = AvgTable.loc[tokeep]  # eliminate zero rows
        if not (GlobalConstants.getVerbose() == VerboseLevel.SILENT) and not self.solveropt.obj.verbose == VerboseLevel.SILENT.value:  # and not is_interactive():
            print(AvgTable)

        return AvgTable

    def getAvgQLenTable(self):
        table = self.obj.getAvgTable()

        # convert to NumPy
        QLen = np.array(list(table.getQLen()))

        cols = ['QLen']
        stations = list(table.getStationNames())
        statnames = []
        for i in range(len(stations)):
            statnames.append(str(stations[i]))
        jobclasses = list(table.getClassNames())
        classnames = []
        for i in range(len(jobclasses)):
            classnames.append(str(jobclasses[i]))
        AvgTable = pd.DataFrame(np.concatenate([[QLen]]).T, columns=cols)
        tokeep = ~(AvgTable <= 0.0).all(axis=1)
        AvgTable.insert(0, "JobClass", classnames)
        AvgTable.insert(0, "Station", statnames)
        AvgTable = AvgTable.loc[tokeep]  # eliminate zero rows
        if not (GlobalConstants.getVerbose() == VerboseLevel.SILENT) and not self.solveropt.obj.verbose == VerboseLevel.SILENT.value:  # and not is_interactive():
            print(AvgTable)

        return AvgTable

    def getAvgTputTable(self):
        table = self.obj.getAvgTable()

        # convert to NumPy
        Tput = np.array(list(table.getTput()))

        cols = ['Tput']
        stations = list(table.getStationNames())
        statnames = []
        for i in range(len(stations)):
            statnames.append(str(stations[i]))
        jobclasses = list(table.getClassNames())
        classnames = []
        for i in range(len(jobclasses)):
            classnames.append(str(jobclasses[i]))
        AvgTable = pd.DataFrame(np.concatenate([[Tput]]).T, columns=cols)
        tokeep = ~(AvgTable <= 0.0).all(axis=1)
        AvgTable.insert(0, "JobClass", classnames)
        AvgTable.insert(0, "Station", statnames)
        AvgTable = AvgTable.loc[tokeep]  # eliminate zero rows
        if not (GlobalConstants.getVerbose() == VerboseLevel.SILENT) and not self.solveropt.obj.verbose == VerboseLevel.SILENT.value:  # and not is_interactive():
            print(AvgTable)

        return AvgTable

    def getAvgArvRTable(self):
        table = self.obj.getAvgTable()

        # convert to NumPy
        ArvR = np.array(list(table.getArvR()))

        cols = ['ArvR']
        stations = list(table.getStationNames())
        statnames = []
        for i in range(len(stations)):
            statnames.append(str(stations[i]))
        jobclasses = list(table.getClassNames())
        classnames = []
        for i in range(len(jobclasses)):
            classnames.append(str(jobclasses[i]))
        AvgTable = pd.DataFrame(np.concatenate([[ArvR]]).T, columns=cols)
        tokeep = ~(AvgTable <= 0.0).all(axis=1)
        AvgTable.insert(0, "JobClass", classnames)
        AvgTable.insert(0, "Station", statnames)
        AvgTable = AvgTable.loc[tokeep]  # eliminate zero rows
        if not (GlobalConstants.getVerbose() == VerboseLevel.SILENT) and not self.solveropt.obj.verbose == VerboseLevel.SILENT.value:  # and not is_interactive():
            print(AvgTable)

        return AvgTable

    def getAvgSysTable(self):
        table = self.obj.getAvgSysTable()

        # convert to NumPy
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
        AvgSysTable = AvgSysTable.loc[tokeep]  # eliminate zero rows
        if not (
                GlobalConstants.getVerbose() == VerboseLevel.SILENT) and not self.solveropt.obj.verbose == VerboseLevel.SILENT.value:
            print(AvgSysTable)
        return AvgSysTable

    def getAvgTput(self):
        Tput = jlineMatrixToArray(self.obj.getAvgTput())
        if not (
                GlobalConstants.getVerbose() == VerboseLevel.SILENT) and not self.solveropt.obj.verbose == VerboseLevel.SILENT.value:
            print(Tput)
        return Tput

    def getAvgResidT(self):
        ResidT = jlineMatrixToArray(self.obj.getAvgResidT())
        if not (
                GlobalConstants.getVerbose() == VerboseLevel.SILENT) and not self.solveropt.obj.verbose == VerboseLevel.SILENT.value:
            print(ResidT)
        return ResidT

    def getAvgArvR(self):
        ArvR = jlineMatrixToArray(self.obj.getAvgArvR())
        if not (
                GlobalConstants.getVerbose() == VerboseLevel.SILENT) and not self.solveropt.obj.verbose == VerboseLevel.SILENT.value:
            print(ArvR)
        return ArvR

    def getAvgUtil(self):
        Util = jlineMatrixToArray(self.obj.getAvgUtil())
        if not (
                GlobalConstants.getVerbose() == VerboseLevel.SILENT) and not self.solveropt.obj.verbose == VerboseLevel.SILENT.value:
            print(Util)
        return Util

    def getAvgQLen(self):
        QLen = jlineMatrixToArray(self.obj.getAvgQLen())
        if not (
                GlobalConstants.getVerbose() == VerboseLevel.SILENT) and not self.solveropt.obj.verbose == VerboseLevel.SILENT.value:
            print(QLen)
        return QLen

    def getAvgRespT(self):
        RespT = jlineMatrixToArray(self.obj.getAvgRespT())
        if not (
                GlobalConstants.getVerbose() == VerboseLevel.SILENT) and not self.solveropt.obj.verbose == VerboseLevel.SILENT.value:
            print(RespT)
        return RespT

    def getAvgSysTput(self):
        SysTput = jlineMatrixToArray(self.obj.getAvgSysTput())
        if not (
                GlobalConstants.getVerbose() == VerboseLevel.SILENT) and not self.solveropt.obj.verbose == VerboseLevel.SILENT.value:
            print(SysTput)
        return SysTput

    def getAvgSysRespT(self):
        SysRespT = jlineMatrixToArray(self.obj.getAvgSysRespT())
        if not (
                GlobalConstants.getVerbose() == VerboseLevel.SILENT) and not self.solveropt.obj.verbose == VerboseLevel.SILENT.value:
            print(SysRespT)
        return SysRespT

    def getCdfRespT(self):
        try:
            table = self.obj.getCdfRespT()
            distribC = self.obj.fluidResult.distribC
            CdfRespT = []
            for i in range(distribC.length):
                for c in range(distribC[i].length):
                    F = jlineMatrixToArray(distribC[i][c])
                    CdfRespT.append(F)
            return np.asarray(CdfRespT)
        except:
            return [[]]


class SolverCTMC(Solver):
    def __init__(self, *args):
        options = SolverOptions(jpype.JPackage('jline').lang.constant.SolverType.CTMC)
        super().__init__(options, args)
        model = args[0]
        self.obj = jpype.JPackage('jline').solvers.ctmc.SolverCTMC(model.obj, self.solveropt.obj)

    def getStateSpace(self):
        return self.obj.getStateSpace()

    def getGenerator(self):
        return self.obj.getGenerator()


class SolverEnv(Solver):
    def __init__(self, *args):
        options = SolverOptions(jpype.JPackage('jline').lang.constant.SolverType.Env)
        super().__init__(options, args)
        model = args[0]
        solvers = jpype.JPackage('jline').solvers.NetworkSolver[len(args[1])]
        for i in range(len(solvers)):
            solvers[i] = args[1][i].obj
        self.obj = jpype.JPackage('jline').solvers.env.SolverEnv(model.obj, solvers, self.solveropt.obj)

    def getEnsembleAvg(self):
        return self.obj.getEnsembleAvg()

    def printAvgTable(self):
        self.obj.printAvgTable()


class SolverFluid(Solver):
    def __init__(self, *args):
        options = SolverOptions(jpype.JPackage('jline').lang.constant.SolverType.Fluid)
        super().__init__(options, args)
        model = args[0]
        self.obj = jpype.JPackage('jline').solvers.fluid.SolverFluid(model.obj, self.solveropt.obj)


class SolverJMT(Solver):
    def __init__(self, *args):
        options = SolverOptions(jpype.JPackage('jline').lang.constant.SolverType.JMT)
        super().__init__(options, args)
        model = args[0]
        self.jmtPath = jpype.JPackage('java').lang.String(os.path.dirname(os.path.abspath(__file__)) + "/JMT.jar")
        self.obj = jpype.JPackage('jline').solvers.jmt.SolverJMT(model.obj, self.solveropt.obj, self.jmtPath)

    def jsimwView(self):
        self.obj.jsimwView(self.jmtPath)

    def jsimgView(self):
        self.obj.jsimgView(self.jmtPath)


class SolverMAM(Solver):
    def __init__(self, *args):
        options = SolverOptions(jpype.JPackage('jline').lang.constant.SolverType.MAM)
        super().__init__(options, args)
        model = args[0]
        self.obj = jpype.JPackage('jline').solvers.mam.SolverMAM(model.obj, self.solveropt.obj)


class SolverMVA(Solver):
    def __init__(self, *args):
        options = SolverOptions(jpype.JPackage('jline').lang.constant.SolverType.MVA)
        super().__init__(options, args)
        model = args[0]
        self.obj = jpype.JPackage('jline').solvers.mva.SolverMVA(model.obj, self.solveropt.obj)

class SolverLQNS(Solver):
    def __init__(self, *args):
        options = SolverOptions(jpype.JPackage('jline').lang.constant.SolverType.LQNS)
        super().__init__(options, args)
        model = args[0]
        self.obj = jpype.JPackage('jline').solvers.lqns.SolverLQNS(model.obj, self.solveropt.obj)


    def getAvgTable(self):
        table = self.obj.getAvgTable()
        # convert to NumPy

        QLen = np.array(list(table.getQLen()))
        Util = np.array(list(table.getUtil()))
        RespT = np.array(list(table.getRespT()))
        ResidT = np.array(list(table.getResidT()))
        Tput = np.array(list(table.getTput()))

        cols = ['QLen', 'Util', 'RespT', 'ResidT', 'Tput']
        nodenames = list(table.getNodeNames())
        mynodenames = []
        for i in range(len(nodenames)):
            mynodenames.append(str(nodenames[i]))
        nodetypes = list(table.getNodeTypes())
        mynodetypes = []
        for i in range(len(nodetypes)):
            mynodetypes.append(str(nodetypes[i]))
        AvgTable = pd.DataFrame(np.concatenate([[QLen, Util, RespT, ResidT, Tput]]).T, columns=cols)
        tokeep = ~(AvgTable <= 0.0).all(axis=1)
        AvgTable.insert(0, "NodeType", mynodetypes)
        AvgTable.insert(0, "Node", mynodenames)
        AvgTable = AvgTable.loc[tokeep]  # eliminate zero rows
        if not (
                GlobalConstants.getVerbose() == VerboseLevel.SILENT) and not self.solveropt.obj.verbose == VerboseLevel.SILENT.value:  # and not is_interactive():
            print(AvgTable)

        return AvgTable


class SolverLN(Solver):
    def __init__(self, *args):
        options = SolverOptions(jpype.JPackage('jline').lang.constant.SolverType.LN)
        super().__init__(options, args)
        model = args[0]
        self.obj = jpype.JPackage('jline').solvers.ln.SolverLN(model.obj, self.solveropt.obj)


class SolverNC(Solver):
    def __init__(self, *args):
        options = SolverOptions(jpype.JPackage('jline').lang.constant.SolverType.NC)
        super().__init__(options, args)
        model = args[0]
        self.obj = jpype.JPackage('jline').solvers.nc.SolverNC(model.obj, self.solveropt.obj)


class SolverSSA(Solver):
    def __init__(self, *args):
        options = SolverOptions(jpype.JPackage('jline').lang.constant.SolverType.SSA)
        super().__init__(options, args)
        model = args[0]
        self.obj = jpype.JPackage('jline').solvers.ssa.SolverSSA(model.obj, self.solveropt.obj)


class SolverOptions():
    def __init__(self, solvertype):
        self.obj = jpype.JPackage('jline').solvers.SolverOptions(solvertype)
