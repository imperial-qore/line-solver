import jpype
import jpype.imports
import numpy as np

from . import jlineMatrixToArray, jlineMapMatrixToArray, jlineMatrixFromArray
from .constants import *


class JobClass:
    def __init__(self):
        pass


class Node:
    def __init__(self):
        pass

    def setRouting(self, jobclass, strategy):
        self.obj.setRouting(jobclass.obj, strategy.value)

    def setProbRouting(self, jobclass, node, prob):
        self.obj.setProbRouting(jobclass.obj, node.obj, prob)


class Station(Node):
    def __init__(self):
        super().__init__()


class RoutingMatrix:
    def __init__(self, rt):
        self.obj = rt

    def set(self, class_source, class_dest, stat_source, stat_dest, prob):
        return self.obj.set(class_source.obj, class_dest.obj, stat_source.obj, stat_dest.obj, prob)

    def setRoutingMatrix(self, jobclass, node, pmatrix):
        if isinstance(jobclass, JobClass):
            for i in range(len(node)):
                for j in range(len(node)):
                    self.set(jobclass, jobclass, node[i], node[j], pmatrix[i][j])
        else:
            for i in range(len(node)):
                for j in range(len(node)):
                    for k in range(len(jobclass)):
                        self.set(jobclass[k], jobclass[k], node[i], node[j], pmatrix[k][i][j])


class Model:
    def __init__(self):
        pass

    def getName(self):
        return self.obj.getName()

    def setName(self, name):
        self.obj.setName(name)

    def getVersion(self):
        return self.obj.getVersion()


class NetworkStruct(dict):
    def fromJline(self, jsn):
        super().__setitem__("nstations", int(jsn.nstations))
        super().__setitem__("nstateful", int(jsn.nstateful))
        super().__setitem__("nnodes", int(jsn.nnodes))
        super().__setitem__("nclasses", int(jsn.nclasses))
        super().__setitem__("nclosedjobs", int(jsn.nclosedjobs))
        super().__setitem__("nchains", int(jsn.nchains))
        super().__setitem__("refstat", jlineMatrixToArray(jsn.refstat))
        super().__setitem__("njobs", jlineMatrixToArray(jsn.njobs))
        super().__setitem__("nservers", jlineMatrixToArray(jsn.nservers))
        super().__setitem__("connmatrix", jlineMatrixToArray(jsn.connmatrix))
        super().__setitem__("scv", jlineMatrixToArray(jsn.scv))
        super().__setitem__("isstation", jlineMatrixToArray(jsn.isstation))
        super().__setitem__("isstateful", jlineMatrixToArray(jsn.isstateful))
        super().__setitem__("isstatedep", jlineMatrixToArray(jsn.isstatedep))
        super().__setitem__("nodeToStateful", jlineMatrixToArray(jsn.nodeToStateful))
        super().__setitem__("nodeToStation", jlineMatrixToArray(jsn.nodeToStation))
        super().__setitem__("stationToNode", jlineMatrixToArray(jsn.stationToNode))
        super().__setitem__("stationToStateful", jlineMatrixToArray(jsn.stationToStateful))
        super().__setitem__("statefulToNode", jlineMatrixToArray(jsn.statefulToNode))
        super().__setitem__("rates", jlineMatrixToArray(jsn.rates))
        super().__setitem__("classprio", jlineMatrixToArray(jsn.classprio))
        super().__setitem__("phases", jlineMatrixToArray(jsn.phases))
        super().__setitem__("phasessz", jlineMatrixToArray(jsn.phasessz))
        super().__setitem__("phaseshift", jlineMatrixToArray(jsn.phaseshift))
        super().__setitem__("schedparam", jlineMatrixToArray(jsn.schedparam))
        super().__setitem__("chains", jlineMatrixToArray(jsn.chains))
        super().__setitem__("rt", jlineMatrixToArray(jsn.rt))
        super().__setitem__("nvars", jlineMatrixToArray(jsn.nvars))
        super().__setitem__("rtnodes", jlineMatrixToArray(jsn.rtnodes))
        super().__setitem__("csmask", jlineMatrixToArray(jsn.csmask))
        super().__setitem__("isslc", jlineMatrixToArray(jsn.isslc))
        super().__setitem__("cap", jlineMatrixToArray(jsn.cap))
        super().__setitem__("refclass", jlineMatrixToArray(jsn.refclass))
        super().__setitem__("lldscaling", jlineMatrixToArray(jsn.lldscaling))
        super().__setitem__("fj", jlineMatrixToArray(jsn.fj))
        super().__setitem__("classcap", jlineMatrixToArray(jsn.classcap))
        super().__setitem__("inchain", jlineMapMatrixToArray(jsn.inchain))
        super().__setitem__("visits", jlineMapMatrixToArray(jsn.visits))
        super().__setitem__("nodevisits", jlineMapMatrixToArray(jsn.nodevisits))
        super().__setitem__("classnames", tuple(jsn.classnames))
        super().__setitem__("nodetypes", tuple(map(lambda x: NodeType.fromJLine(x), jsn.nodetypes)))
        super().__setitem__("nodenames", tuple(jsn.nodenames))

        sched = np.empty(int(jsn.nstations), dtype=object)
        space = np.empty(int(jsn.nstations), dtype=object)
        mu = np.empty(shape=(int(jsn.nstations), int(jsn.nclasses)), dtype=object)
        phi = np.empty(shape=(int(jsn.nstations), int(jsn.nclasses)), dtype=object)
        pie = np.empty(shape=(int(jsn.nstations), int(jsn.nclasses)), dtype=object)
        proctype = np.empty(shape=(int(jsn.nstations), int(jsn.nclasses)), dtype=object)
        routing = np.empty(shape=(int(jsn.nstations), int(jsn.nclasses)), dtype=object)
        droprule = np.empty(shape=(int(jsn.nstations), int(jsn.nclasses)), dtype=object)
        proc = np.empty(shape=(int(jsn.nstations), int(jsn.nclasses), 2), dtype=object)
        # TODO: missing in Jline, rtorig always set to None?
        #rtorig = np.empty(shape=(int(jsn.nstations), int(jsn.nclasses)), dtype=object)
        for ist in range(int(jsn.nstations)):
            sched[ist] = SchedStrategy(jsn.sched.get(jsn.stations[ist]))
            space[ist] = jlineMatrixToArray(jsn.space.get(jsn.stations[ist]))
            for jcl in range(int(jsn.nclasses)):
                mu[ist, jcl] = jlineMatrixToArray(jsn.mu.get(jsn.stations[ist]).get(jsn.jobclasses[jcl]))
                phi[ist, jcl] = jlineMatrixToArray(jsn.phi.get(jsn.stations[ist]).get(jsn.jobclasses[jcl]))
                pie[ist, jcl] = jlineMatrixToArray(jsn.pie.get(jsn.stations[ist]).get(jsn.jobclasses[jcl]))
                #rtorig[ist, jcl] = jlineMatrixToArray(jsn.rtorig.get(jsn.stations[ist]).get(jsn.jobclasses[jcl]))
                proctype[ist, jcl] = ProcessType(jsn.proctype.get(jsn.stations[ist]).get(jsn.jobclasses[jcl]))
                routing[ist, jcl] = RoutingStrategy(jsn.routing.get(jsn.stations[ist]).get(jsn.jobclasses[jcl]))
                droprule[ist, jcl] = DropStrategy(jsn.droprule.get(jsn.stations[ist]).get(jsn.jobclasses[jcl]))
                proc[ist, jcl, 0] = jlineMatrixToArray(jsn.proc.get(jsn.stations[ist]).get(jsn.jobclasses[jcl]).get(0))
                proc[ist, jcl, 1] = jlineMatrixToArray(jsn.proc.get(jsn.stations[ist]).get(jsn.jobclasses[jcl]).get(1))
        super().__setitem__("sched", sched)
        super().__setitem__("space", space)
        super().__setitem__("mu", mu)
        super().__setitem__("phi", phi)
        super().__setitem__("pie", pie)
        super().__setitem__("proctype", proctype)
        super().__setitem__("routing", routing)
        super().__setitem__("droprule", droprule)
        #super().__setitem__("rtorig", rtorig)
        super().__setitem__("proc", proc)

        # TODO: NodeParam class to be ported in Python
        # nodeparam = np.empty(int(jsn.nstations), dtype=object)
        # for ind in range(int(jsn.nnodes)):
        #    nodeparam[ind] = NodeParam(jsn.sched.get(jsn.nodes[ind]))

        # TODO: fields missing in JLINE
        #state = np.empty(int(jsn.nstateful), dtype=object)
        #stateprior = np.empty(int(jsn.nstateful), dtype=object)
        #for isf in range(int(jsn.nstateful)):
        #    state[isf] = jlineMatrixToArray(jsn.state.get(jsn.stateful[isf]))
        #    stateprior[isf] = jlineMatrixToArray(jsn.state.get(jsn.stateprior[isf]))
        #super().__setitem__("state", state)
        #super().__setitem__("stateprior", stateprior)

        # TODO: fields not parsed yet
        # SerializableFunction<Pair<Map<Node, Matrix>, Map<Node, Matrix>>, Matrix> rtfun;
        # public Map<Station, Map<JobClass, SerializableFunction<Double, Double>>> lst;
        # public Map<Station, SerializableFunction<Matrix, Double>> cdscaling;
        # public Map<Integer, Sync> sync;

class Network(Model):
    def __init__(self, *argv):
        super().__init__()
        if isinstance(argv[0], jpype.JPackage('jline').lang.Network):
            self.obj = argv[0]
        else:
            name = argv[0]
            self.obj = jpype.JPackage('jline').lang.Network(name)

    def serialRouting(*argv):
        ctr = 0
        if len(argv) == 1:
            rtlist = jpype.JPackage('jline').lang.nodes.Node[len(argv[0])]
            for arg in argv[0]:
                rtlist[ctr] = jpype.JObject(arg.obj, 'jline.lang.nodes.Node')
                ctr += 1
        else:
            rtlist = jpype.JPackage('jline').lang.nodes.Node[len(argv)]
            for arg in argv:
                rtlist[ctr] = jpype.JObject(arg.obj, 'jline.lang.nodes.Node')
                ctr += 1

        return RoutingMatrix(jpype.JPackage('jline').lang.Network.serialRouting(rtlist))

    def reset(self, hard=True):
        self.obj.reset(hard)

    def link(self, routing):
        self.obj.link(routing.obj)

    def addLink(self, source, dest):
        self.obj.addLink(source.obj, dest.obj)

    def initRoutingMatrix(self):
        rt = self.obj.initRoutingMatrix()
        return RoutingMatrix(rt)

    def getNumberOfNodes(self):
        return self.obj.getNumberOfNodes()

    def getNumberOfStations(self):
        return self.obj.getNumberOfStations()

    def getNumberOfClasses(self):
        return self.obj.getNumberOfClasses()

    def getTranHandles(self):
        Qt, Ut, Tt = self.obj.getTranHandles()
        return Qt, Ut, Tt

    def jsimgView(self):
        from line_solver import SolverJMT
        SolverJMT(self).jsimgView()

    def jsimwView(self):
        from line_solver import SolverJMT
        SolverJMT(self).jsimgView()

    def addLinks(self, linkPairs):
        for i in range(len(linkPairs)):
            self.obj.addLink(linkPairs[i][0].obj, linkPairs[i][1].obj)

    def getStruct(self, force=True):
        jsn = self.obj.getStruct(force)
        sn = NetworkStruct()
        sn.fromJline(jsn)
        return sn


class Cache(Node):
    def __init__(self, model, name, nitems, itemLevelCap, replPolicy, graph=()):
        super().__init__()
        if isinstance(itemLevelCap, int):
            if len(graph) == 0:
                self.obj = jpype.JPackage('jline').lang.nodes.Cache(model.obj, name, nitems,
                                                                    jpype.JPackage('jline').util.Matrix(itemLevelCap),
                                                                    replPolicy.value)
            else:
                self.obj = jpype.JPackage('jline').lang.nodes.Cache(model.obj, name, nitems,
                                                                    jpype.JPackage('jline').util.Matrix(itemLevelCap),
                                                                    replPolicy.value, graph)
        else:
            if len(graph) == 0:
                self.obj = jpype.JPackage('jline').lang.nodes.Cache(model.obj, name, nitems,
                                                                    jpype.JPackage('jline').util.Matrix(itemLevelCap),
                                                                    replPolicy.value)
            else:
                self.obj = jpype.JPackage('jline').lang.nodes.Cache(model.obj, name, nitems,
                                                                    jpype.JPackage('jline').util.Matrix(itemLevelCap),
                                                                    replPolicy.value, graph)

    def setRead(self, jobclass, distrib):
        self.obj.setRead(jobclass.obj, distrib.obj)

    def setHitClass(self, jobclass1, jobclass2):
        self.obj.setHitClass(jobclass1.obj, jobclass2.obj)

    def setMissClass(self, jobclass1, jobclass2):
        self.obj.setMissClass(jobclass1.obj, jobclass2.obj)


class Ensemble:
    def __init__(self):
        pass

    def getModel(self, stagenum):
        return Network(self.obj.getModel(stagenum))

    def getEnsemble(self):
        jensemble = self.obj.getEnsemble()
        ensemble = np.empty(jensemble.size(), dtype=object)
        for i in range(len(ensemble)):
            ensemble[i] = Network(jensemble.get(i))
        return ensemble


class Env(Ensemble):
    def __init__(self, name, nstages):
        super().__init__()
        self.obj = jpype.JPackage('jline').lang.Env(name, nstages)

    def addStage(self, stage, envname, envtype, envmodel):
        self.obj.addStage(stage, envname, envtype, envmodel.obj)

    def addTransition(self, envname0, envname1, rate):
        self.obj.addTransition(envname0, envname1, rate.obj)

    def getStageTable(self):
        return self.obj.getStageTable()


class Source(Station):
    def __init__(self, model, name):
        super().__init__()
        self.obj = jpype.JPackage('jline').lang.nodes.Source(model.obj, name)

    def setArrival(self, jobclass, distribution):
        self.obj.setArrival(jobclass.obj, distribution.obj)


class ClassSwitch(Node):
    def __init__(self, model, name):
        super().__init__()
        self.obj = jpype.JPackage('jline').lang.nodes.ClassSwitch(model.obj, name)

    def initClassSwitchMatrix(self):
        return jlineMatrixToArray(self.obj.initClassSwitchMatrix())

    def setClassSwitchingMatrix(self, csmatrix):
        self.obj.setClassSwitchingMatrix(jlineMatrixFromArray(csmatrix))


class Sink(Node):
    def __init__(self, model, name):
        super().__init__()
        self.obj = jpype.JPackage('jline').lang.nodes.Sink(model.obj, name)


class Fork(Node):
    def __init__(self, model, name):
        super().__init__()
        self.obj = jpype.JPackage('jline').lang.nodes.Fork(model.obj, name)


class Join(Station):
    def __init__(self, model, name, forknode):
        super().__init__()
        self.obj = jpype.JPackage('jline').lang.nodes.Join(model.obj, name, forknode.obj)


class Queue(Station):
    def __init__(self, model, name, strategy):
        super().__init__()
        self.obj = jpype.JPackage('jline').lang.nodes.Queue(model.obj, name, strategy.value)

    def setService(self, jobclass, distribution):
        self.obj.setService(jobclass.obj, distribution.obj)

    def setNumberOfServers(self, nservers):
        self.obj.setNumberOfServers(nservers)

    def setLoadDependence(self, ldscaling):
        self.obj.setLoadDependence(jlineMatrixFromArray(ldscaling))


class Delay(Station):
    def __init__(self, model, name):
        super().__init__()
        self.obj = jpype.JPackage('jline').lang.nodes.Delay(model.obj, name)

    def setService(self, jobclass, distribution):
        self.obj.setService(jobclass.obj, distribution.obj)


class Router(Node):
    def __init__(self, model, name):
        super().__init__()
        self.obj = jpype.JPackage('jline').lang.nodes.Router(model.obj, name)


class OpenClass(JobClass):
    def __init__(self, model, name, prio=0):
        super().__init__()
        self.obj = jpype.JPackage('jline').lang.OpenClass(model.obj, name, prio)
        self.completes = False

    def __setattr__(self, name, value):
        object.__setattr__(self, name, value)
        if name == 'completes' and hasattr(self, 'completes'):
            if self.completes:
                self.obj.setCompletes(True)
            else:
                self.obj.setCompletes(False)


class ClosedClass(JobClass):
    def __init__(self, model, name, njobs, refstat, prio=0):
        super().__init__()
        self.obj = jpype.JPackage('jline').lang.ClosedClass(model.obj, name, njobs, refstat.obj, prio)
        self.completes = False

    def __setattr__(self, name, value):
        object.__setattr__(self, name, value)
        if name == 'completes' and hasattr(self, 'completes'):
            if self.completes:
                self.obj.setCompletes(True)
            else:
                self.obj.setCompletes(False)


class SelfLoopingClass(JobClass):
    def __init__(self, model, name, njobs, refstat, prio=0):
        super().__init__()
        self.obj = jpype.JPackage('jline').lang.SelfLoopingClass(model.obj, name, njobs, refstat.obj, prio)
        self.completes = False

    def __setattr__(self, name, value):
        object.__setattr__(self, name, value)
        if name == 'completes' and hasattr(self, 'completes'):
            if self.completes:
                self.obj.setCompletes(True)
            else:
                self.obj.setCompletes(False)
