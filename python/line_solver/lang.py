import jpype
import jpype.imports

from line_solver import jlineMatrixToArray, jlineMapMatrixToArray, jlineArrayToMatrix
from .constants import NodeType


class RoutingMatrix:
    def __init__(self, rt):
        self.obj = rt

    def set(self, class_source, class_dest, stat_source, stat_dest, prob):
        return self.obj.set(class_source.obj, class_dest.obj, stat_source.obj, stat_dest.obj, prob)

    def setRoutingMatrix(self, jobclass, node, pmatrix):
        for i in range(len(node)):
            for j in range(len(node)):
                for k in range(len(jobclass)):
                    self.set(jobclass[k], jobclass[k], node[i], node[j], pmatrix[k][i][j])
class Network:
    def __init__(self, name):
        self.obj = jpype.JPackage('jline').lang.Network(name)

    def serialRouting(*argv):
        rtlist = jpype.JPackage('jline').lang.nodes.Node[len(argv)]
        ctr = 0
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

    def jsimgView(self):
        self.obj.jsimgView()

    def addLinks(self, linkPairs):
        for i in range(len(linkPairs)):
            self.obj.addLink(linkPairs[i][0].obj,linkPairs[i][1].obj)

    def getStruct(self, force=True):
        jsn = self.obj.getStruct(force)
        sn = dict(nstations=int(jsn.nstations),
                  nstateful=int(jsn.nstateful),
                  nnodes=int(jsn.nnodes),
                  nclasses=int(jsn.nclasses),
                  nclosedjobs=int(jsn.nclosedjobs),
                  nchains=int(jsn.nchains),
                  refstat=jlineMatrixToArray(jsn.refstat),
                  njobs=jlineMatrixToArray(jsn.njobs),
                  nservers=jlineMatrixToArray(jsn.nservers),
                  connmatrix=jlineMatrixToArray(jsn.connmatrix),
                  scv=jlineMatrixToArray(jsn.scv),
                  isstation=jlineMatrixToArray(jsn.isstation),
                  isstateful=jlineMatrixToArray(jsn.isstateful),
                  isstatedep=jlineMatrixToArray(jsn.isstatedep),
                  nodeToStateful=jlineMatrixToArray(jsn.nodeToStateful),
                  nodeToStation=jlineMatrixToArray(jsn.nodeToStation),
                  stationToNode=jlineMatrixToArray(jsn.stationToNode),
                  stationToStateful=jlineMatrixToArray(jsn.stationToStateful),
                  statefulToNode=jlineMatrixToArray(jsn.statefulToNode),
                  rates=jlineMatrixToArray(jsn.rates),
                  classprio=jlineMatrixToArray(jsn.classprio),
                  phases=jlineMatrixToArray(jsn.phases),
                  phasessz=jlineMatrixToArray(jsn.phasessz),
                  phaseshift=jlineMatrixToArray(jsn.phaseshift),
                  schedparam=jlineMatrixToArray(jsn.schedparam),
                  chains=jlineMatrixToArray(jsn.chains),
                  rt=jlineMatrixToArray(jsn.rt),
                  nvars=jlineMatrixToArray(jsn.nvars),
                  rtnodes=jlineMatrixToArray(jsn.rtnodes),
                  csmask=jlineMatrixToArray(jsn.csmask),
                  isslc=jlineMatrixToArray(jsn.isslc),
                  cap=jlineMatrixToArray(jsn.cap),
                  refclass=jlineMatrixToArray(jsn.refclass),
                  lldscaling=jlineMatrixToArray(jsn.lldscaling),
                  fj=jlineMatrixToArray(jsn.fj),
                  classcap=jlineMatrixToArray(jsn.classcap),
                  inchain = jlineMapMatrixToArray(jsn.inchain),
                  visits = jlineMapMatrixToArray(jsn.visits),
                  nodevisits = jlineMapMatrixToArray(jsn.nodevisits),
                  classnames=tuple(jsn.classnames),
                  nodetypes=tuple(map(lambda x: NodeType.fromJava(x), jsn.nodetypes)),
                  nodenames=tuple(jsn.nodenames))

        # SerializableFunction<Pair<Map<Node, Matrix>, Map<Node, Matrix>>, Matrix> rtfun;
        # public Map<JobClass, Map<JobClass, Matrix>> rtorig;
        # public Map<Station, Map<JobClass, SerializableFunction<Double, Double>>> lst;
        # public Map<StatefulNode, Matrix> state;
        # public Map<StatefulNode, Matrix> stateprior;
        # public Map<Station, Matrix> space;
        # public Map<Node, Map<JobClass, RoutingStrategy>> routing;
        # public Map<Station, Map<JobClass, ProcessType>> proctype;
        # public Map<Station, Map<JobClass, Matrix>> mu;
        # public Map<Station, Map<JobClass, Matrix>> phi;
        # public Map<Station, Map<JobClass, Map<Integer, Matrix>>> proc;
        # public Map<Station, Map<JobClass, Matrix>> pie;
        # public Map<Station, SchedStrategy> sched;
        # public Map<Station, Map<JobClass, DropStrategy>> droprule;	//This represents dropid in LINE
        # public Map<Node, NodeParam> nodeparam;
        # public Map<Integer, Sync> sync;
        # public Map<Station, SerializableFunction<Matrix, Double>> cdscaling;

        # public List<Station> stations;
        # public List<StatefulNode> stateful;
        # public List<JobClass> jobclasses;
        # public List<Node> nodes;

        return sn


class Cache:
    def __init__(self, model, name, nitems, itemLevelCap, replPolicy, graph=()):
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

class Env:
    def __init__(self, name, nstages):
        self.obj = jpype.JPackage('jline').lang.Env(name, nstages)

    def addStage(self, stage, envname, envtype, envmodel):
        self.obj.addStage(stage, envname, envtype, envmodel.obj)

    def addTransition(self, envname0, envname1, rate):
        self.obj.addTransition(envname0, envname1, rate.obj)

    def getStageTable(self):
        return self.obj.getStageTable()

class Source:
    def __init__(self, model, name):
        self.obj = jpype.JPackage('jline').lang.nodes.Source(model.obj, name)

    def setArrival(self, jobclass, distribution):
        self.obj.setArrival(jobclass.obj, distribution.obj)

class ClassSwitch:
    def __init__(self, model, name):
        self.obj = jpype.JPackage('jline').lang.nodes.ClassSwitch(model.obj, name)

    def initClassSwitchMatrix(self):
        return jlineMatrixToArray(self.obj.initClassSwitchMatrix())
    def setClassSwitchingMatrix(self, csmatrix):
        self.obj.setClassSwitchingMatrix(jlineArrayToMatrix(csmatrix))

class Sink:
    def __init__(self, model, name):
        self.obj = jpype.JPackage('jline').lang.nodes.Sink(model.obj, name)

class Fork:
    def __init__(self, model, name):
        self.obj = jpype.JPackage('jline').lang.nodes.Fork(model.obj, name)

class Join:
    def __init__(self, model, name, forknode):
        self.obj = jpype.JPackage('jline').lang.nodes.Join(model.obj, name, forknode.obj)

class Queue:
    def __init__(self, model, name, strategy):
        self.obj = jpype.JPackage('jline').lang.nodes.Queue(model.obj, name, strategy.value)

    def setService(self, jobclass, distribution):
        self.obj.setService(jobclass.obj, distribution.obj)

    def setNumberOfServers(self, nservers):
        self.obj.setNumberOfServers(nservers)


class Delay:
    def __init__(self, model, name):
        self.obj = jpype.JPackage('jline').lang.nodes.Delay(model.obj, name)

    def setService(self, jobclass, distribution):
        self.obj.setService(jobclass.obj, distribution.obj)


class Router:
    def __init__(self, model, name):
        self.obj = jpype.JPackage('jline').lang.nodes.Router(model.obj, name)
    def setRouting(self, jobclass, strategy):
        self.obj.setRouting(jobclass.obj, strategy.value)

class OpenClass:
    def __init__(self, model, name, prio=0):
        self.obj = jpype.JPackage('jline').lang.OpenClass(model.obj, name, prio)
        self.completes = False

    def __setattr__(self, name, value):
        object.__setattr__(self, name, value)
        if name == 'completes' and hasattr(self, 'completes'):
            if self.completes:
                self.obj.setCompletes(True)
            else:
                self.obj.setCompletes(False)


class ClosedClass:
    def __init__(self, model, name, njobs, refstat, prio=0):
        self.obj = jpype.JPackage('jline').lang.ClosedClass(model.obj, name, njobs, refstat.obj, prio)
        self.completes = False

    def __setattr__(self, name, value):
        object.__setattr__(self, name, value)
        if name == 'completes' and hasattr(self, 'completes'):
            if self.completes:
                self.obj.setCompletes(True)
            else:
                self.obj.setCompletes(False)


class SelfLoopingClass:
    def __init__(self, model, name, njobs, refstat, prio=0):
        self.obj = jpype.JPackage('jline').lang.SelfLoopingClass(model.obj, name, njobs, refstat.obj, prio)
        self.completes = False

    def __setattr__(self, name, value):
        object.__setattr__(self, name, value)
        if name == 'completes' and hasattr(self, 'completes'):
            if self.completes:
                self.obj.setCompletes(True)
            else:
                self.obj.setCompletes(False)
