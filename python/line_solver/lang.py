import jpype
import jpype.imports

from line_solver import jlineToArray
from .constants import NodeType


class RoutingMatrix:
    def __init__(self, rt):
        self.obj = rt

    def set(self, class_source, class_dest, stat_source, stat_dest, prob):
        return self.obj.set(class_source.obj, class_dest.obj, stat_source.obj, stat_dest.obj, prob)


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

    def link(self, routing):
        self.obj.link(routing.obj)

    def addLink(self, source, dest):
        self.obj.addLink(source.obj, dest.obj)

    def initRoutingMatrix(self):
        rt = self.obj.initRoutingMatrix()
        return RoutingMatrix(rt)

    def jsimgView(self):
        self.obj.jsimgView()

    def getStruct(self, force=True):
        jsn = self.obj.getStruct(force)
        sn = dict(nstations=int(jsn.nstations),
                  nstateful=int(jsn.nstateful),
                  nnodes=int(jsn.nnodes),
                  nclasses=int(jsn.nclasses),
                  nclosedjobs=int(jsn.nclosedjobs),
                  nchains=int(jsn.nchains),
                  refstat=jlineToArray(jsn.refstat),
                  njobs=jlineToArray(jsn.njobs),
                  nservers=jlineToArray(jsn.nservers),
                  connmatrix=jlineToArray(jsn.connmatrix),
                  scv=jlineToArray(jsn.scv),
                  isstation=jlineToArray(jsn.isstation),
                  isstateful=jlineToArray(jsn.isstateful),
                  isstatedep=jlineToArray(jsn.isstatedep),
                  nodeToStateful=jlineToArray(jsn.nodeToStateful),
                  nodeToStation=jlineToArray(jsn.nodeToStation),
                  stationToNode=jlineToArray(jsn.stationToNode),
                  stationToStateful=jlineToArray(jsn.stationToStateful),
                  statefulToNode=jlineToArray(jsn.statefulToNode),
                  rates=jlineToArray(jsn.rates),
                  classprio=jlineToArray(jsn.classprio),
                  phases=jlineToArray(jsn.phases),
                  phasessz=jlineToArray(jsn.phasessz),
                  phaseshift=jlineToArray(jsn.phaseshift),
                  schedparam=jlineToArray(jsn.schedparam),
                  chains=jlineToArray(jsn.chains),
                  rt=jlineToArray(jsn.rt),
                  nvars=jlineToArray(jsn.nvars),
                  rtnodes=jlineToArray(jsn.rtnodes),
                  csmask=jlineToArray(jsn.csmask),
                  isslc=jlineToArray(jsn.isslc),
                  cap=jlineToArray(jsn.cap),
                  refclass=jlineToArray(jsn.refclass),
                  lldscaling=jlineToArray(jsn.lldscaling),
                  fj=jlineToArray(jsn.fj),
                  classcap=jlineToArray(jsn.classcap),
                  classnames=tuple(jsn.classnames),
                  nodetypes=tuple(map(lambda x: NodeType.fromJava(x), jsn.nodetypes)),
                  nodenames=tuple(jsn.nodenames))

        # SerializableFunction<Pair<Map<Node, Matrix>, Map<Node, Matrix>>, Matrix> rtfun;
        #
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
        # public Map<Integer, Matrix> inchain;
        # public Map<Integer, Matrix> visits;	//The integer represents the chain's ID (inchain)
        # public Map<Integer, Matrix> nodevisits; //The integer represents the chain's ID (inchain)
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


class Source:
    def __init__(self, model, name):
        self.obj = jpype.JPackage('jline').lang.nodes.Source(model.obj, name)

    def setArrival(self, jobclass, distribution):
        self.obj.setArrival(jobclass.obj, distribution.obj)


class Sink:
    def __init__(self, model, name):
        self.obj = jpype.JPackage('jline').lang.nodes.Sink(model.obj, name)


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


class OpenClass:
    def __init__(self, model, name):
        self.obj = jpype.JPackage('jline').lang.OpenClass(model.obj, name)
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
