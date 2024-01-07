import jpype
import jpype.imports
import numpy as np
from pprint import pformat
from . import SchedStrategy, jlineFromDistribution, SchedStrategyType, CallType
from .lang import jlineMatrixToArray


class LayeredNetwork:
    def __init__(self, name):
        self.obj = jpype.JPackage('jline').lang.layered.LayeredNetwork(name)

    def parseXML(self, filename, verbose=False):
        self.obj.parseXML(filename, verbose)

    def writeXML(self, filename, abstractNames=False):
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
        return self.obj.getNumberOfLayers()

    def getNumberOfModels(self):
        return self.obj.getNumberOfModels()

    def summary(self):
        return self.obj.summary()

    def parseXML(self, filename, verbose):
        return self.obj.parseXML(filename, verbose)

    def getStruct(self):
        jsn = self.obj.getStruct()
        lsn = LayeredNetworkStruct()
        lsn.fromJline(jsn)
        return lsn


class LayeredNetworkStruct():
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
    def __init__(self, model, name, mult, schedStrategy):
        self.obj = jpype.JPackage('jline').lang.layered.Processor(model.obj, name, mult, schedStrategy.value)


class Task:
    def __init__(self, model, name, mult, schedStrategy):
        self.obj = jpype.JPackage('jline').lang.layered.Task(model.obj, name, mult, schedStrategy.value)

    def on(self, proc):
        self.obj.on(proc.obj)
        return self

    def setThinkTime(self, distrib):
        self.obj.setThinkTime(distrib.obj)
        return self

    def addPrecedence(self, prec):
        self.obj.addPrecedence(prec)
        return self


class Entry:
    def __init__(self, model, name):
        self.obj = jpype.JPackage('jline').lang.layered.Entry(model.obj, name)

    def on(self, proc):
        self.obj.on(proc.obj)
        return self


class Activity:
    def __init__(self, model, name, distrib):
        self.obj = jpype.JPackage('jline').lang.layered.Activity(model.obj, name, distrib.obj)

    def on(self, proc):
        self.obj.on(proc.obj)
        return self

    def boundTo(self, proc):
        self.obj.boundTo(proc.obj)
        return self

    def repliesTo(self, entry):
        self.obj.repliesTo(entry.obj)
        return self

    def synchCall(self, entry, callmult):
        self.obj.synchCall(entry.obj, callmult)
        return self


class ActivityPrecedence:
    def __init__(self, name):
        self.obj = jpype.JPackage('jline').lang.layered.ActivityPrecedence(name)

    @staticmethod
    def Serial(act0, act1):
        return jpype.JPackage('jline').lang.layered.ActivityPrecedence.Serial(
            jpype.java.util.ArrayList((act0.obj, act1.obj)))

    @staticmethod
    def Sequence(act0, act1):
        return jpype.JPackage('jline').lang.layered.ActivityPrecedence.Sequence(
            act0.obj.getName(), act1.obj.getName())
