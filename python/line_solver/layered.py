import jpype
import jpype.imports


class LayeredNetwork:
    def __init__(self, name):
        self.obj = jpype.JPackage('jline').lang.layered.LayeredNetwork(name)

    def writeXML(self, filename, abstractNames=False):
        self.obj.writeXML(filename, abstractNames)


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
