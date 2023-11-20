import jpype
import jpype.imports

from line_solver import jlineArrayToMatrix

class APH:
    def __init__(self, *args):
        if len(args) == 1:
            self.obj = args[0]
        else:
            alpha = args[0]
            T = args[1]
            self.obj = jpype.JPackage('jline').lang.distributions.APH(jlineArrayToMatrix(alpha).toList1D(), jlineArrayToMatrix(T))

class Cox2:
    def __init__(self, *args):
        if len(args) == 1:
            self.obj = args[0]
        else:
            mu1 = args[0]
            mu2 = args[1]
            phi1 = args[2]
            self.obj = jpype.JPackage('jline').lang.distributions.Cox2(mu1, mu2, phi1)

    def fitMeanAndSCV(mean, scv):
        return Cox2(jpype.JPackage('jline').lang.distributions.Cox2.fitMeanAndSCV(mean, scv))


class Det:
    def __init__(self, value):
        self.obj = jpype.JPackage('jline').lang.distributions.Det(value)


class Disabled:
    def __init__(self, value):
        self.obj = jpype.JPackage('jline').lang.distributions.Disabled()


class Exp:
    def __init__(self, rate):
        self.obj = jpype.JPackage('jline').lang.distributions.Exp(rate)

    def fitMean(mean):
        return Erlang(jpype.JPackage('jline').lang.distributions.Exp.fitMean(mean))


class Erlang:
    def __init__(self, *args):
        if len(args) == 1:
            self.obj = args[0]
        else:
            rate = args[0]
            nphases = args[1]
            self.obj = jpype.JPackage('jline').lang.distributions.Erlang(rate, nphases)

    def fitMeanAndSCV(mean, scv):
        return Erlang(jpype.JPackage('jline').lang.distributions.Erlang.fitMeanAndSCV(mean, scv))
    def fitMeanAndOrder(mean, order):
        return Erlang(jpype.JPackage('jline').lang.distributions.Erlang.fitMeanAndOrder(mean, order))


class HyperExp:
    def __init__(self, *args):
        if len(args) == 1:
            self.obj = args[0]
        else:
            p = args[0]
            lambda1 = args[1]
            lambda2 = args[2]
            self.obj = jpype.JPackage('jline').lang.distributions.HyperExp(p, lambda1, lambda2)

    def fitMeanAndSCV(mean, scv):
        return HyperExp(jpype.JPackage('jline').lang.distributions.HyperExp.fitMeanAndSCV(mean, scv))


class Immediate:
    def __init__(self):
        self.obj = jpype.JPackage('jline').lang.distributions.Immediate()

class MAP:
    def __init__(self, *args):
        if len(args) == 1:
            self.obj = args[0]
        else:
            D0 = args[0]
            D1 = args[1]
            self.obj = jpype.JPackage('jline').lang.processes.MAP(jlineArrayToMatrix(D0), jlineArrayToMatrix(D1))

    def toPH(self):
        self.obj.toPH()

class PH:
    def __init__(self, *args):
        if len(args) == 1:
            self.obj = args[0]
        else:
            alpha = args[0]
            T = args[1]
            self.obj = jpype.JPackage('jline').lang.distributions.PH(jlineArrayToMatrix(alpha).toList1D(), jlineArrayToMatrix(T))


class Zipf:
    def __init__(self, *args):
        if len(args) == 1:
            self.obj = args[0]
        else:
            s = args[0]
            n = args[1]
            self.obj = jpype.JPackage('jline').lang.distributions.Zipf(s, n)
