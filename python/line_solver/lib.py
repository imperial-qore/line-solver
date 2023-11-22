import jpype
import jpype.imports

class SMCSolver:
    def __init__(self):
        pass
    def qbd_caudal(self, A0, A1, A2):
        return jpype.JPackage('jline').lib.thirdparty.SMCSolver(A0,A1,A2)