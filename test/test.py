'''
Created on Jan 30, 2012

@author: Moises Osorio [WCoder]
'''
import subprocess
import time

class MOCDETest():
    
    __EXE__ = "./bin/mocderunner"
    __VAR_OUT__ = "var.out"
    __OBJ_OUT__ = "front.out"
    __RESULTS__ = "results/history"
    __TESTS__ = [
                 ["deb2", 2, 2],
                 ["deb3", 2, 2],
                 ["fonseca2", 5, 2],
                 ["kursawe", 3, 2],
                 ["zdt1", 30, 2],
                 ["zdt2", 30, 2],
                 ["zdt3", 30, 2],
                 #["wfg1", 30, 2],
                 #["wfg2", 30, 2],
                 #["wfg3", 30, 2],
                 ["uf1", 30, 2],
                 ["uf2", 30, 2],
                 ["uf3", 30, 2],
                 ["uf4", 30, 2],
                 ["uf5", 30, 2],
                 ["uf6", 30, 2],
                 ["uf7", 30, 2],
                 ["dtlz1", 12, 3],
                 ["dtlz2", 12, 3],
                 #["r_dtlz2", 12, 3],
                 ["dtlz3", 12, 3],
                 #["dtlz5im", 12, 3],
                 ["dtlz7", 12, 3],
                 ["uf8", 30, 3],
                 ["uf9", 30, 3],
                 ["uf10", 30, 3],
                 ]
    
    def __init__(self, name):
        self.name = name
        self.path = "%s%s.%s" % (MOCDETest.__RESULTS__, self.name, time.strftime("%Y%m%d-%H%M%S"))
        
    def testDimension(self, dim):
        for test in MOCDETest.__TESTS__:
            if test[2] == dim:
                self.testFunction(test[0], test[1], test[2])
                
    def testFunction(self, function, nreal, nobj):
        varFile = "%s/%s_%s" % (self.path, function, MOCDETest.__VAR_OUT__)
        objFile = "%s/%s_%s" % (self.path, function, MOCDETest.__OBJ_OUT__)
        cmd = [MOCDETest.__EXE__, function, "%s" % nreal, varFile, objFile, "--silent"]
        
        print "----------------------------------"
        print "Running '%s'" % (cmd[:3].join(" "))
        start = time.time()
        returnCode = subprocess.call(cmd, shell=True)
        if returnCode != 0:
            print "    ERROR: Received return code %d!" % (returnCode)
            raise "Received return code %d!" % (returnCode)
        end = time.time()
        print "Test took %.2f seconds" % (end - start)
        print "----------------------------------"


if __name__ == '__main__':
    testName = input("Test name: ")
    dim = input("Dimensions to test: ").split()
    dim = [int(d) for d in dim]
    
    test = MOCDETest(testName)
    for d in xrange(2, 100):
        if len(dim) == 0 or d in dim:
            test.testDimension(dim)
    