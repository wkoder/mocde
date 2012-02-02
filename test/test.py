#! /usr/bin/python

'''
Created on Jan 30, 2012

@author: Moises Osorio [WCoder]
'''
from subprocess import Popen, PIPE
import os
import random
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
        dirName = "" if self.name is None else self.name + "."
        self.path = "%s/%s%s" % (MOCDETest.__RESULTS__, dirName, time.strftime("%Y%m%d-%H%M%S"))
        if not os.path.exists(self.path):
            os.makedirs(self.path)
        random.seed()
        
    def testDimension(self, dim, times=1):
        for test in MOCDETest.__TESTS__:
            if test[2] == dim:
                self.testFunction(test[0], test[1], test[2], times)
                
    def testFunctionName(self, function, times=1):
        function = function.lower()
        tests = [test for test in MOCDETest.__TESTS__ if test[0] == function or (function.endswith("*") and test[0].startswith(function[:-1]))]
        start = time.time()
        i = 0
        for test in tests:
            i += 1
            print "Testing %s (%d/%d)" % (test[0], i, len(tests))
            self.testFunction(test[0], test[1], test[2], times)
            
        end = time.time()
        timeTook = end - start
        print "Tests took %.2f seconds" % (timeTook)
                
    def testFunction(self, function, nreal, nobj, times=1):
        timeSum = 0.0
        print "----------------------------------"
        inputFile = open('test/test.in', 'r')
        inputData = inputFile.read()
        inputFile.close()
        for i in xrange(times):
            iStr = "_%d" % (i)
            if times == 1:
                iStr = "" 
            varFile = "%s/%s%s_%s" % (self.path, function, iStr, MOCDETest.__VAR_OUT__)
            objFile = "%s/%s%s_%s" % (self.path, function, iStr, MOCDETest.__OBJ_OUT__)
            cmd = [MOCDETest.__EXE__, function, "%s" % nreal, varFile, objFile, "--silent"]
            
            print "Running '%s' (%d/%d)" % (" ".join(cmd[:3]), i+1, times)
            ran = random.random()
            testInput = "%s\n%f\n" % (inputData, ran)
            start = time.time()
            
            run = Popen(cmd, stdin=PIPE)
            run.communicate(input=testInput)
            returnCode = run.wait()
            if returnCode != 0:
                error = "Received return code %d!" % (returnCode)
                print "    ERROR: %s" % (error)
                raise Exception(error)
            
            end = time.time()
            timeTook = end - start
            timeSum += timeTook
            print "    Test took %.2f seconds" % (timeTook)
            print
            
        print "Stats for %s" % (function)
        print "    Average time: %.2f seconds" % (timeSum / times)
        print "----------------------------------"


if __name__ == '__main__':
    name = raw_input("Test name: ")
    times = raw_input("Times to run: ")
    if len(times) == 0:
        times = 1
    else:
        times = int(times)
    
    functions = raw_input("Functions to test: ").split()
    test = MOCDETest(name)
    if len(functions) > 0:
        for function in functions:
            test.testFunctionName(function, times)
    else:
        dim = raw_input("Dimensions to test: ").split()
        dim = [int(d) for d in dim]
        for d in xrange(2, 101):
            if len(dim) == 0 or d in dim:
                test.testDimension(d, times)
    