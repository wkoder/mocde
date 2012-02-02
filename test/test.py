#! /usr/bin/python

'''
Created on Jan 30, 2012

@author: Moises Osorio [WCoder]
'''
from subprocess import Popen, PIPE
import os
import random
import sqlite3
import sys
import time

class MOCDETest():
    
    __EXE__ = "./bin/mocderunner"
    __VAR_OUT__ = "var.out"
    __OBJ_OUT__ = "front.out"
    __RESULTS__ = "results/history"
    __DB__ = "test_runs.db"
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
        self._initDB()
        
    def _initDB(self):
        self.db = sqlite3.connect(MOCDETest.__DB__)
        cur = self.db.cursor()
        try:
            cur.execute("CREATE TABLE runs (id INT, function TEXT, running_time REAL, start_time DATETIME)")
            self.db.commit()
        except:
            pass
        
    def _addDBRun(self, function, runningTime, startTime):
        cur = self.db.cursor()
        cur.execute("INSERT INTO runs (function, running_time, start_time) VALUES ('%s', %f, '%s')" % (function, runningTime, startTime))
        self.db.commit()
        
    def _clearDBRuns(self, function):
        cur = self.db.cursor()
        cur.execute("DELETE FROM runs WHERE function LIKE '%s'" % (function))
        self.db.commit()
        
    def _calculateETC(self, function):
        cur = self.db.cursor()
        cur.execute("SELECT AVG(running_time) FROM runs WHERE function LIKE '%s' GROUP BY function" % (function))
        data = cur.fetchone()
        return None if data is None else data[0]
        
    def testDimension(self, dim, times=1):
        for test in MOCDETest.__TESTS__:
            if test[2] == dim:
                self.testFunction(test[0], test[1], test[2], times)
                
    def _calculateTestsETC(self, functions):
        return map(self._calculateETC, functions)
        
    def testFunctionName(self, function, times=1):
        function = function.lower()
        tests = [test for test in MOCDETest.__TESTS__ if test[0] == function or (function.endswith("*") and test[0].startswith(function[:-1]))]
        startSecond = time.time()
        etcAll = self._calculateTestsETC([test[0] for test in tests])
        i = 0
        for test in tests:
            self._printBars()
            etcRemaining = etcAll[i:]
            if None in etcRemaining:
                etcRemainingSum = None
            else:
                etcRemainingSum = sum(etcRemaining)
            print "Testing %s (%d/%d), test ETC: %s s, remaining tests ETC: %s s" % (test[0], i+1, len(tests), \
                         self._strETC(etcAll[i], times), self._strETC(etcRemainingSum, times))
            self.testFunction(test[0], test[1], test[2], times)
            self._printBars()
            i += 1
            
        endSecond = time.time()
        runningTime = endSecond - startSecond
        print "Tests took %.2f s" % (runningTime)
    
    def _printBars(self, bars=1):
        for _ in xrange(bars):
            print "-" * 89
            
    def _strETC(self, etc, times=1):
        return "N/A" if etc is None else "%.2f" % (etc * times)
        
    def testFunction(self, function, nreal, nobj, times=1):
        timeSum = 0.0
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
            
            etc = self._calculateETC(function)
            print
            print "    Running '%s' (%d/%d), run ETC: %s s, remaining runs ETC: %s s" % (" ".join(cmd[:3]), i+1, times, self._strETC(etc), self._strETC(etc, times-i))
            ran = random.random()
            testInput = "%s\n%f\n" % (inputData, ran)
            startTime = time.strftime("%Y-%m-%d %H:%M:%S")
            startSecond = time.time()
            
            run = Popen(cmd, stdin=PIPE)
            run.communicate(input=testInput)
            returnCode = run.wait()
            if returnCode != 0:
                error = "Received return code %d!" % (returnCode)
                print "    ERROR: %s" % (error)
                raise Exception(error)
            
            endSecond = time.time()
            runningTime = endSecond - startSecond
            timeSum += runningTime
            if i == 0:
                self._clearDBRuns(function)
            self._addDBRun(function, runningTime, startTime)
            print "        Test took %.2f s" % (runningTime)

def _getDimension(arr):
    if len(arr) == 1 and arr[0].endswith("d"):
        try:
            return int(arr[0][:-1])
        except:
            pass
    return 0

if __name__ == '__main__':
    if len(sys.argv) > 1:
        name = sys.argv[1]
    else:
        name = raw_input("Test name: ")
    if len(sys.argv) > 2:
        times = int(sys.argv[2])
    else:
        times = int(raw_input("Times to run: "))
    if len(sys.argv) > 3:
        functions = sys.argv[3].split()
    else:
        functions = raw_input("Functions to test: ").split()
    test = MOCDETest(name)
    dimension = _getDimension(functions)
    if not dimension:
        for function in functions:
            test.testFunctionName(function, times)
    else:
        test.testDimension(dimension, times)
    