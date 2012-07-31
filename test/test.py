#! /usr/bin/python

'''
Created on Jan 30, 2012

@author: Moises Osorio [WCoder]
'''
from subprocess import Popen, PIPE
import os
import random
import sqlite3
import time
import argparse

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
                 ["zdt4", 10, 2],
                 ["zdt6", 10, 2],
                 ["uf1", 30, 2],
                 ["uf2", 30, 2],
                 ["uf3", 30, 2],
                 ["uf4", 30, 2],
                 ["uf5", 30, 2],
                 ["uf6", 30, 2],
                 ["uf7", 30, 2],
                 ["uf8", 30, 3],
                 ["uf9", 30, 3],
                 ["uf10", 30, 3],
                 ["dtlz1", 12, 3],
                 ["dtlz2", 12, 3],
                 ["dtlz3", 12, 3],
                 ["dtlz4", 12, 3],
                 ["dtlz5", 12, 3],
                 ["dtlz6", 12, 3],
                 ["dtlz7", 22, 3],
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
        
    def _calculateTestsETC(self, functions):
        return map(self._calculateETC, functions)
        
    def testFunctions(self, exe, functions, times, populationSize, maxEvaluations, mutationRate, crossoverProbability, elitism):
        tests = []
        for test in MOCDETest.__TESTS__:
            if True in [self.testMatches(function, test[0], test[2]) for function in functions]:
                tests.append(test)
                
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
            print "Testing %s (%d/%d), test ETC: %s, remaining tests ETC: %s" % (test[0], i+1, len(tests), \
                         self._strETC(etcAll[i], times), self._strETC(etcRemainingSum, times))
            self.testFunction(exe, test[0], test[1], test[2], times, populationSize, maxEvaluations, mutationRate, crossoverProbability, elitism)
            self._printBars()
            i += 1
            
        endSecond = time.time()
        runningTime = endSecond - startSecond
        print "Tests took %s" % (self._strETC(runningTime))
    
    def _printBars(self, bars=1):
        for _ in xrange(bars):
            print "-" * 89
            
    def _strETC(self, etc, times=1):
        if etc is None:
            return "N/A"
        units = "s"
        etc = etc * times
        if etc >= 60:
            etc = etc / 60
            units = "m"
        if etc >= 60:
            etc = etc / 60
            units = "h"
        return "%.2f %s" % (etc, units)
        
    def testFunction(self, exe, function, nreal, nobj, times, populationSize, maxEvaluations, mutationRate, crossoverProbability, elitism):
        timeSum = 0.0
        for i in xrange(times):
            iStr = "_%d" % (i)
            if times == 1:
                iStr = "" 
            filePrefix = "%s/%s%s" % (self.path, function, iStr)
            cmd = [exe, function, "%s" % nreal, filePrefix, "--silent"]
            
            etc = self._calculateETC(function)
            print
            print "    Running '%s' (%d/%d), run ETC: %s, remaining runs ETC: %s" % (" ".join(cmd[:3]), i+1, times, self._strETC(etc), self._strETC(etc, times-i))
            ran = random.random()
            testInput = "%d\n%d\n%f\n%f\n" % (populationSize, maxEvaluations, mutationRate, crossoverProbability)
            if elitism:
                testInput = "%s%d\n" % (testInput, elitism)
            testInput = "%s%f\n" % (testInput, ran)
            
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
            print "        Test took %s" % (self._strETC(runningTime))

    def testMatches(self, desc, testName, testDim):
        desc = desc.lower()
        testName = testName.lower()
        if desc == testName:
            return True
        if desc.endswith("*") and testName.startswith(desc[:-1]):
            return True
        try:
            if desc.endswith("d") and testDim == int(desc[:-1]):
                return True
        except:
            None
        return False
    
if __name__ == '__main__':
    parser = argparse.ArgumentParser(description="Runs mocDE against several test functions.")
    parser.add_argument("functions", metavar="FUNCTION", nargs="+", \
                        help="function to test (matches '*' as wildcard and 'd' as postfix for dimension")
    parser.add_argument("--name", "-n", help="run name")
    parser.add_argument("--runs", "-r", type=int, help="number of times to run the tests")
    parser.add_argument("--population", "-p", type=int, help="population size")
    parser.add_argument("--evaluations", "-e", type=int, help="maximum number of evaluations")
    parser.add_argument("--mutation", "-m", "--differential", "-d", type=float,  \
                        help="mutation rate / differential variation")
    parser.add_argument("--crossover", "-c", type=float, help="crossover probability")
    parser.add_argument("--elitism", "-l", type=int, help="maximum number of iterations to be survived by the elite")
    parser.add_argument("--exe", nargs="?", default=MOCDETest.__EXE__, help="executable to test")
    args = parser.parse_args()
    
    test = MOCDETest(args.name)
    test.testFunctions(args.exe, args.functions, args.runs, args.population, args.evaluations, args.mutation, args.crossover, 
                       args.elitism)
    