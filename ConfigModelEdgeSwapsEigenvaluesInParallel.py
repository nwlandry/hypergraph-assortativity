from GenerativeModels import *
from HyperEigenvalues import *
from MeanFieldTheory import *
import numpy as np
import Hypergraph
import copy
import os
import shelve
import multiprocessing as mp

def parallelRun(hyperedgeList, a, m, assortativityTolerance, maxShufflingIterations, temperature, maxEigenvalueIterations, eigenvalueTolerance):
    hypergraph = Hypergraph.HypergraphGenerator(hyperedgeList, type="hyperedge-list")
    hypergraph.shuffleHyperedges(a, m, assortativityTolerance, maxShufflingIterations, temperature)
    hyperedgeList = hypergraph.getHyperedgeList()

    assortativity = getAssortativity(hyperedgeList, m)
    meanFieldCECEigenvalue = getCECEigenvalue(hyperedgeList, m)

    weights = np.ones(len(hyperedgeList))
    T = SparseTensor(hyperedgeList, weights, n)
    cec = T.getCEC(maxEigenvalueIterations, eigenvalueTolerance)[0]
    trueCECEigenvalue = cec
    return assortativity, meanFieldCECEigenvalue, trueCECEigenvalue

numProcesses = len(os.sched_getaffinity(0))

mainFolder = os.getcwd()


dataFolder = "Data"
datasetFolder = "Power-Law"

# Hypergraph parameters
n = 10000
m = 3
isDegreeCorrelated = True
parameters = [{"degree-distribution":"power-law","exponent":3,"min-degree":10,"max-degree":100,"hyperedge-size":m,"size":n,"is-correlated":isDegreeCorrelated}]

h = Hypergraph.HypergraphGenerator(parameters)

# Shuffle parameters
spacing = 0.05
originalAssortativity = h.getAssortativity(m)
changeAssortativityList = np.unique(np.concatenate([np.arange(-1 - originalAssortativity, 0, spacing), np.arange(0, 1 - originalAssortativity, spacing)]))
targetAssortativityList = originalAssortativity + changeAssortativityList
assortativityTolerance = 0.01
temperature = 0.00001
maxShufflingIterations = 1e6

maxEigenvalueIterations = 1000
eigenvalueTolerance = 1e-5

assortativities = np.zeros(len(targetAssortativityList))
meanFieldCECEigenvalues = np.zeros(len(targetAssortativityList))
trueCECEigenvalues = np.zeros(len(targetAssortativityList))

argList = list()

for a in targetAssortativityList:
    argList.append((list(h.getHyperedgeList()), a, m, assortativityTolerance, maxShufflingIterations, temperature, maxEigenvalueIterations, eigenvalueTolerance))

with mp.Pool(processes=numProcesses) as pool:
    data = pool.starmap(parallelRun, argList)

for i in range(len(data)):
    assortativities[i] = data[i][0]
    meanFieldCECEigenvalues[i] = data[i][1]
    trueCECEigenvalues[i] = data[i][2]

    if abs(originalAssortativity - assortativities[i]) < 0.0001:
        print("There is an assortativity which matches the original assortativity")
        originalEigenvalue = trueCECEigenvalues[i]

with shelve.open(os.path.join(mainFolder, dataFolder, datasetFolder, datasetFolder + "_eigenvalues")) as data:
    data["rho"] = assortativities
    data["mean-field-eigenvalues"] = meanFieldCECEigenvalues
    data["true-eigenvalues"] = trueCECEigenvalues
    data["original-assortativity"] = originalAssortativity
    data["original-eigenvalue"] = originalEigenvalue