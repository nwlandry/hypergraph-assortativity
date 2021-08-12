import Hypergraph
import numpy as np
import HypergraphContagion
import MeanFieldTheory
import os
import shelve
from HyperEigenvalues import SparseTensor
import multiprocessing as mp


def parallelRun(hyperedgeList, gamma, beta, initialFractionInfected, tmax, fractionToAverage, numRuns, m, maxEigenvalueIterations, eigenvalueTolerance):
    H = Hypergraph.Hypergraph(hyperedgeList)
    assortativity = MeanFieldTheory.getAssortativity(hyperedgeList, m)

    equilibrium = HypergraphContagion.getEquilibriumPoint(H, gamma, beta, initialFractionInfected, tmax, fractionToAverage, numRuns, True)
    weights = np.ones(len(hyperedgeList))
    T = SparseTensor(hyperedgeList, weights, n)
    eigenvalue = T.getCEC(maxEigenvalueIterations, eigenvalueTolerance)[0]
    return equilibrium, assortativity, eigenvalue




# filename = sys.argv[1]
mainFolder = os.getcwd()
# Import the hypergraph
dataFolder = "Data"


datasetFolder = "Eu-Emails"
datasetFolder = "congress-bills"
# datasetFolder = "tags-ask-ubuntu"
# datasetFolder = "email-Enron"



filename = os.path.join(mainFolder, dataFolder, datasetFolder, datasetFolder + "_epidemics")

# Epidemic parameters
gamma = 1
tmax = 300
fractionToAverage = 0.1
isVerbose = True
numProcesses = len(os.sched_getaffinity(0))
numRuns = 5
initialFractionInfected = 1.0
m = 3

maxEigenvalueIterations = 1000
eigenvalueTolerance = 1e-5

with shelve.open(filename) as data:
    listOfHyperedgeLists = data["list-of-hyperedge-lists"]

hyperedgeList = listOfHyperedgeLists[-1]

n = max([max(index) for index in hyperedgeList]) + 1

weights = np.ones(len(hyperedgeList))
T = SparseTensor(hyperedgeList, weights, n)
spectralRadius = T.getCEC(maxEigenvalueIterations, eigenvalueTolerance)[0]
print(spectralRadius)
beta3Crit = gamma/spectralRadius

argList = list()

for i in range(len(listOfHyperedgeLists)):
    hyperedgeList = listOfHyperedgeLists[i]
    beta3 = 2*beta3Crit
    beta = {3 : beta3}
    argList.append((hyperedgeList, gamma, beta, initialFractionInfected, tmax, fractionToAverage, numRuns, m, maxEigenvalueIterations, eigenvalueTolerance))

with mp.Pool(processes=numProcesses) as pool:
    data = pool.starmap(parallelRun, argList)

equilibria, assortativities, eigenvalues = map(list, zip(*data))

with shelve.open(filename) as data:
    data["gamma"] = gamma
    data["beta3"] = beta3
    data["equilibria"] = equilibria
    data["assortativities"] = assortativities
    data["eigenvalues"] = eigenvalues

