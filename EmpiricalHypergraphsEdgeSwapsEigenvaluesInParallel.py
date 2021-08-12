from GenerativeModels import *
from HyperEigenvalues import *
from MeanFieldTheory import *
import numpy as np
from math import factorial
import Hypergraph
import copy
import os
import shelve
import utilities
import os
import multiprocessing as mp

def parallelRun(hyperedgeList, a, m, assortativityTolerance, maxShufflingIterations, temperature, maxEigenvalueIterations, eigenvalueTolerance):
    hypergraph = Hypergraph.HypergraphGenerator(hyperedgeList, type="hyperedge-list")
    hypergraph.shuffleHyperedges(a, m, assortativityTolerance, maxShufflingIterations, temperature)
    hyperedgeList = hypergraph.getHyperedgeList()

    assortativity = getAssortativity(hyperedgeList, m)
    meanFieldEigenvalue = getCliqueExpansionEigenvalue(hyperedgeList, m)

    weights = np.ones(len(hyperedgeList))
    T = SparseTensor(hyperedgeList, weights, n)
    cec = T.getCEC(maxEigenvalueIterations, eigenvalueTolerance)[0]
    trueEigenvalue = cec
    print(a, flush=True)
    return assortativity, meanFieldEigenvalue, trueEigenvalue

numProcesses = len(os.sched_getaffinity(0))

m = 3

mainFolder = os.getcwd()
# Import the hypergraph
dataFolder = "Data"

# datasetFolder = "Eu-Emails"
# sizeFile = "email-Eu-nverts.txt"
# memberFile = "email-Eu-simplices.txt"

# datasetFolder = "congress-bills"
# sizeFile = "congress-bills-nverts.txt"
# memberFile = "congress-bills-simplices.txt"
# #
datasetFolder = "tags-ask-ubuntu"
sizeFile = "tags-ask-ubuntu-nverts.txt"
memberFile = "tags-ask-ubuntu-simplices.txt"

# datasetFolder = "email-Enron"
# sizeFile = "email-Enron-nverts.txt"
# memberFile = "email-Enron-simplices.txt"

# datasetFolder = "tags-math-sx"
# sizeFile = "tags-math-sx-nverts.txt"
# memberFile = "tags-math-sx-simplices.txt"


hyperedgeSizeFile = os.path.join(dataFolder, datasetFolder, sizeFile)
memberIDFile = os.path.join(dataFolder, datasetFolder, memberFile)

hyperedgeSizes = [m]
hyperedgeList = utilities.readScHoLPData(hyperedgeSizeFile, memberIDFile)
hyperedgeList = utilities.filterHyperedgesBySize(hyperedgeList, hyperedgeSizes)
h = Hypergraph.HypergraphGenerator(hyperedgeList, type="hyperedge-list")

n = max([max(index) for index in hyperedgeList]) + 1

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
meanFieldEigenvalues = np.zeros(len(targetAssortativityList))
trueEigenvalues = np.zeros(len(targetAssortativityList))

argList = list()

for a in targetAssortativityList:
    argList.append((list(h.getHyperedgeList()), a, m, assortativityTolerance, maxShufflingIterations, temperature, maxEigenvalueIterations, eigenvalueTolerance))

with mp.Pool(processes=numProcesses) as pool:
    data = pool.starmap(parallelRun, argList)

for i in range(len(data)):
    assortativities[i] = data[i][0]
    meanFieldEigenvalues[i] = data[i][1]
    trueEigenvalues[i] = data[i][2]

    if abs(originalAssortativity - assortativities[i]) < 0.0001:
        originalEigenvalue = trueEigenvalues[i]


with shelve.open(os.path.join(mainFolder, dataFolder, datasetFolder, datasetFolder + "_eigenvalues")) as data:
    data["rho"] = assortativities
    data["mean-field-eigenvalues"] = meanFieldEigenvalues
    data["true-eigenvalues"] = trueEigenvalues
    data["original-assortativity"] = originalAssortativity
    data["original-eigenvalue"] = originalEigenvalue
