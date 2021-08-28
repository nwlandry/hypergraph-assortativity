from HyperEigenvalues import *
from MeanFieldTheory import *
import numpy as np
import Hypergraph
from HyperEigenvalues import SparseTensor
import os
import shelve
import utilities
import multiprocessing as mp

def parallelRun(hyperedgeList, a, m, assortativityTolerance, maxShufflingIterations, temperature):
    hypergraph = Hypergraph.HypergraphGenerator(hyperedgeList, type="hyperedge-list")
    hypergraph.shuffleHyperedges(a, m, assortativityTolerance, maxShufflingIterations, temperature)
    hyperedgeList = hypergraph.getHyperedgeList()
    return hyperedgeList

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
# datasetFolder = "tags-ask-ubuntu"
# sizeFile = "tags-ask-ubuntu-nverts.txt"
# memberFile = "tags-ask-ubuntu-simplices.txt"

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

# datasetFolder = "Power-Law"
# # Hypergraph parameters
# n = 1000
# m = 3
# isDegreeCorrelated = True
# parameters = [{"degree-distribution":"power-law","exponent":3,"min-degree":10,"max-degree":100,"hyperedge-size":m,"size":n,"is-correlated":isDegreeCorrelated}]
# h = Hypergraph.HypergraphGenerator(parameters)

# Shuffle parameters
spacing = 0.05
originalAssortativity = h.getAssortativity(m)
changeAssortativityList = np.unique(np.concatenate([np.arange(-1 - originalAssortativity, 0, spacing), [0]]))
# changeAssortativityList = np.arange(0, 1 - originalAssortativity, spacing)

targetAssortativityList = originalAssortativity + changeAssortativityList
assortativityTolerance = 0.01
temperature = 0.00001
maxShufflingIterations = 1e6

argList = list()

for a in targetAssortativityList:
    argList.append((list(h.getHyperedgeList()), a, m, assortativityTolerance, maxShufflingIterations, temperature))

with mp.Pool(processes=numProcesses) as pool:
    listOfHyperedgeLists = pool.starmap(parallelRun, argList)

with shelve.open(os.path.join(mainFolder, dataFolder, datasetFolder, datasetFolder + "_epidemics")) as data:
    data["list-of-hyperedge-lists"] = listOfHyperedgeLists