from HyperEigenvalues import *
from MeanFieldTheory import *
import numpy as np
import Hypergraph
import copy
import os
import shelve
import utilities
import os

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
assortativityList = originalAssortativity + changeAssortativityList

assortativityTolerance = 0.001
temperature = 0.00001
maxShufflingIterations = 1e5

maxEigenvalueIterations = 1000
eigenvalueTolerance = 1e-5

assortativities = np.zeros(len(assortativityList))
meanFieldEigenvalues = np.zeros(len(assortativityList))
trueEigenvalues = np.zeros(len(assortativityList))

for i in range(len(assortativityList)):
    a = assortativityList[i]
    hNew = copy.deepcopy(h)
    hNew.shuffleHyperedges(a, m, type, tolerance=assortativityTolerance, maxIterations=maxShufflingIterations, temperature=temperature)

    hyperedgeList = hNew.getHyperedgeList()

    assortativities[i] = getAssortativity(hyperedgeList, m)
    meanFieldEigenvalues[i] = getExpansionEigenvalue(hyperedgeList, m)

    weights = np.ones(len(hyperedgeList))
    T = SparseTensor(hyperedgeList, weights, n)
    cec = T.getCEC(maxEigenvalueIterations, eigenvalueTolerance)[0]
    trueEigenvalues[i] = cec
    if a == originalAssortativity:
        originalEigenvalue = cec
    print(i)


with shelve.open(os.path.join(mainFolder, dataFolder, datasetFolder, datasetFolder + "_eigenvalues")) as data:
    data["rho"] = assortativities
    data["mean-field-eigenvalues"] = meanFieldEigenvalues
    data["true-eigenvalues"] = trueEigenvalues
    data["original-assortativity"] = originalAssortativity
    data["original-eigenvalue"] = originalEigenvalue
