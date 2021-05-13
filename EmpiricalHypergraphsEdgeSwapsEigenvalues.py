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

m = 3

mainFolder = os.getcwd()
# Import the hypergraph
dataFolder = "Data"

# datasetFolder = "Eu-Emails"
# sizeFile = "email-Eu-nverts.txt"
# memberFile = "email-Eu-simplices.txt"

datasetFolder = "congress-bills"
sizeFile = "congress-bills-nverts.txt"
memberFile = "congress-bills-simplices.txt"
# #
# datasetFolder = "tags-ask-ubuntu"
# sizeFile = "tags-ask-ubuntu-nverts.txt"
# memberFile = "tags-ask-ubuntu-simplices.txt"

# datasetFolder = "email-Enron"
# sizeFile = "email-Enron-nverts.txt"
# memberFile = "email-Enron-simplices.txt"
#
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
type = "large-degrees"
originalAssortativity = h.getAssortativity(m, type)

changeAssortativityList = np.unique(np.concatenate([np.linspace(-0.5, 0, 11), np.linspace(0, 0.5, 11)]))
assortativityList = originalAssortativity + changeAssortativityList
assortativityTolerance = 0.01
temperature = 0.00001
maxShufflingIterations = 1e5

maxEigenvalueIterations = 1000
eigenvalueTolerance = 1e-5

assortativities = np.zeros(len(assortativityList))
meanFieldCECEigenvalues = np.zeros(len(assortativityList))
trueCECEigenvalues = np.zeros(len(assortativityList))

for i in range(len(assortativityList)):
    a = assortativityList[i]
    hNew = copy.deepcopy(h)
    hNew.shuffleHyperedges(a, m, type, tolerance=assortativityTolerance, maxIterations=maxShufflingIterations, temperature=temperature)

    hyperedgeList = hNew.getHyperedgeList()

    assortativities[i] = getAssortativity(hyperedgeList, m)
    meanFieldCECEigenvalues[i] = getCECEigenvalue(hyperedgeList, m)

    weights = np.ones(len(hyperedgeList))
    T = SparseTensor(hyperedgeList, weights, n)
    cec = T.getCEC(maxEigenvalueIterations, eigenvalueTolerance)[0]
    trueCECEigenvalues[i] = cec
    print(i)


with shelve.open(os.path.join(mainFolder, dataFolder, datasetFolder, datasetFolder + "_eigenvalues")) as data:
    data["rho"] = assortativities
    data["mean-field-eigenvalues"] = meanFieldCECEigenvalues
    data["true-eigenvalues"] = trueCECEigenvalues
