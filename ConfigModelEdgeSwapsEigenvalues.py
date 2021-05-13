from GenerativeModels import *
from HyperEigenvalues import *
from MeanFieldTheory import *
import numpy as np
from math import factorial
import Hypergraph
import copy
import os
import shelve

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
assortativityList = np.linspace(-0.5, 0.5, 11)
assortativityTolerance = 0.01
temperature = 0.00001
maxShufflingIterations = 1e6

type = "large-degrees"
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
