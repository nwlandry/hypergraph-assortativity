import Hypergraph
from HyperEigenvalues import *
from MeanFieldTheory import *
import numpy as np
import os
import shelve

mainFolder = os.getcwd()
dataFolder = "Data"
datasetFolder = "Eigenvector"

# Hypergraph parameters
n = 1000
m = 3
isDegreeCorrelated = True
parameters = [{"degree-distribution":"power-law","exponent":3,"min-degree":10,"max-degree":100,"hyperedge-size":m,"size":n,"is-correlated":isDegreeCorrelated}]

h = Hypergraph.HypergraphGenerator(parameters)

a = 0.2
assortativityTolerance = 0.00001
temperature = 0.00001
maxShufflingIterations = 1e5

maxEigenvalueIterations = 1000
eigenvalueTolerance = 1e-5

h.shuffleHyperedges(a, m, tolerance=assortativityTolerance, maxIterations=maxShufflingIterations, temperature=temperature)
hyperedgeList = h.getHyperedgeList()

if len(hyperedgeList) > 0:

    weights = np.ones(len(hyperedgeList))
    T = SparseTensor(hyperedgeList, weights, n)
    trueEigenvector = T.getCEC(maxEigenvalueIterations, eigenvalueTolerance)[1]

    epsilon = getEpsilon(hyperedgeList, m)

    print(epsilon)

    zerothOrder = zerothOrderEigenvector(hyperedgeList)
    firstOrder = firstOrderEigenvector(hyperedgeList, epsilon, m)

with shelve.open(os.path.join(mainFolder, dataFolder, datasetFolder, datasetFolder)) as data:
    data["zeroth-order-eigenvector"] = zerothOrder
    data["first-order-eigenvector"] = firstOrder
    data["actual-eigenvector"] = trueEigenvector
