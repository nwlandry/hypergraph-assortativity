from GenerativeModels import *
from HyperEigenvalues import *
from MeanFieldTheory import *
import numpy as np
from math import factorial
import time
import shelve

mainFolder = os.getcwd()
dataFolder = "Data"
datasetFolder = "Large-Degrees"

n = 100
kmin = 10
kmax = 100
r = 3
#k = np.random.randint(low=10, high=31, size=(n))
m = 3

numSims = 100
epsilonList = np.linspace(-0.1, 0.6, 20)
maxIterations = 20
tolerance = 1e-5
numEdgesSorted = list()
numEdges = list()

trueCECEigenvalues = list()
meanFieldCECEigenvalues = list()
assortativities = list()
epsilons = list()

for i in range(len(epsilonList)):
    for sim in range(numSims):
        epsilon = epsilonList[i]
        try:
            k = generatePowerLawDegreeSequence(n, kmin, kmax, r)
            hyperedgeList = createAssortativeProd(k, m, epsilon, type="large-degrees")
            if len(hyperedgeList) > 0:
                assortativities.append(getAssortativity(hyperedgeList, m))
                meanFieldCECEigenvalues.append(getCECEigenvalue(hyperedgeList, m))

                weights = np.ones(len(hyperedgeList))
                T = SparseTensor(hyperedgeList, weights, n)
                cec = T.getCEC(maxIterations, tolerance)[0]
                trueCECEigenvalues.append(cec)
                epsilons.append(epsilon)
        except:
            pass
    print(i)

with shelve.open(os.path.join(mainFolder, dataFolder, datasetFolder, datasetFolder + "_eigenvalues")) as data:
    data["rho"] = assortativities
    data["epsilon"] = epsilons
    data["mean-field-eigenvalues"] = meanFieldCECEigenvalues
    data["true-eigenvalues"] = trueCECEigenvalues