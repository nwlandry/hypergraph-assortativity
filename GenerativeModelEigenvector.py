from GenerativeModels import *
from HyperEigenvalues import *
from MeanFieldTheory import *
import numpy as np
import matplotlib.pyplot as plt

mainFolder = os.getcwd()
dataFolder = "Data"
datasetFolder = "Large-Degrees"

n = 100
epsilon = 0.7
kmin = 10
kmax = 30
r = 3
m = 3

maxEigenvalueIterations = 1000
eigenvalueTolerance = 1e-5

k = np.random.randint(kmin, kmax, n)#generatePowerLawDegreeSequence(n, kmin, kmax, r)

hyperedgeList = createAssortativeProd(k, m, epsilon, type="large-degrees")

if len(hyperedgeList) > 0:

    weights = np.ones(len(hyperedgeList))
    T = SparseTensor(hyperedgeList, weights, n)
    trueCECEigenvector = T.getCEC(maxEigenvalueIterations, eigenvalueTolerance)[1]

    epsilon = getEpsilon(hyperedgeList, m)

    print(epsilon)

    zerothOrderCECEigenvector = zerothOrderEigenvector(hyperedgeList)
    firstOrderCECEigenvector = firstOrderEigenvector(hyperedgeList, epsilon, m)

    dataMin = min(np.min(trueCECEigenvector), np.min(zerothOrderCECEigenvector))
    dataMax = max(np.max(trueCECEigenvector), np.max(zerothOrderCECEigenvector))

    plt.figure()
    plt.scatter(zerothOrderCECEigenvector, trueCECEigenvector, color="black", marker="x", label="0th order approximation")
    plt.scatter(firstOrderCECEigenvector, trueCECEigenvector, color="black", marker="o", label="1st order approximation")
    plt.plot([dataMin, dataMax], [dataMin, dataMax], 'k--')
    plt.xlabel(r"$u^{(MF)}_i$", fontsize=18)
    plt.ylabel(r"$u_i$", fontsize=18)
    plt.legend(loc="upper left")
    plt.savefig("Figures/eigenvector-approximation.png", dpi=1000)
    plt.savefig("Figures/eigenvector-approximation.pdf", dpi=1000)
    plt.show()