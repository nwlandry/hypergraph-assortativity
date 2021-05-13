from GenerativeModels import *
from HyperEigenvalues import *
from MeanFieldTheory import *
import numpy as np
from math import factorial
import time
import matplotlib.pyplot as plt

n = 100
k = generatePowerLawDegreeSequence(n, 10, 100, 3)
#k = np.random.randint(low=10, high=31, size=(n))
m = 3
epsilonList = np.linspace(0, 0.3, 10)
maxIterations = 20
tolerance = 1e-5
numEdgesSorted = list()
numEdges = list()

trueCECEigenvalue = np.zeros(len(epsilonList))
meanFieldCECEigenvalue = np.zeros(len(epsilonList))
assortativity = np.zeros(len(epsilonList))

for i in range(len(epsilonList)):
    epsilon = epsilonList[i]
    hyperedgeList = createAssortativeProd(k, m, epsilon, type="large-degrees")

    assortativity[i] = getAssortativity(hyperedgeList, m)
    meanFieldCECEigenvalue[i] = getCECEigenvalue(hyperedgeList, m)

    weights = np.ones(len(hyperedgeList))
    T = SparseTensor(hyperedgeList, weights, n)
    cec = T.getCEC(maxIterations, tolerance)[0]
    trueCECEigenvalue[i] = cec
    print(i)

plt.figure()
plt.plot(assortativity, trueCECEigenvalue, 'ko-', label="CEC Eigenvalue (True)")
plt.plot(assortativity, meanFieldCECEigenvalue, 'ro-', label="CEC Eigenvalue (MF)")
plt.legend()
plt.xlabel(r"$\rho$")
plt.ylabel(r"$\lambda$")
plt.show()
