from GenerativeModels import *
from HyperEigenvalues import *
from MeanFieldTheory import *
import numpy as np
import shelve

mainFolder = os.getcwd()
dataFolder = "Data"
datasetFolder = "EpsilonRho"

n = 100
kmin = 10
kmax = 100
r = 3
m = 3

numSims = 100
epsilonList = np.linspace(-0.1, 0.6, 20)
numEdgesSorted = list()
numEdges = list()

assortativities = list()
predictedAssortativities = list()
epsilons = list()

for i in range(len(epsilonList)):
    for sim in range(numSims):
        epsilon = epsilonList[i]
        try:
            k = generatePowerLawDegreeSequence(n, kmin, kmax, r)
            hyperedgeList = createAssortativeProd(k, m, epsilon, type="large-degrees")
            if len(hyperedgeList) > 0:
                assortativities.append(getAssortativity(hyperedgeList, m))
                predictedAssortativities.append(epsilonToRho(epsilon, hyperedgeList, m))
                epsilons.append(epsilon)
        except:
            pass
    print(i)

with shelve.open(os.path.join(mainFolder, dataFolder, datasetFolder, datasetFolder + "_eigenvalues")) as data:
    data["rho"] = assortativities
    data["predicted-rho"] = predictedAssortativities
    data["epsilon"] = epsilons