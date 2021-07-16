from GenerativeModels import *
from HyperEigenvalues import *
from MeanFieldTheory import *
import numpy as np
import shelve
import os
import multiprocessing as mp

def parallelRun(k, m, epsilon):
    hyperedgeList = createAssortativeProd(k, m, epsilon, type="large-degrees")
    if len(hyperedgeList) > 0:
        assortativity = getAssortativity(hyperedgeList, m)
        predictedAssortativity = epsilonToRho(epsilon, hyperedgeList, m)
        return epsilon, assortativity, predictedAssortativity
    else:
        return np.NaN, np.NaN, np.NaN

numProcesses = len(os.sched_getaffinity(0))

mainFolder = os.getcwd()
dataFolder = "Data"
datasetFolder = "EpsilonRho"

n = 100
kmin = 10
kmax = 100
r = 3
m = 3

numSims = 10
epsilonList = np.linspace(-0.1, 0.6, 100)
numEdgesSorted = list()
numEdges = list()

argList = list()

for i in range(len(epsilonList)):
    for sim in range(numSims):
        epsilon = epsilonList[i]
        k = generatePowerLawDegreeSequence(n, kmin, kmax, r)
        argList.append((k, m, epsilon))

with mp.Pool(processes=numProcesses) as pool:
    data = pool.starmap(parallelRun, argList)

epsilons =  np.zeros(len(epsilonList)*numSims) 
assortativities = np.zeros(len(epsilonList)*numSims)
predictedAssortativities =  np.zeros(len(epsilonList)*numSims)


for i in range(len(data)):
    epsilons[i] = data[i][0]
    assortativities[i] = data[i][1]
    predictedAssortativities[i] = data[i][2]

with shelve.open(os.path.join(mainFolder, dataFolder, datasetFolder, datasetFolder + "_generative")) as data:
    data["rho"] = assortativities
    data["predicted-rho"] = predictedAssortativities
    data["epsilon"] = epsilons