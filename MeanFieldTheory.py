import numpy as np
from collections import defaultdict
from math import factorial
from scipy.special import comb
from itertools import combinations

def getCECEigenvalue(edgeList, m):
    k = getDegreeSequence(edgeList)
    rho = getAssortativity(edgeList, m)
    k1 = np.mean(k)
    k2 = np.mean(np.power(k, 2))
    l0 = (m-1)*k2/k1
    return l0*(1 + rho)

def getAssortativity(edgeList, m):
    k = getDegreeSequence(edgeList)
    k1 = np.mean(k)
    k2 = np.mean(np.power(k, 2))
    kk1 = meanTriangles(edgeList, k, m)
    return k1**2*kk1/k2**2 - 1

def meanTriangles(edgeList, k, m):
    numEdges = len(edgeList)
    numCombos = comb(m, 2)
    sumTriangles = sum([np.sum(np.prod(k[list(indices)]) for indices in combinations(edge, 2)) for edge in edgeList])
    return sumTriangles/(numEdges*numCombos)

def getDegreeSequence(edgeList):
    # edgeList must be the same size
    d = defaultdict(lambda : 0)
    for edge in edgeList:
        for node in edge:
            d[node] += 1
    degreeSequence = np.zeros(max(d.keys()) + 1)
    for node, degree in d.items():
        degreeSequence[node] = degree
    return degreeSequence

def epsilonToRho(epsilon, edgeList, m):
    k = getDegreeSequence(edgeList)
    k1 = np.mean(k)
    k2 = np.mean(np.power(k, 2))
    k3 = np.mean(np.power(k, 3))
    return (k2/k1**2)**m*(k1**2*k3**2/k2**4 - 1)*epsilon