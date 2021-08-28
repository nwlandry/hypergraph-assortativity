import numpy as np
from collections import defaultdict
from math import factorial
from scipy.special import comb
from itertools import combinations
from numpy.linalg import norm


def getExpansionEigenvalue(edgeList, m):
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
    return (m - 1)*(k1/k2)**2*kk1/2 - 1

def meanTriangles(edgeList, k, m):
    numEdges = len(edgeList)
    numCombos = comb(m, 2)
    sumTriangles = sum([np.sum(np.prod(k[list(indices)])/(numEdges*numCombos) for indices in combinations(edge, 2)) for edge in edgeList])
    return sumTriangles

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

def zerothOrderEigenvector(edgeList):
    k = getDegreeSequence(edgeList)
    u = k
    return u/norm(u)

def firstOrderEigenvector(edgeList, epsilon, m, type="large-degrees"):
    k = getDegreeSequence(edgeList)

    k1 = np.mean(k)
    k2 = np.mean(np.power(k, 2))
    k3 = np.mean(np.power(k, 3))

    u0 = k
    if type == "large-degrees":
        u1 = (k2/k1**2)**m*(k1**2*k3/k2**3*np.power(k, 2) - k)

    u = u0 + epsilon*u1

    return u/norm(u)

def getEpsilon(edgeList, m, type="large-degrees"):
    
    k = getDegreeSequence(edgeList)
    gAvg = getAverageG(edgeList, m, type=type)

    if type == "large-degrees":

        k1 = np.mean(k)
        k2 = np.mean(np.power(k, 2))
        k3 = np.mean(np.power(k, 3))
        return gAvg/((k3/k1**3)**m - (k2/k1**2)**(2*m))


def getAverageG(edgeList, m, type="large-degrees"):
    k = getDegreeSequence(edgeList)
    if type == "large-degrees":
        k1 = np.mean(k)
        k2 = np.mean(np.power(k, 2))
        return np.mean([np.prod(k[list(edge)]/k1) - (k2/k1**2)**m for edge in edgeList])