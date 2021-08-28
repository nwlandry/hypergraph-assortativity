import math
import random
import numpy as np
import copy
import itertools
from itertools import combinations, product
from scipy.special import comb
import os
import Hypergraph
import shelve
import time
import matplotlib.pyplot as plt

def createChungLuProd(k, m):
    sortedIndices = np.argsort(k)
    k = k[sortedIndices]
    n = len(k)
    S = np.sum(k)

    p = min(S/m*(k[-1]/S)**m, 1)

    maxM = n**m
    majorIndex = n-1
    edgeList = list()
    index = maxM - np.random.geometric(p)
    while index >= 0:
        edge = getEdgeProd(index, n, m)
        if edge[0] < majorIndex:
            majorIndex = edge[0]
            p = min(k[majorIndex]/m*(k[-1]/S)**(m-1), 1)

        q = min(np.prod(k[edge]/S)*S/m, 1)

        if q >= p:
            edgeList.append(tuple(sortedIndices[edge]))
        elif random.random() <= q/p:
            edgeList.append(tuple(sortedIndices[edge]))
        index -= np.random.geometric(p)
    return edgeList

def createAssortativeProd(k, m, epsilon, type="large-degrees"):
    sortedIndices = np.argsort(k)
    k = k[sortedIndices]
    n = len(k)
    volume = np.sum(k)
    k1 = np.mean(k)
    k2 = np.mean(np.power(k, 2))

    if type == "large-degrees":
        h = lambda k: np.prod(k/k1) - (k2/k1**2)**m
    elif type == "aligned-degrees":
        h = lambda k: np.sum([(combo[0] - k1)*(combo[1] - k1) for combo in combinations(k, 2)])/(comb(len(k), 2, exact=True)*k1**2) - (k2 - k1**2)**2/k1**4

    if getMaxH(k[0], k[-1], m, epsilon, h, type) < 0 or getMinH(k[0], k[-1], m, epsilon, h, type) < 0:
        print("Invalid assortativity value")
        return list()

    p = getMaxAssortativeP(k[-1], k[0], k[-1], m, volume, epsilon, h, type)

    maxM = n**m
    majorIndex = n-1
    edgeList = list()
    index = maxM - np.random.geometric(p)# - 1 # -1 b/c zero indexing
    while index >= 0:
        edge = getEdgeProd(index, n, m)
        # reduce max probability
        if edge[0] < majorIndex:
            majorIndex = edge[0]
            p = getMaxAssortativeP(k[majorIndex], k[0], k[-1], m, volume, epsilon, h, type)
        q = getAssortativeP(k[edge], volume, epsilon, h)
        if q > p:
            print("hi")
            edgeList.append(tuple(sortedIndices[edge]))
        elif random.random() <= q/p:
            edgeList.append(tuple(sortedIndices[edge]))
        index -= np.random.geometric(p) # index backwards
    return edgeList

def createAssortativeProdExhaustive(k, m, epsilon, type="large-degrees"):
    n = len(k)
    volume = np.sum(k)
    k1 = np.mean(k)
    k2 = np.mean(np.power(k, 2))

    if type == "large-degrees":
        g = lambda k: np.prod(k/k1) - (k2/k1**2)**m
    elif type == "aligned-degrees":
        g = lambda k: np.sum([(combo[0] - k1)*(combo[1] - k1) for combo in combinations(k, 2)])/(comb(len(k), 2, exact=True)*k1**2) - (k2 - k1**2)**2/k1**4

    edgeList = list()
    for edge in list(product(range(n), repeat=m)):
        p = getAssortativeP(k[list(edge)], volume, epsilon, g)
        if random.random() < p:
            edgeList.append(edge)
    return edgeList

def getAssortativeP(k, volume, epsilon, h):
    m = len(k)
    return min(np.prod(k/volume)*volume/m*(1 + epsilon*h(k)), 1)

def getMaxAssortativeP(kMajor, kMin, kMax, m, volume, epsilon, h, type):
    maxConfig = kMajor*(kMax/volume)**(m-1)/m
    maxH = getMaxH(kMin, kMax, m, epsilon, h, type)
    return min(maxConfig*maxH, 1)

def getMaxH(kMin, kMax, m, epsilon, h, type):
    if type == "large-degrees":
        return max(1 + epsilon*h([kMin]*m), 1 + epsilon*h([kMax]*m))
    elif type == "aligned-degrees":
        return max(1 + epsilon*h([kMax]*m), 1 + epsilon*h([kMin]*int(m/2) + [kMax]*(m-int(m/2))))

def getMinH(kMin, kMax, m, epsilon, h, type):
    if type == "large-degrees":
        return min(1 + epsilon*h([kMin]*m), 1 + epsilon*h([kMax]*m))
    elif type == "aligned-degrees":
        return min(1 + epsilon*h([kMax]*m), 1 + epsilon*h([kMin]*int(m/2) + [kMax]*(m-int(m/2))))

#https://stackoverflow.com/questions/53834707/element-at-index-in-itertools-product
def getEdgeProd(index, n, hyperedgeSize):
    return [(index // (n**r) % n) for r in range(hyperedgeSize-1, -1, -1)]

def generatePowerLawDegreeSequence(n, minDegree, maxDegree, exponent):
    return np.round(invCDFPowerLaw(np.random.rand(n), minDegree, maxDegree, exponent))

def invCDFPowerLaw(u, minDegree, maxDegree, exponent):
    return (minDegree**(1-exponent) + u*(maxDegree**(1-exponent) - minDegree**(1-exponent)))**(1/(1-exponent))
