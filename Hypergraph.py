import math
import random
import time
from collections import defaultdict
import numpy as np
import copy
import itertools
from numpy.linalg import norm

class Hypergraph:
    def __init__(self, hyperedges, weightedEdges=False):
        self.addEdges(hyperedges, weightedEdges)
        self.deleteDegenerateHyperedges()
        self.findHyperedgeSizes()
        self.generateNeighbors()
        self.nodeLabels = list(self.nodes.keys())

    def __iter__(self):
        """Iterate over the nodes. Use: 'for n in G'.
        Returns
        -------
        niter : iterator
            An iterator over all nodes in the graph.
        Examples
        --------
        >>> G = nx.path_graph(4)  # or DiGraph, MultiGraph, MultiDiGraph, etc
        >>> [n for n in G]
        [0, 1, 2, 3]
        >>> list(G)
        [0, 1, 2, 3]
        """
        return iter(self.nodes)

    def __contains__(self, n):
        """Returns True if n is a node, False otherwise. Use: 'n in G'.
        Examples
        --------
        >>> G = nx.path_graph(4)  # or DiGraph, MultiGraph, MultiDiGraph, etc
        >>> 1 in G
        True
        """
        try:
            return n in self.nodes
        except TypeError:
            return False

    def __len__(self):
        """Returns the number of nodes in the graph. Use: 'len(G)'.
        Returns
        -------
        nnodes : int
            The number of nodes in the graph.
        See Also
        --------
        number_of_nodes, order  which are identical
        Examples
        --------
        >>> G = nx.path_graph(4)  # or DiGraph, MultiGraph, MultiDiGraph, etc
        >>> len(G)
        4
        """
        return len(self.nodes)

    def addEdges(self, hyperedges, weightedEdges):
        # unweighted format for hyperedges: {"id0":{"members":(1,2,3)}, "id1":{"members":(1,2)},...}
        # weighted format for hyperedges: {"id0":{"members":(1,2,3),"weight":1.1}, "id1":{"members":(1,2),"weight":0.5},...}
        self.weightedEdges = weightedEdges
        self.nodes = dict()
        nodes = set()
        # if list of tuples
        if isinstance(hyperedges, list):
            self.hyperedges = dict()
            uid = 0
            for hyperedge in hyperedges:
                if self.weightedEdges:
                    self.hyperedges[uid] = {"members":hyperedge[:-1],"weight":hyperedge[-1]}
                else:
                    self.hyperedges[uid] = {"members":hyperedge}
                    nodes.update(hyperedge)
                uid += 1

        elif isinstance(hyperedges, dict):
            self.hyperedges = hyperedges.copy()
            for edgeData in self.hyperedges.values():
                nodes.update(edgeData["members"])

        for nodeLabel in list(nodes):
            self.nodes[nodeLabel] = dict()
        # need a better way to check whether the format is correct

    def addNodeAttributes(self, nodeAttributes):
        # find unique nodes in the hyperedges
        for label, attribute in nodeAttributes.items():
            try:
                self.nodes[label] = attribute
            except:
                print("invalid label")

    def deleteDegenerateHyperedges(self):
        cleanedHyperedges = dict()
        for uid, hyperedge in self.hyperedges.items():
            if len(hyperedge["members"]) >= 2:
                cleanedHyperedges[uid] = hyperedge
        self.hyperedges = cleanedHyperedges

    def number_of_nodes(self):
        return len(self.nodes)

    def has_node(self, n):
        try:
            return n in self.nodes
        except TypeError:
            return False

    def findHyperedgeSizes(self):
        hyperedgeSizes = set()
        for edgeData in list(self.hyperedges.values()):
            hyperedgeSizes.add(len(edgeData["members"]))
        self.hyperedgeSizes = list(hyperedgeSizes)

    def getHyperedgeSizes(self):
        return self.hyperedgeSizes

    def generateNeighbors(self):
        self.neighbors = defaultdict(dict)
        if self.weightedEdges:
            self.generateWeightedNeighbors()
        else:
            self.generateUnweightedNeighbors()

    def generateUnweightedNeighbors(self):
        for uid, edgeData in self.hyperedges.items():
            try:
                members = edgeData["members"]
            except:
                print("Incorrect input format for hyperedge list")
                break
            for index in range(len(members)):
                self.neighbors[members[index]][uid] = {"neighbors":tuple([members[i] for i in range(len(members)) if i != index])}

    def generateWeightedNeighbors(self):
        for uid, edgeData in self.hyperedges.items():
            try:
                members = edgeData["members"]
                weight = edgeData["weight"]
            except:
                print("Incorrect input format for weighted hyperedge list")
                break
            for index in range(len(members)):
                self.neighbors[members[index]][uid] = {"neighbors":tuple([members[i] for i in range(len(members)) if i != index]), "weight":weight}

    def getDegreeSequenceBySize(self, hyperedgeSize):
        degreeSequence = defaultdict(lambda: 0)
        for uid, hyperedge in self.hyperedges.items():
            if len(hyperedge["members"]) == hyperedgeSize:
                for node in hyperedge["members"]:
                    degreeSequence[node] += 1
        return degreeSequence

    def getHyperedgesBySize(self, hyperedgeSize):
        hyperedges = list()
        for uid, data in self.hyperedges.items():
            if len(data["members"]) == hyperedgeSize:
                hyperedges.append(data["members"])
        return hyperedges


class HypergraphGenerator:
    def __init__(self, data, type="structure"):
        if type == "structure":
            self.generateHyperdegreeSequence(data)
            self.generateHyperedges()
        elif type == "hyperedge-list":
            self.hyperedgeListToDictionary(data)
            self.getHyperdegreeSequence()
        elif type == "hyperedge-dictionary":
            self.hyperedges = copy.deepcopy(data)
            self.getHyperdegreeSequence()
        else:
            print("invalid option")

    def getHyperedges(self):
        return self.hyperedges

    def getHyperedgeList(self):
        return [e["members"] for e in list(self.hyperedges.values())]

    def hyperedgeListToDictionary(self, hyperedges):
        self.hyperedges = dict()
        uid = 0
        for hyperedge in hyperedges:
            self.hyperedges[uid] = {"members":hyperedge}
            uid += 1

    def getHyperdegreeSequence(self):
        try:
            return self.hyperdegreeSequence
        except:
            self.hyperdegreeSequence = defaultdict(lambda: defaultdict(lambda: 0))
            for uid, data in self.hyperedges.items():
                members = data["members"]
                hyperedgeSize = len(members)
                for node in members:
                    self.hyperdegreeSequence[node][hyperedgeSize] += 1

    def generateHyperdegreeSequence(self, structure):
        self.hyperdegreeSequence = defaultdict(lambda: defaultdict(lambda: 0))
        self.hyperedgeSizes = list()
        correlatedSequence = list()
        for info in structure:
            try:
                hyperedgeSize = info["hyperedge-size"]
                self.hyperedgeSizes.append(hyperedgeSize)
            except:
                print("Error in specified distribution parameters")

            if info["degree-distribution"] == "power-law":
                try:
                    numNodes = info["size"]
                    minDegree = info["min-degree"]
                    maxDegree = info["max-degree"]
                    exponent = info["exponent"]
                except:
                    print("Error in specified distribution parameters")
                    break
                sequence = self.generatePowerLawDegreeSequence(numNodes, minDegree, maxDegree, exponent)
            elif info["degree-distribution"] == "uniform":
                try:
                    numNodes = info["size"]
                    minDegree = info["min-degree"]
                    maxDegree = info["max-degree"]
                except:
                    print("Error in specified distribution parameters")
                    break
                sequence = self.generateUniformDegreeSequence(numNodes, minDegree, maxDegree)
            elif info["degree-distribution"] == "poisson":
                try:
                    numNodes = info["size"]
                    meanDegree = info["mean-degree"]
                except:
                    print("Error in specified distribution parameters")
                    break
                sequence = self.generatePoissonDegreeSequence(numNodes, meanDegree)
            else:
                print("Invalid selection")
                break
            try:
                isCorrelated = info["is-correlated"]
            except:
                print("Specify whether this hyperedge size is correlated or not")
                break
            if isCorrelated:
                if correlatedSequence == []:
                    correlatedSequence = sequence
                self.updateHyperdegreeSequence(correlatedSequence, hyperedgeSize)
                #self.hyperdegreeSequence[hyperedgeSize] = correlatedSequence
            else:
                self.updateHyperdegreeSequence(sequence, hyperedgeSize)
                #self.hyperdegreeSequence[hyperedgeSize] = sequence

    def updateHyperdegreeSequence(self, sequence, hyperedgeSize):
        try:
            for node in range(len(sequence)):
                self.hyperdegreeSequence[node][hyperedgeSize] = sequence[node]
        except:
            self.hyperdegreeSequence = defaultdict(lambda: defaultdict(lambda: 0))
            for node in range(len(sequence)):
                self.hyperdegreeSequence[node][hyperedgeSize] = sequence[node]

    def generateHyperedges(self):
        self.hyperedges = dict()
        for hyperedgeSize in self.hyperedgeSizes:
            self.hyperedges.update(self.generateHyperedgesBySize(hyperedgeSize))

    def generateHyperedgesBySize(self, hyperedgeSize):
        import string
        k = dict()
        for node, hyperdegree in self.hyperdegreeSequence.items():
            k[node] = hyperdegree[hyperedgeSize]

        # Making sure we have the right number of stubs
        if (sum(k.values()) % hyperedgeSize) != 0:
            remainder = sum(k.values()) % hyperedgeSize
            for i in range(int(round(hyperedgeSize - remainder))):
                j = random.randrange(len(k))
                k[j] = k[j] + 1

        stubs = list()
        hyperedges = dict()
        # Creating the list to index through
        for index in range(len(k)):
            stubs.extend([index]*int(k[index]))

        while len(stubs) != 0:
            uid = ''.join(random.choice(string.ascii_lowercase) for i in range(8))
            u = random.sample(range(len(stubs)), hyperedgeSize)
            hyperedge = list()
            for index in u:
                hyperedge.append(stubs[index])

            hyperedges[uid] = {"members":tuple(hyperedge)}

            for index in sorted(u, reverse=True):
                del stubs[index]
        return hyperedges

    def generatePowerLawDegreeSequence(self, numNodes, minDegree, maxDegree, exponent):
        degreeSequence = list()
        for i in range(numNodes):
            u = random.uniform(0, 1)
            degreeSequence.append(round(self.invCDFPowerLaw(u, minDegree, maxDegree, exponent)))
        return degreeSequence # originally this was sorted but I'm worried about between-size correlations

    def invCDFPowerLaw(self, u, minDegree, maxDegree, exponent):
        return (minDegree**(1-exponent) + u*(maxDegree**(1-exponent) - minDegree**(1-exponent)))**(1/(1-exponent))

    def generateUniformDegreeSequence(self, numNodes, minDegree, maxDegree):
        degreeSequence = list()
        for i in range(numNodes):
            u = random.randrange(round(minDegree), round(maxDegree))
            degreeSequence.append(round(u))
        return degreeSequence

    def generatePoissonDegreeSequence(self, numNodes, meanDegree):
        return np.random.poisson(lam=meanDegree, size=numNodes).tolist()

    def getHyperedgesBySize(self, hyperedgeSize, ids=True):
        if ids:
            hyperedgesBySize = dict()
            for uid, hyperedgeData in self.hyperedges.items():
                if len(hyperedgeData["members"]) == hyperedgeSize:
                    hyperedgesBySize[uid] = dict(hyperedgeData)
            return hyperedgesBySize
        else:
            hyperedgesBySize = list()
            for hyperedgeData in list(self.hyperedges.values()):
                if len(hyperedgeData["members"]) == hyperedgeSize:
                    hyperedgesBySize.append(hyperedgeData["members"])
            return hyperedgesBySize

    def updateHyperedgesAfterShuffle(self, hyperedgeListBySize):
        for key, hyperedge in hyperedgeListBySize.items():
            self.hyperedges[key] = dict(hyperedge)

    # assortativity function and averaging
    def degreeAssortativityFunction(self, hyperdegrees, type, **kwargs):
        if type == "aligned-degrees":
            value = 0
            meanDegree = kwargs["meandegree"]
            isVariance = kwargs["isvariance"]
            hyperedgeSize = len(hyperdegrees)
            if isVariance:
                for i in range(hyperedgeSize):
                        value += (hyperdegrees[i] - meanDegree)**2/meanDegree**2
                return value/math.comb(hyperedgeSize, 2)
            else:
                for i in range(hyperedgeSize):
                    for j in range(i):
                        value += (hyperdegrees[i] - meanDegree)*(hyperdegrees[j] - meanDegree)/meanDegree**2
                return value/math.comb(hyperedgeSize, 2)

        elif type == "large-degrees":
            isVariance = kwargs["isvariance"]
            meanDegree = kwargs["meandegree"]
            hyperedgeSize = len(hyperdegrees)
            if isVariance:
                value = 0
                for degree in hyperdegrees:
                    value += (degree/meanDegree)**hyperedgeSize
                return value/hyperedgeSize
            else:
                value = 1
                for degree in hyperdegrees:
                    value *= degree/meanDegree
                return value

        elif type == "top-bottom":
            isVariance = kwargs["isvariance"]
            meanMinDegree = kwargs["meanmindegree"]
            meanMaxDegree = kwargs["meanmaxdegree"]
            if isVariance:
                return 0.5*((max(hyperdegrees) - meanMaxDegree)**2 + (min(hyperdegrees) - meanMinDegree)**2)
            else:
                return (max(hyperdegrees) - meanMaxDegree)*(min(hyperdegrees) - meanMinDegree)

        elif type == "top-2":
            isVariance = kwargs["isvariance"]
            meanSecondMaxDegree = kwargs["meansecondmaxdegree"]
            meanMaxDegree = kwargs["meanmaxdegree"]
            if isVariance:
                hyperdegrees = sorted(hyperdegrees)
                return 0.5*((hyperdegrees[-1] - meanMaxDegree)**2 + (hyperdegrees[-2] - meanSecondMaxDegree)**2)
            else:
                hyperdegrees = sorted(hyperdegrees)
                return (hyperdegrees[-1] - meanMaxDegree)*(hyperdegrees[-2] - meanSecondMaxDegree)

        else:
            print("Not a valid selection")

    def averageHValue(self, hyperedgeSize, type="aligned-degrees"):
        if type == "aligned-degrees":
            meanDegree = self.getDegreeMoment(hyperedgeSize, power=1)
            kwargs = {"meandegree":meanDegree}
        elif type == "large-degrees":
            meanDegree = self.getDegreeMoment(hyperedgeSize, power=1)
            kwargs = {"meandegree":meanDegree}
        elif type == "top-bottom":
            meanMinDegree = self.getOrderedDegreeMoment(hyperedgeSize, power=1, type="smallest")
            meanMaxDegree = self.getOrderedDegreeMoment(hyperedgeSize, power=1, type="largest")
            kwargs = {"meanmindegree":meanMinDegree, "meanmaxdegree":meanMaxDegree}
        elif type == "top-2":
            meanSecondMaxDegree = self.getOrderedDegreeMoment(hyperedgeSize, power=1, type="second-largest")
            meanMaxDegree = self.getOrderedDegreeMoment(hyperedgeSize, power=1, type="largest")
            kwargs = {"meansecondmaxdegree":meanSecondMaxDegree, "meanmaxdegree":meanMaxDegree}

        hyperedges = self.getHyperedgesBySize(hyperedgeSize, ids=False)
        numEdges = len(hyperedges)
        hAvg = 0
        for edge in hyperedges:
            degrees = [self.hyperdegreeSequence[index][hyperedgeSize] for index in edge]
            hAvg += self.degreeAssortativityFunction(degrees, type, isvariance=False, **kwargs)/numEdges
        return hAvg

    def updateAverageHValue(self, hAvg, degrees1, degrees2, newDegrees1, newDegrees2, hyperedgeSize, numEdges, type="aligned-degrees", **kwargs):
        return hAvg + self.degreeAssortativityFunction(newDegrees1, type, isvariance=False, **kwargs)/numEdges + self.degreeAssortativityFunction(newDegrees2, type, isvariance=False, **kwargs)/numEdges - self.degreeAssortativityFunction(degrees1, type, isvariance=False, **kwargs)/numEdges - self.degreeAssortativityFunction(degrees2, type, isvariance=False, **kwargs)/numEdges

    def maxHValue(self, hyperedgeSize, type="aligned-degrees"):
        if type == "aligned-degrees":
            meanDegree = self.getDegreeMoment(hyperedgeSize, power=1)
            kwargs = {"meandegree":meanDegree}
        elif type == "large-degrees":
            meanDegree = self.getDegreeMoment(hyperedgeSize, power=1)
            kwargs = {"meandegree":meanDegree}
        elif type == "top-bottom":
            meanMinDegree = self.getOrderedDegreeMoment(hyperedgeSize, power=1, type="smallest")
            meanMaxDegree = self.getOrderedDegreeMoment(hyperedgeSize, power=1, type="largest")
            kwargs = {"meanmindegree":meanMinDegree, "meanmaxdegree":meanMaxDegree}
        elif type == "top-2":
            meanSecondMaxDegree = self.getOrderedDegreeMoment(hyperedgeSize, power=1, type="second-largest")
            meanMaxDegree = self.getOrderedDegreeMoment(hyperedgeSize, power=1, type="largest")
            kwargs = {"meansecondmaxdegree":meanSecondMaxDegree, "meanmaxdegree":meanMaxDegree}

        hyperedges = self.getHyperedgesBySize(hyperedgeSize, ids=False)
        numEdges = len(hyperedges)
        hMax = 0
        for edge in hyperedges:
            degrees = [self.hyperdegreeSequence[index][hyperedgeSize] for index in edge]
            hMax += self.degreeAssortativityFunction(degrees, type, isvariance=True, **kwargs)/numEdges
        return hMax

    def updateMaxHValue(self, hMax, degrees1, degrees2, newDegrees1, newDegrees2, hyperedgeSize, numEdges, type="aligned-degrees", **kwargs):
        if type == "aligned-degrees":
            return hMax
        elif type == "large-degrees":
            return hMax
        elif type == "top-bottom":
            return hMax + self.degreeAssortativityFunction(newDegrees1, type, isvariance=True, **kwargs)/numEdges + self.degreeAssortativityFunction(newDegrees2, type, isvariance=True, **kwargs)/numEdges - self.degreeAssortativityFunction(degrees1, type, isvariance=True, **kwargs)/numEdges - self.degreeAssortativityFunction(degrees2, type, isvariance=True, **kwargs)/numEdges
        elif type == "top-2":
            return hMax + self.degreeAssortativityFunction(newDegrees1, type, isvariance=True, **kwargs)/numEdges + self.degreeAssortativityFunction(newDegrees2, type, isvariance=True, **kwargs)/numEdges - self.degreeAssortativityFunction(degrees1, type, isvariance=True, **kwargs)/numEdges - self.degreeAssortativityFunction(degrees2, type, isvariance=True, **kwargs)/numEdges

    def getDegreeMoment(self, hyperedgeSize, power=1):
        hyperdegreeSequence = [hyperdegree[hyperedgeSize] for hyperdegree in list(self.hyperdegreeSequence.values())]
        return sum([hyperdegree**power for hyperdegree in hyperdegreeSequence])/float(len(hyperdegreeSequence))

    def getOrderedDegreeMoment(self, hyperedgeSize, power=1, type="largest"):
        hyperedges = self.getHyperedgesBySize(hyperedgeSize, ids=False)
        if type == "largest":
            return np.mean([max([self.hyperdegreeSequence[index][hyperedgeSize] for index in edge]) for edge in hyperedges])
        elif type == "second-largest":
            return np.mean([sorted([self.hyperdegreeSequence[index][hyperedgeSize] for index in edge])[-2] for edge in hyperedges])
        elif type == "smallest":
            return np.mean([min([self.hyperdegreeSequence[index][hyperedgeSize] for index in edge]) for edge in hyperedges])

    def updateOrderedDegreeMoment(self, currentValue, degrees1, degrees2, newDegrees1, newDegrees2, numEdges, type="largest"):
            if type == "largest":
                return currentValue + (max(newDegrees1) + max(newDegrees2) - max(degrees1) - max(degrees2))/numEdges
            elif type == "second-largest":
                return currentValue + (sorted(newDegrees1)[-2] + sorted(newDegrees2)[-2] - sorted(degrees1)[-2] - sorted(degrees1)[-2])/numEdges
            elif type == "smallest":
                return currentValue + (min(newDegrees1) + min(newDegrees2) - min(degrees1) - min(degrees2))/numEdges

    def nullHValue(self, hyperedgeSize, type="aligned-degrees"):
        if type == "aligned-degrees":
            meanDegree = self.getDegreeMoment(hyperedgeSize, power=1)
            meanSquaredDegree = self.getDegreeMoment(hyperedgeSize, power=2)
            return (meanSquaredDegree - meanDegree**2)**2/meanDegree**4
        elif type == "large-degrees":
            meanDegree = self.getDegreeMoment(hyperedgeSize, power=1)
            meanSquaredDegree = self.getDegreeMoment(hyperedgeSize, power=2)
            return (meanSquaredDegree/meanDegree**2)**hyperedgeSize
        elif type == "top-bottom":
            return self.topBottomNullValue(hyperedgeSize)
        elif type == "top-2":
            return self.top2NullValue(hyperedgeSize)

    def topBottomNullValue(self, hyperedgeSize):
        k = [hyperdegree[hyperedgeSize] for hyperdegree in list(self.hyperdegreeSequence.values())]
        degrees, counts = np.unique(k, return_counts=True)
        p = counts/sum(counts)

        meanDegree = np.sum(np.multiply(degrees, p))
        combos = [list(edge) for  edge in itertools.product(range(len(degrees)), repeat=hyperedgeSize)]
        meanMinDegree = np.sum([np.prod(p[indices])*np.min(degrees[indices])*np.prod(degrees[indices]) for indices in combos])/meanDegree**hyperedgeSize
        meanMaxDegree = np.sum([np.prod(p[indices])*np.max(degrees[indices])*np.prod(degrees[indices]) for indices in combos])/meanDegree**hyperedgeSize
        product = np.sum([np.prod(p[indices])*np.min(degrees[indices])*np.max(degrees[indices])*np.prod(degrees[indices]) for indices in combos])/meanDegree**hyperedgeSize
        return product - meanMinDegree*meanMaxDegree

    def top2NullValue(self, hyperedgeSize):
        k = [hyperdegree[hyperedgeSize] for hyperdegree in list(self.hyperdegreeSequence.values())]
        degrees, counts = np.unique(k, return_counts=True)
        p = counts/sum(counts)

        meanDegree = np.sum(np.multiply(degrees, p))
        combos = [list(edge) for  edge in itertools.product(range(len(degrees)), repeat=hyperedgeSize)]
        meanSecondMaxDegree = np.sum([np.prod(p[indices])*np.sort(degrees[indices])[-2]*np.prod(degrees[indices]) for indices in combos])/meanDegree**hyperedgeSize
        meanMaxDegree = np.sum([np.prod(p[indices])*np.max(degrees[indices])*np.prod(degrees[indices]) for indices in combos])/meanDegree**hyperedgeSize
        product = np.sum([np.prod(p[indices])*np.sort(degrees[indices])[-2]*np.max(degrees[indices])*np.prod(degrees[indices]) for indices in combos])/meanDegree**hyperedgeSize
        return product - meanSecondMaxDegree*meanMaxDegree

    def shuffleHyperedges(self, targetAssortativity, hyperedgeSize, type, tolerance=0.01, maxIterations=10000, temperature=0.001, isVerbose=False):
        shuffledEdges = self.getHyperedgesBySize(hyperedgeSize, ids=True)
        numEdges = len(shuffledEdges)

        if type == "aligned-degrees" or type == "large-degrees":
            meanDegree = self.getDegreeMoment(hyperedgeSize, power=1)
        elif type == "top-bottom":
            meanMinDegree = self.getOrderedDegreeMoment(hyperedgeSize, power=1, type="smallest")
            meanMaxDegree = self.getOrderedDegreeMoment(hyperedgeSize, power=1, type="largest")
        elif type == "top-2":
            meanSecondMaxDegree = self.getOrderedDegreeMoment(hyperedgeSize, power=1, type="second-largest")
            meanMaxDegree = self.getOrderedDegreeMoment(hyperedgeSize, power=1, type="largest")

        hAvg = self.averageHValue(hyperedgeSize, type=type)
         # these values don't change for aligned-degrees
        hNull = self.nullHValue(hyperedgeSize, type=type)
        hMax = self.maxHValue(hyperedgeSize, type=type)

        assortativity = (hAvg - hNull)/(hMax - hNull)
        if isVerbose:
            assortativityList = [assortativity]
        iteration = 0
        while (abs(targetAssortativity - assortativity) > tolerance) and (iteration < maxIterations):
            uid1, uid2 = random.sample(shuffledEdges.keys(), 2)
            member1, member2 = random.sample(range(hyperedgeSize), 2)
            edge1 = shuffledEdges[uid1]["members"]
            edge2 = shuffledEdges[uid2]["members"]
            newEdge1 = list(edge1)
            newEdge2 = list(edge2)
            newEdge1[member1] = edge2[member2]
            newEdge2[member2] = edge1[member1]
            # ensure that we don't create/destroy multiedges
            if len(set(edge1)) == len(set(newEdge1)) and len(set(edge2)) == len(set(newEdge2)):
                degrees1 = [self.hyperdegreeSequence[index][hyperedgeSize] for index in edge1]
                degrees2 = [self.hyperdegreeSequence[index][hyperedgeSize] for index in edge2]
                newDegrees1 = [self.hyperdegreeSequence[index][hyperedgeSize] for index in newEdge1]
                newDegrees2 = [self.hyperdegreeSequence[index][hyperedgeSize] for index in newEdge2]
                # handle different cases
                if type == "aligned-degrees":
                    avgKwargs = {"meandegree":meanDegree}
                    maxKwargs = {}
                elif type == "large-degrees":
                    avgKwargs = {"meandegree":meanDegree}
                    maxKwargs = {}
                elif type == "top-bottom":
                    newMeanMinDegree = self.updateOrderedDegreeMoment(meanMinDegree, degrees1, degrees2, newDegrees1, newDegrees2, numEdges, type="smallest")
                    newMeanMaxDegree = self.updateOrderedDegreeMoment(meanMaxDegree, degrees1, degrees2, newDegrees1, newDegrees2, numEdges, type="largest")
                    avgKwargs = {"meanmindegree":newMeanMinDegree, "meanmaxdegree":newMeanMaxDegree}
                    maxKwargs = {"meanmindegree":newMeanMinDegree, "meanmaxdegree":newMeanMaxDegree}
                elif type == "top-2":
                    newMeanSecondMaxDegree = self.updateOrderedDegreeMoment(meanSecondMaxDegree, degrees1, degrees2, newDegrees1, newDegrees2, numEdges, type="second-largest")
                    newMeanMaxDegree = self.updateOrderedDegreeMoment(meanMaxDegree, degrees1, degrees2, newDegrees1, newDegrees2, numEdges, type="largest")
                    avgKwargs = {"meansecondmaxdegree":newMeanSecondMaxDegree, "meanmaxdegree":newMeanMaxDegree}
                    maxKwargs = {"meansecondmaxdegree":newMeanSecondMaxDegree, "meanmaxdegree":newMeanMaxDegree}

                newHAvg = self.updateAverageHValue(hAvg, degrees1, degrees2, newDegrees1, newDegrees2, hyperedgeSize, numEdges, type=type, **avgKwargs)

                newHMax = self.updateMaxHValue(hMax, degrees1, degrees2, newDegrees1, newDegrees2, hyperedgeSize, numEdges, type=type, **maxKwargs)

                newAssortativity = (newHAvg - hNull)/(newHMax - hNull)

                difference = (assortativity - targetAssortativity)**2 - (newAssortativity - targetAssortativity)**2
                if random.random() <= self.boltzmannProbability(difference, temperature):
                    shuffledEdges[uid1]["members"] = tuple(newEdge1)
                    shuffledEdges[uid2]["members"] = tuple(newEdge2)
                    hAvg = newHAvg
                    hMax = newHMax
                    assortativity = newAssortativity
                    if isVerbose:
                        assortativityList.append(assortativity)
                    if type == "top-bottom":
                        meanMinDegree = newMeanMinDegree
                        meanMaxDegree = newMeanMaxDegree
                    elif type == "top-2":
                        meanSecondMaxDegree = newMeanSecondMaxDegree
                        meanMaxDegree = newMeanMaxDegree
                    iteration += 1
        print("Number of double-edge swaps: " + str(iteration), flush=True)
        self.updateHyperedgesAfterShuffle(shuffledEdges)
        if isVerbose:
            return assortativityList

    def boltzmannProbability(self, difference, temperature):
        if difference >= 0:
            return 1
        if temperature == 0:
            return 0
        else:
            return math.exp(difference/temperature)

    def getAssortativity(self, hyperedgeSize, type):
        if type == "aligned-degrees":
            meanDegree = self.getDegreeMoment(hyperedgeSize, power=1)
        elif type == "top-bottom":
            meanMinDegree = self.getOrderedDegreeMoment(hyperedgeSize, power=1, type="smallest")
            meanMaxDegree = self.getOrderedDegreeMoment(hyperedgeSize, power=1, type="largest")
        elif type == "top-2":
            meanSecondMaxDegree = self.getOrderedDegreeMoment(hyperedgeSize, power=1, type="second-largest")
            meanMaxDegree = self.getOrderedDegreeMoment(hyperedgeSize, power=1, type="largest")

        hAvg = self.averageHValue(hyperedgeSize, type=type)
        hNull = self.nullHValue(hyperedgeSize, type=type)
        hMax = self.maxHValue(hyperedgeSize, type=type)
        return (hAvg - hNull)/(hMax - hNull)

    def shuffleHyperedgesByCommunity(self, targetFractionOfEdges, hyperedgeSize, communities, tolerance=0.01, maxIterations=10000, temperature=0.001):
        if abs(np.sum(targetFractionOfEdges) - 1) > 0.00001:
            print("Fraction of edges does not add up to 1")
            return
        shuffledEdges = self.getHyperedgesBySize(hyperedgeSize, ids=True)
        numEdges = len(shuffledEdges)
        fracNumCommunities = np.zeros(min(len(np.unique(communities)), hyperedgeSize))
        for edgeData in list(shuffledEdges.values()):
            fracNumCommunities[self.getNumCommunities(edgeData["members"], communities)-1] += 1/numEdges

        print(fracNumCommunities)

        iteration = 0
        while (norm(fracNumCommunities - targetFractionOfEdges) > tolerance) and (iteration < maxIterations):
            uid1, uid2 = random.sample(shuffledEdges.keys(), 2)
            member1, member2 = random.sample(range(hyperedgeSize), 2)
            edge1 = shuffledEdges[uid1]["members"]
            edge2 = shuffledEdges[uid2]["members"]
            newEdge1 = list(edge1)
            newEdge2 = list(edge2)
            newEdge1[member1] = edge2[member2]
            newEdge2[member2] = edge1[member1]
            # ensure that we don't create/destroy multiedges
            if len(set(edge1)) == len(set(newEdge1)) and len(set(edge2)) == len(set(newEdge2)):
                numComm1 = self.getNumCommunities(edge1, communities)
                numComm2 = self.getNumCommunities(edge2, communities)
                newNumComm1 = self.getNumCommunities(newEdge1, communities)
                newNumComm2 = self.getNumCommunities(newEdge2, communities)
                newFracNumCommunities = fracNumCommunities.copy()
                newFracNumCommunities[numComm1-1] -= 1.0/numEdges
                newFracNumCommunities[numComm2-1] -= 1.0/numEdges
                newFracNumCommunities[newNumComm1-1] += 1.0/numEdges
                newFracNumCommunities[newNumComm2-1] += 1.0/numEdges

                if norm(newFracNumCommunities - targetFractionOfEdges) <= norm(fracNumCommunities - targetFractionOfEdges):
                    shuffledEdges[uid1]["members"] = tuple(newEdge1)
                    shuffledEdges[uid2]["members"] = tuple(newEdge2)
                    fracNumCommunities = newFracNumCommunities.copy()
                    iteration += 1
                elif random.random() <= math.exp((norm(fracNumCommunities - targetFractionOfEdges) - norm(newFracNumCommunities - targetFractionOfEdges))/temperature):
                    shuffledEdges[uid1]["members"] = tuple(newEdge1)
                    shuffledEdges[uid2]["members"] = tuple(newEdge2)
                    fracNumCommunities = newFracNumCommunities.copy()
                    iteration += 1
        print(fracNumCommunities)

        print("Number of double-edge swaps: " + str(iteration), flush=True)
        self.updateHyperedgesAfterShuffle(shuffledEdges)

    def getNumCommunities(self, hyperedge, communities):
        return len(np.unique(communities[list(hyperedge)]))
