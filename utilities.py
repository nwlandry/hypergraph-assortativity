import csv
from collections import defaultdict

def readScHoLPData(hyperedgeSizeFile, memberIDFile, hyperedgeIDFile=None):
    hyperedgeList = list()
    with open(hyperedgeSizeFile) as sizeFile, open(memberIDFile) as idFile:
        sizes = sizeFile.read().splitlines()
        members = idFile.read().splitlines()
        memberIndex = 0
        for index in range(len(sizes)):
            edge = list()
            hyperedgeSize = int(sizes[index])
            for i in range(memberIndex, memberIndex+hyperedgeSize):
                member = int(members[i])
                edge.append(member)
            hyperedgeList.append(tuple(edge))
            memberIndex += hyperedgeSize
    return hyperedgeList

def readScHoLPDataInTimeRange(hyperedgeSizeFile, memberIDFile, times, minTime, maxTime):
    hyperedgeList = list()
    ids = set()
    with open(hyperedgeSizeFile) as sizeFile, open(memberIDFile) as idFile:
        sizes = sizeFile.read().splitlines()
        members = idFile.read().splitlines()
        memberIndex = 0
        for index in range(len(sizes)):
            edge = list()
            hyperedgeSize = int(sizes[index])
            if times[index] >= minTime and times[index] <= maxTime:
                for i in range(memberIndex, memberIndex+hyperedgeSize):
                    member = int(members[i])
                    edge.append(member)
                    ids.add(member)
                hyperedgeList.append(tuple(edge))
            memberIndex += hyperedgeSize
        return hyperedgeList, ids

def filterHyperedgesBySize(hyperedgeList, sizeList):
    filteredHyperedgeList = list()
    sizes = set(sizeList)
    for edge in hyperedgeList:
        if len(edge) in sizes:
            filteredHyperedgeList.append(edge)
    return filteredHyperedgeList

def getModularity(edgeList, communities, m):
    # group labels must be 0-indexed and sequential
    filteredEdges = filterHyperedgesBySize(edgeList, [m])
    edgesInSingleGroup = 0
    numEdges = len(edgeList)
    volume = numEdges*m
    groupIDs = set(communities.values())
    numberOfGroups = len(groupIDs)

    k = defaultdict(lambda : 0)

    for edge in filteredEdges:
        for node in edge:
            k[node] += 1
        if len(set([communities[node] for node in edge])) == 1:
            edgesInSingleGroup += 1

    volumeOfGroups = dict()
    nullValue = 0
    for node in list(communities.keys()):
        if communities[node] not in volumeOfGroups:
            volumeOfGroups[communities[node]] = k[node]
        else:
            volumeOfGroups[communities[node]] += k[node]

    for group in list(groupIDs):
        nullValue += (volumeOfGroups[group]/volume)**m

    return edgesInSingleGroup/numEdges - nullValue
