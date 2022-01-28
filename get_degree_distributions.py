from HyperEigenvalues import *
from MeanFieldTheory import *
import numpy as np
import Hypergraph
import os
import shelve
import utilities
import os

m = 3

mainFolder = os.getcwd()
# Import the hypergraph
dataFolder = "Data"


#CONFIGURATION-MODEL
n = 10000
m = 3
isDegreeCorrelated = True
parameters = [{"degree-distribution":"power-law","exponent":3,"min-degree":10,"max-degree":100,"hyperedge-size":m,"size":n,"is-correlated":isDegreeCorrelated}]

h1 = Hypergraph.HypergraphGenerator(parameters)
k1 = h1.getHyperdegreeSequenceBySize(m)
print(np.mean(k1))

# TAGS-ASK-UBUNTU
datasetFolder = "tags-ask-ubuntu"
sizeFile = "tags-ask-ubuntu-nverts.txt"
memberFile = "tags-ask-ubuntu-simplices.txt"

hyperedgeSizeFile = os.path.join(dataFolder, datasetFolder, sizeFile)
memberIDFile = os.path.join(dataFolder, datasetFolder, memberFile)

hyperedgeSizes = [m]
hyperedgeList = utilities.readScHoLPData(hyperedgeSizeFile, memberIDFile)
hyperedgeList = utilities.filterHyperedgesBySize(hyperedgeList, hyperedgeSizes)
h2 = Hypergraph.HypergraphGenerator(hyperedgeList, type="hyperedge-list")
k2 = h2.getHyperdegreeSequenceBySize(m)
print(np.mean(k2))

# CONGRESS BILLS
datasetFolder = "congress-bills"
sizeFile = "congress-bills-nverts.txt"
memberFile = "congress-bills-simplices.txt"

hyperedgeSizeFile = os.path.join(dataFolder, datasetFolder, sizeFile)
memberIDFile = os.path.join(dataFolder, datasetFolder, memberFile)

hyperedgeSizes = [m]
hyperedgeList = utilities.readScHoLPData(hyperedgeSizeFile, memberIDFile)
hyperedgeList = utilities.filterHyperedgesBySize(hyperedgeList, hyperedgeSizes)
h3 = Hypergraph.HypergraphGenerator(hyperedgeList, type="hyperedge-list")
k3 = h3.getHyperdegreeSequenceBySize(m)
print(np.mean(k3))

#EU-EMAILS
datasetFolder = "Eu-Emails"
sizeFile = "email-Eu-nverts.txt"
memberFile = "email-Eu-simplices.txt"

hyperedgeSizeFile = os.path.join(dataFolder, datasetFolder, sizeFile)
memberIDFile = os.path.join(dataFolder, datasetFolder, memberFile)

hyperedgeSizes = [m]
hyperedgeList = utilities.readScHoLPData(hyperedgeSizeFile, memberIDFile)
hyperedgeList = utilities.filterHyperedgesBySize(hyperedgeList, hyperedgeSizes)
h4 = Hypergraph.HypergraphGenerator(hyperedgeList, type="hyperedge-list")
k4 = h4.getHyperdegreeSequenceBySize(m)
print(np.mean(k4))

with shelve.open(os.path.join(mainFolder, dataFolder,"degree_distribution")) as data:
    data["CM"] = k1
    data["TAU"] = k2
    data["CB"] = k3
    data["EE"] = k4