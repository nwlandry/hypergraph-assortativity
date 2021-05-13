import Hypergraph
import numpy as np
import HypergraphContagion
import MeanFieldTheory
import shelve
import copy
import os

mainFolder = os.getcwd()

dataFolder = "Data"
datasetFolder = "Power-Law"
# Hypergraph parameters
n = 10000
isDegreeCorrelated = True
parameters = [{"degree-distribution":"power-law","exponent":3,"min-degree":10,"max-degree":1000,"hyperedge-size":2,"size":n,"is-correlated":isDegreeCorrelated},{"degree-distribution":"power-law","exponent":3,"min-degree":10,"max-degree":1000,"hyperedge-size":3,"size":n,"is-correlated":isDegreeCorrelated}]

h = Hypergraph.HypergraphGenerator(parameters)

# Shuffle parameters
sizeToShuffle = 3
changeAssortativity = 0.25
tolerance = 0.01
temperature = 0.00001
maxIterations = 1e5

types = ["aligned-degrees", "large-degrees", "top-bottom", "top-2"]

for type in types:
    hLower = copy.deepcopy(h)
    hNeutral = copy.deepcopy(h)
    hHigher = copy.deepcopy(h)

    originalAssortativity = hNeutral.getAssortativity(sizeToShuffle, type)
    print(originalAssortativity)
    hLower.shuffleHyperedges(originalAssortativity-changeAssortativity, sizeToShuffle, type, tolerance=tolerance, maxIterations=maxIterations, temperature=temperature)
    lowerAssortativity = hLower.getAssortativity(sizeToShuffle, type)
    print(lowerAssortativity)
    hHigher.shuffleHyperedges(originalAssortativity+changeAssortativity, sizeToShuffle, type, tolerance=tolerance, maxIterations=maxIterations, temperature=temperature)
    higherAssortativity = hHigher.getAssortativity(sizeToShuffle, type)
    print(higherAssortativity)

    with shelve.open(os.path.join(mainFolder, dataFolder, datasetFolder, datasetFolder + "_" + type)) as data:
        data["lower-assortativity-hypergraph"] = hLower.getHyperedgeList()
        data["zero-assortativity-hypergraph"] = hNeutral.getHyperedgeList()
        data["higher-assortativity-hypergraph"] = hHigher.getHyperedgeList()
