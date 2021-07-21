import numpy as np
from math import factorial
from numpy.linalg import norm

class SparseTensor:
    def __init__(self, indices, weights, numNodes):
        self.indices = indices
        self.weights = weights
        self.numNodes = numNodes
        self.dimension = len(indices[0])

    def apply(self, vector, type="multiply"):
        newVector = np.zeros(self.numNodes)
        for index in self.indices:
            # ordered permutations
            for shift in range(self.dimension):
                newVector[index[shift]] += self.g(vector, index[shift+1:] + index[:shift], type=type)
        return newVector
    def g(self, vector, indices, type="multiply"):
        if type == "multiply":
            return np.prod(vector[list(indices)])
        if type == "add":
            return np.sum(vector[list(indices)])

    def f(self, vector, dimension=2, type="id"):
        if type == "id":
            return vector
        if type == "power":
            return np.power(vector, 1.0/(dimension - 1))

    def getCEC(self, maxIter, tolerance):
        l = 0
        x = np.random.uniform(size=(self.numNodes))
        for i in range(maxIter):
            newX = self.apply(x, type="add")
            newX = self.f(newX, dimension=self.dimension, type="id")
            newL = norm(newX)/norm(x)
            newX = newX/norm(newX)
            if abs(l - newL) <= tolerance:
                break
            x = newX.copy()
            l = newL
        return newL, newX

    def getZEC(self, maxIter, tolerance):
        l = 0
        x = np.random.uniform(size=(self.numNodes))
        for i in range(maxIter):
            newX = self.apply(x, type="multiply")
            newX = self.f(newX, dimension=self.dimension, type="id")
            newL = norm(newX)/norm(x)
            newX = newX/norm(newX)
            if abs(l - newL) <= tolerance:
                break
            x = newX.copy()
            l = newL
        return newL, newX

    def getHEC(self, maxIter, tolerance):
        l = 0
        x = np.random.uniform(size=(self.numNodes))
        for i in range(maxIter):
            newX = self.apply(x, type="multiply")
            newX = self.f(newX, dimension=self.dimension, type="power")
            newL = norm(newX)/norm(x)
            newX = newX/norm(newX)
            if abs(l - newL) <= tolerance:
                break
            x = newX.copy()
            l = newL
        return newL, newX
