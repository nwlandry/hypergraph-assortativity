#%%
from GenerativeModels import *
from HyperEigenvalues import *
from MeanFieldTheory import *
import numpy as np
from math import factorial
import matplotlib.pyplot as plt
import Hypergraph
import copy
import os
import shelve
#%%
mainFolder = os.getcwd()

dataFolder = "Data"
datasetFolder = "congress-bills"

with shelve.open(os.path.join(mainFolder, dataFolder, datasetFolder, datasetFolder + "_eigenvalues")) as data:
    assortativities = data["rho"]
    meanFieldCECEigenvalues = data["mean-field-eigenvalues"]
    trueCECEigenvalues = data["true-eigenvalues"]

plt.figure()
plt.plot(assortativities, trueCECEigenvalues, 'ko-', label="CEC Eigenvalue (True)")
plt.plot(assortativities, meanFieldCECEigenvalues, 'ro-', label="CEC Eigenvalue (MF)")
plt.legend()
plt.xlabel(r"$\rho$")
plt.ylabel(r"$\lambda$")
plt.ylim([0, 1.1*max(max(trueCECEigenvalues), max(meanFieldCECEigenvalues))])
plt.tight_layout()
plt.show()
#%%
plt.figure()
plt.plot(assortativities, np.divide(np.abs(trueCECEigenvalues - meanFieldCECEigenvalues), trueCECEigenvalues), 'ko-', label="Error")
plt.legend()
plt.xlabel(r"$\rho$")
plt.ylabel(r"$|\lambda_{MF}-\lambda_{True}|/\lambda_{True}$")
plt.tight_layout()
plt.show()

# %%
