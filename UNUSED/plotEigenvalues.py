#%%
from GenerativeModels import *
from HyperEigenvalues import *
from MeanFieldTheory import *
import numpy as np
from math import factorial
import matplotlib.pyplot as plt
import os
import shelve

mainFolder = os.getcwd()

dataFolder = "Data"
datasetFolder = "Power-Law"
saveFilename = "Figures/Power-Law.pdf"

with shelve.open(os.path.join(mainFolder, dataFolder, datasetFolder, datasetFolder + "_eigenvalues")) as data:
    assortativities = data["rho"]
    meanFieldEigenvalues = data["mean-field-eigenvalues"]
    trueEigenvalues = data["true-eigenvalues"]
    originalAssortativity = data["original-assortativity"]
    originalEigenvalue = data["original-eigenvalue"]

plt.figure()
plt.plot(assortativities, trueEigenvalues, 'ko-', label="CEC Eigenvalue (True)")
plt.plot(assortativities, meanFieldEigenvalues, 'ro-', label="CEC Eigenvalue (Mean Field)")
plt.scatter(originalAssortativity, originalEigenvalue, s=500, marker="x", linewidth=2, color="blue")
plt.legend()
plt.xlabel(r"$\rho$", fontsize=18)
plt.ylabel(r"$\lambda$", fontsize=18)
plt.ylim([0, 1.1*max(max(trueEigenvalues), max(meanFieldEigenvalues))])
plt.tight_layout()
plt.savefig(saveFilename, dpi=1000)
plt.show()
#%%
