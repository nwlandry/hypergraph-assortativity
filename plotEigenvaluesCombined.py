#%%
from matplotlib.collections import PolyCollection
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


datasetFolderList = ["Power-Law", "tags-ask-ubuntu", "congress-bills", "Eu-Emails"]
labels = ['(a)', '(b)', '(c)', '(d)']

saveFilename = "Figures/eigenvalue.pdf"


plt.figure(figsize=(6.4, 4.8))

for i in range(4):
    datasetFolder = datasetFolderList[i]
    with shelve.open(os.path.join(mainFolder, dataFolder, datasetFolder, datasetFolder + "_eigenvalues")) as data:
        assortativities = data["rho"]
        meanFieldEigenvalues = data["mean-field-eigenvalues"]
        trueEigenvalues = data["true-eigenvalues"]
        originalAssortativity = data["original-assortativity"]
        originalEigenvalue = data["original-eigenvalue"]

    plt.subplot(2, 2, i + 1)

    plt.plot(assortativities, trueEigenvalues, 'k^-', label=r"$\lambda$", markersize=3, linewidth=0.5)
    plt.plot(assortativities, meanFieldEigenvalues, 'ko-', label=r"$\lambda^{(0)} + \epsilon \lambda^{(1)}$", markersize=3, linewidth=0.5)
    # plt.scatter(originalAssortativity, originalEigenvalue, s=100, c="k", marker="X", linewidth=0.1, markeredgewidth=0.1)
    # plt.scatter(originalAssortativity, originalEigenvalue, s=100, c="r", marker="x", linewidth=1)
    plt.plot(originalAssortativity, originalEigenvalue, markersize=5, marker="s", markerfacecolor="red", markeredgewidth=0.5, markeredgecolor="white")

    if i == 0:
        plt.legend(fontsize=8, loc="lower left")
        arrowdict = dict(facecolor='black', shrink=0.1, width=0.1, headwidth=3, headlength=5)
        plt.annotate(r"Original $(\rho, \lambda)$", xy=(originalAssortativity, originalEigenvalue), xytext=(originalAssortativity + 0.05, originalEigenvalue + 40), arrowprops=arrowdict)
        
    if i == 2:
        plt.xlabel(r"$\rho$", fontsize=12)
        plt.ylabel(r"$\lambda$", fontsize=12)

    plt.ylim([0, 1.1*max(max(trueEigenvalues), max(meanFieldEigenvalues))])

    plt.gca().tick_params(axis='both', which='major', labelsize=6)
    plt.gca().tick_params(axis='both', which='minor', labelsize=6)

    xPos = min(assortativities) + 0.1*(max(assortativities) - min(assortativities))
    yPos = 0.9*max(max(meanFieldEigenvalues), max(trueEigenvalues))

    plt.text(xPos, yPos, labels[i], fontsize=12)

plt.tight_layout()
plt.savefig("Figures/eigenvalue_vs_rho.png", dpi=1000)
plt.savefig("Figures/eigenvalue_vs_rho.pdf", dpi=1000)
plt.show()
#%%
