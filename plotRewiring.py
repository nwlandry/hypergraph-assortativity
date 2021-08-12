import matplotlib.pyplot as plt
import os
import shelve
import numpy as np
from HyperEigenvalues import SparseTensor

mainFolder = os.getcwd()
dataFolder = "Data"
datasetFolder = "congress-bills"
# datasetFolder = "tags-ask-ubuntu"
# datasetFolder = "email-Enron"
# datasetFolder = "Eu-Emails"

with shelve.open(os.path.join(mainFolder, dataFolder, datasetFolder, datasetFolder + "_epidemics")) as data:
    rho = data["assortativities"]
    equilibria = data["equilibria"]
    eigenvalues = data["eigenvalues"]


upperBound = np.power(eigenvalues, -1)
fixedBeta3 = 2*upperBound[-1]*np.ones(len(rho))

plt.figure()
plt.subplot(211)
plt.plot(rho, upperBound, 'ko-', label="Upper bound for extinction")
plt.plot(rho, fixedBeta3, 'k--', label=r"$\beta_3^{fixed}$")
plt.ylim([0, 1.1*max(max(upperBound), fixedBeta3[0])])
plt.ylabel(r"$\beta_3$", fontsize=18)
plt.legend(loc="lower left")

plt.subplot(212)
plt.plot(rho, equilibria, 'ko-', label=r"Epidemic extent for $\beta_3^{fixed}$")
plt.ylim([0, 1.1*max(equilibria)])
plt.xlabel(r"$\rho$", fontsize=18)
plt.ylabel(r"$U$", fontsize=18)
plt.legend(loc="upper left")
plt.tight_layout()
plt.savefig("Figures/epidemic_rewiring.png", dpi=1000)
plt.savefig("Figures/epidemic_rewiring.pdf", dpi=1000)
plt.show()