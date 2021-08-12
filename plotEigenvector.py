import matplotlib.pyplot as plt
import os
import shelve
import numpy as np

mainFolder = os.getcwd()
dataFolder = "Data"
datasetFolder = "Eigenvector"

with shelve.open(os.path.join(mainFolder, dataFolder, datasetFolder, datasetFolder)) as data:
    zerothOrder = data["zeroth-order-eigenvector"]
    firstOrder = data["first-order-eigenvector"]
    trueEigenvector = data["actual-eigenvector"]


dataMin = min(np.min(trueEigenvector), np.min(zerothOrder))
dataMax = max(np.max(trueEigenvector), np.max(zerothOrder))

plt.figure()
plt.scatter(zerothOrder, trueEigenvector, color="black", marker="x", label="0th order approximation", s=20)
plt.scatter(firstOrder, trueEigenvector, color="black", marker="o", label="1st order approximation", s=20)
plt.plot([dataMin, dataMax], [dataMin, dataMax], 'k--')
plt.xlabel(r"$u^{(MF)}_i$", fontsize=18)
plt.ylabel(r"$u_i$", fontsize=18)
plt.legend(loc="upper left")
plt.tight_layout()
plt.savefig("Figures/eigenvector_approximation.png", dpi=1000)
plt.savefig("Figures/eigenvector_approximation.pdf", dpi=1000)
plt.show()