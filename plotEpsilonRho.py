#%%
import shelve
import matplotlib.pyplot as plt
import os
import numpy as np

mainFolder = os.getcwd()

dataFolder = "Data"
datasetFolder = "EpsilonRho"

with shelve.open(os.path.join(mainFolder, dataFolder, datasetFolder, datasetFolder + "_generative")) as data:
    epsilons = data["epsilon"]
    assortativities = data["rho"]
    predictedAssortativities = data["predicted-rho"]

dataMin = min(np.nanmin(assortativities), np.nanmin(predictedAssortativities))
dataMax = max(np.nanmax(assortativities), np.nanmax(predictedAssortativities))
xMin = np.nanmin(assortativities)
xMax = np.nanmax(assortativities)

plt.figure()
ax = plt.gca()
pts = plt.scatter(assortativities, predictedAssortativities, c=epsilons, cmap="viridis")
plt.plot([dataMin, dataMax], [dataMin, dataMax], 'k--')
cb = plt.colorbar(pts)
cb.set_label(label=r'$\epsilon$', weight='bold', fontsize=18)
plt.xlabel(r"$\rho$", fontsize=18)
plt.ylabel(r"$\rho_{MF}$", fontsize=18)
plt.savefig("Figures/epsilon-rho.png", dpi=1000)
plt.savefig("Figures/epsilon-rho.pdf", dpi=1000)
plt.show()
