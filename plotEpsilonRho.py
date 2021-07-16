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
yMin = 0
yMax = 25
#%%
print(dataMin)
#%%

plt.figure()
ax = plt.gca()
pts = plt.scatter(assortativities, predictedAssortativities, c=epsilons, cmap="viridis")
plt.plot([dataMin, dataMax], [dataMin, dataMax], 'k--')
plt.xlim([xMin, xMax])
# plt.ylim([yMin, yMax])

axins = ax.inset_axes([0.5, 0.5, 0.47, 0.47])
axins.scatter(assortativities, predictedAssortativities, c=epsilons, cmap="viridis")
axins.plot([dataMin, dataMax], [dataMin, dataMax], 'k--')
x1, x2, y1, y2 = -0.06, 0.15, -0.2, 4
axins.set_xlim(x1, x2)
axins.set_ylim(y1, y2)

ax.indicate_inset_zoom(axins, edgecolor="black")

# plt.ylim([-5, 27])
plt.colorbar(pts)
plt.xlabel(r"$\rho$", fontsize=18)
plt.ylabel(r"$\rho_{MF}$", fontsize=18)
plt.savefig("Figures/epsilon-rho.png", dpi=1000)
plt.savefig("Figures/epsilon-rho.pdf", dpi=1000)
plt.show()
