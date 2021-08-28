import matplotlib.pyplot as plt
import os
import shelve
import numpy as np
from HyperEigenvalues import SparseTensor

mainFolder = os.getcwd()
dataFolder = "Data"

# # congress-bills
# datasetFolder = "congress-bills"
# qualifier = ""
# height = 0.6
# ylim = [0, 2.5]
# hasLegend = True

# # Power-Law
# datasetFolder = "Power-Law"
# qualifier = ""
# height = 0.6
# ylim = [0, 3]
# hasLegend = False

# # Eu-Emails small beta
# datasetFolder = "Eu-Emails"
# qualifier = "_small_beta"
# height = 0.65
# ylim = [0, 1]
# hasLegend = False

# # Eu-Emails medium beta
# datasetFolder = "Eu-Emails"
# qualifier = "_medium_beta"
# height = 0.65
# ylim = [0, 1]
# hasLegend = False

# # Eu-Emails large beta
# datasetFolder = "Eu-Emails"
# qualifier = "_large_beta"
# height = 0.5
# ylim = [0, 8]
# hasLegend = False

# # tags-ask-ubuntu small beta
# datasetFolder = "tags-ask-ubuntu"
# qualifier = "_small_beta"
# height = 0.25
# ylim = [0, 1]
# hasLegend = False

# # tags-ask-ubuntu large beta
# datasetFolder = "tags-ask-ubuntu"
# qualifier = "_large_beta"
# height = 0.75
# ylim = [0, 2]
# hasLegend = False

with shelve.open(os.path.join(mainFolder, dataFolder, datasetFolder, datasetFolder + "_epidemics" + qualifier)) as data:
    rho = data["assortativities"]
    equilibria = data["equilibria"]
    eigenvalues = data["eigenvalues"]
    gamma = data["gamma"]
    beta3 = data["beta3"]

upperBound = gamma*np.power(eigenvalues, -1)
beta3cFraction = np.divide(beta3, upperBound)
print("beta3 = " + str(beta3))

plt.figure(figsize=(6.4, 4.8))
plt.subplot(211)

plotHeight =  1.1*max(max(beta3cFraction), 1)

plt.ticklabel_format(axis="y", style="sci", scilimits=(-2, 0)) 

plt.plot(rho, np.divide(beta3, upperBound), 'ko-', label=r"$\beta_3/\beta_3^c$")
plt.plot(rho, np.ones(len(rho)), 'k--', label=r"$\beta_3=\beta_3^c$")
plt.fill_between(rho, np.zeros(len(rho)), np.ones(len(rho)), facecolor="lightgrey")

plt.xlim([min(rho), max(rho)])
plt.ylim([0, plotHeight])
plt.ylabel(r"$\beta_3/\beta_3^c$", fontsize=18)
xPos = min(rho) + 0.05*(max(rho) - min(rho))
yPos = height*plotHeight
plt.text(xPos, yPos, "(a)", fontsize=20)
if hasLegend:
    plt.legend(loc="lower left", facecolor="none", edgecolor="black")

plt.subplot(212)
plt.ticklabel_format(axis="y", style="sci", scilimits=(-2, 0)) 

plt.plot(rho, 100*np.array(equilibria), 'ko-', label=r"Epidemic equilibrium given $\beta_3$")
yPos = height*ylim[1]
plt.xlim([min(rho), max(rho)])
plt.ylim(ylim)
plt.xlabel(r"$\rho$", fontsize=18)
plt.ylabel("% infected", fontsize=18)
plt.text(xPos, yPos, "(b)", fontsize=20)
if hasLegend:
    plt.legend(loc="upper left", facecolor="none", edgecolor="black")
plt.tight_layout()
plt.savefig("Figures/epidemic_rewiring_" + datasetFolder  + qualifier + ".png", dpi=1000)
plt.savefig("Figures/epidemic_rewiring_" + datasetFolder  + qualifier + ".pdf", dpi=1000)
plt.show()