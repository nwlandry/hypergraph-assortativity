import matplotlib.pyplot as plt
import os
import shelve
import numpy as np

mainFolder = os.getcwd()
dataFolder = "Data"


with shelve.open(os.path.join(mainFolder, dataFolder, "degree_distribution")) as data:
    k_CM = data["CM"]
    k_TAU = data["TAU"]
    k_CB = data["CB"]
    k_EE = data["EE"]

# label fontsize=24
fs = 24
tfs = 16

#CM
plt.figure(figsize=(6, 3))
bins = 10**(np.linspace(0, 4, 100))
plt.hist(k_CM, density=True, log=True, color="black", bins=bins)
plt.xscale('log')
plt.xlabel(r"$k^{(3)}$", fontsize=fs)
plt.ylabel(r"$P(k^{(3)})$", fontsize=fs)
plt.xticks(fontsize=tfs)
plt.yticks(fontsize=tfs)
plt.tight_layout()
plt.savefig("Figures/degree_distribution_CM.png", dpi=1000)
plt.savefig("Figures/degree_distribution_CM.pdf", dpi=1000)
plt.show()


#TAU
plt.figure(figsize=(6, 3))
plt.hist(k_TAU, density=True, log=True, color="black", bins=bins)
plt.xscale('log')
plt.xlabel(r"$k^{(3)}$", fontsize=fs)
plt.ylabel(r"$P(k^{(3)})$", fontsize=fs)
plt.xticks(fontsize=tfs)
plt.yticks(fontsize=tfs)
plt.tight_layout()
plt.savefig("Figures/degree_distribution_TAU.png", dpi=1000)
plt.savefig("Figures/degree_distribution_TAU.pdf", dpi=1000)
plt.show()


#CB
plt.figure(figsize=(6, 3))
plt.hist(k_CB, density=True, log=True, color="black", bins=bins)
plt.xscale('log')
plt.xlabel(r"$k^{(3)}$", fontsize=fs)
plt.ylabel(r"$P(k^{(3)})$", fontsize=fs)
plt.xticks(fontsize=tfs)
plt.yticks(fontsize=tfs)
plt.tight_layout()
plt.savefig("Figures/degree_distribution_CB.png", dpi=1000)
plt.savefig("Figures/degree_distribution_CB.pdf", dpi=1000)
plt.show()


#EE
plt.figure(figsize=(6, 3))
plt.hist(k_EE, density=True, log=True, color="black", bins=bins)
plt.xscale('log')
plt.xlabel(r"$k^{(3)}$", fontsize=fs)
plt.ylabel(r"$P(k^{(3)})$", fontsize=fs)
plt.xticks(fontsize=tfs)
plt.yticks(fontsize=tfs)
plt.tight_layout()
plt.savefig("Figures/degree_distribution_EE.png", dpi=1000)
plt.savefig("Figures/degree_distribution_EE.pdf", dpi=1000)
plt.show()