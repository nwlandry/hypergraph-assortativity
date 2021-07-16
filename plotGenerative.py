import shelve
import matplotlib.pyplot as plt

mainFolder = os.getcwd()

dataFolder = "Data"
datasetFolder = "Large-Degrees"

with shelve.open(os.path.join(mainFolder, dataFolder, datasetFolder, datasetFolder + "_eigenvalues")) as data:
    assortativities = data["rho"]
    epsilons = data["epsilon"]
    meanFieldCECEigenvalues = data["mean-field-eigenvalues"]
    trueCECEigenvalues = data["true-eigenvalues"]

plt.figure()
pts = plt.scatter(meanFieldCECEigenvalues, trueCECEigenvalues, c=epsilons, cmap="viridis")
plt.plot([min(trueCECEigenvalues), max(trueCECEigenvalues)], [min(trueCECEigenvalues), max(trueCECEigenvalues)], 'k--', label="Perfect")
plt.legend()
plt.colorbar(pts)
plt.xlabel(r"$\rho$")
plt.ylabel(r"$\lambda$")
plt.show()
