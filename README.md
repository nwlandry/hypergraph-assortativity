This repository contains the code for the forthcoming paper *Hypergraph assortativity: the effect on the largest eigenvalue and a strategy for suppressing epidemics* by Nicholas W. Landry and Juan G. Restrepo.

The *Figures* folder contains all the figures used in the paper, the *UNUSED* folder contains relevant but unused code, and the (not-included for space reasons) *Data* folder contains the datasets on which our analysis runs.

A description of the files in the main directory:
- ConfigModelEdgeSwapsEigenvaluesInParallel.py: This file generates the $$\lambda$$ vs. $$\rho$$ plot for Fig. 2 in the main text for a configuration null model.
- EmpiricalHypergraphsEdgeSwapsEigenvaluesInParallel.py: This file generates the $$\lambda$$ vs. $$\rho$$ plot for Fig. 2 in the main text for the datasets from Austin Benson's [data repository](https://www.cs.cornell.edu/~arb/data/)
- GetAssortativeHypergraphsInParallel.py: This file preferentially rewires a given hypergraph for various target $$\rho$$ values and stores each resulting hypergraph as an item in a list. These hypergraphs are used in generating Fig. 3 in the main text and the figures in the Supplemental Material.
- HyperEigenvalues.py: This file contains the SparseTensor class which is used to compute the true eigenvalue that is plotted in Fig. 2 and indirectly used in Fig. 3.
- Hypergraph.py: This file contains the Hypergraph class (used to represent hypergraphs for the contagion simulations) and the HypergraphGenerator class (used to generate configuration model hypergraphs as well as perform preferential edge swaps).
- MeanFieldTheory.py: This file is used to calculate the approximate eigenvalue used in Fig. 2 and calculate the assortativity as well as other (unused) miscellaneous functions.
- PlotEigenvaluesCombined.py: Generates Fig. 2 from the data generated with ConfigModelEdgeSwapsEigenvaluesInParallel.py and EmpiricalHypergraphsEdgeSwapsEigenvaluesInParallel.py
- PlotRewiring.py: Generates Fig. 3 from the data generated with GetAssortativeHypergraphsInParallel.py and RunContagionInParallel.py.
- utilities.py: Imports the datasets from Austin Benson's [data repository](https://www.cs.cornell.edu/~arb/data/) and filters by hyperedge size.