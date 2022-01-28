This repository contains the code for the paper *Hypergraph assortativity: a dynamical systems perspective* by Nicholas W. Landry and Juan G. Restrepo.

The *Figures* folder contains all the figures used in the paper, the *UNUSED* folder contains relevant but unused code, and the (not-included because of Github storage limitations) *Data* folder contains the datasets on which our analysis runs.

A description of the scripts in the main directory:
- configuration_model_eigenvalues_in_parallel.py: This script generates the data for the lambda vs. rho plot for Fig. 2 in the main text for a configuration null model.
- empirical_hypergraph_eigenvalues_in_parallel.py: This script generates the data for the lambda vs. rho plot for Fig. 2 in the main text for the datasets from Austin Benson's [data repository](https://www.cs.cornell.edu/~arb/data/)
- get_degree_distributions.py: This script stores the 3rd order degree sequence, the mean degree, and size of the hypergraph datasets used and stores them as a shelve file.
- HyperEigenvalues.py: This file contains the SparseTensor class which is used to compute the true eigenvalue that is plotted in Fig. 2 and indirectly used in Fig. 3.
- Hypergraph.py: This file contains the Hypergraph class (used to represent hypergraphs for the contagion simulations) and the HypergraphGenerator class (used to generate configuration model hypergraphs as well as perform preferential edge swaps).
- HypergraphContagion.py: This file contains the functions for running the hypergraph SIS model and returning the equilibrium value of an epidemic simulation.
- MeanFieldTheory.py: This file is used to calculate the approximate eigenvalue used in Fig. 2.
- plot_degree_distributions: This script generates the degree distribution plots in Table 1.
- plot_eigenvalues_vs_rho.py: Generates Fig. 2 from the data generated with configuration_model_eigenvalues_in_parallel.py and empirical_hypergraph_eigenvalues_in_parallel.py
- plot_epidemic_rewiring.py: Generates Fig. 3 (and Fig. 4) from the data generated with RewireHypergraphsInParallel.py and RunContagionInParallel.py.
- rewire_hypergraphs_in_parallel.py: This script generates hypergraphs from an original hypergraph dataset with different assortative structure constructed with double-edge swaps on which to run SIS contagion with the run_contagion_in_parallel.py script.
- run_contagion_in_parallel.py: This script runs SIS contagion on each assortative hypergraph generated with the rewire_hypergraphs_in_parallel.py script.
- utilities.py: Imports the datasets from Austin Benson's [data repository](https://www.cs.cornell.edu/~arb/data/) and filters a hyperedge list by hyperedge size.