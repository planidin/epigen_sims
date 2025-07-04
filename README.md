# Epigenetic RI Simulations

The purpose of this code is to simulate reproductive isolation between two demes at a neutral genetic locus, which is in linkage with an epigenetic locus under divergent selection. The simulation code used to perform this analysis is a simulation for the interaction between 3 loci however, the parameters of the third locus are set such that only the neutral locus and epigenetic locus influence RI.

The first step to running the simulation is to install the required R packages.
The packages used for simulations are:
	doParallel
	getopt
	Rcpp
The packages used for figures are:
	ggplot2
	cowplot
	viridis
	dplyr
	tidyr
	EigenR
	Matrix

These can be installed with the following line of code

```
install.packages(c(“doParallel”,”getopt”,”Rcpp”,”ggplot2”,”cowplot”,”viridis”,”dplyr”,”tidyr”,”EigenR”,”Matrix”))
```

Once Rcpp is installed, you will have to compile the epigen_funcs1.4.cpp to an R package called EpiSims, which can be done using the following code

```
library(Rcpp)
Rcpp.package.skeleton("EpiSims", example_code = F, cpp_files = "epigen_funcs1.4.cpp")
install.packages("./EpiSims",repos=NULL, type = "source")
```

To set up the simulation script you have to change the working directory in epigen_caller_rev2_v4.sh to your working directory
Then, in the terminal, you can run all of the required simulations with

```
sh epigen_caller_rev2_v4.sh
```

To plot the simulation results, you have to change the .csv files that are loaded in analysis_revision.R to those that were just generated.
The order that they are loaded in analysis_rev2_v4.R is the same order that the simulations are ran in epigen_caller_rev2_v4.sh

analysis_rev2_v4.R can then be run interactively or as an Rscript
