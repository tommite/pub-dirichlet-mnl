General
-------
This repository contains scripts to execute the tests reported in
T. Tervonen, F. Pignatti, D. Postmus: "From individual to population
preferences: Comparison of discrete choice and Dirichlet models for
treatment benefit-risk trade-offs"; published in Medical Decision Making, 2019.

Requirements
------------
R3.4 or higher. Tested with R3.4.1 on Windows.

The following packages available at CRAN are required:
-plyr
-smaa
-MCMCprecision
-mlogit
-support.CEs
-ggplot2
-reshape2
-gridExtra
-ggthemes
-devEMF
-hitandrun
-evd

Files
-----
err.f.R, load.dce.R, dirichlet-cvm.R, simulate-cbm.R, pinkfloyd-plot.R: Helper files
analysis.mnlscale.R: MNL scale simulation tests (Figure S1)
analysis.convergence.R: Main tests (Figures 1 and 2)
analysis.weightplot.R: Figure 3
analysis.fullsample.R: Table 1
run.all.analyses.R: a main script to execute all analyses

Execution
---------
Working directory should be root of the repository. Then, start R and just

source('run.all.analyses.R')

This will create output in the console for Table 1 and the figures as PDFs.

===
In London, 14 March 2019
Tommi Tervonen
Evidera

