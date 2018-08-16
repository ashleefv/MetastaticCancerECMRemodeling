# MetastaticCancerECMRemodeling
PDE-based model of ECM remodeling and its impact on metastatic cancer migration

[![DOI](https://zenodo.org/badge/125926100.svg)](https://zenodo.org/badge/latestdoi/125926100)

## Overview
This mathematical model describes the dynamic process of extracellular matrix remodeling in the vicinity of a primary tumor at the onset 
of metastasis. The biological network under investigation consists of cancer cells, two populations of collagen fibers,
and two enzymes that react to remodel the microenvironment, impacting cancer cell migration away from the primary tumor.

## Mathematical Model for Metastatic Cancer Migration through a Remodeling Extracellular Matrix
### Code Authors
Yen T. Nguyen Edalgo and Ashlee N. Ford Versypt, 
School of Chemical Engineering,
Oklahoma State University.
Corresponding author: A. N. Ford Versypt, ashleefv@okstate.edu

### Related publication for model details
Yen T. Nguyen Edalgo and Ashlee N. Ford Versypt, Mathematical Modeling of Metastatic Cancer Migration through a Remodeling Extracellular Matrix, Processes, 6, 58, 2018. https://doi.org/10.3390/pr6050058

### Scripts

* solve_pdepd_CancerECM.m
This file includes the model definition and the parameters. The model is solved using MATLAB's built in pdepe function.
It can be run as a stand alone script or to be called as a subroutine by localsensitivity_CancerECM.m

* localsensitivity_CancerECM.m
This file conducts a local sensitivity analysis of the model.

### Recommended Supplementary Package
export_fig for exporting MATLAB figures: https://github.com/altmany/export_fig
