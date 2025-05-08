
# Table of content
- [Reconstruction Framework](#reconstruction-framework)
- [Getting Started](#getting-started)
  - [Quick start: run the reconstruction procedure](#quick-start-run-the-reconstruction-procedure)
  - [Exponential Scaling](#exponential-scaling)


# Reconstruction Framework

The reconstruction procedure is a inversion-problem designed to reconstruct the effective heterogeneity and asymmetric structural connectivity from observed neural activity. In brief, the reconstruction procedure consists of two part: Temporal Reconstruction: we use existed reconstruction framework: Dynamical Differential Covariance to infer the Jacobian matrix from neural activity under the assumption of stable stochastic process. Spatial Reconstruction: by incorporating the symmetric structural connectivity, we further seperate the directionality contributed from effective heterogeneity and asymmetric structural connectivity. 


# Getting Started

## Quick start: run the reconstruction procedure

The following provides a step-by-step instruction to simulate the neural mass model with heterogeneity and asymmetric connectivity and run the reconstruction procedure.
Refer to the Manuscript, the steps below can get the results from Fig. 1 and 2.
1. Open the main performing file, **HetergeneousMainTestScript1.m**.
2. Run the first 4 blocks in the main performing file to create the model parameters, run the simulation and calculate the ground-truth Jacobian matrix.
3. Run the 5th block for Temporal Reconstruction. This call function, **LinearReconst.m** to estimate the Jacobian Matrix.
5. Run the 6th block for Spatial Reconstruction, **RevealHHetero1.m** is called to further separate the Jacobian to effective heterogeneity and asymmetric structural connectivity.
6. Run the following blocks until **Line193** for validation and evaluation in **HetergeneousMainTestScript1.m**.

## Exponential Scaling
7. 

