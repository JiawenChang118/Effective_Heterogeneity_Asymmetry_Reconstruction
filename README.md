
# Table of content
- [Reconstruction Framework](#reconstruction-framework)
- [Getting Started](#getting-started)
  - [Quick start: run the reconstruction procedure](#quick-start-run-the-reconstruction-procedure)



# Reconstruction Framework

The reconstruction procedure is a inversion-problem designed to reconstruct the effective heterogeneity and asymmetric structural connectivity from observed neural activity. In brief, the reconstruction procedure consists of two part: Temporal Reconstruction: we use existed reconstruction framework: Dynamical Differential Covariance to infer the Jacobian matrix from neural activity under the assumption of stable stochastic process. Spatial Reconstruction: by incorporating the symmetric structural connectivity, we further seperate the directionality contributed from effective heterogeneity and asymmetric structural connectivity. 


# Getting Started

## Quick start: run the reconstruction procedure

The following provides a step-by-step instruction to simulate the neural mass model with heterogeneity and asymmetric connectivity and run the reconstruction procedure.

1. Open the main performing file, **HetergeneousMainTestScript1.m**.
2. Run the first 4 blocks in the main performing file to create the model parameters, run the simulation and calculate the ground-truth Jacobian matrix.
3. Run the 5th block for Temporal Reconstruction. This call function, **LinearReconst.m**.
4. Create two directories, **./checkpoint/** (for saving trained model results) and **./data/** (for downloading the MNIST dataset).
5. Run the following command to call the script named `mnist.py` with the config file specified through the option:

   ```
   python mnist.py --config=./mnist_config.yaml
   ```

After training is finished, you should find four files in the **./checkpoint/mnist/** folder：

- Two '.ph' files which contain the trained model parameters.
- One '.yaml' file which is a copy of the config file used for running the training the model.
- One '.txt' log file that prints the standard output during training (such as model performance).
- One directroy called `mnn_net_snn_result` that stores the simulation result of the SNN reconstructed from the trained MNN (if enabled).


