
# Table of content
- [Reconstruction Framework](#Reconstruction Framework)
- [Getting Started](#getting-started)
  - [Quick start: run the reconstruction procedure](#quick-start-run-the-reconstruction-procedure)



# Reconstruction Framework

The moment neural network is a type of second-order artificial neural network model designed to capture the nonlinear coupling of correlated activity of spiking neurons. In brief, the moment neural networks extend conventional rate-based artificial neural network models by incorporating the covariance of fluctuating neural activity. This repository provides a comprehensive framework for simulating and training moment neural networks based on the standard workflow of Pytorch. 


## The architecture of this repository

* `mnn_core`: core modules implementing the moment activation and other building blocks of MNN.
* `models`: a module containging various network architectures for fast and convenient model construction
* `snn`: modules for reconstructing SNN from MNN and for simulating the corresponding SNN in a flexible manner.
* `utils`: a collection of useful utilities for training MNN (ANN compatible).

# Getting Started

## Quick start: run the reconstruction procedure

The following provides a step-by-step instruction to train an MNN to learn MNIST image classification task with a multi-layer perceptron structure.

1. Clone the repository to your local drive.
2. Copy the demo files, **./example/mnist/mnist.py** and **./example/mnist/mnist_config.yaml** to the root directory.
3. Create two directories, **./checkpoint/** (for saving trained model results) and **./data/** (for downloading the MNIST dataset).
4. Run the following command to call the script named `mnist.py` with the config file specified through the option:

   ```
   python mnist.py --config=./mnist_config.yaml
   ```

After training is finished, you should find four files in the **./checkpoint/mnist/** folder：

- Two '.ph' files which contain the trained model parameters.
- One '.yaml' file which is a copy of the config file used for running the training the model.
- One '.txt' log file that prints the standard output during training (such as model performance).
- One directroy called `mnn_net_snn_result` that stores the simulation result of the SNN reconstructed from the trained MNN (if enabled).


