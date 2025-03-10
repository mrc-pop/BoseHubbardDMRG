# BoseHubbardDMRG

## Introduction

This repository contains the files needed to simulate the 1D Bose-Hubbard model, using the finite-size DMRG algorithm in
order to extract the ground-state energy and other observables optimizing over MPS states. The project was written using
the Julia programming language and the ITensors package. The work was related to the  university course "Numerical Methods
for Physics", held in AY 2023/2024 by Professor Davide Rossini at the University of Pisa.

This work was done by Alessandro Gori ([@nepero27178](https://github.com/nepero27178)) and me.

## Outline of the repository

- The files under `/modules` contain the routines (DMRG algorithms, calculations of expectation values over MPS states, 
plot functions etc.) which are used throughout the other codes.

- The file `convergence.jl` is used to find the optimal values of the algorithm (maximum bond dimension, cutoff to the
truncation, etc.) before performing the subsequent simulations.

- The file `simulations.jl` contains different options to analyze several aspects of the phase diagram of the model. 
Its settings are stored in the `setup/simulations_setup.jl` script.

- The file `analysis.jl` is used to read the data generated in the previous step and make fits and plots for the
quantities of interest, such as the boundaries between the Mott insulating phase and the superfluid phase of the model,
as well as the two-point correlation function.