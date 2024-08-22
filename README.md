# ODE Model Fitting Using Genetic Algorithm with Parallel Computing

This repository contains MATLAB code designed to fit an ODE model to cellular fluorescence data using a Genetic Algorithm (GA). The GA is optimized through parallel computing, significantly reducing computation time.

## Overview

The primary objective of this project is to model calcium ion concentration (`Ca^{+2}`) within cells using experimental fluorescence data. The fitting process leverages MATLAB's Genetic Algorithm toolbox, enhanced by parallel computing to optimize multiple cells' parameters simultaneously.

## Process

1. **Data Importing**: Experimental data is imported and organized into time series, cell area, and fluorescence means.
2. **Distance Matrix Calculation**: The script calculates the pairwise Euclidean distances between cell vertices and normalizes them.
3. **Heatmaps**: Visualizations of the distance matrix (`Beta`) and its adjusted version (`Beta dash`) are created to understand cell interactions.
4. **ODE Model Pre-Fitting**: Initial simulation of the calcium concentration using predefined parameters.
5. **Genetic Algorithm Optimization**: Parallel computing is utilized to optimize model parameters for each cell using a genetic algorithm.
6. **ODE Model Post-Fitting**: The optimized parameters are used for final ODE simulation, and results are plotted against experimental data.
7. **Performance Boost**: Parallel computing reduces computational time significantly by dividing tasks across available GPUs.