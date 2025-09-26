# Efficient_Response_to_Realtime_Price

_This work presents a computationally efficient method for the flexible response of thermostatically controlled loads (TCLs) to uncertain real-time prices. We formulate the demand response problem using stochastic dynamic programming with analytically piecewise linear value functions. This achieves extreme computation performance by simplifying the multi-stage, stochastic, nonlinear optimization problem to a single-stage, deterministic, mixed-integer linear type._

Codes for Submitted Paper "Efficient Response of Thermostatically Controlled Loads to Uncertain Real-time Prices".

Authors: Xueyuan Cui, Liudong Chen, Yi Wang, and Bolun Xu.

## Experiments

### Reproduction
To create the required data, please run
```
cd Codes/
Data_create.py
```
This code will generate both the post-processing real-time price and the probabilistic distributions of prices. Please note that the generated data has been saved in ```Data```for direct use. The raw data of prices can be collected at [NYISO]{https://www.nyiso.com/energy-market-operational-data}.

To reproduce the comparisons between the proposed method and SDDP, please run the matlab codes in ```Codes```: ```Com_main.m``` and ```matlab Com_SDDP.m```. We acknowledge the open-source toolbox named [FAST]{https://stanford.edu/~lcambier/cgi-bin/fast/tuto.php} to support our experiment on SDDP.

To reproduce the experiments of real-time demand response with different comfort functions, please run the matlab codes in ```Codes```: ```Response.m```.

The figures are generated with the results from matlab codes and the python codes in ```Codes```: ```Figures.py```.

## Citation
```
```
