# Examples of Typical Tasks for Handling ODE Systems Biology Models

## Directory: LV_param_estimation
This directory contains scripts for estimating parameters of an ODE model using Lotka-Volterra model as an example.
Script `simulate_experimental_data.R` can be used to simulate the data. 

![Untransformed simulated data.](./LV_param_estimation/lv_sim1_nontransformed.png)

The simulated data saved in `lv_nontransformed.txt` file.
Script `parameter_optimization_nonequally_spaced_timepoints.py` performs parameter optimization.

![par_opt_1](./LV_param_estimation/par_opt_1.png)

---

The following steps in the script perform a few transforms to make the data look more like real experimental biological data.
