# Code to Reproduce Simulations

This directory contains code to reproduce the latent class model simulations 
presented in Appendix C, using the `R` package `LCMCR`. The simulations are 
broken up into 7 `R` scripts, each labelled `example_*.R`. The `R` script 
`helper_functions.R` contains helper functions used to run the simulations. 
While your working directory in `R` is the current directory, simply run a given 
`example_*.R` script in order to reproduce the results for that example from the 
paper. However, please note that the package `LCMCR` used in the simulations 
may not play well with your `R` session, and it may look as though your session
is frozen even though the simulation is still running (you can verify this by
checking the simulation results which are saved every 10 iterations in the file 
`example_*_simulation.RData`).


