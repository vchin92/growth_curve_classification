# Multiclass classification of growth curves using random change points and heterogeneous random effects

Author: Vincent Chin, Jarod Lee, Louise M. Ryan, Robert Kohn, and Scott A. Sisson
 
## Installation of software required
1. Install R.

   For Windows: <br/>Download the binary setup file for R [here](https://cran.r-project.org/bin/windows/base/) and open the downloaded .exe file.

   For MacOS: <br/>Download the appropriate version of .pkg file [here](https://cran.r-project.org/bin/macosx/) and open the downloaded .pkg file.

2. Install RStudio. Choose the appropriate installer file for your operating system [here](https://posit.co/products/open-source/rstudio/), download it and then run it to install RStudio.

3. Install Matlab from [here](https://au.mathworks.com/campaigns/products/trials.html).

## Instructions
### Results for simulation experiment
1. Start a Matlab session from the `code` folder and run `simulation_data` in the command window to generate simulated data and Figure 1. A copy of the simulated data is saved in the `data` folder.

2. Run `simulation_results` in the command window to generate results from MCMC algorithm on simulated data. `data_type` on Line 17 can be set as 'df_fixed' to import simulated data generated from a model assuming fixed knot locations or 'df_random' to import simulated data generated from a model assuming random knot locations. `knot_type` on Line 18 can be set as 'fixed' for an estimation model assuming fixed knot locations or 'random' for an estimation model assuming random knot locations. Intermediate results are compiled and saved in `~/results/simulation` folder as `.Rdata` format.


   Intermediate result naming format: `df_mr_4.Rdata` includes all results based on a random knot locations estimation model on 'df_fixed' with `N=4*50=200`.

3. Start an R session from the `code` folder and run `source("simulation_results_plot.R")` in the console to generate Figures 2 and 3.

### Results for application data
1. Start a Matlab session from the `code` folder and run `application_results` in the command window. A copy of the intermediate results is saved in `~/results/application` folder.

2. Start an R session from the `code` folder and run `source("application_results_plot.R")` in the console to generate Figure 6.

3. Run `application_results` in the console of the Matlab session to generate Figures 4 and 5.

## Contact
For questions, comments or remarks about the code please contact Vincent Chin (vincent.chin@unsw.edu.au).
