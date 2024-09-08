# covariateselection
Evaluation of forward selection procedures for selecting covariates to include in a GLM describing neuronal tuning, used in Nevjen, F. and Dunn, B. (2024) Insights in neuronal tuning: Navigating the statistical challenges of autocorrelation and missing variables.
Methods are tested on the calcium data used by Zong W, Obenhaus HA, Skytøen ER, Eneqvist H, de Jong NL, Vale R, Jorge MR, Moser MB, Moser EI. Large-scale two-photon calcium imaging in freely moving mice. Cell. 2022; 185(7):1240–1256, and on simulated data.

The main script is main.R, intended to be run in chunks, where the calcium data is processed to be used in a GLM framework before the forward selection procedures are run. Core_figures.R contains the code for plotting all the figures contained in Nevjen, F. and Dunn, B. (2023), including those that rely on results from the main script, and also further simulations illustrating various challenges discussed in the article. Major (and minor) functions for running the other two scripts are included in functions_model_selection.R.

Feel free to contact fredrik.nevjen@ntnu.no if you have any questions!
