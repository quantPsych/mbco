# Improved Inference in Mediation Analysis: Introducing the Model-based Constrained Optimization Procedure 

This repository contains the supplemental materials for the manuscript entitled *Improved Inference in Mediation Analysis: Introducing the Model-based Constrained Optimization Procedure*. Here is a list of the documents:

- `SuppMaterials-MBCO-rev1.Rmd`: The document explains the R code used in the manuscript. We discuss in detail the steps required to analyze our empirical example and re-produce the results.
- `memory_study.csv`: The data file that contains the empirical data example. When using the data, please cite the relevant study by MacKinnon et al. (2018).
- `rsq.R`: This file contain the `rsq()` function that produces R^2^. 
- `simulation` folder contains the simulation functions used in the study.
  - The main two files for the single-mediator and two-mediator model are `run_sim_1m.R` and `run_sim_2m.R`, respectively. The other files in the folder are read into these two main files during the simulation runs. 

## References
MacKinnon, D. P., Valente, M. J., & Wurpts, I. C. (2018). Benchmark validation of statistical models: Application to mediation analysis of imagery and memory. _Psychological Methods_, 23, 654–671. [https://doi.org/10.1037/met0000174](https://doi.org/10.1037/met0000174)
