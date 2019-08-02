Work in progress version of code to run closed-loop simulations of Kukskokwim "like" system to evaluate the performance of alternative management procedures in the face of environmental stochasticity and uncertainty in drivers of recruitment:

>Connors B. M., B. Staton, L. Coggins, C. Walters, M. Jones, D. Gwinn, M. Catalano and S. Fleischman. In preparation. Incorporating harvest-population diversity trade-offs into salmon management in large river basins: a case study of Kuskokwim River Chinook.

Posterior samples required to parameterize the closed loop simulations are available upon request from brendan.connors@dfo-mpo.gc.ca and are from analyses described in:

>Staton B., M. Catalano, B.M. Connors, L. Coggins, M. Jones, C. Walters, S. Fleischman and D. Gwinn. In preparation. Evaluation of methods for quantifying population diversity in mixed-stock Pacific salmon fisheries.

Project made possible by grant 1503 from the [Arctic-Yukon-Kuskokwim Sustainable Salmon Initiative](https://www.aykssi.org/)

## Files
- `make.R`: All code required to reproduce current set of simulations and figures.

- `load.R`: Loads packages and scripts necessary for analysis. This file should be sourced prior to running other scripts.

- `functions.R`: All functions written for the analysis should be placed in this file.
  
- `close_loop_sims.R`: Run closed loop forwad simulations.

- `figures.R`: Generate figures not associated with simulations (e.g., equilibrium trade-offs, populaiton diversity )
  
- `simulation_summary.Rmd`: R Markdown doc that summarizes simulations, and generates figures for manuscript (in `images` folder)

- `appendix_A.R`: R Markdown doc that details of the stationary Ricker to time-varying Beverton-Holt formulation we use and  simulations to justify its parameterization in our closed-loop simulations. 


