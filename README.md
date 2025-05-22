# Universal-Inference
- Application experiments using Universal Inference

## LASSO
### $H_0$ Null Hypothesis
- testing the significance of residual
- $\beta_0 \neq 0$, $\beta_1 = \beta_2 = ... = \beta_k = 0$


### Experiment Implementation
- `./run-lasso-exp.sh`: to run the simulation in parallel (38 cores setting)
- the results will be saved to 'results.csv'
- run `python3 csv2boxplot.py --csv_file 'results.csv'` to plot the results stored in csv format as boxplot visualization.
- 