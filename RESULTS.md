# Results

| Horizon | Matched origins | No change | Maturity AR(1) | Dynamic Nelson–Siegel |
|---|---:|---:|---:|---:|
| 1 month | 312 | 25.51 bp | 25.74 bp | 27.26 bp |
| 3 months | 310 | 47.57 bp | 48.48 bp | 53.15 bp |
| 12 months | 301 | 109.55 bp | 113.44 bp | 129.31 bp |

No change beat both alternatives in 17 of the 18 horizon/maturity combinations. The 12-month forecast for the one-year yield was the exception. Maturity AR(1) scored 139.86 bp with a relative RMSE of 0.994. Its 95% paired MSE-gap interval was about −1,183 to +679 bp², which includes zero.

At the 2-year maturity and 3-month horizon, no change beat the other models in the selected financial crisis and 2022 windows. In the March–December 2020 realization window, DNS scored 57.14 bp versus 59.43 for no change (ten realizations). The [complete table](report/accuracy.csv) gives common counts and paired circular-block intervals for the full sample. Stress-window losses are descriptive.

The DNS level, slope and curvature come from projecting the Fed's fitted Svensson curve onto a fixed basis. At very short maturities level plus slope approximates the short yield, and curvature tends to zero. At very long maturities the level remains. This experiment compares three forecast rules with the same revised end-of-month data and expanding monthly windows. Results can change with the estimation window, different factor dynamics, different target vintages, or quoted-instrument inputs.

## Reconstruct one forecast origin

On Dec. 31, 1999, the 2-year yield was 6.183212%. Fit a regression of next month's rate on a constant and current month, using pairs entirely within Jan. 1985–Dec. 1999. The intercept is 0.216762 percentage points and the slope is 0.9634495. Multiply the December rate by that slope and add the intercept to obtain 6.173975%, versus Jan. 31's realized 6.545794%. No change predicts 6.183212%. At that realization the baseline missed by 36.26 bp and AR(1) missed by 37.18 bp. [Origin audit](report/one-origin.json)
