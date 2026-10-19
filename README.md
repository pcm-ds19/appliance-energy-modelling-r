# Appliance Energy Modelling in R

Analysis of household appliance energy consumption using environmental measurements, data transformations and weighted aggregation models.

Developed from my SIG718 Real World Analytics midterm assessment at Deakin University, with an additional exploratory holdout evaluation for this portfolio.

## Objective

Compare four aggregation models for estimating appliance energy consumption and assess whether they improve on a simple mean baseline.

## Data

The analysis uses a reproducible sample of 650 observations from the course-provided, modified Appliances Energy Prediction dataset.

Variables examined:
- Living-room temperature and relative humidity
- Office-room temperature and relative humidity
- Atmospheric pressure
- Appliance energy consumption, measured in Wh

Four temperature and humidity variables are used as model inputs; pressure is excluded from the fitted models.

## Methods

- Exploratory analysis using summary statistics, histograms and scatter plots
- Log transformations of living-room temperature and energy consumption
- Min–max scaling
- Weighted Arithmetic Mean (WAM)
- Weighted Power Mean with p = 0.5 and p = 2
- Ordered Weighted Average (OWA)
- Model comparison using RMSE and MAE

The aggregation models use non-negative weights constrained to sum to one. OWA weights apply to sorted input positions rather than named predictors.

## Exploratory Holdout Evaluation

The 650 observations were split into 520 training and 130 validation rows using a fixed random seed. Scaling parameters and aggregation weights were fitted using training data only.

Predictions were converted back to Wh before evaluation. The baseline predicts the training-set mean energy consumption for every validation observation.

| Model | Validation RMSE (Wh) | Validation MAE (Wh) |
| --- | ---: | ---: |
| Training-mean baseline | 92.46 | 56.74 |
| WAM | 110.83 | 66.35 |
| WPM (p = 0.5) | 110.25 | 65.78 |
| WPM (p = 2) | 115.95 | 68.85 |
| OWA | 114.64 | 66.20 |

Lower values indicate better performance.

**Finding:** WPM (p = 0.5) performed best among the four aggregation models, but none outperformed the training-mean baseline on this validation split.

These results do not support replacing the baseline with an aggregation model for this setup.

## Original Assessment Analysis

The original analysis fitted and compared models using all 650 sampled observations.

WPM (p = 0.5) achieved a transformed-scale fitting RMSE of 0.16930. For the assessment's supplied prediction case, it estimated 125.97 Wh against a measured value of 120 Wh.

This single-case result is separate from the holdout evaluation and does not establish general predictive accuracy.

## Limitations

- Validation uses one random split and does not assess future-period forecasting.
- Exploratory analysis and transformation choices preceded the holdout extension; this is not an untouched final test.
- Data represent one household, limiting generalisation.
- One validation input slightly exceeded the training scaling range and was retained without clipping.
- Model fitting objectives differ from the reported original-scale evaluation metrics.

## Tools

R, base R graphics and lpSolve, with course-provided aggregation functions.

## Project Status

Documentation is being added first. The reviewed R script and supporting reproducibility instructions will follow.

## Acknowledgements

Original dataset research: Candanedo, Feldheim and Deramaix (2017), “Data driven prediction models of energy use of appliances in a low-energy house.”

Course context: SIG718 Real World Analytics, Deakin University.

Author: Priyanka Chandrahar Mane
