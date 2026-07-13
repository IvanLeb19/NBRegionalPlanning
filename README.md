# NBRegionalPlanning
Simulation tool for regional sample size planning and treatment effect consistency assessment in MRCTs under negative binomial models.
# RegionSizeR (negative binomial, rate ratio)
English | [Русский]([README.ru.md](https://github.com/IvanLeb19/NBRegionalPlanning/blob/main/README_RU.md))

Simulation tool for planning the regional sample size of a multi-regional
clinical trial (MRCT) and estimating the probability of demonstrating
consistency of the treatment effect in a region of interest. The endpoint is an
event count modelled with a negative binomial distribution, the target of
inference is the rate ratio (treatment vs control), the hypothesis is
non-inferiority, and the design is repeated measures (one group measured twice,
before and after the intervention) or parallel.
## The problem
When a drug is developed through a single global trial, regulators in some
regions (Japan in particular) ask for evidence that the effect seen in their
region is consistent with the overall effect, rather than an accident of a small
local sample. ICH E17 frames this as an evaluation of consistency across regions;
the Japanese MHLW guidance operationalises it with the preservation-of-effect
criterion: the region should retain a fixed fraction of the overall effect, most
conservatively at least 50%.
The practical planning question is then: how many subjects does the region need
so that, if the drug truly works, the trial has a high probability of meeting the
consistency criterion? This tool answers that by Monte Carlo simulation. It draws
trials under the assumed effect, fits the model, and counts how often both the
overall non-inferiority test succeeds and the region retains at least the target
fraction of the effect. That frequency is the consistency probability, the
quantity used to size the region.
## Method
The design follows the preservation-of-effect approach described in ICH E17 and
the MHLW guidance, using the consistency-probability framework of Quan et al.
(2010). It is a negative binomial, repeated-measures extension in the spirit of
the RegionSizeR application of Sun et al. (2024).
Consistency is measured as the fraction of the overall effect retained in the
region:
- log scale: `(logRR_cn - log(nim)) / (logRR_main - log(nim))`
- rate-ratio scale: `(1 - RR_cn/nim) / (1 - RR_main/nim)`
where `_cn` denotes the region of interest, `_main` the overall population, and
`nim` the non-inferiority margin. A simulated trial is consistent when this
fraction exceeds the threshold `thd` (0.5 for the Japanese rule).
## Two analysis modes
- **`correlated`** selects the model. `FALSE` treats the two measurements as
  independent arms and fits `glm.nb` with a profile-likelihood confidence
  interval. `TRUE` treats them as a within-subject pair and fits a negative
  binomial GEE (`geeM`) with an unstructured working correlation and robust
  (sandwich) standard errors on a Wald interval.
- **`conservative`** selects the denominator. `FALSE` divides by the number of
  converged fits. `TRUE` counts every non-converged fit as a study failure,
  divides by the full number of simulations, and reports the failure count.
## Requirements
R (>= 4.1) and the following packages:
```r
install.packages(c("MASS", "simstudy", "geeM", "tidyverse", "flextable"))
```
`MASS` and base `stats` cover the independent mode. `simstudy` and `geeM` are only
needed for the paired (`correlated = TRUE`) mode. `tidyverse` and `flextable` are
used for the results table.
## Usage
Open `RegionSizeR_negbin_RR.Rmd` and set the two switches in the `options` chunk:
```r
correlated   <- FALSE   # FALSE = independent arms, TRUE = paired GEE
conservative <- FALSE   # FALSE = converged fits only, TRUE = failures count
```
Then knit to HTML. One run produces one results table.
To run without editing the file, knit with parameters:
```r
rmarkdown::render(
  "RegionSizeR_negbin_RR.Rmd",
  params = list(correlated = TRUE, conservative = TRUE)
)
```
To sweep all four combinations:
```r
grid <- expand.grid(correlated = c(FALSE, TRUE), conservative = c(FALSE, TRUE))
for (i in seq_len(nrow(grid))) {
  rmarkdown::render(
    "RegionSizeR_negbin_RR.Rmd",
    params = list(correlated = grid$correlated[i],
                  conservative = grid$conservative[i]),
    output_file = sprintf("result_%d.html", i)
  )
}
```
## Parameters
Set in the `params` chunk.
| Parameter   | Meaning                                                            |
|-------------|-------------------------------------------------------------------|
| `n`         | Subjects. 58 for the independent mode, 29 for the paired mode (both are 58 observations). Set automatically from `correlated`. |
| `ratio`     | Allocation ratio between arms (independent mode only).            |
| `ratio_cn`  | Allocation ratio inside the region of interest.                   |
| `pct`       | Region-of-interest share of the total sample.                     |
| `pct_in`    | Share of region subjects enrolled in the global trial.            |
| `mean_trt`  | Mean event rate under treatment.                                  |
| `mean_ctrl` | Mean event rate under control.                                    |
| `nsim`      | Number of simulations.                                            |
| `SigL`      | Significance level.                                               |
| `thd`       | Consistency threshold (fraction of effect retained; 0.5 for Japan).|
| `od_ctrl`   | Overdispersion Var(X)/mu, control. Must be greater than 1.        |
| `od_trt`    | Overdispersion Var(X)/mu, treatment. Must be greater than 1.      |
| `nim`       | Non-inferiority margin on the rate ratio scale.                   |
| `hma`       | Hypothesis direction: 1 = higher is worse, 0 = better.            |
| `rho`       | Within-subject correlation (paired mode only).                    |
| `seed`      | Random seed.                                                      |
## Output
| Field             | Meaning                                                       |
|-------------------|---------------------------------------------------------------|
| `n_global`        | Global sample size.                                           |
| `n_chinese`       | Region-of-interest sample size.                               |
| `pct`, `pct_in`   | Region shares, in percent.                                    |
| `conditional_p`   | Consistency probability (log scale), conditional on a significant global test. |
| `unconditional_p` | Consistency probability (log scale), unconditional.           |
| `cond_p_red`      | As `conditional_p`, on the rate-ratio scale.                  |
| `uncond_p_red`    | As `unconditional_p`, on the rate-ratio scale.                |
| `power_global`    | Probability the global non-inferiority test is significant.   |
| `power_cn`        | The same, restricted to the region.                           |
| `threshold`       | Consistency threshold used.                                   |
| `n_failed`        | Non-converged simulations (conservative mode only).           |
## Reproducibility
The seed is fixed inside the simulation driver, so results are reproducible across
runs. The `conservative` switch changes only
post-processing (the denominator), so it does not alter the simulated data, and
the conservative and non-conservative modes share identical draws.
## References
1. ICH E17 (2017). General Principles for Planning and Design of Multi-Regional
   Clinical Trials. International Council for Harmonisation. Available from the ICH
   website (ich.org).
2. ICH E5(R1) (1998). Ethnic Factors in the Acceptability of Foreign Clinical
   Data. International Council for Harmonisation. Available from the ICH website
   (ich.org).
3. Ministry of Health, Labour and Welfare of Japan (2007). Basic Concepts for
   Joint International Clinical Trials (Basic Principles on Global Clinical
   Trials), 28 September 2007.
   https://www.pmda.go.jp/files/000153265.pdf
4. Quan H, Li M, Chen J, et al. (2010). Assessment of Consistency of Treatment
   Effects in Multiregional Clinical Trials. Therapeutic Innovation & Regulatory
   Science (Drug Information Journal) 44(5): 617-632.
   https://link.springer.com/article/10.1177/009286151004400509
5. Sun G, Sun X, Su H, et al. (2024). RegionSizeR - A Novel App for Regional
   Sample Size Planning in MRCTs. Therapeutic Innovation & Regulatory Science.
   https://pmc.ncbi.nlm.nih.gov/articles/PMC11530476/
   App: https://github.com/rsr-ss/RegionSizeR
## Files
- `RegionSizeR_negbin_RR.Rmd` - the merged report.
- `README.md` - this file (English).
- `README.ru.md` - Russian version.
