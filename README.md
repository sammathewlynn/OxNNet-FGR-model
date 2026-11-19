# First-trimester FGR prediction model (OxNNet Toolkit – placenta)

A first-trimester competing-risks prediction model for fetal growth restriction
(FGR), combining established maternal/biochemical risk factors with automated
placental ultrasound biomarkers. This repository contains the fitted model only,
as an R object, together with everything needed to compute an individual's
predicted risk. **No individual-level (patient) data are included.**

> **Research use only.** This model is provided to support reproducibility and
> external evaluation. It has not been approved or validated for clinical
> decision-making. Before any use in a new population it should be externally
> validated and, as shown below, will benefit from recalibration to the local
> setting (see *Performance* and *Recalibration*).

---

## Contents

| File | Description |
|------|-------------|
| `FGR_complete_model.rds` | The fitted model (R, serialised with `saveRDS`). |
| `README.md` | This file. |

The `.rds` is a `rms::cph` Fine–Gray (subdistribution hazard) fit, saved
**data-free**. The design matrix, response, per-row linear predictors,
residuals, Fine–Gray weights and missing-data indices have all been removed,
and the formula, terms and call carry no captured R environments (they point
only to the global environment, so no workspace objects are serialised with
the model). It retains the coefficients and their variance–covariance matrix,
the baseline (cumulative incidence) function, the variable definitions, the
restricted-cubic-spline knots and summary ranges for each predictor, which is
everything required to predict on new data.

## Requirements

* R (developed under 4.3.x; the distributed file was saved and tested under
  R 4.4.3 with `rms` 8.0.0)
* R packages: **`rms`**, **`survival`**

## Quick start

```r
library(rms)
library(survival)

model <- readRDS("FGR_complete_model.rds")

# two small helpers
fg_lp        <- function(model, newdata) as.numeric(predict(model, newdata, type = "lp"))
fg_pred_risk <- function(model, newdata, day) as.numeric(1 - rms::survest(model, newdata, times = day)$surv)

# prediction horizons, expressed as gestational age in DAYS
H <- c(`<32+0` = 223, `<37+0` = 258, `<=40+0` = 280)

# one example individual (column names and codings must match the table below)
newdata <- data.frame(
  ASPIRINuse          = 0,
  Race                = factor("White",
                          levels = c("White","Black","East Asian","South Asian","Other or mixed")),
  Smoking             = 0,
  CH                  = 0,
  BMI                 = 24.5,
  SLE.APS             = FALSE,
  prevPExlastBWz      = 0,
  prevFGRxlastBWz     = 0,
  prevBOTHxlastBWz    = 0,
  prevNEITHERxlastBWz = 0,
  map.mom.log10       = 0,
  p.mom.log10         = 0,
  plgf.mom.log10      = 0,
  vol_z               = 0,
  fmbv3_z             = 0,
  ut.mom.log10        = 0
)

fg_lp(model, newdata)                       # linear predictor
fg_pred_risk(model, newdata, H[["<=40+0"]]) # predicted cumulative incidence of FGR by term
```

`predict(type = "lp")` and `survest()` need only the model object. `rms::Predict()`,
`summary()` and `nomogram()` additionally require a `datadist` object, which is
not distributed here; the reference values needed to interpret the linear
predictor are summarised below and in the source publication.

You may see a warning from `survest` that standard errors are *"approximate
except at predictor means"*. This is expected for a data-free object and affects
only the confidence interval of the survival estimate, **not the point risk**,
which is exact.

## Inputs

The model expects a data frame with the following columns. Continuous markers
must be prepared in exactly the same way as in the development study (see
*Citation*); the multiples-of-the-median (MoM) values are entered on the
**log10** scale.

| Variable | Description | Type / coding |
|----------|-------------|---------------|
| `ASPIRINuse` | Aspirin prophylaxis (forced into the model) | 0 = no, 1 = yes |
| `Race` | Maternal race / ethnicity | factor: `White` (reference), `Black`, `East Asian`, `South Asian`, `Other or mixed` |
| `Smoking` | Current smoker | 0 = no, 1 = yes |
| `CH` | Chronic hypertension | 0 = no, 1 = yes |
| `BMI` | Body mass index | continuous (kg/m²) |
| `SLE.APS` | Systemic lupus erythematosus or antiphospholipid syndrome | `FALSE` / `TRUE` (or 0 / 1) |
| `prevPExlastBWz` | Obstetric-history interaction term (see note) | continuous; 0 for nulliparous |
| `prevFGRxlastBWz` | Obstetric-history interaction term (see note) | continuous; 0 for nulliparous |
| `prevBOTHxlastBWz` | Obstetric-history interaction term (see note) | continuous; 0 for nulliparous |
| `prevNEITHERxlastBWz` | Obstetric-history interaction term (see note); restricted cubic spline (3 knots) | continuous; 0 for nulliparous |
| `map.mom.log10` | Mean arterial pressure, log10 MoM | continuous |
| `p.mom.log10` | PAPP-A, log10 MoM | continuous |
| `plgf.mom.log10` | Placental growth factor (PlGF), log10 MoM | continuous |
| `vol_z` | Placental volume z-score (automated, OxNNet) | continuous |
| `fmbv3_z` | svFMBV-M z-score — single-vessel fractional moving blood volume at the myometrium (automated, OxNNet) | continuous |
| `ut.mom.log10` | Uterine artery pulsatility index, log10 MoM; restricted cubic spline (3 knots) | continuous |

**Obstetric-history interaction terms.** The four `prev*xlastBWz` variables
encode the most recent birth's birthweight z-score, split by the woman's
previous-pregnancy history of pre-eclampsia (PE) and/or FGR:

* `prevPExlastBWz` — previous PE only (no previous FGR): set to the last-birth
  birthweight z-score, otherwise 0.
* `prevFGRxlastBWz` — previous FGR only (no previous PE): last-birth z-score, else 0.
* `prevBOTHxlastBWz` — previous PE **and** FGR: last-birth z-score, else 0.
* `prevNEITHERxlastBWz` — a previous pregnancy with neither PE nor FGR:
  last-birth z-score, else 0.

For any parous woman, exactly one of these four is non-zero; **for a nulliparous
woman all four are 0.**

The two spline terms (`prevNEITHERxlastBWz`, `ut.mom.log10`) are supplied as
single raw values — the spline basis is generated automatically from the knots
stored in the model, so no manual transformation is needed. For reference, the
knots are:

* `prevNEITHERxlastBWz`: −0.96466362, 0.07363606, 1.35096910
* `ut.mom.log10`: −0.176652260, 0.005106235, 0.170129250

## Prediction horizons and baseline

The model is fitted on a gestational-age timescale in **days**, and produces the
predicted cumulative incidence of FGR by a chosen horizon. The three horizons
used in the paper:

| Horizon | Time (days) | Baseline CIF₀ (centred LP = 0) | log cumulative baseline subdistribution hazard |
|---------|-------------|-------------------------------|------------------------------------------------|
| < 32⁺⁰ weeks | 223 | 0.14% | −6.57 |
| < 37⁺⁰ weeks | 258 | 0.62% | −5.07 |
| ≤ 40⁺⁰ weeks | 280 | 3.75% | −3.26 |

An individual's predicted risk is

```
CIF(t) = 1 − (1 − CIF₀(t))^exp(LP)
```

where `CIF₀(t)` is the baseline cumulative incidence above and `LP` is the
linear predictor (centred at the development-cohort reference values, so LP = 0
corresponds to the average development profile). `fg_pred_risk()` performs this
calculation for you; the baseline table is provided so the model can also be
applied by hand or reimplemented elsewhere. The 21 coefficients are available
via `coef(model)`.

## Model summary

* **Type:** Fine–Gray competing-risks (subdistribution hazard) model
  (`rms::cph` on a `survival::finegray`-expanded dataset; ridge penalty 0.1),
  with non-FGR delivery as the competing event.
* **Outcome:** a composite definition of FGR (severe smallness, placental
  insufficiency or associated neonatal morbidity), as defined in the paper.
* **Development cohort (FirstPLUS):** n = 3479 (251 FGR events, 3228 competing
  deliveries).
* **External validation cohort (OxPLUS):** n = 2717 (132 FGR events). Aspirin
  was not recorded in this cohort and was set to 0.

### Performance (from the source study)

* **Discrimination (external, ≤ 40⁺⁰ horizon):** AUC 0.75 (95% CI 0.71–0.79).
* **Calibration (external, overall time-integrated slope):** 0.776
  (95% CI 0.63–0.93), i.e. predicted risks are somewhat too extreme in the
  external cohort.

Full development, internal-validation (optimism-corrected) and external results,
including horizon-specific calibration, are reported in the source publication
and its appendix.

### Recalibration

Because external calibration was imperfect, in line with TRIPOD+AI guidance the
model should be recalibrated (intercept and slope) to a new target population
before its absolute risk estimates are used there. Recalibration adjusts the
absolute risks without changing the model's discrimination (which depends only
on the ranking of predictions). The source publication reports a worked
intercept–slope recalibration in the external cohort.

## Data availability and privacy

The distributed model object contains **no individual-level data**: only the
model coefficients and their variance–covariance matrix, the baseline
function, variable metadata, spline knots and summary ranges for each
predictor. The underlying FirstPLUS and OxPLUS datasets are not included and
are not available here; enquiries about data access should be directed to the
study authors.

**File integrity.** MD5 of `FGR_complete_model.rds`:
`257c09708db89948402f26ff9379c596`

## Citation

If you use this model, please cite the source publication

## Contact

Sam Mathewlynn, University of Oxford
sam.mathewlynn@wrh.ox.ac.uk