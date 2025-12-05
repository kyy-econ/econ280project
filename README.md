# Replication and Extension of Lee (2017): Extensive and Intensive Margin Effects in Sample Selection Models

This repository contains Stata code and data to replicate and extend  
Myoung-jae Lee (2017), *“Extensive and intensive margin effects in sample selection models: Racial effects on wages.”*  

This repository exactly replicates the linear (independence), normality (Heckit), and quadratic control-function specifications reported in Lee (2017), and then implements an additional cubic control-function model that is not estimated in the original paper. Lee (2017) notes that higher-order polynomial approximations of the selection index are theoretically admissible, but the published article reports results only up to the quadratic case. Building on that observation, the ado file here augments the outcome equation with the corresponding third-order term in the selection index and computes the associated intensive margin effect (IME), extensive margin effect (EME), and total effect (TE). The rows labeled IME_cub, EME_cub, and TE_cub in the output table are therefore new estimates produced by this replication/extension, using the original data and covariates, and do not appear in Lee (2017).

The main components are:

- Stata ado-files implementing the decomposition (`MJleeIMEEME_BASE2.ado`, `MJleeIMEEMEs2.ado`)
- Example dataset (`MJLEEdata.dta`)
- A driver do-file (e.g. `replication.do`) that reproduces the treatment-effect table with bootstrap confidence intervals

---

## 1. Overview

The repository is organized around a single example using the wage data from Lee (2017). The workflow is:

1. Load the example dataset `MJLEEdata.dta`.
2. Run `MJleeIMEEMEs2`, which:
   - Fits the probit, Tobit, independence, normality (Heckit), quadratic, and cubic models
   - Computes the intensive margin effect (IME), extensive margin effect (EME), and total effect (TE) for each model
   - Performs a nonparametric bootstrap
   - Returns a matrix of estimates, standard errors, and percentile 95% confidence intervals in `r(results)`

Expected runtime on a standard laptop (Stata 14–18) with 200 bootstrap replications is well under 10 minutes.

---

## 2. Data Availability and Provenance

### Data used

The example uses a single dataset:

| Data file       | Description                                                 | Format | Provided | Notes                                |
|-----------------|-------------------------------------------------------------|--------|----------|--------------------------------------|
| `MJLEEdata.dta` | Individual-level wage sample used in Lee (2017)             | `.dta` | Yes*     | Original replication / example data  |


---

## 3. Software Requirements

The replication and extension are entirely in Stata.

- **Stata**: tested with Stata 14; should work with Stata 15+ as well

No other external software or languages are required.

Approximate compute requirements:

- OS: any platform supported by Stata (Windows, macOS, Linux)
- CPU/RAM: standard modern laptop/desktop is sufficient
- Runtime: `< 10` minutes for 200 bootstrap replications on the example data

Randomness is controlled by a fixed seed passed via the `seed()` option of `MJleeIMEEMEs2` (e.g., `seed(12345)`), ensuring reproducible bootstrap results.

---

## 4. Code Description and File Structure

A typical layout for this repository is:

- `MJleeIMEEME_BASE2.ado`  
  Core decomposition routine. Given an outcome, treatment, and selection indicator, plus covariate lists, it:
  - Estimates the probit selection equation
  - Constructs the necessary control-function terms
  - Estimates outcome equations under:
    - Tobit
    - Independence (OLS on selected sample)
    - Normality (Heckit)
    - Quadratic control function
    - Cubic control function
  - Computes IME, EME, and TE for each model
  - Returns these as scalars in `r()` and as a matrix `r(vif)`

- `MJleeIMEEMEs2.ado`  
  Wrapper program that:
  - Runs and displays the baseline regressions (probit, Tobit, independence, normality, quadratic, cubic)
  - Calls `MJleeIMEEME_BASE2` once on the full sample
  - Performs a nonparametric bootstrap with a user-specified number of replications
  - Aggregates the bootstrap draws into:
    - Point estimates (`full_est`)
    - Bootstrap standard errors
    - Percentile 95% confidence intervals
  - Returns the final `results` matrix in `r(results)` and prints it in a readable format

- `MJLEEdata.dta`  
  Example dataset used to replicate the decomposition for Mexican-American vs white male workers, as in Lee (2017).

- `replication.do` 
  Script showing how to run the ado programs, generate replications and extentions, and format the output. 


---

## 5. Instructions to Replicators

To reproduce the example decomposition and bootstrap table:

1. **Install the ado-files**

   Place `MJleeIMEEME_BASE2.ado` and `MJleeIMEEMEs2.ado` either:
   - In Stata’s `PERSONAL` ado directory (see `sysdir` in Stata), or
   - In the working directory from which you will run the replication

2. **Obtain the data**

   - Ensure `MJLEEdata.dta` is in the working directory  
   - If you cannot redistribute the data, obtain it from the original replication archive or author and save it under that filename

3. **Run the driver script**

   From within Stata:

       cd "path/to/this/repo"
       do replication.do

   This will:
   - Run all models
   - Perform the bootstrap
   - Print a summary table in the Stata Results window
   - Return the matrix `r(results)`


---

## 6. Output: Mapping to Quantities in Lee (2017)

The final `results` matrix has 13 rows and 4 columns.

**Rows (models × effects):**

1. `Tobit`   – Tobit benchmark effect  
2. `IME_ind` – IME under independence model  
3. `EME_ind` – EME under independence model  
4. `TE_ind`  – Total (IME + EME) under independence  
5. `IME_nor` – IME under normality (Heckit)  
6. `EME_nor` – EME under normality  
7. `TE_nor`  – Total under normality  
8. `IME_qua` – IME under quadratic control function  
9. `EME_qua` – EME under quadratic control function  
10. `TE_qua` – Total under quadratic control function  
11. `IME_cub` – IME under cubic control function  
12. `EME_cub` – EME under cubic control function  
13. `TE_cub` – Total under cubic control function  

**Columns:**

1. Estimate  
2. Bootstrap standard error  
3. 2.5th percentile (lower 95% CI)  
4. 97.5th percentile (upper 95% CI)

These correspond to the decomposed racial wage gap components discussed in Lee (2017), with the cubic specification providing an additional robustness check suggested by the author.

## Extensive and Intensive Margin Effects (Lee, 2017)

| Model      | Estimate  | Std. Err. | 95% CI lower | 95% CI upper |
|-----------|-----------|-----------|--------------|--------------|
| Tobit     | 0.1435579 | 0.0133751 | 0.1195571    | 0.1675532    |
| `IME_ind` | -0.0626832| 0.0049002 | -0.0715186   | -0.0518146   |
| `EME_ind` | 0.2094379 | 0.0167560 | 0.1792414    | 0.2401934    |
| `TE_ind`  | 0.1467547 | 0.0171911 | 0.1167288    | 0.1809648    |
| `IME_nor` | -0.1487247| 0.0098656 | -0.1661907   | -0.1297252   |
| `EME_nor` | 0.3045402 | 0.0245510 | 0.2587204    | 0.3485208    |
| `TE_nor`  | 0.1558155 | 0.0181942 | 0.1238322    | 0.1922616    |
| `IME_qua` | -0.1572057| 0.0109561 | -0.1774103   | -0.1345250   |
| `EME_qua` | 0.3071671 | 0.0246718 | 0.2600560    | 0.3525467    |
| `TE_qua`  | 0.1499613 | 0.0175364 | 0.1193335    | 0.1847584    |
| `IME_cub` | -0.1577669| 0.0110575 | -0.1782499   | -0.1355234   |
| `EME_cub` | 0.3070517 | 0.0246589 | 0.2600691    | 0.3525870    |
| `TE_cub`  | 0.1492848 | 0.0175237 | 0.1188180    | 0.1839596    |


---

## 7. Reference

Myoung-jae Lee (2017).  
“Extensive and intensive margin effects in sample selection models: Racial effects on wages.”  
*Journal of the Royal Statistical Society. Series A (Statistics in Society)*, Vol. 180, No. 3, pp. 817-839.
