# new_msna_analysis — `fastmsna`

A performance-optimised replacement for the calculation layer of
[`analysistools`](https://github.com/impact-initiatives/analysistools)
(`create_analysis()` and the `create_analysis_*()` family) for large MSNA datasets.

* **Same results.** Estimates, confidence intervals, `n`, `n_total`, `n_w`,
  `n_w_total`, analysis keys and row order match `analysistools` 0.0.0.903
  (srvyr 1.3.1 / survey 4.5). Validated row by row, see [Validation](#9-validation).
* **Same inputs and outputs.** A srvyr design plus a list of analysis (LOA) in,
  the analysistools long results table out. `Step 3 - merge analysis.R` and
  `Step 4 - formatted analysis.R` work unchanged.
* **Much faster.** **100–180× on the real NGA 2026 household data**
  (MSNA-wide: 11 min → 6.5 s; State: 23 min → 7.8 s) and **100–200× on
  synthetic data**. The whole household sheet (34 disaggregations) runs in
  **3.8 minutes**; the existing outputs took about 29 hours of wall-clock time.
  See [Benchmarks](#10-benchmark-results).

---

## Contents

1. [What the tool does](#1-what-the-tool-does)
2. [Why it is faster](#2-why-it-is-faster)
3. [Project structure](#3-project-structure)
4. [Installation / setup](#4-installation--setup)
5. [Input data requirements](#5-input-data-requirements)
6. [Functions](#6-functions)
7. [Examples](#7-examples)
8. [Methodology](#8-methodology): weighting, variance, CIs, medians, missing values, sample sizes
9. [Validation](#9-validation)
10. [Benchmark results](#10-benchmark-results)
11. [Methodological differences with analysistools](#11-methodological-differences-with-analysistools)
12. [Migrating the existing analysis scripts](#12-migrating-the-existing-analysis-scripts)
13. [Limitations and notes](#13-limitations-and-notes)

---

## 1. What the tool does

`fastmsna` computes, for every row of a list of analysis (LOA) and every group of
its disaggregation (`group_var`):

| `analysis_type` | statistic | CI |
|---|---|---|
| `mean` | weighted mean | Taylor-linearised, *t* distribution |
| `median` | weighted median ("school" rule, as analysistools) | Woodruff interval |
| `prop_select_one` | weighted proportion of each answer | linearised (mean of indicator) |
| `prop_select_multiple` | weighted proportion selecting each choice | linearised |
| `ratio` | `sum(w * numerator) / sum(w * denominator)` | linearised ratio variance |

It supports weights, strata, clusters (first stage), finite population
correction, disaggregation by one or several variables (`"admin1, pop_group"`),
domain (subgroup) estimation, missing values, sample sizes (weighted and
unweighted), and the `survey.lonely.psu` / `survey.adjust.domain.lonely` options.

`analysistools` has no suppression or minimum sample-size rule, so `fastmsna`
does not apply one either. Every output row carries `n_total` / `n_w_total`, so
thresholds can be applied downstream (for example
`filter(results, n_total >= 30)`).

## 2. Why it is faster

### Where the old tool spends its time

Profiling `analysistools::create_analysis()` with the NGA scripts showed that the
time is spent re-building survey designs, not computing statistics:

1. **One full dataset copy per group, per indicator.** `srvyr` evaluates a grouped
   `summarise()` by splitting the design into one survey-design object per
   group (`[.survey.design2`). Each split copies **all columns** of
   `design$variables`. The NGA household sheet has 1,035 columns, so every
   `State` group of every indicator copies about 20 MB.
   - `prop_select_one` groups by `group_var × answer`, so the number of
     copies is groups × answers.
   - `filter(!is.na(x))` copies the whole dataset again for every indicator.
2. **Row-by-row LOA loop.** `purrr::map()` runs one full srvyr pipeline per LOA
   row. Each pipeline repeats `group_by`, `filter`, `summarise`, joins,
   `pivot_longer` / `pivot_wider`, `unite` and key building, and nothing is
   shared between indicators that use the same disaggregation.
3. **Statistics that are computed and then dropped.**
   - `n_w` is computed with `survey_total(vartype = "ci")`, i.e. a full variance
     calculation, and then the CI is discarded.
   - For select-multiple questions, 5 survey statistics per choice are computed
     (including two totals with CIs), and most are discarded.
4. **Per-stratum R loops** in `survey::svyrecvar()` (via `tapply`) for every group
   of every indicator. The median is computed twice: once for the statistic and
   once in a separate totals pipeline joined back on.
5. **`do.call(rbind, ...)`** of hundreds of tibbles at the end, plus
   `run_disaggregated_analysis()` re-doing all of the above for every
   disaggregation.

### What fastmsna does instead

```
LOA ──► check once ──► split by analysis type ──► chunks of indicators (n × k matrices)
                                                        │
design ─► as_fast_design() (once):                      ▼
  weights, stratum codes, PSU codes,     for each disaggregation (group spec, cached):
  nPSU per stratum, fpc                    one kernel call = all groups × all indicators
                                             rowsum() over group / stratum / PSU cells
                                                        │
                                                        ▼
                                  long analysistools table + analysis keys (vectorised)
```

* **Prepare once.** `as_fast_design()` turns the design into a few integer/double
  vectors (weights, stratum and PSU codes, full-design PSU counts, fpc). The
  dataset is never copied; only the columns an indicator needs are read, and
  each converted column is cached across disaggregations.
* **Group once.** Each distinct `group_var` is turned into an integer group id,
  plus group × stratum and group × PSU cells, a single time, and is shared by
  every indicator.
* **All groups and all indicators in one pass.** Indicators of the same type are
  packed into an `n × k` matrix: numeric values, 0/1 answers or 0/1 choices.
  Weighted totals, domain means/ratios, influence values and variance
  components for **all groups and all k indicators** are then computed with
  a few `rowsum()` calls (grouped column sums in C). The variance formula is
  `survey`'s closed form (stratum-centred PSU totals of the linearised
  values), not a per-stratum R loop.
* **Nothing computed twice or thrown away.** No CIs on `n_w`, one pass for the
  median and its totals, `data.table` for assembling the results, and
  analysis keys built with vectorised `paste()`.
* **Bounded memory.** Matrices are processed in chunks of at most
  `getOption("fastmsna.max_cells", 4e6)` cells (~32 MB per matrix).

## 3. Project structure

```
new_msna_analysis/
├── R/                          # package source (fastmsna); file prefix = layer
│   ├── functions-*.R           #   public API: create_analysis_fast(), run_disaggregated_analysis_fast(), ...
│   ├── statistics-*.R          #   design layer, variance engine, estimator kernels, type runners
│   ├── utils-*.R               #   grouping, column preparation/caching, synthetic data
│   └── fastmsna-package.R
├── man/                        # help pages (roxygen2)
├── tests/testthat/             # unit + equivalence tests (vs analysistools and survey)
├── scripts/
│   ├── 00_setup.R              # install / load fastmsna
│   ├── 01_example_usage.R      # how to use
│   ├── Step 2 - main analysis - fast.R      # migrated NGA Step 2
│   └── validate_against_analysistools.R     # full validation report
├── benchmarks/
│   ├── run_benchmarks.R        # synthetic benchmarks (old vs new)
│   ├── benchmark_real_data.R   # NGA real data, old vs new + validation
│   ├── benchmark_real_data_full_sheet.R     # NGA real data, all disaggregations (new)
│   └── bench_functions.R
├── data/                       # local inputs (not part of the package)
├── outputs/
│   ├── validation/             # VALIDATION_REPORT.md, row-by-row comparisons
│   └── benchmarks/             # BENCHMARKS.md, benchmark_results.csv, real-data results
├── logs/
├── renv/  + renv.lock          # recorded package versions (renv not activated, see renv/README.md)
├── DESCRIPTION, NAMESPACE, LICENSE
└── new_msna_analysis.Rproj     # RStudio project (Build pane = package)
```

`msna_analysis_main` keeps its helpers in `src/` and sources them. Here the code
is an **R package** so that it can be installed once and loaded with
`library(fastmsna)` from `msna_analysis_main` or any other project, and so that it
has tests and help pages. R only builds packages from files directly inside
`R/`, so the requested `R/functions`, `R/statistics` and `R/utils` sub-folders
are file-name prefixes instead (`functions-`, `statistics-`, `utils-`).

## 4. Installation / setup

Requirements: R ≥ 4.1 and **data.table**. That is the only hard dependency,
and it is already installed in your environment.

**From GitHub** (any computer):

```r
# install.packages("remotes")   # once, if needed
remotes::install_github("Tomas201513/fastmsna", upgrade = "never")
library(fastmsna)
```

`upgrade = "never"` installs only fastmsna (plus data.table if it is missing)
and leaves your other packages untouched. Re-run the same line to get updates.

**From the local folder:**

```r
source("C:/Users/User/Music/new_msna_analysis/scripts/00_setup.R")
library(fastmsna)
```

Without installing (development): `pkgload::load_all("C:/Users/User/Music/new_msna_analysis")`.
In RStudio: open `new_msna_analysis.Rproj` → *Build* → *Install*.

| package | needed for | status in your library |
|---|---|---|
| data.table | everything (Imports) | installed (1.18.6.1) |
| srvyr, survey | passing srvyr designs (optional: a data frame works too), validation | installed |
| analysistools | validation / benchmarks only | installed (0.0.0.903) |
| writexl | `save_format = "xlsx"` | installed |
| tibble | results returned as tibbles (else data.frame) | installed |
| callr, pkgload, testthat, readxl | benchmarks / tests | installed |

**No additional package is needed.** `collapse` was considered, but base
`rowsum()` plus `data.table` were enough, so nothing new had to be installed.
renv: `msna_analysis_main` does not use renv, so it is not activated here (that
would re-install ~80 packages into a private library). `renv.lock` records the
exact versions used; see `renv/README.md` to switch to renv later.

## 5. Input data requirements

* **Design**: a `srvyr::as_survey_design()` / `survey::svydesign()` object (weights,
  strata, optional cluster `ids`, optional `fpc`). Alternatively, pass a data
  frame plus `weights = "weight", strata = "Strata", ids = NULL, fpc = NULL`
  (column names). Missing weights or strata are an error, as in survey.
  Replicate-weight, calibrated/post-stratified, PPS and multi-stage-with-fpc
  designs are rejected with an explicit error.
* **LOA**: same columns as analysistools: `analysis_type`, `analysis_var`,
  `group_var` (`NA` = no disaggregation; `"a, b"` = a × b), `level` (default
  0.95). Ratios also use `analysis_var_numerator`, `analysis_var_denominator`,
  `numerator_NA_to_0`, `filter_denominator_0` (default `TRUE`).
* **Variables**:
  * mean / median / ratio: numeric columns. Text columns that contain numbers
    (data read with `col_types = "text"`) are converted automatically;
    non-numeric text gives an error for that indicator.
  * select one: any type (character, numeric, logical, factor).
  * select multiple: the parent column (its `NA` = not asked) plus choice columns
    `<parent><sm_separator><choice>` holding 0/1 (numeric, logical, or text
    "0"/"1"/"TRUE"/"FALSE").
* Column names may contain spaces and symbols (`"State - Gender of HoH"`,
  `"Age & Gender HoH"`).

## 6. Functions

| analysistools (old) | fastmsna (new) | notes |
|---|---|---|
| `create_analysis(design, loa, group_var, sm_separator)` | `create_analysis_fast(design, loa, group_var, sm_separator, on_error = "stop", ratio_legacy = FALSE, ...)` | returns `results_table`, `dataset`, `loa` **+ `skipped`** |
| `create_analysis_mean()` | `create_analysis_mean_fast()` | same arguments |
| `create_analysis_median()` | `create_analysis_median_fast(..., qrule = "school")` | same arguments |
| `create_analysis_prop_select_one()` | `create_analysis_prop_select_one_fast()` | same arguments |
| `create_analysis_prop_select_multiple()` | `create_analysis_prop_select_multiple_fast()` | same arguments |
| `create_analysis_ratio()` | `create_analysis_ratio_fast(..., ratio_legacy = FALSE)` | same arguments |
| `create_loa()` / `check_loa()` | `create_loa_fast()` / `check_loa_fast()` | same rules / messages |
| `add_weights()` | `add_weights_fast()` | same arguments and formula |
| `create_analysis_safe()` *(NGA src)* | `create_analysis_fast(..., on_error = "skip")` | `$skipped` has the same columns |
| `run_disaggregated_analysis()` *(NGA src)* | `run_disaggregated_analysis_fast()` | same arguments and file names |
| | `as_fast_design()` | prepare a design once (re-use for many calls) |
| | `compare_analysis_results(old, new)` | row-by-row validation of two results tables |
| | `make_synthetic_msna()` | synthetic MSNA data for tests / benchmarks |
| | `clear_design_cache()` | free the converted-column cache |

The results table has the analysistools columns, in the same order:
`analysis_type, analysis_var, analysis_var_value, group_var, group_var_value, stat,
stat_low, stat_upp, n, n_total, n_w, n_w_total, analysis_key`.
`create_analysis_key_table()`, `create_loa_from_results()`, `review_analysis()` and
all of `presentresults` work on it unchanged; keep using them from analysistools.

## 7. Examples

```r
library(fastmsna)
library(srvyr)
options(survey.lonely.psu = "adjust", survey.adjust.domain.lonely = TRUE)  # as src/init.R

my_design <- as_survey_design(data_main, weights = "weight", strata = "Strata")

# 1. exactly like analysistools
res <- create_analysis_fast(my_design, loa = my_loa_sh)
res$results_table
res$skipped

# 2. many disaggregations of the same data: prepare the design once
fd <- as_fast_design(my_design)
loa_state <- transform(my_loa_sh, group_var = "State")
res_state <- create_analysis_fast(fd, loa = loa_state, on_error = "skip")

# 3. one file per disaggregation (drop-in for run_disaggregated_analysis)
run_disaggregated_analysis_fast(
  sh = "hh data", loa_sheet = my_loa_sh, disagg_vars = disagg,
  design = fd, output_dir = "./outputs/analysis/2scs/results table/hh data/",
  save_format = "xlsx"
)

# 4. a single indicator
create_analysis_ratio_fast(fd, group_var = "State",
                           analysis_var_numerator = "exclusive_breastfeeding",
                           analysis_var_denominator = "under6_month_age", level = 0.9)

# 5. validate against the old tool
old <- analysistools::create_analysis(my_design, loa = loa_state)
compare_analysis_results(old, res_state)
```

More in `scripts/01_example_usage.R`; full migrated workflow in
`scripts/Step 2 - main analysis - fast.R`.

## 8. Methodology

Everything below is what analysistools/srvyr/survey do. fastmsna implements the
same formulas, including their edge cases.

### Weighting
* Weights are the design weights `1 / prob` of the survey design (`weight` in the
  NGA 2SCS data; all 1 for the quota design). `add_weights_fast()` reproduces
  `analysistools::add_weights()`:
  `w_h = (N_h / ΣN) / (n_h / Σn)` (normalised weights, Σw = number of interviews).
* Point estimates: weighted mean `Σ w y / Σ w`, proportion `Σ w 1[y = l] / Σ w`,
  ratio `Σ w y / Σ w x`, over the valid rows of the group.

### Variance and confidence intervals (Taylor linearisation)
* Linearised values: mean / proportion `z = w (y − ŷ) / Σw`; ratio
  `z = w (y − R̂ x) / Σ w x`.
* Domain estimation: for a group (domain) the variance is computed on the
  domain subset but with the **full-design number of PSUs per stratum**
  (`fpc$sampsize`). PSUs of the stratum that are not in the domain count as
  zeros, exactly as survey's `onestrat()` does. Per stratum h:
  `v_h = f_h · n_h/(n_h−1) · Σ_j (z_hj − z̄_h)²`, where z_hj are PSU totals,
  z̄_h is their mean over the n_h PSUs, and `f_h = (N_h − n_h)/N_h` when an
  fpc is given (else 1). `V = Σ_h v_h`.
* Lonely PSUs follow the `survey.lonely.psu` option (`fail`, `remove`,
  `certainty`, `adjust`: centred on the grand mean of PSU totals, `average`)
  and `survey.adjust.domain.lonely`. Options are read at run time, like survey.
  The NGA `init.R` uses `adjust` + `TRUE`.
* CI: `est ± qt(1 − (1 − level)/2, df) · SE` with `df = degf(design)` =
  number of PSUs − number of strata **of the filtered design** (rows with a
  non-missing analysis variable / parent question / valid ratio, all groups
  together), as `srvyr` does.
* Without `ids` (as in the NGA Step 2 design) every household is its own PSU.

### Medians
* Weighted median with the **"school"** rule (analysistools default): with
  `cumw` the cumulative weights in sorted order, the median is the average of
  the two middle values when `cumw` hits exactly 50% of the weight, otherwise
  the next value (`qrule = "math"` is available for the lower value).
* CI: **Woodruff** (`interval.type = "mean"`): p̂ = weighted share of values ≤
  median, CI on p̂ with `df = degf(group design)`, transformed back with the same
  quantile rule; a bound outside [0, 1] gives `NaN`, as in survey. survey's exact
  floating-point arithmetic is reproduced for p̂, because the "> 1" test is
  decided on the last bits when the median is the group's largest value.

### Missing values
| type | rows used | `n` | `n_total` | `n_w` / `n_w_total` |
|---|---|---|---|---|
| mean / median | analysis variable not NA | valid rows in group | = n | Σw valid / = n_w |
| prop_select_one | analysis variable not NA | rows with this answer | valid rows in group | Σw with this answer / Σw valid |
| prop_select_one, row `NA` | rows with missing answer | number of NA | NaN | NaN |
| prop_select_multiple | parent not NA and choice not NA | Σ choice (number selecting) | valid rows | Σ w·choice / Σw valid |
| prop_select_multiple, row `"NA"` | rows with missing parent | number of NA (row only if > 0) | NA | NA |
| ratio | denominator not NA (and ≠ 0 if `filter_denominator_0`), numerator NA → 0 if `numerator_NA_to_0`, else removed | rows used | = n | Σw used / = n_w |

* Groups where nothing is valid are kept (`n = 0`). As in analysistools, `stat`,
  `stat_low`, `stat_upp`, `n_total`, `n_w` and `n_w_total` are `NaN` when
  `n_total = 0`.
* Missing values of the **grouping** variable form their own group (`"NA"`), as in
  dplyr. `Step 4` filters them out.

### Ordering and keys
Groups are ordered like `dplyr::group_by()` (each variable ascending, C locale,
`NA` last); answers likewise; select-multiple choices in dataset column order.
`analysis_key` follows the analysistools format
(`type @/@ var %/% value @/@ group %/% value -/- ...`).

## 9. Validation

Validated with `scripts/validate_against_analysistools.R`. The report is in
`outputs/validation/VALIDATION_REPORT.md`, with row-by-row CSVs in
`outputs/validation/`. Rows are matched on `analysis_key`. A value matches when
both sides are missing or `|old − new| ≤ 1e-8 · max(1, |old|)`.

**analysistools test data and synthetic MSNA data** (1,000 households with
edge cases: a single-PSU stratum, missing group values, a group where a
variable is entirely missing, all-missing questions, single-observation
groups and heavily tied values; 171 LOA rows over 5 disaggregations including
`"admin1, pop_group"`; levels 0.9 and 0.95):

| scenario | design | result rows | value mismatches | same row order | failed rows old / new | max rel. diff |
|---|---|---|---|---|---|---|
| template_loa | analysistools template data + LOA | 143 | 0 | yes | 0 / 0 | 2e-14 |
| template_loa_with_ratio | + ratios | 147 | 0 | yes | 0 / 0 | 2e-14 |
| template_no_loa | automatic LOA, `group_var = c("admin1", "admin1, admin2")` | 399 | 0 | yes | 0 / 0 | 2e-14 |
| A_2scs_strata_weights | weights + strata (= NGA Step 2), adjust + domain lonely | 2,830 | 0 | yes | 1 / 1 | 4e-13 |
| B_clustered | + cluster ids | 2,830 | 0 | yes | 1 / 1 | 1e-13 |
| C_quota_unweighted | no weights / strata (= NGA quota) | 2,830 | 0 | yes | 1 / 1 | 5e-14 |
| D_lonely_fail | survey defaults (`lonely.psu = "fail"`) | 93 | 0 | yes | 156 / 156 | 0 |
| E_lonely_remove / certainty / average | clustered, each lonely option | 2,830 each | 0 | yes | 1 / 1 | ≤ 6e-13 |
| E_adjust_no_domain | `adjust.domain.lonely = FALSE` | 2,830 | 0 | yes | 1 / 1 | 1e-13 |
| F_fpc | + finite population correction | 2,830 | 0 | yes | 1 / 1 | 1e-13 |
| G_ratio_documented_filters | A with the ratio fix on | 2,802 | 0 | yes | 1 / 1 | 4e-13 |

The single failing row in A–F is a select-one grouped by itself, which both
tools reject. In G, all rows are identical except the 8 ratio rows with
`numerator_NA_to_0 = FALSE`, which differ by design (84 values, listed in
`outputs/validation/ratio_legacy_differences.csv`; see
[section 11](#11-methodological-differences-with-analysistools)). Scenarios A–F
run with `ratio_legacy = TRUE`, i.e. they reproduce the old ratio behaviour
exactly.

**Real data** (`benchmarks/benchmark_real_data.R`, NGA 2026 R1 2SCS household
sheet, 20,858 households × 1,035 columns, DAP `final` sheet, 253 indicators,
design `weights = weight, strata = Strata`):

| disaggregation | result rows | value mismatches | same row order | max rel. diff | old | new | speed-up | peak memory old → new |
|---|---|---|---|---|---|---|---|---|
| MSNA-wide | 2,251 | 0 | yes | 3e-13 | 654 s | 6.5 s | 101× | 1,009 → 670 MB |
| State | 12,441 | 0 | yes | 6e-12 | 1,394 s | 7.8 s | 179× | 1,031 → 707 MB |

Both tools produce identical results on the real data. The nutrition sheet
(9,096 children, 7 disaggregations, including the DAP ratio
`exclusive_breastfeeding / under6_month_age`) runs in 3.1 s; its ratio matches
`analysistools::create_analysis_ratio()`.

Unit tests (`devtools::test()`, 16 test blocks, all passing) also check:
* the variance engine against `survey::svymean`/`svyratio` directly, for all
  lonely-PSU options × domain-lonely settings × (rows / clusters / fpc) designs;
* the median and its Woodruff CI against `survey::svyquantile` on heavily tied
  data (a stress run with 3,600 groups had 0 differences);
* `create_loa_fast` vs `create_loa`, `add_weights_fast` vs `add_weights`, data
  frame vs srvyr input, file names of `run_disaggregated_analysis_fast`.

Differences are floating-point noise only: the largest relative difference is
below 1e-12. The only intentional differences are listed in [section 11](#11-methodological-differences-with-analysistools).

## 10. Benchmark results

Setup: Windows 11, R 4.5.1, 8 cores; analysistools 0.0.0.903 vs fastmsna 0.1.0.
Each case runs in a fresh R session. Memory = peak R heap during the analysis.
"Row-indicators / s" = dataset rows × LOA rows / seconds. Results were compared
for every case and are identical in all of them. Full tables:
`outputs/benchmarks/BENCHMARKS.md`, `benchmark_results.csv`,
`real_data_*.csv`.

### Real NGA data (household sheet: 20,858 rows × 1,035 columns, 253 indicators)

| run | old (analysistools) | new (fastmsna) | speed-up |
|---|---|---|---|
| MSNA-wide (1 group) | 654 s, 1,009 MB | 6.5 s, 670 MB | **101×** |
| State (7 groups) | 1,394 s, 1,031 MB | 7.8 s, 707 MB | **179×** |
| all 34 disaggregations, incl. saving the files | ≈ 29 h of wall-clock time (timestamps of the existing outputs, 4 Oct 04:23 → 5 Oct 09:27, several sessions) | **3.8 min** (one session, peak 1.6 GB) | — |

The largest disaggregation, `Strata - Age & Gender HoH` (937,689 result rows),
takes 17 s; `State - Setting` takes 4.3 s (3 h 42 min in the existing run).
Per-disaggregation times are in `outputs/benchmarks/real_data_full_sheet_hh_data.csv`.

### Synthetic MSNA data (364 columns, 31 indicators × 5 disaggregations = 154 LOA rows)

| rows | old | new | speed-up | peak memory old → new | row-indicators / s old → new |
|---|---|---|---|---|---|
| 1,000 | 158 s | 0.92 s | **172×** | 326 → 93 MB | 976 → 167,391 |
| 5,000 | 152 s | 1.05 s | **145×** | 253 → 92 MB | 5,052 → 733,333 |
| 20,000 | 212 s | 1.97 s | **108×** | 366 → 179 MB | 14,545 → 1,563,452 |
| 100,000 | 645 s | 4.84 s | **133×** | 1,177 → 531 MB | 23,881 → 3,181,818 |

### Which functions gain the most (100,000 rows, disaggregation `admin1, pop_group`, 18 groups)

| analysis type | LOA rows | old | new | speed-up |
|---|---|---|---|---|
| prop_select_multiple | 4 | 38.4 s | 0.22 s | **174×** |
| mean | 6 | 17.8 s | 0.11 s | **162×** |
| prop_select_one | 13 | 121.9 s | 0.81 s | **151×** |
| median | 6 | 28.0 s | 0.40 s | **70×** |
| ratio | 2 | 7.0 s | 0.10 s | **70×** |

* The largest gains are for **select-one and select-multiple** questions.
  - These are 80% of the NGA DAP (159 + 47 of 253 hh indicators).
  - In the old tool they cost the most: one design copy per group × answer,
    plus 5 survey statistics per choice for select-multiple.
* **Medians** gain less: sorting per group stays an R-level step, and the old
  median already skipped some of the overhead.
* **The speed-up grows with dataset width.** The old cost is proportional to the
  number of columns copied per group, so the real 1,035-column sheet gains more
  (101–179×) than the 364-column synthetic data.
* **The new tool scales roughly linearly**: about 1.6 → 3.2 million
  row-indicators per second from 20,000 to 100,000 rows. The old tool processes
  1,000–24,000.
* At n = 20,000 the per-function numbers vary between runs (0.06–0.3 s for the
  new tool, because the cases are very short). See `BENCHMARKS.md` for the full table.

Reproduce with `Rscript benchmarks/run_benchmarks.R` (≈ 45 min, mostly the old
tool), `benchmarks/benchmark_real_data.R` and `benchmarks/benchmark_real_data_full_sheet.R`.

## 11. Methodological differences with analysistools

By default, results are identical to analysistools, with the exceptions below.
Each one fixes a clear problem in the old tool, and where it changes numbers it
can be switched back.

1. **Ratio filters when `numerator_NA_to_0 = FALSE` (bug fix).** In
   `create_analysis_ratio()`, the numerator filter overwrites the filtered
   design, so the documented denominator filters (drop missing / zero
   denominators) are silently lost. Zero denominators then stay in the ratio and
   in `n`/`n_w`, and rows with a missing denominator stay in `n`/`n_w`.
   fastmsna applies the documented filters; `ratio_legacy = TRUE` reproduces
   the old numbers. Validation scenario G: all other rows are identical, and
   only the 8 affected ratio rows differ (84 values).
2. **Missing `numerator_NA_to_0` / `filter_denominator_0` on a ratio row** (as in the
   NGA DAP: exclusive breastfeeding): analysistools stops with
   "missing value where TRUE/FALSE needed". fastmsna uses the default `TRUE`
   with a warning.
3. **Ratio rows in `run_disaggregated_analysis()` (NGA helper).**
   `filter(analysis_var != group_var)` drops every row with
   `analysis_var = NA`, so the DAP's ratio (nut data, exclusive breastfeeding)
   was never computed by Step 2. `run_disaggregated_analysis_fast()` keeps
   these rows; `keep_ratio_rows = FALSE` restores the old behaviour.
4. **Select-multiple choices named `*_low` / `*_upp` (bug fix).** analysistools
   parses its own `stat_low` / `stat_upp` column names, so a choice such as
   `very_low` is split into fake rows (`very`, `very_low`) with `NA` statistics.
   fastmsna reports it correctly. This does not occur in the NGA household data.
5. **Text `"TRUE"`/`"FALSE"` select-multiple dummies** (logical columns written to
   Excel and read back as text) become `NA` with `as.numeric()` in analysistools,
   which silently gives `NaN` proportions. fastmsna maps them to 1/0. 0/1
   dummies, as in the NGA data, are unaffected.
6. **Text columns for mean / median / ratio** are converted to numbers when they
   contain numbers (srvyr would stop), and are an explicit per-indicator error
   otherwise.
7. **Errors.** `create_analysis()` stops the whole run at the first failing
   indicator. With `on_error = "skip"` (default in
   `run_disaggregated_analysis_fast()`), the indicator is logged in `$skipped`
   / `_logs/*.csv` and the rest is computed. Validation scenario D
   (`lonely.psu = "fail"`): the same 156 rows fail in both tools.
8. Output column types: `analysis_var_value` is always character and the `n*`
   columns are always double (analysistools mixes integer and double, depending
   on the type). For a *logical* select-one variable, analysistools'
   final `rbind()` can turn `analysis_var_value` into `"0"`/`"1"` while its own
   `analysis_key` says `FALSE`/`TRUE`. fastmsna writes `"FALSE"`/`"TRUE"` in
   both. Keys and numbers are identical. This does not apply to text data such
   as the NGA datasets.

Not reproduced (unsupported designs, error raised): replicate-weight,
calibrated/post-stratified, PPS and multi-stage-with-fpc designs.

**Observation, not changed:** the NGA 2SCS design in Step 2
(`as_survey_design(weights = "weight", strata = "Strata")`) has no cluster ids,
so CIs treat households as independent within strata and are likely too narrow
for a two-stage cluster sample. If the sampling team agrees, add
`ids = <cluster id column>`; fastmsna supports it (validated in scenario B).

## 12. Migrating the existing analysis scripts

`Step 1`, `Step 3`, `Step 4` and `Step 5` do not change. In `Step 2`:

1. Install fastmsna once: `remotes::install_github("Tomas201513/fastmsna", upgrade = "never")`
   (or, from the local folder, `source("C:/Users/User/Music/new_msna_analysis/scripts/00_setup.R")`).
2. After `source("src/init.R")`, add `library(fastmsna)`. Keep analysistools loaded
   for its other helpers; the function names do not clash.
3. After building `my_design`, add `my_fast_design <- as_fast_design(my_design)`.
4. Replace `run_disaggregated_analysis(` with `run_disaggregated_analysis_fast(`
   and `design = my_design` with `design = my_fast_design`. All other arguments
   are the same.
5. Anywhere else, rename:
   * `create_analysis(` → `create_analysis_fast(`
   * `create_analysis_safe(` → `create_analysis_fast(..., on_error = "skip")`
   * `create_analysis_mean(` etc. → `create_analysis_mean_fast(` etc.
6. Optional: since a full sheet now takes minutes, loop over all sheets and all
   disaggregations in one run.

A ready-to-use version is `scripts/Step 2 - main analysis - fast.R`. Copy it into
`msna_analysis_main` next to the other Step scripts and run it from there. It
writes the same files to the same folders, so Step 3 merges them as before.
Skipped indicators are logged in `<output_dir>/_logs/` (a sub-folder, so Step 3
ignores it).

To check a migration on your own data, run both tools on one disaggregation and
compare:

```r
old <- analysistools::create_analysis(my_design, loa = loa_sub)
new <- create_analysis_fast(my_design, loa = loa_sub)
compare_analysis_results(old, new)   # => EQUIVALENT
```

## 13. Limitations and notes

* Memory: matrices are chunked; `options(fastmsna.max_cells = 2e6)` lowers peak
  memory (slightly slower), and `clear_design_cache(fd)` frees the cache of
  converted columns.
* Medians are computed one variable at a time (sorting); everything else is
  batched.
* `level` is per LOA row. Several levels for the same indicator are supported.
* Run the tests with `devtools::test()` (≈ 1–2 minutes, mostly analysistools).
