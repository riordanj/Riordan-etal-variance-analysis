# Mean and variance effects across MMPC two-cohort comparisons

Analysis code and inputs for Riordan et al. Running `variance_analysis.R`
reproduces the per-trait statistics reported in the paper from the released
MMPC per-animal data.

## Contents

| File | Description |
|---|---|
| `variance_analysis.R` | The analysis. Reads `trait_index.csv`, selects the two cohorts for each trait from the files in `data/`, and writes `results.csv`. |
| `trait_index.csv` | One row per trait (369). Gives the source file, the assay column, the definition of each cohort, any factor held fixed, and any excluded animal. This is what the script runs from. |
| `cohort_definitions.csv` | The same information organised by comparison (41 rows) and in prose, for reading rather than execution. Includes the reason for each exclusion and notes on comparisons that deviate from the general rule. |
| `data/` | The 28 MMPC per-animal files the analysis reads, unmodified except as noted below. |
| `results.csv` | Produced by running the script; see Running below. |
| `sessionInfo.txt` | R version and platform used to produce `results.csv`. |

## Requirements

R >= 3.6.0. **No packages are required** — the script uses only base R.

The minimum version is set by `sample()`, whose algorithm changed in R 3.6.0.
The script calls `RNGkind(sample.kind = "Rejection")` and `set.seed()` so that
the permutation results are identical on any R at or above that version.

## Running

    Rscript variance_analysis.R [index] [data_dir] [output] [nperm] [seed]

All arguments are optional:

    Rscript variance_analysis.R                      # defaults, reproduces results.csv
    Rscript variance_analysis.R trait_index.csv data out.csv 9999 none
    Rscript variance_analysis.R my_index.csv  my_data out.csv

Defaults are `trait_index.csv`, `data`, `results.csv`, 9999 permutations and
seed 20261023. Runtime is roughly 20 minutes for 369 traits at 9999
permutations. The script reports the permutation count and seed it used.

### On the permutation p-values

The six permutation p-values are Monte Carlo estimates of exact permutation
tests, so they are not fixed quantities. Re-running with a different seed
changes them slightly: near p = 0.05 the standard error is about 0.002 at
9999 permutations, and about 0.0007 at 99999. A result that moves from 0.048
to 0.052 on reseeding has not failed to reproduce; it is a trait sitting close
enough to the threshold that the Monte Carlo draw decides which side it falls
on, and that is worth knowing about such a trait.

The default seed is supplied so that `results.csv` can be checked line for
line. Pass `none` as the seed to draw a fresh one each run, which is the more
informative choice when judging how stable a particular result is, and the
natural one when applying the script to new data.

## Applying this to other cohorts

The analysis is parameterised over `trait_index.csv` rather than written
against this dataset, so running it on a different study means writing a new
index rather than editing the script. Each row needs a source file, the column
holding the measurements, and a definition of the two cohorts in terms of
columns present in that file:

    source_file,assay_column,cohort_A,cohort_B,held_fixed,excluded_animals
    mystudy.csv,glucose (mg/dL),Genotype=WT,Genotype=KO,Diet=chow,

`cohort_A` and `cohort_B` select the two groups; `held_fixed` constrains both
to a shared level; `excluded_animals` lists values of an `Animal` column to
drop. Terms within a field are separated by ` | `. Any file readable as
delimited text with one row per animal will work.

## Method

For each trait, animals are divided into cohorts A and B as specified in
`trait_index.csv`, and six effect metrics are computed:

| Metric | Definition |
|---|---|
| `meanDif` | absolute difference in cohort means |
| `meanDif2` | that difference divided by the sum of the means |
| `meanRat` | larger mean / smaller mean |
| `madDif` | absolute difference in cohort MADs |
| `madDif2` | that difference divided by the sum of the MADs |
| `madRat` | larger MAD / smaller MAD |

MAD is the median absolute deviation, unscaled (`mad(x, constant = 1)`).
`meanRat` and `madRat` take the larger of the two possible ratios, so both are
at least 1 and neither records which cohort is the more variable one.

Each metric is tested by permuting the cohort labels 9999 times and comparing
the observed value with the permuted distribution. The reported p-value is
`(number of permuted values >= observed + 1) / 10000`, one-sided, which suits
these statistics because each is a non-negative magnitude.

Three further tests are reported per trait:

- `p_welch_t` — Welch's t-test on the measured values.
- `p_mad_t` — Welch's t-test on each animal's absolute deviation from its own
  cohort median.
- `p_brown_forsythe` — the pooled-variance form of the same test on those
  deviations.

## Notes on the data

The files in `data/` are as released by MMPC, with one exception: the study
8104 traits are recorded in a combined 8103-8104 workbook rather than a
per-experiment file, and are supplied here as `MMPC_Experiment_8104.csv`, an
unmodified export of that worksheet.

Three points that are not recoverable from the released files by inspection
and are therefore specified in the index:

- Study 7714's hormone assays are recorded in the study 7725 file, and are
  compared within the LabDiet arm only. Its other traits come from its own
  file with no diet constraint.
- One animal is excluded from the control cohort of study 8003's glucose
  infusion rate, which records 41,900 mg/kg/min against a cohort range of
  4 to 52. The study's other traits use the full cohorts.
- One control animal is excluded from both traits of study 9207. The reason is
  not recorded in the source data; see `cohort_definitions.csv`.

The MMPC files are latin1-encoded and carry a few stray bytes before the first
column name, so the script reads them as raw bytes and converts explicitly.
Within the cohort specifications, `column=value` terms are separated by ` | `
rather than `;`, because one strain name (`B6;129S4-Apoa4<tm1Bres>/J`)
contains a semicolon.

## Data source and attribution

The files in `data/` are per-animal data from publicly released studies in the
Mouse Metabolic Phenotyping Centers database (MMPC-Live, https://www.mmpc.org),
redistributed here under the MMPC-Live terms of use, which permit reuse for any
purpose without restriction subject to scholarly citation norms.

If you use these data, acknowledge MMPC-Live with the text they specify:

> The MMPC-Live data used in this manuscript was supported by the NIDDK National
> Mouse Metabolic Phenotyping Centers (MMPC-Live, RRID:SCR_008997,
> [www.mmpc.org](https://www.mmpc.org)).

The database itself is described in Laughlin MR, Lloyd KCK, Cline GW, Wasserman
DH. NIH Mouse Metabolic Phenotyping Centers: the power of centralized
phenotyping. *Mamm Genome*. 2012;23:623-631.

MMPC-Live's own disclaimer applies to the data in `data/`:

> THE MMPC-Live MAKES NO REPRESENTATION ABOUT THE SUITABILITY OR ACCURACY OF THE
> SOFTWARE OR DATA FOR ANY PURPOSE, AND MAKES NO WARRANTIES, EITHER EXPRESS OR
> IMPLIED, INCLUDING MERCHANTABILITY AND FITNESS FOR A PARTICULAR PURPOSE OR
> THAT THE USE OF THE SOFTWARE OR DATA WILL NOT INFRINGE ANY THIRD PARTY
> PATENTS, COPYRIGHTS, TRADEMARKS, OR OTHER RIGHTS. THE SOFTWARE AND DATA ARE
> PROVIDED "AS IS".

MMPC-Live note that their datasets are not necessarily static and that they may
post corrections and updates. The copies in `data/` are the versions the
published analysis was run on, and are included so that the analysis reproduces
from this repository alone.

## License

The code in this repository — `variance_analysis.R` and the index files that
drive it — is released under the MIT License; see `LICENSE`.

The license covers the code only. The contents of `data/` are MMPC-Live's and are
redistributed under the MMPC-Live terms of use described above, not relicensed.
