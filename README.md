# cstats

**Simple statistical analysis and data visualization in R for biomedical researchers.**

`cstats` is a lightweight R toolkit designed to make statistical analysis and data visualization accessible to biomedical researchers with minimal programming experience.

It provides a straightforward workflow for importing experimental data, defining factors and variables, aggregating replicates, performing statistical tests, and generating plots with statistical significance annotations.

## Features

- **Data import:** Read tabular datasets and handle comma decimal separators.
- **Experimental design:** Define experimental factors, quantitative variables, and replicate identifiers.
- **Data aggregation:** Aggregate technical replicates using group means.
- **Statistical analysis:** Two-way ANOVA, Tukey post hoc comparisons, Welch's t-tests, and Wilcoxon rank-sum tests.
- **Data visualization:** Generate bar plots with individual observations, mean ± SD, and statistical significance annotations.

## Repository structure

```text
cstats/
├── README.md
├── cstats.R
└── training/
    ├── tutorial.xlsx
    └── tutorial.R
```

**Files:**

- `cstats.R` — Main statistical analysis and visualization toolkit.
- `training/tutorial.xlsx` — Synthetic experimental dataset for training, including the experimental design and example data.
- `training/tutorial.R` — Step-by-step tutorial demonstrating data manipulation, statistical analysis, and visualization using `cstats`.

## Requirements

Install [R](https://www.r-project.org/) and optionally [RStudio](https://posit.co/download/rstudio-desktop/).

Install the required R packages:

```r
install.packages(c(
  "dplyr", "tidyr", "ggplot2", "factoextra",
  "patchwork", "rstatix", "ggpubr", "emmeans"
))
```

## Getting started

Clone or download this repository.

Open RStudio and set the working directory to the repository root, where `cstats.R` is located.

Load the toolkit:

```r
source("cstats.R")
```

The toolkit functions are now available in your R session.

### Import experimental data

`cstats` reads tabular text files (tab-separated values).

For example, after exporting the `R_DATA` worksheet from `training/tutorial.xlsx` as a CSV file:

```r
data <- cstats.read("training/tutorial.txt")
```

The CSV file must first be exported from Excel; it is not included separately in the repository.

### Define experimental factors and variables

```r
data <- cstats.set_factors(data, c("Time", "Treatment"))

data <- cstats.set_vars(data, c("Viability", "Marker", "Proliferation"))

data <- cstats.set_reps(data, "BioRep")
```

### Aggregate technical replicates

```r
data_agg <- cstats.aggregate(data)
```

This computes mean values within the groups defined by experimental factors and replicate identifiers.

### Perform statistical analysis

Two-way ANOVA:

```r
anova_results <- cstats.2wanova(data_agg)
```

Tukey post hoc comparisons:

```r
tukey_results <- cstats.tukey(data_agg, "Viability")
```

Welch's t-tests:

```r
t_results <- cstats.t(data_agg, "Viability", base = "48")
```

Wilcoxon rank-sum tests:

```r
wilcox_results <- cstats.wilcox(data_agg, "Viability", base = "48")
```

### Generate plots

```r
plots <- cstats.barplot.tukey(data_agg, "Viability")

plots[[1]]
```

Additional plotting functions:

```r
cstats.barplot.t(data_agg, "Viability")
cstats.barplot.wilcox(data_agg, "Viability")
```

These functions return lists of `ggplot2` objects that can be further customized.

## Training tutorial

A practical training example is available in the `training/` directory.

The tutorial demonstrates a complete workflow:

1. Preparing experimental data in Excel.
2. Importing and manipulating data in R.
3. Loading `cstats`.
4. Defining experimental design.
5. Aggregating technical replicates.
6. Performing statistical tests.
7. Generating plots with statistical annotations.

**Training files:**

- [Synthetic experimental dataset](training/tutorial.xlsx)
- [R training tutorial](training/tutorial.R)

The synthetic dataset represents a cell culture experiment with three experimental conditions (Control, Drug0, Drug1), nine time points, biological and technical replicates, and three quantitative measurements.

The training data are entirely synthetic and intended for educational purposes.

## Statistical considerations

Users are responsible for selecting appropriate statistical tests and verifying their assumptions.

Important considerations:

- Technical replicates should not be treated as independent biological observations.
- The current t-test and Wilcoxon implementations assume independent samples.
- Pairwise t-tests and Wilcoxon tests do not currently adjust for multiple comparisons.
- Two-way ANOVA assumes an appropriate experimental design and model structure.
- Aggregation must preserve the correct experimental unit.

## License

No license has been specified yet.
