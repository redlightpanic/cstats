# cstats

**Simple statistical analysis and data visualization in R for biomedical researchers.**

`cstats` is a lightweight R toolkit designed to make statistical analysis and data visualization accessible to researchers with minimal programming experience.

It provides a straightforward workflow for importing experimental data, defining variables and factors, aggregating replicates, performing statistical tests, and generating publication-oriented plots.

## Features

- **Data import:** Read tabular datasets and handle comma decimal separators.
- **Experimental design:** Define factors, quantitative variables, and replicate identifiers.
- **Data aggregation:** Aggregate technical replicates using group means.
- **Statistical analysis:** Two-way ANOVA, Tukey post hoc comparisons, Welch's t-tests, and Wilcoxon rank-sum tests.
- **Visualization:** Generate bar plots displaying individual observations, mean ± SD, and statistical significance annotations.

## Requirements

Install R and RStudio, then install the required packages:

```r
install.packages(c(
  "dplyr", "tidyr", "ggplot2", "factoextra",
  "patchwork", "rstatix", "ggpubr", "emmeans"
))
```

## Getting started

Download `cstats.R` and the example dataset from this repository.

Load the toolkit:

```r
source("cstats.R")
```

Import the data:

```r
df <- cstats.read("data.csv", separator = ",")
```

Define experimental factors and variables:

```r
df <- cstats.set_factors(df, c("STAGE", "Treatment"))
df <- cstats.set_vars(df, c("Viability", "Marker", "Proliferation"))
df <- cstats.set_reps(df, "BioRep")
```

Aggregate technical replicates:

```r
df_agg <- cstats.aggregate(df)
```

Perform statistical analysis:

```r
anova_results <- cstats.2wanova(df_agg)
tukey_results <- cstats.tukey(df_agg, "Viability")
```

Generate plots:

```r
plots <- cstats.barplot.tukey(df_agg, "Viability")
plots[[1]]
```

**Note:** The example assumes that the dataset already contains a `STAGE` column. This can be created from the original time points using `dplyr::mutate()`.

## Statistical considerations

The toolkit currently supports independent-sample comparisons and standard two-way ANOVA.

Users should verify experimental independence, model assumptions, and the appropriate unit of replication before performing statistical tests.

Pairwise t-tests and Wilcoxon tests do not currently apply multiple-comparison adjustments.

## Example dataset

The repository includes a fully synthetic cell culture dataset for training purposes, featuring:

- Three experimental conditions: Control, Drug0, and Drug1.
- Nine time points.
- Four biological replicates and three technical replicates.
- Three quantitative measurements: cell viability, marker expression, and proliferation.

The accompanying Excel workbook illustrates how experimental data can be organized and consolidated into a format suitable for R.

## License

License to be defined.
