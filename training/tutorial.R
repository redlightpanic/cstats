################################################################################################################
# cstats | HANDS-ON TUTORIAL (30 - 45 min)
# From Excel to R: data wrangling, statistics and plots
# Training dataset: synthetic cell-culture experiment (NOT real research data)
################################################################################################################

# BEFORE STARTING
# 1. Open tutorial.xlsx and review INTRO / Drug0 / Drug1 / Control / R_DATA.
# 2. Export the R_DATA sheet as a tab-separated TSV named training.tsv.
# 3. Save training_data.csv in the same folder as this script and cstats.R.
# 3. In RStudio, open this script and run commands with Ctrl+Enter.
# 4. Install required packages ONLY ONCE:
#    install.packages(c("dplyr", "tidyr", "ggplot2", "factoextra", "patchwork", "rstatix", "ggpubr", "emmeans"))

################################################################################################################
# HOWTO 1 | LOADING THE TOOLKIT AND THE DATA (1 min)
################################################################################################################
# Load toolkit
source("cstats.R")

# Read data
df <- cstats.read("training.tsv", separator = ",")

################################################################################################################
# HOWTO 2 | INSPECTING THE DATA (5 min)
################################################################################################################
# Check / Verify data
head(df)
str(df)
names(df)
dim(df)

# Inspect data - BASIC
df[1,1]         # Row 1, Column 1 (Time)
df[1:5,1]       # Rows 1 to 5, Column 1 (Time)
df[1,]          # Row 1, All Columns
df[1, c(1,3,5)] # Row 1, Columns 1 (Time) + 3 () + 5 ()

################################################################################################################
# 3 | MANIPULATE THE DATA
################################################################################################################
# Inspect the groups
count(df, Treatment)
count(df, Time)

# Filter a subset (the original df is unchanged)
control1 <- df[df$Treatment == "control",]
control2 <- df %>% filter(Treatment == "control")
head(control)

# Select specific columns
measurements <- df %>% select(Time, Treatment, BioRep, Viability)
head(measurements)

# Create a new column with the nine time points grouped into three stages
# 0/6/12 h = EARLY; 24/36/48 h = MID; 60/72/96 h = LATE
df <- df %>%
  mutate(Stage = case_when(
    Time %in% c(0, 6, 12)   ~ "EARLY",
    Time %in% c(24, 36, 48) ~ "MID",
    Time %in% c(60, 72, 96) ~ "LATE",
    TRUE ~ NA_character_
  ))
stopifnot(!anyNA(df$Stage))
df$Stage <- factor(df$Stage, levels = c("EARLY", "MID", "LATE"))

# Optional: inspect summary values without changing the original dataset
df %>%
  group_by(Stage, Treatment) %>%
  summarise(mean_viability = mean(Viability), .groups = "drop")

################################################################################################################
# 3 | DEFINE FACTORS, VARIABLES AND REPLICATES (5 min)
################################################################################################################
# Define Factors
df <- cstats.set_factors(df, c("Stage", "Treatment"))

# Define variables
df <- cstats.set_vars(df, c("Viability", "Marker", "Proliferation"))

# Define biological replicates
df <- cstats.set_reps(df, "BioRep")

################################################################################################################
# 4 | AGGREGATE TECHNICAL REPLICATES (5 min)
################################################################################################################

# Here it averages the 3 technical replicates AND the 3 time points per STAGE for each biological replicate.
df_agg <- cstats.aggregate(df)
head(df_agg)
dim(df_agg)   # Expected: 3 stages x 3 treatments x 4 BioRep = 36 rows
count(df_agg, STAGE, Treatment)

################################################################################################################
# 5 | STATISTICAL TESTS (10 min)
################################################################################################################
# Two-way ANOVA with interaction: STAGE * Treatment
# Requires independent experimental units and model assumptions.
anova_results <- cstats.2wanova(df_agg)
anova_results$Viability

# Tukey-adjusted treatment comparisons WITHIN each stage
# (not comparisons of all stage-treatment combinations)
tukey_results <- cstats.tukey(df_agg, "Viability")
tukey_results

# Welch t-tests and Mann-Whitney/Wilcoxon rank-sum tests:
# Compare treatments WITHIN the EARLY stage (first factor = STAGE).
t_results <- cstats.t(df_agg, "Viability", base = "EARLY")
t_results
wilcox_results <- cstats.wilcox(df_agg, "Viability", base = "EARLY")
wilcox_results

# Note: cstats.t and cstats.wilcox use unpaired tests, and their pairwise
# p-values are NOT corrected for multiple comparisons.

###############################################################################
# 6 | PLOT DATA AND TEST RESULTS (10 min)
###############################################################################
# Plot using tukey statistical test
plots_tukey <- cstats.barplot.tukey(df_agg, "Viability")
names(plots_tukey)
plots_tukey$EARLY

# Plot using Welsch t-test statistical test
plots_t <- cstats.barplot.t(df_agg, "Viability")
plots_t$MID

# Plot using Wilcoxon Ranked Sum statistical test
plots_wilcox <- cstats.barplot.wilcox(df_agg, "Viability")
plots_wilcox$LATE

# ggplot objects can be customized after they are created.
my_plot <- plots_tukey$EARLY + ggplot2::labs(title = "EARLY stage: cell viability", y = "Viability (%)")
print(my_plot)

# Export the figure
# ggplot2::ggsave("viability_early.png", plot = my_plot, width = 7, height = 5, dpi = 300)
