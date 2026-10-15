################################################################################################################
# cstats | HANDS-ON TUTORIAL
# From Excel to R: data wrangling, statistics and plots
# Training dataset: synthetic cell-culture experiment (NOT real research data)
################################################################################################################

# BEFORE STARTING
# 1. Open tutorial.xlsx and review INTRO / Drug0 / Drug1 / Control / tutorial.
# 2. Export the "tutorial" sheet as a tab-separated TXT named tutorial.txt.
# 3. Save tutorial.txt in the training folder.
# 3. In RStudio, open this script
# 4. Install required packages ONLY ONCE:
#    install.packages(c("dplyr", "tidyr", "ggplot2", "factoextra", "patchwork", "rstatix", "ggpubr", "emmeans"))

################################################################################################################
# HOWTO 1 | LOADING THE TOOLKIT AND THE DATA (1 min)
################################################################################################################

# Load toolkit
source("cstats.R")

# Read data
df <- cstats.read("training/tutorial.txt")

################################################################################################################
# HOWTO 2 | INSPECTING THE DATA (10 min)
################################################################################################################

# Check / Verify data
head(df)
head(df, n = 10)
tail(df)
colnames(df)
dim(df)

# Inspect data - BASIC
df[1,5]           # Row 1, Column 5 (Viability)
df[1,"Viability"] # Row 1, Column Viability
df[1,1:4]         # Row 1, Columns 1 to 4 (Time, Treatment, BioRep, TechRep)
df[4,]            # Row 4, All Columns
df[1, c(1,3,5)]   # Row 1, Columns 1, 3 and 5 (Time, BioRep, Viability)
df[4,c("Time","Proliferation")]

#!! EXERCISE : Inspect columns Viability, Marker and Proliferation of rows 10 to 20

# Inspect data - ADVANCED
df$Viability
df$Viability[5:10]
df$Viability[df$Treatment == "Drug0"]
df[df$Treatment == "Drug0",]

# Inspect data - STATS
mean(df$Viability)
sd(df$Viability)
summary(df$Viability)
table(df$Treatment)
table(df$Treatment, df$Time)

#!! EXERCISE : find the average (mean) of the Viability of the Control samples


################################################################################################################
# HOWTO 3 | MANIPULATING THE DATA (10 min)
################################################################################################################

# Filter a subset
control <- df[df$Treatment == "Control",]
dim(control)

# Select specific columns
measurements <- df[,c("Time", "Treatment", "BioRep", "Viability")]
head(measurements)

#!! EXERCISE : create a subset comprised of Time and Proliferation for Control samples only

# Create a new "Effect" column equal to Viability * Marker
df$Effect <- df$Viability * df$Marker
head(df)
df$Effect <- NULL
head(df)

################################################################################################################
# HOWTO 4 | DEFINING FACTORS, VARIABLES AND REPLICATES (2 min)
################################################################################################################

# Define Factors
df <- cstats.set_factors(df, c("Time", "Treatment"))

# Define variables
df <- cstats.set_vars(df, c("Viability", "Marker", "Proliferation"))

# Define biological replicates
df <- cstats.set_reps(df, "BioRep")

################################################################################################################
# HOWTO 5 | AGGREGATING TECHNICAL REPLICATES (2 min)
################################################################################################################

# Here it averages the 3 technical replicates per Time and Treatment, for each biological replicate.
df_agg <- cstats.aggregate(df)
head(df_agg)

################################################################################################################
# HOWTO 6 | PERFORMING STATISTICAL TESTS (10 min)
################################################################################################################

# Two-way ANOVA with interaction: Time * Treatment
anova_results <- cstats.2wanova(df_agg)
anova_results$Viability

# Tukey-adjusted treatment comparisons
tukey_results <- cstats.tukey(df_agg, "Viability")
tukey_results

# Welch t-tests and Mann-Whitney/Wilcoxon rank-sum tests:
t_results <- cstats.t(df_agg, "Viability", base = "T12")
t_results
wilcox_results <- cstats.wilcox(df_agg, "Viability", base = "T72")
wilcox_results

###############################################################################
# HOWTO 7 | PLOTTING DATA AND TEST RESULTS (10 min)
###############################################################################

# Plot using tukey statistical test
plots_tukey <- cstats.barplot.tukey(df_agg, "Viability")
names(plots_tukey)
plots_tukey$T12
p1 <- plots_tukey$T12
p2 <- plots_tukey$T48
p3 <- plots_tukey$T72
p1 | p2 | p3

# Plot using Welch t-test statistical test
plots_t <- cstats.barplot.t(df_agg, "Marker")
plots_t$T48

# Plot using Wilcoxon Ranked Sum statistical test
plots_wilcox <- cstats.barplot.wilcox(df_agg, "Proliferation")
plots_wilcox$T72

