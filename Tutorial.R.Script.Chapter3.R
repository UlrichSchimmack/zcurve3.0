

###############################################################################
# Chapter 3 – Analyzing the OSC Reproducibility Project - Replications
###############################################################################

### 0.0 – SETUP

# Install required packages (run once only)
install.packages(c("zcurve", "KernSmooth", "stringr", "parallel"))

# Load z-curve function from a local file
setwd("C:/Users/ulric/Dropbox/PHPCurves/DOCS/z-curve/Tutorial")
zcurve3 <- "Zing.25.07.11.test.R"
source(zcurve3)

# Alternatively, load z-curve function directly from GitHub
zcurve3 <- "https://raw.githubusercontent.com/UlrichSchimmack/zcurve3.0/refs/heads/main/Zing.25.07.11.test.R"
source(zcurve3)

### 1.0 – DATA PREPARATION (skip if already done in Chapter 2)

# Download OSC reproducibility data from OSF
url <- "https://osf.io/download/fgjvw/"
osc <- read.csv(url)
osc <- data.frame(osc)
dim(osc)  # Expected dimensions: 168 x 138

# Extract p-values from original and replication studies
pval <- data.frame(
  authors = osc$"Authors..O.",
  ori     = osc$"T_pval.recalc..O.",
  rep     = osc$"T_pval..R."
)

# Keep rows with complete p-values only
pval <- pval[complete.cases(pval), ]

# Check distribution of replication p-values
table(cut(pval$rep, c(0, .01, .05, .10, .20, 1)))
# Interpretation:
#   < .01 = "clearly significant"
#   .01 – .05 = "just significant"
#   .05 – .10 = "marginal"
#   .10 – .20 = "borderline"
#   > .20 = "not significant"

# Create a contingency table comparing significance in originals vs replications
tab <- table(pval$ori < .05, pval$rep < .05)
tab
tab / sum(tab)             # Proportion of each cell
tab[2, ] / sum(tab[2, ])   # Replication success rate for significant originals

# Note: Only about 33% of originally significant results replicated

### 1.1 – Convert p-values to z-values
z.ori <- qnorm(1 - pval$ori / 2)
z.rep <- qnorm(1 - pval$rep / 2)

# Truncate extreme z-values for display
z.ori[z.ori > 10] <- 10
z.rep[z.rep > 10] <- 10

# Plot: z-scores for original vs. replication studies
plot(z.ori, z.rep, xlim = c(0, 6), ylim = c(0, 6))
abline(a = 0, b = 1)  # 45° line
# Note: Replication z-values are generally lower

# Correlation between original and replication z-values
cor(z.ori, z.rep)  # Expected: ~.40

# Check replication success by bins of z-ori
tab <- table(
  cut(z.ori, c(0, 1.5, 2, 2.6, 3, 4, 6, 11)),
  cut(z.rep, c(0, 1.5, 2, 2.6, 3, 4, 6, 11))
)
rownames(tab) <- paste0("ori_", rownames(tab))
colnames(tab) <- paste0("rep_", colnames(tab))
tab

### 2.0 – Z-CURVE ANALYSIS

# Run basic z-curve model
source(zcurve3)
Zing(z.rep)

# Adjust plot appearance
ymax <- 0.8
ymin <- -0.05
hist.bar.width <- 0.4
Zing(z.rep)
# Note: ODR ≈ EDR suggests no publication bias

### 3.0 – TESTING FOR PUBLICATION BIAS

Int.Beg     <- 0       # Include all z > 0
TEST4BIAS   <- TRUE    # Enable test
just        <- 0.6     # Define “just significant” region: z = 2.0–2.6
Zing(z.rep)

# Result: No evidence of publication bias or "reverse p-hacking"
# Reverse p-hacking = altering results to appear less significant

### 4.0 – TEST FOR P-HACKING

# Restrict analysis to z > 2.6 (above just-significant range)
Int.Beg <- z.crit + just
Zing(z.rep)

# Result: No evidence of p-hacking in replications

### 5.0 – TEST FOR HETEROGENEITY (may take time)

Int.Beg             <- 0        # Use all z-values
TEST4HETEROGENEITY  <- 500      # Enable test with 500 bootstrap samples
z.res <- Zing(z.rep)
round(z.res$fit.comp,3)

# Interpretation:
# SD ~1.80, 90% CI = [1.20, 2.20] → supports heterogeneity
# Note: z > 6 excluded from estimation (assumed power ≈ 1)

### Illustrating Heterogeneity Test

### 5.1 – Fit Model with One Component, Free Mean, Fixed SD

source(zcurve3)
ymax <- 0.8
ymin <- -0.05
hist.bar.width <- 0.4
Int.Beg <- 0
Est.Method <- "EXT"
ncz <- 2         # Starting value for 1 component
z.res <- Zing(z.rep)
z.res 
#Result:  - model does not predict significant results


### 5.2 – Fit Model with One Component, Free Mean, Free SD

source(zcurve3)
ymax <- 0.8
ymin <- -0.05
hist.bar.width <- 0.4
Int.Beg <- 0
Est.Method <- "EXT"
ncz <- 2         # Starting value for 1 component
SD.FIXED = FALSE
zsd = 4
z.res <- Zing(z.rep)
z.res 
#Result:  - model fits, but the ERR estimate is low 


### 5.3 – Fit Model with Two Free Components

source(zcurve3)
ymax <- 0.8
ymin <- -0.05
hist.bar.width <- 0.4
Int.Beg <- 0
Est.Method <- "EXT"
ncz <- c(0, 2)    # Start values for 2 components
zsd <- 1
z.res <- Zing(z.rep)
#Result:  - model fits, and shows higher ERR 
#This result is preferred because it does not make assumptions about 
#distribution of true power of studies. 
#This model allows for one set of low powered and one set of high powered studies
z.res
#Results
# - One component mean close to zero with high weight (.73)
# - One component mean high (3.27) with low weight (.27)



### 6.0 – Fit EM Algorithm Model

source(zcurve3)
ymax <- 0.8
ymin <- -0.04
hist.bar.width <- 0.4
Int.Beg <- 0
Est.Method <- "EM"   # Expectation-Maximization
boot.iter <- 500
z.res <- Zing(z.rep)


### 7.0 – FALSE POSITIVE RISK FOR NON-SIGNIFICANT REPLICATIONS

# Estimating FDR for non-significant replication results

z.rep.ns <- z.rep[z.rep < 1.96] # subset of only non-significant results

source(zcurve3)
ymax <- 1
ymin <- -0.05
hist.bar.width <- 0.4
Int.Beg <- 0
Int.End <= z.crit
Est.Method <- "EXT"   # To get confidence intervals, not yet possible with "EM"
ncz = c(0,3)
boot.iter <- 500
z.res <- Zing(z.rep.ns)

# Interpretation:
# ODR underestimates power because of selection for non-significance
# EDR corrects for this. 
# EDR can be converted into FDR 
# FDR is high; 95% CI upper bound includes 100%...
# but lower value is 27%
# Absence of Evidence, but not Evidence of Absence of an Effect



###############################################################################
# Summary Notes:
#
# 1. Use EXT method to estimate freely varying means and SDs (Est.Method = "EXT").
# 2. Use all tests when there is no evidence of bias
# 3. Deeper understanding of heterogeneity test
# 4. !!! Use non-significant replication results to estimate FDR of original findings.
###############################################################################



