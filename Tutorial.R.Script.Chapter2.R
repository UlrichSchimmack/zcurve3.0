### Install required packages (only needs to be done once)
install.packages(c("zcurve", "KernSmooth", "stringr", "parallel"))

### Set working directory to the folder containing the zcurve3.0 function
setwd("C:\\Users\\ulric\\Dropbox\\PHPCurves\\DOCS\\z-curve\\Tutorial")

### Load the zcurve function fron github (renamed as zcurve3)
zcurve3 = "https://raw.githubusercontent.com/UlrichSchimmack/zcurve3.0/refs/heads/main/Zing.25.06.09.test.R"
source(zcurve3)

###############################################################
### Chapter 2 – Analyzing the OSC Reproducibility Project
###############################################################

### 1.0 – DATA PREPARATION

### Download data directly from OSF
url <- "https://osf.io/download/fgjvw/"
osc <- read.csv(url)
osc <- data.frame(osc)
dim(osc)  # Should return 168 x 138

### Extract p-values from original and replication studies
pval <- data.frame(
  osc$"Authors..O.",
  osc$"T_pval.recalc..O.",
  osc$"T_pval..R."
)

### Keep only rows with complete p-values
pval <- pval[complete.cases(pval), ]
colnames(pval) <- c("authors", "ori", "rep")

### Check distribution of original p-values
table(cut(pval$ori, c(0, .01, .05, .10, .20, 1)))
# Interpretation: how many are 
#    "really" significant, p < .01
#    "just" significant, p < .05 & p > .01
#    "marginally" significant, p > .05 & p < .10
#    "marginally marginally significant, p < .20 & p > .20
#    "clearly not" significant, p > .20


### Contingency table of significance (ori vs. rep)
tab <- table(pval$ori < .05, pval$rep < .05)
tab
tab / sum(tab)            # Proportions
tab[2, ] / sum(tab[2, ])  # Replication success for p < .05 originals

### Only ~33% of original p < .05 results replicated

### Convert two-sided p-values to z-values
z.ori <- qnorm(1 - pval$ori / 2)
z.rep <- qnorm(1 - pval$rep / 2)

### alternative formula: -qnorm(pval$ori/2)

### Scatterplot of original vs. replication z-values
plot(z.ori, z.rep, xlim = c(0, 6), ylim = c(0, 6))
abline(a = 0, b = 1)  # 45-degree line
# Replication z-values are generally lower

#######################################################################
### 2.0 – Z-CURVE ANALYSIS OF ORIGINAL STUDIES
#######################################################################

### 2.1 Run basic z-curve model
source(zcurve3)
Zing(z.ori)

### 2.2 Adjust plot parameters for better visualization
source(zcurve3)
ymax <- 1.0
ymin <- -0.05
hist.bar.width <- 0.4
Zing(z.ori)

### 3.0 – TESTING FOR PUBLICATION BIAS

### Set parameters for bias detection
Int.Beg <- 0        # Include all z > 0 (significant and non-significant)
TEST4BIAS <- TRUE   # Enable bias test
just <- 0.6         # "Just significant" region: z = 2.0 to 2.6 (default)
Zing(z.ori)

### Result: Significant excess of just-significant results (p = .0025)
### Interpretation: Evidence of bias
### Could be selection bias or p-hacking

### 4.0 – TEST FOR P-HACKING

### P-hacking inflates just-significant results but not high z-values
### Restrict model to z-values > 2.6 (just above just-significant region)
Int.Beg <- z.crit + just
Zing(z.ori)

### p > .05, fallacy to infer no p-hacking from non-significant result
### Try more sensitive threshold (just = 0.4)
just <- 0.4
Int.Beg <- z.crit + just
Zing(z.ori)

### Result: p = .0009
### Did we just p-hack a p-hacking test? :)
### Interpretation: Even after correction for multiple comparisons (α = .025),
### the result is significant → consistent with p-hacking

### 5.0 – TEST FOR HETEROGENEITY (slow!)
source(zcurve3)
het.res = TEST4HETEROGENEITY(z.ori,500)
round(het.res,3)

### Interpretation: SD CI = [0.50, 1.71] → no strong evidence of heterogeneity
### Note: z > 6 studies excluded from this test (assumed power ≈ 1)

### Summary of Results
# (a) Evidence of bias
# (b) Evidence of p-hacking
# (c) Surprisingly little heterogeneity, despite varied study topics and methods

### 6.0 – FINAL MODEL (EM algorithm + bootstrapping)

### When bias is present, ignoring it overestimates average power.
### When p-hacking is present, zcurve underestimates EDR. 
### Best practice: Use selection model, treat EDR as lower-bound.

source(zcurve3)
ymax <- 1.0
ymin <- -0.05
hist.bar.width <- 0.4
Est.Method <- "EM"      # Use EM algorithm for better estimates
boot.iter <- 100        # Bootstrap CIs
z.res <- Zing(z.ori)    # Store results

### Final Results Summary:
# 1. Bias: 91% significant results, EDR = 5%-73%, confirms bias 
# 2. EDR → False Discovery Risk ≈ 20% (CI: 2% to 94%), inconclusive about false positives
# 3. ERR = 60% (CI: 41% to 77%) → expected replication rate for hypothetical exact replications
# 4. Predicted actual replication rate to be between EDR and ERR (~41%) 
#    Actual replication rate was 33%, between EDR and ERR. 
# 5. Local power estimates suggest ...
#      ~70% replicability for z = 3-4
#      ~80% replicability for z = 4-5
#      -85% replicability for z = 5-6
#      z > 6 implies 100% replicability

### → Chapter 3 will examine the replication results in detail
	  and test the z-curve predictions based on these results. 


### 7.0 – ADJUSTING ALPHA 
### The false discovery risk (FDR) depends on alpha
### Lower alpha to see the FDR
source(zcurve3)
ymax <- 1.0
ymin <- -0.05
hist.bar.width <- 0.4
Est.Method <- "EM"      # Use EM algorithm for better estimates
boot.iter <- 500        # Bootstrap CIs
alpha <- .01            # LOWER ALPHA TO .01, .005 
z.res <- Zing(z.ori)    # Store results

### With alpha = .01, the FDR decreases to 4%, and
### the upper limit of the 95%CI shrinks from 94% to 14%
### Only 56% significant results at p < .01

### With alpha = .005, the FDR decreases to 2%, and
### the upper limit of the 95%CI shrinks to 9%
### Only 45% significant results at p < .005

### Trade off: Lower false positive risk, but higher false negative risk with alpha = .005


### Summary
### Chapter 2 introduced you to ... 
### 1. hist.bar.width #make histograms look better with small datasets
### 2. a Int.Beg <- 0;   # Use all z-values > 0 and test for bias
###    b TEST4BIAS <- TRUE   
### 3. a Int.Beg <- z.crit + just;   # Use only "really" significant results 
###    b TEST4BIAS <- TRUE			  # to fit model and test for p-hacking   
### 4. TEST4HETEROGENEITY <- bootstraps of test #test heterogeneity
### 5. Est.Method <- "EM"            # Use EM algorithm for better estimates
### 6. boot.iter <- 100              # Get confidence intervals
### 7. alpha =                       # change alpha level, default alpha = .05


