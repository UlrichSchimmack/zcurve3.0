### Install required packages (only needs to be done once)
install.packages(c("zcurve", "KernSmooth", "stringr", "parallel"))

### download zcurve functions and run locally
### Set working directory to the folder containing the zcurve3.0 function
setwd("C:\\Users\\ulric\\Dropbox\\PHPCurves\\DOCS\\z-curve\\Tutorial")
zcurve3 = "Zing.25.07.11.test.R"
source(zcurve3)

### load zcurve function from github (no download needed)
### Load the zcurve function fron github (renamed as zcurve3)
zcurve3 = "https://raw.githubusercontent.com/UlrichSchimmack/zcurve3.0/refs/heads/main/Zing.25.07.11.test.R"
source(zcurve3)

#########################################################################
### Chapter 3 – Analyzing the OSC Reproducibility Project - Replications
#########################################################################

### 1.0 – DATA PREPARATION (not required, if you just worked on Chapter 2)

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
table(cut(pval$rep, c(0, .01, .05, .10, .20, 1)))
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

### limit maximum to 10 (has no influence on z-curve estimates)
z.ori[z.ori > 10] = 10
z.rep[z.rep > 10] = 10

### test correlation
cor(z.ori,z.rep)
# result: r = .4, studies with higher power produce stronger z-scores in ori and rep

### check replication rates for different levels of z-values in ori
tab = table(cut(z.ori,c(0,1.5,2,2.6,3,4,6,11)),cut(z.rep,c(0,1.5,2,2.6,3,4,6,11)))
rownames(tab) = paste0("ori",rownames(tab))
colnames(tab) = paste0("rep",colnames(tab))
tab
# result: z-values up to 4 are more likely to end with failure than success

#######################################################################
### 2.0 – Z-CURVE ANALYSIS OF REPLICATION STUDIES
#######################################################################

### 2.1 Run basic z-curve model
source(zcurve3)
Zing(z.rep)

### 2.2 Adjust plot parameters for better visualization
source(zcurve3)
ymax <- 0.8
ymin <- -0.05
hist.bar.width <- 0.4		# adjust the bar.width for a smooth histogram
Zing(z.rep)
# ODR == EDR suggests no publication bias


### 3.0 – TESTING FOR PUBLICATION BIAS

### Set parameters for bias detection
Int.Beg <- 0        # Include all z > 0 (significant and non-significant)
TEST4BIAS <- TRUE   # Enable bias test
just <- 0.6         # "Just significant" region: z = 2.0 to 2.6 (default)
Zing(z.rep)

### Result: No evidence of bias
### Also means, no evidence of reverse p-hacking
### reverse p-hacking: turning significant results into non-signifiant results
### to make replication studies more interesting

### 4.0 – TEST FOR P-HACKING

### P-hacking inflates just-significant results but not high z-values
### Restrict model to z-values > 2.6 (just above just-significant region)
Int.Beg <- z.crit + just
Zing(z.rep)

### no evidence of p-hacking

### 5.0 – TEST FOR HETEROGENEITY (slow!)
Int.Beg = 0                 # reset to 0 to use all z-values 
TEST4HETEROGENEITY = 50     # 0 means FALSE, > 0 is # of bootstrap runs
het.res = Zing(z.rep)
round(het.res,3)

### Interpretation: SD M = 180, 90CI = [1.20 - 2.20] evidence of heterogeneity
### Note: z > 6 studies excluded from this test (assumed power ≈ 1)
### comparison of model 2 (2 components vs. model 3 (SD > 1) shows no difference
### a model with a normal distribution of z-values fits the data. 
### SD of non-central z-values is SD - 1.  

### 5.1a - Fit Model with One Free Component
source(zcurve3)
ymax <- 0.8
ymin <- -0.05
hist.bar.width <- 0.4		# adjust the bar.width for a smooth histogram
Int.Beg = 0               # reset to 0 to use all z-values 
Est.Method = "EXT"
ncz = 2                   # specifies 1 component, value is starting value
zsd = 1                   # fix at 1
boot.iter = 500           # get confidence intervals for all parameters
z.res = Zing(z.rep)

### 5.2a – Fit Model with Two Free Components
source(zcurve3)
ymax <- 0.8
ymin <- -0.05
hist.bar.width <- 0.4		# adjust the bar.width for a smooth histogram
Int.Beg = 0                 # reset to 0 to use all z-values 
Est.Method = "EXT"
ncz = c(0,2)              # specifies 1 component, value is starting value
zsd = 1                   # fix at 1
boot.iter = 500           # get confidence intervals for all parameters
z.res = Zing(z.rep)


### 6.1 – Fit Model with EM algorithm
source(zcurve3)
ymax = .8
ymin = -.04
hist.bar.width <- 0.4     # adjust the bar.width for a smooth histogram
Int.Beg = 0               # reset to 0 to use all z-values 
boot.iter = 500           # get confidence intervals for all parameters
Est.Method = "EM"         # Use the EM algorithm of zcurve.2.0
z.res = Zing(z.rep)


### 7 - Estimating False Positive Risk of Non-Significant Replications

### FDR in the previous figure applies to significant results in replication studies
### The more interesting question is the false positive risk in the original studies
### We can estimate it, by computing the EDR of the non-significant results
### The ODR is, of course, 0, but the EDR estimates power taking the selection
### for non-significance into account. 

z.rep.ns = z.rep[z.rep < 1.96]  # select only non-significant results
z.res = Zing(z.rep.ns)          # fit z-curve 
 
#Result: The EDR is low. Even the upper limit of the 95%CI is only 22%. 
#        The low EDR implies a high FDR, but the 95%CI is wide. 
#        We cannot say that there are many false positive results,
#        but we can say that many could be fales positive results. 

### Summary
### Chapter 3 introduced you to ... 
### 1. Est.Method <- "EXT"           # Use EXT for free means and SDs of ncz
### 2. Use only non-significant results of replication studies to estimate 
       FDR of original studies. 


