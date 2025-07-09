

### Install required packages (only needs to be done once)
install.packages(c("zcurve","KernSmooth","stringr","parallel"))

### Set working directory to the folder containing the zcurve.3.0 function
setwd("C:\\Users\\ulric\\Dropbox\\PHPCurves\\R-Code\\")

### Rename the zcurve function -> zcurve3 and run the code
### I use = instead of ->
zcurve3 = "Zing.25.06.09.test.R"
source(zcurve3)

### 

### set the working directory (folder) with the the zcurve.3.0 function
setwd("C:\\Users\\ulric\\Dropbox\\PHPCurves\\R-Code\\")

### make the function available with the source function
zcurve3 = "Zing.25.06.09.test.R"

###############################################################
### Chapter 1 - Chapter 1 - Chapter 1 - Chapter 1 - Chapter 1 
###############################################################

### 1.0 Simulate z-scores with 40% power (for Figure 1)
### We simulate z-values from a normal distribution with mean = 1.7, sd = 1
### This gives about 40% power for detecting z > 1.96
zval = rnorm(10000, 1.7)

### 1.1 Select only the statistically significant results
fig1.zval = zval[zval > 1.96]

### 1.2 Compute power using one of three equivalent approaches
pnorm(1.96, 1.7, lower.tail = FALSE)  # Method 1: recommended
1 - pnorm(1.96, 1.7)                  # Method 2: same as above
pnorm(1.7, 1.96)                      # Method 3: valid only for z-values

### 1.3 check power, version 3 (the simplfied version that changes observed value and mean of distribution)
### only works with z-scores, do not do this with t-values
pnorm(1.7,1.96)

#######################################################################

### 2. RUNNING Z-CURVE

### 2.1 Load the Zing function again to reset parameters
source(zcurve3)

### 2.2 Set model to assume a single true non-central z (ncz = 1.7)
ncz = 1.7

### 2.3. Generate z-curve plot (note: axis scaling may need adjustments)
Zing(fig1.zval)

### 2.4 Improve plot appearance by adjusting y-axis limits
ymax = 1.2   # upper limit of histogram bars
ymin = -0.06 # allows local power estimates to display below x-axis

### 2.5 Regenerate the figure with new y-axis
Zing(fig1.zval)

###############################################################
### Extension of Chapter 1 example
###############################################################

### 2.6 Try a mismatched model with ncz = 2.8 (80% power, p ≈ .005)

source(zcurve3)  # Optional to reset to default
ymax = 1.2       # Optinal to respecify model
ymin = -0.06     # Optinal to respecify model
ncz = 2.8        # set the false parameter (ncz = 2.8 == power = 80%)
Zing(fig1.zval)  # Visual inspection shows poor model fit


### 3. LET ZING ESTIMATE THE MODEL

### Now let Zing estimate the model using default 7 components (ncz = 0:6)
source(zcurve3)
ymax = 1.2
ymin = -0.06
zing.res = Zing(fig1.zval)  # Save output

### Print result
zing.res

### Commentary:
# - Estimated ERR is slightly too high (~5 pp)
# - Estimated EDR is slightly too low (~5 pp)

### Model weights interpretation:
# - w.all: weights for full model (used to compute EDR)
# - w.sig: weights for significant values only (used to compute ERR)
# - Example: mixture of ncz=1 and ncz=2 may average to ~1.5
#   (i.e., 0.5*1 + 0.5*2 = 1.5), power = pnorm(1.5,1.96) = 32%
# - the w.sig estimates overestimate the weight for component ncz = 2
    (i.e., 0.2*1 + 0.8*2) = 1.8, power = pnorm(1.8,1.96) = 44%
# - This is an artifact of homogeneous data falling between fixed ncz values.
# - Not a concern with real data where power is more variable.


### 3. BOOTSTRAPPED CONFIDENCE INTERVAL

### Add bootstrapped 95% CIs for the estimates
### For testing, use a small number like 50
### For final results, use at least 500 (results converge with increasing iterations)

### 4. Final Notes

### Summary:
# - You now know how to: 
#   (1) install required packages
#   (2) simulate z-values with known power
#   (3) run zcurve3.0 to fit and visualize the model
#   (4) estimate confidence intervals


### If you run into problems, try using an AI or reach out to:
### ulrich.schimmack@utoronto.ca


### Bonus Exercise


### Let's have some fun with z-curve

# Simulate two sets of z-values:
# - 10,000 moderately powered ~50% power (mean = 2, SD = 1)
# - 5,000 artificially narrow (mean = 2, SD = 0.1) to simulate p-hacking
fun.dat = c(rnorm(10000, 2, 1), rnorm(5000, 2, 0.1))

# Limit the z-scores to the range 0.75 to 3.25
fun.dat = fun.dat[fun.dat > 0.75 & fun.dat < 3.25]

# Load z-curve function
source(zcurve3)

# Set plotting parameters
hist.bar.width = 0.5  # wide bins to highlight central peak
ymax = 2
ymin = -0.08
Int.Beg = 0
Int.End = 6
x.lim.min = 0.25
x.lim.max = 3.75

# Run z-curve model and plot
Zing(fun.dat)

### This is sometimes called a "finger plot" (you’ll see why...)
# - We restrict z-scores to between 0.75 and 3.25
# - Histogram bars are 0.5 wide, so we get 5 bars
# - Because of the spike near z = 2.0, the middle bar is overrepresented
# - This creates a visual "middle finger" — a playful way to highlight
#   how marginally significant results are overrepresented
# - Despite a mean of 2 (suggesting ~50% power), the pattern indicates possible bias



