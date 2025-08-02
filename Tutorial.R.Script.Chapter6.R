



###############################################################################
# Chapter 6 – Comparing default z-curve and NDIST-zcurve 
###############################################################################

### 0.0 – SETUP

# disable scientific notation
options(scipen = 999)

# Set work directory
setwd("C:/Users/ulric/Dropbox/PHPCurves/DOCS/z-curve/Tutorial")

# Install required packages (run once only)
#install.packages(c("zcurve", "KernSmooth", "stringr", "parallel"))

# Load z-curve and P-curve functions from a local file
setwd("C:/Users/ulric/Dropbox/PHPCurves/DOCS/z-curve/Tutorial")
zcurve3 <- "Zing.25.07.11.test.R"
source(zcurve3)

# Alternatively, load z-curve function directly from GitHub
zcurve3 <- "https://raw.githubusercontent.com/UlrichSchimmack/zcurve3.0/refs/heads/main/Zing.25.07.11.test.R"
source(zcurve3)

### 1.0 – ILLUSTRATION WITH EXAMPLES

### 1.1 - SIMULATING A NORMAL DISTRIBUTION OF NON-CENTRAL Z-VALUES (POWER) 

set.seed(20250731)
sim.ncz = rnorm(200000,2,.5)                   # simulate the distribution of the non-central z-values
table(sim.ncz > 1 & sim.ncz < 3)              # ~95% of ncz-values are within a range from 1 to 3
z.crit = qnorm(1-.05/2)                       # the critical z-value for alpha = .05, two-sided
power = 1-pnorm(z.crit,sim.ncz)               # convert ncz to power 
hist(power)
table(power > .17 & power < .85)              # 95% of studies have power between 17% and 85%
true.edr = mean(power);true.edr               # average power / expected discovery rate
true.err = mean(power^2)/mean(power);true.err # average power of significant studies / expected replication rate

### now we add sampling error 
sim.z = unlist(lapply(sim.ncz,function(x) rnorm(1,x) )) 
sim.z = abs(sim.z)                            # sign is not important; z-curve uses absolute values


### Analyze with default z-curve
source(zcurve3)
Title = "Default z-curve: True EDR = 51%, TRUE EDR = 58%"
res.1 = Zing(sim.z)
res.1$res                # ERR estimate is good, EDR is underestimated by 5 percentage points
###
 

### Analyze with extended z-curve: 1 component with free SD
source(zcurve3)
Est.Method = "EXT"
ncz = 2  # starting value, free to fit the data
zsds = 5 # starting value, free to fit the data 
res.2 = Zing(sim.z)
res.2$res                # ERR and EDR estimates are good
res.2$ncz                # ncz estimate is 1.97, close to the simulated value of 2
res.2$zsds               # zsds estimate is 1.13, close to the simulated value 
                         #	 true.zsds = sqrt(1 + .5^2); true.zsds


### The model that assumes a normal distribution outperforms default z-curve
### because the simulated distribution matches the model's assumption 
### The estimates of the default model are still useful 


### 1.2 - SIMULATING A NON-NORMAL DISTRIBUTION OF NON-CENTRAL Z-VALUES (POWER) 

set.seed(20250731)
weights = c(100,10,100)
weights/sum(weights)
weights
sim.ncz = c(rep(1,weights[1]*1000),rep(2,weights[2]*1000),rep(3,weights[3]*1000))
table(sim.ncz >= 1 & sim.ncz <= 3)            # 100% of ncz-values are within a range from 1 to 3
z.crit = qnorm(1-.05/2)                       # the critical z-value for alpha = .05, two-sided
power = 1-pnorm(z.crit,sim.ncz)               # convert ncz to power 
table(power)
hist(power)
table(power >= .16 & power <= .86)            # 100% of studies have power between 17% and 85%
true.edr = mean(power);true.edr               # average power / expected discovery rate
true.err = mean(power^2)/mean(power);true.err # average power of significant studies / expected replication rate

### now we add sampling error 
sim.z = unlist(lapply(sim.ncz,function(x) rnorm(1,x) )) 
sim.z = abs(sim.z)                            # sign is not important; z-curve uses absolute values

### Analyze with default z-curve
source(zcurve3)
Title = "Default z-curve: True EDR = 51%, TRUE EDR = 73%"
res.1 = Zing(sim.z)
res.1$res                # ERR estimate is good, EDR estimate is good
###
 

### Analyze with extended z-curve: 1 component with free SD
source(zcurve3)
Est.Method = "EXT"
ncz = 2  # starting value, free to fit the data
zsds = 5 # starting value, free to fit the data 
res.2 = Zing(sim.z)
res.2$res                # ERR estimate is good; EDR estimates is 14 percentage points too high
res.2$ncz                # ncz estimate is 1.97, close to the simulated value of 2
res.2$zsds               # zsds estimate is 1.24, too high
                         #	 true.zsds: sd(sim.ncz)


### The normal model "fails" to predict the distribution of the
### non-significant results because it assumes a normal distribution
### see z-curve plot



### 2.0 – SIMULATION 

sim.k.sig = 10000
sim.ncz = c(0,2,4) # change or add more components
sim.zsd = c(1,1,1) #
                   # has to stay at 1 to compute true power            
                   # other values require estimation with large sample simulations and no bias

# make sure that there is an equal number of means and sds of components
if (length(sim.ncz) != length(sim.zsd)) sim.zsd = rep(1,length(sim.zsd))


source(zcurve3)    # load zcurve3 and default settings

# create weights for the components 
# code needs to be changed if there are more than three components
weights = c()
for (i1 in seq(1,11) ) {
	for (i2 in seq(1,11) ) {
		for (i3 in seq(1,11) ) {
			add = c((i1-1)/10,(i2-1)/10,(i3-1)/10) 
			weights = rbind(weights,add)
}}}

weights = weights[rowSums(weights) == 1,]
dim(weights)

# 66 different combinations of weights ranging from 0 to 1 in .1 steps


# run the simulation and collect the z-curve estimates for 
# z-curve and p-curve 

res = c()

run.i = 60
weights[run.i,]

b = 1                 # begin 
e = 5                 # end for testing
e = nrow(weights);e   # end for all runs

for (run.i in b:e ) { # for loop to run the simulations

print(paste("Run: ",run.i))

sim.w.all = weights[run.i,];sim.w.all

sim.pow.dir = pnorm(abs(sim.ncz),z.crit);sim.pow.dir 
sim.sign.error = 1-pnorm(sim.ncz,-z.crit);round(sim.sign.error,3)
sim.pow = sim.pow.dir + sim.sign.error;sim.pow

sim.ext = pnorm(6,sim.ncz,lower.tail=FALSE)
sim.ext = sum(sim.ext*sim.w.all)
sim.ext

sim.w.all.ext = c(sim.w.all*(1-sim.ext),sim.ext)
sim.w.all.ext 

sim.pow.ext = c(sim.pow,1)

sim.w.sig.ext = sim.w.all.ext * sim.pow.ext
sim.w.sig.ext = sim.w.sig.ext / sum(sim.w.sig.ext)
sum(sim.w.sig.ext)

sim.edr = sum(sim.w.all.ext * sim.pow.ext);sim.edr
sim.err = sum(sim.w.sig.ext * c(sim.pow.dir,1));sim.err

sim.components = length(sim.ncz)

sim.z = c()
i = 3
for (i in 1:sim.components) sim.z = c(sim.z,
	rnorm(sim.k.sig*21*sim.w.all[i],sim.ncz[i],sim.zsd[i]) 
)

sim.z = abs(sim.z)
sim.z = sim.z[sim.z > z.crit]
sim.z = sim.z[order(runif(length(sim.z),0,1))]
sim.z = sim.z[1:sim.k]
table(sim.z > z.crit)
summary(sim.z)

# Run basic z-curve model with "OF" method
source(zcurve3)    # load zcurve3 every time to get default settings
Title = paste0("True Average Power of Significant Results (ERR): ",round(sim.err*100));Title
ymax = 1
ymin = -.05
boot.iter = 500
res.1 = Zing(sim.z);res.1
res.1 = c(res.1$res[2,],res.1$res[3,])
res.1

# Run Extended z-curve model with 1 component and free SD
# A one-component model with fixed SD is equivalent to p-curve
# The free SD parameter allows z-curve to model heterogeneity, ...
# assuming a normal distribution of non-central z-values
Est.Method = "EXT"    # specify extended method with free parameters (zcurve.3.0)
ncz = 2				# specify only one component with any starting value
zsds = 5				# specify SD with maximum value to allow for free estimation
res.2 = Zing(sim.z);res.2
res.2 = c(res.2$res[2,],res.2$res[3,])

res.run = c(run.i,weights[run.i,],sim.edr,sim.err,res.1,res.2)
res.run

res = rbind(res,res.run)

write.table(res,"sim10k.ZcurveVSZcurveNorm.dat") # write results each trial
                                   # can resume if stopped 
                                   # by loading the completed results

} # End of for loops

dim(res)  # 66 rows and 18 columns

#write the completed results / overwrites the data from the for loop
#write.table(res,"sim.ZcurveVSZcurveNorm.dat")


########################################################
### GET STORED RESULTS
########################################################

# load the saved data
setwd("C:/Users/ulric/Dropbox/PHPCurves/DOCS/z-curve/Tutorial")
check = read.table("sim10k.ZcurveVSZcurveNorm.dat",row.names=NULL)
check = data.frame(check)
check = check[,2:dim(check)[2]]
dim(check)

# careful, if you run this code any results in the res file will be replaced 
# with the saved data. 
#res = check  # use only if sure to replace res; run analysis with "res" matrix 

dim(res)

res[1,]

########################################################################
### Evaluation of ERR estimates: z-curve (OF) VS z-curve (EXT free ZSD)
########################################################################

# columns
# -  6: true ERR
# -  7: res.1
# - 13: res.2 

round(res[,c(1:4,6,7,13)],3)

res = res[order(res[,6]),]

lwd.ci = .5

plot(res[,6],res[,7],xlim=c(0,1),ylim=c(0,1),pch=16,cex=1,col="forestgreen",ylab="Estimated ERR",
	xlab = "Simulated True EDR")  # plot OF and ER estimates
par(new=TRUE)
plot(res[,6],res[,8],type="l",lty=2,lwd=lwd.ci,xlim=c(0,1),ylim=c(0,1),col="forestgreen",ylab="",xlab="")  # plot OF and ER estimates
par(new=TRUE)
plot(res[,6],res[,9],type="l",lty=2,lwd=lwd.ci,xlim=c(0,1),ylim=c(0,1),col="forestgreen",ylab="",xlab="")  # plot OF and ER estimates

par(new=TRUE)
plot(res[,6],res[,13],pch=15,cex=1,xlim=c(0,1),ylim=c(0,1),col="purple3",ylab="",xlab="")  # plot OF and ER estimates
par(new=TRUE)
plot(res[,6],res[,14],type="l",lty=2,lwd=lwd.ci,xlim=c(0,1),ylim=c(0,1),col="purple3",ylab="",xlab="")  # plot OF and ER estimates
par(new=TRUE)
plot(res[,6],res[,15],type="l",lty=2,lwd=lwd.ci,xlim=c(0,1),ylim=c(0,1),col="purple3",ylab="",xlab="")  # plot OF and ER estimates

abline(a = 0, b = 1,lty=2)
legend(.1,.9,legend=c("Default Z-Curve", "NDIST Z-Curve"),lwd=2,col=c("forestgreen","purple3"))


# compute root mean square error for Z-Curve
rmse.err.1 = sqrt(mean((res[,6]-res[,7])^2))
rmse.err.1

# compute root mean square error for P-Curve
rmse.err.2 = sqrt(mean((res[,6]-res[,13])^2))
rmse.err.2

# compute directional bias
dir.err.1 = res[,7] - res[,6];summary(dir.err.1)
dir.err.2 = res[,13] - res[,6];summary(dir.err.2)

# confidence interval coverage 
tab.r1 = table(res[,8] < res[,6] & res[,9] > res[,6])
tab.r1 / sum(tab.r1)

tab.r2 = table(res[,14] < res[,6] & res[,15] > res[,6])
tab.r2 / sum(tab.r2)


########################################################################
### Evaluation of ERR estimates: z-curve (OF) VS z-curve (EXT free ZSD)
########################################################################

# columns
# -  5: true EDR
# - 10: res.1
# - 16: res.2 


round(res[,c(1:4,5,10,16)],3)

res = res[order(res[,5]),]

lwd.ci = .5

plot(res[,5],res[,10],xlim=c(0,1),ylim=c(0,1),pch=16,cex=1,col="forestgreen",ylab="Estimated ERR",
	xlab = "Simulated True EDR")  # plot OF and ER estimates
par(new=TRUE)
plot(res[,5],res[,11],type="l",lty=2,lwd=lwd.ci,xlim=c(0,1),ylim=c(0,1),col="forestgreen",ylab="",xlab="")  # plot OF and ER estimates
par(new=TRUE)
plot(res[,5],res[,12],type="l",lty=2,lwd=lwd.ci,xlim=c(0,1),ylim=c(0,1),col="forestgreen",ylab="",xlab="")  # plot OF and ER estimates

abline(a = 0, b = 1,lty=2)
legend(.1,.9,legend=c("Default Z-curve"),lwd=2,col=c("forestgreen"))


plot(res[,5],res[,16],pch=15,cex=1,xlim=c(0,1),ylim=c(0,1),col="purple3",ylab="",xlab="")  # plot OF and ER estimates
par(new=TRUE)
plot(res[,5],res[,17],type="l",lty=2,lwd=lwd.ci,xlim=c(0,1),ylim=c(0,1),col="purple3",ylab="",xlab="")  # plot OF and ER estimates
par(new=TRUE)
plot(res[,5],res[,18],type="l",lty=2,lwd=lwd.ci,xlim=c(0,1),ylim=c(0,1),col="purple3",ylab="",xlab="")  # plot OF and ER estimates

abline(a = 0, b = 1,lty=2)
legend(.5,.2,legend=c("Free SD Z-curve"),lwd=2,col=c("purple3"))


# compute root mean square error for Z-Curve
rmse.err.1 = sqrt(mean((res[,5]-res[,10])^2))
rmse.err.1

# compute root mean square error for P-Curve
rmse.err.2 = sqrt(mean((res[,5]-res[,16])^2))
rmse.err.2

# compute directional bias
dir.edr.1 = res[,10] - res[,5];summary(dir.edr.1)
dir.edr.2 = res[,16] - res[,5];summary(dir.edr.2)

# confidence interval coverage 
ci.r1 = as.numeric(res[,11] < res[,5] & res[,12] > res[,5])
tab.r1 = table(ci.r1)
tab.r1 / sum(tab.r1)

ci.r2 = as.numeric(res[,17] < res[,5] & res[,18] > res[,5])
tab.r2 = table(ci.r2)
tab.r2 / sum(tab.r2)


edr.bias = cbind(res[,1:4],dir.edr.1,dir.edr.2)
edr.bias = edr.bias[order(edr.bias[,1]),]
edr.bias

### z-curve 1SD bad
round(edr.bias[abs(edr.bias[,6]) > .50,],3)

### default z-curve CI does not include true value
round(res[ci.r1 == 0,c(1:4,5,10:12)],3)

#63 is the worst because the true value is low and implies a high FDR
#   however, the boundary of the CI is close 
#   follow up with EM method 

run.i = 63
sim.w.all = weights[run.i,];sim.w.all
sim.k = 10000

sim.components = length(sim.ncz)

sim.z = c()
i = 3
for (i in 1:sim.components) sim.z = c(sim.z,
	rnorm(sim.k*21*sim.w.all[i],sim.ncz[i],sim.zsd[i]) 
)

sim.z = abs(sim.z)
sim.z = sim.z[sim.z > z.crit]
sim.z = sim.z[order(runif(length(sim.z),0,1))]
sim.z = sim.z[1:sim.k]
table(sim.z > z.crit)
summary(sim.z)

source(sim.z)
ymax = 1
ymin = -.05
Est.Method = "EM"
boot.iter = 500
res.1 = Zing(sim.z);res.1





