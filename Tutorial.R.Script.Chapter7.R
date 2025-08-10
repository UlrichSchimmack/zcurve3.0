


###############################################################################
# Chapter 7– Comparing z-curve and BACON 
###############################################################################

### 0.0 – SETUP

# disable scientific notation
options(scipen = 999)

# Set work directory
setwd("C:/Users/ulric/Dropbox/PHPCurves/DOCS/z-curve/Tutorial")

# Install required packages (run once only)
#install.packages(c("zcurve", "KernSmooth", "stringr", "parallel"))
#install.packages("bacon","mixtools")

library(bacon) # load bacon package

# Load z-curve and P-curve functions from a local file
setwd("C:/Users/ulric/Dropbox/PHPCurves/DOCS/z-curve/Tutorial")
zcurve3 <- "Zing.25.07.11.test.R"
source(zcurve3)

# Alternatively, load z-curve function directly from GitHub
zcurve3b <- "https://raw.githubusercontent.com/UlrichSchimmack/zcurve3.0/refs/heads/main/Zing.25.07.11.test.R"
source(zcurve3b)

### 1.0 – ILLUSTRATION WITH EXAMPLES

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

### Convert two-sided p-values to z-values
z.rep <- qnorm(1 - pval$rep / 2)


### non.significant / significant / extreme (z > 6)
z.rep[z.rep == Inf] = 10
summary(z.rep)

table(cut(z.rep,c(0,1.96,6,Inf)))
### only 5 extreme z-values

### Observed Discovery Rate: 33% significant results
tab = table(z.rep > 1.96);tab/sum(tab)


### 1.1 Run default z-curve model
source(zcurve3)
ymax <- 1.0
ymin <- -0.05
hist.bar.width <- 0.4
Int.Beg = 0
res.1 = Zing(z.rep)

### 1.2 Run model with three free components
source(zcurve3)
ymax <- 1.0
ymin <- -0.05
hist.bar.width <- 0.4
Est.Method = "EXT"
NCZ.FIXED = FALSE
ncz = c(0,2,4)
zsds = c(1,1,1)
Int.Beg = 0
res.2 = Zing(z.rep)
res.2

### 1.3 Run model with two components with free means and SDs 
source(zcurve3)
ymax <- 1.0
ymin <- -0.05
hist.bar.width <- 0.4
Est.Method = "EXT"
NCZ.FIXED = FALSE
ZSDS.FIXED = FALSE
ncz = c(0,2)
zsds = c(3,3)
Int.Beg = 0
res.3 = Zing(z.rep)
res.3


### run bacon model 
bc.res = bacon(teststatistics=z.rep)

### check the results 
bc.res@estimates

### visualize with zcurve, using fixed parameters from bacon
source(zcurve3)
ymax <- 1.0
ymin <- -0.05
hist.bar.width <- 0.4
Est.Method = "EXT"
NCZ.FIXED = TRUE
ZSDS.FIXED = TRUE
W.FIXED = TRUE
w.fix = bc.res@estimates[1:2]/sum(bc.res@estimates[1:2])
ncz = bc.res@estimates[4:5]
zsds = bc.res@estimates[7:8]
Augment.Factor = 5
Int.Beg = 0
bw.draw = .5
res.4 = Zing(z.rep)
res.4

### alternative method to compute EDR and ERR from bacon results


### visualize with zcurve, using fixed parameters from bacon
source(zcurve3)
ymax <- 1.0
ymin <- -0.05
hist.bar.width <- 0.4
Est.Method = "EXT"
NCZ.FIXED = TRUE
ZSDS.FIXED = TRUE
W.FIXED = TRUE
w.fix = bc.res@estimates[1:2]/sum(bc.res@estimates[1:2])
ncz = bc.res@estimates[4:5]
zsds = bc.res@estimates[7:8]
Augment.Factor = 5
Int.Beg = 0
Int.End = 11
bw.draw = .5
res.5 = Zing(z.rep)
res.5

### alternative method to compute EDR and ERR from bacon results








zx = seq(0,10,.01)
pow.zx.dir = pnorm(abs(zx),1.96) 
pow.zx.sign.error = + pnorm(-1.96,abs(zx))

est.cw.all = bc.res.1@estimates[1:3];est.cw.all
est.ncz = bc.res.1@estimates[4:6];est.ncz
est.zsds = bc.res.1@estimates[7:9];est.zsds
est.zsds[est.zsds < 1] = 1;est.zsds          ### sd < 1 need to be set to 1
est.components = length(est.ncz)
est.zncz.sd = sqrt(est.zsds^2-1);est.zncz.sd

use = which(est.ncz > 0)

i = 1
est.wd.all = c()
for (i in use) {
	if (est.zncz.sd[i] == 0) {
		wd = rep(0,length(zx))
		wd[which(round(zx,2) == round(est.ncz[i],2))] = 1
	} else {
		wd = dnorm(zx,est.ncz[i],est.zncz.sd[i])
	}
	sum(wd)
	wd = wd/sum(wd)
	wd = wd*est.cw.all[i]
	sum(wd)
	est.wd.all = rbind(est.wd.all,wd)
}

dim(est.wd.all)
table(is.na(est.wd.all))

est.wd.all = colSums(est.wd.all)
sum(est.wd.all)

#plot(round(zx,2),round(est.wd.all,2),type="l")

est.wd.sig = est.wd.all*(pow.zx.dir+pow.zx.sign.error)
est.wd.sig = est.wd.sig/sum(est.wd.sig)

bacon1.edr = sum((pow.zx.dir+pow.zx.sign.error)*est.wd.all);bacon1.edr
bacon1.err = sum(pow.zx.dir*est.wd.sig)/sum(est.wd.sig);bacon1.err

res.5 = c(bacon1.err,bacon1.edr)
res.5

### even lower estimates because extremes are not fitted well 

restimates = rbind(
res.1$res[2:3],
res.2$res[2:3],
res.3$res[2:3],
res.4$res[2:3],
res.5
)

round(restimates,2)

### Conclusion: 
### The true EDR is known because there is no selection bias in the OSC:Rep Data
### The true EDR is 33% 
### Bacon underestimates the EDR



### 2.0 – SIMULATION STUDY 

### START 

### 2.0 – SIMULATION 

sim.runs = 1       # repeat each simulation 100 times
sim.k.sig = 10000  # simulate 100 significant results per "study"
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

weights <- weights[rep(1:nrow(weights), times = sim.runs), ]
dim(weights)

check = weights[,1]*100 + weights[,2]*10 + weights[,3]
table(check)

# run the simulation and collect the z-curve estimates for 
# z-curve and p-curve 

res = c()

run.i = 34
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

set.seed(run.i)
sim.z = c()
i = 3
for (i in 1:sim.components) sim.z = c(sim.z,
	rnorm(sim.k*sim.w.all[i],sim.ncz[i],sim.zsd[i]) 
)

sim.z = sim.z[order(runif(length(sim.z),0,1))]
table(sim.z > z.crit)
summary(sim.z)
length(sim.z)
#hist(sim.z)

z.val.input = sim.z


################ Get RESULTS #########################

# Run basic z-curve model with "OF" method
source(zcurve3)    # load zcurve3 every time to get default settings
Title = paste(round(sim.edr*100),round(sim.err*100));Title
ymax = 1
ymin = -.05
Int.Beg = 0
res.1 = Zing(abs(z.val.input));res.1
res.1 = c(res.1$res[3],res.1$res[2])
res.1

# Run with 3 components: free ncz and zsds = 1
source(zcurve3)    # load zcurve3 every time to get default settings
Title = paste(round(sim.edr*100),round(sim.err*100));Title
Int.Beg = 0
Est.Method = "EXT"
NCZ.FIXED = FALSE
ncz = c(0,2,4)
zsds = c(1,1,1)
res.2 = Zing(abs(z.val.input));res.3
res.2 = c(res.2$res[3],res.2$res[2])
res.2

# Run with 2 components: free ncz and free zsds
source(zcurve3)    # load zcurve3 every time to get default settings
Title = paste(round(sim.edr*100),round(sim.err*100));Title
Est.Method = "EXT"
NCZ.FIXED = FALSE
ZSD.FIXED = FALSE
ncz = c(0,4)
zsds = c(2,2)
Int.Beg = 0
res.3 = Zing(abs(z.val.input));res.4
res.3 = c(res.3$res[3],res.3$res[2])
res.3


# Run BACON Analysis with default priors

bc.res = bacon(teststatistics=z.val.input)

# Get estimates with parameters as fixed input to z-curve 

source(zcurve3)
Title = paste(round(sim.edr*100),round(sim.err*100));Title
ymax <- 1.0
ymin <- -0.05
hist.bar.width <- 0.4
Est.Method = "EXT"
NCZ.FIXED = TRUE
ZSDS.FIXED = TRUE
W.FIXED = TRUE
w.fix = bc.res@estimates[1:2]/sum(bc.res@estimates[1:2])
ncz = bc.res@estimates[4:5]
zsds = bc.res@estimates[7:8]
Augment.Factor = 5
Int.Beg = 0
bw.draw = .5
res.4 = Zing(abs(z.val.input))
res.4 = c(res.4$res[3],res.4$res[2])
res.4

### use full distribution 

zx = seq(0,10,.01)
pow.zx.dir = pnorm(abs(zx),1.96) 
pow.zx.sign.error = + pnorm(-1.96,abs(zx))

est.cw.all = bc.res@estimates[1:3];est.cw.all
est.ncz = bc.res@estimates[4:6];est.ncz
est.zsds = bc.res@estimates[7:9];est.zsds
est.zsds[est.zsds < 1] = 1;est.zsds
est.components = length(est.ncz)
est.zncz.sd = sqrt(est.zsds^2-1);est.zncz.sd

use = which(est.ncz > 0)

i = 1
est.wd.all = c()
for (i in use) {
	if (est.zncz.sd[i] == 0) {
		wd = rep(0,length(zx))
		wd[which(round(zx,2) == round(est.ncz[i],2))] = 1
	} else {
		wd = dnorm(zx,est.ncz[i],est.zncz.sd[i])
	}
	sum(wd)
	wd = wd/sum(wd)
	wd = wd*est.cw.all[i]
	sum(wd)
	est.wd.all = rbind(est.wd.all,wd)
}

dim(est.wd.all)
table(is.na(est.wd.all))

est.wd.all = colSums(est.wd.all)
sum(est.wd.all)

est.wd.sig = est.wd.all*(pow.zx.dir+pow.zx.sign.error)
est.wd.sig = est.wd.sig/sum(est.wd.sig)

bacon1.edr = sum((pow.zx.dir+pow.zx.sign.error)*est.wd.all);bacon1.edr
bacon1.err = sum(pow.zx.dir*est.wd.sig)/sum(est.wd.sig);bacon1.err
res.5 = c(bacon1.edr,bacon1.err)
res.5

res.run = rbind(
res.1,
res.2,
res.3,
res.4,
res.5
)

res.run = rbind(c(sim.edr,sim.err),res.run)
print(round(res.run,2))

res.run = c(run.i,weights[run.i,],
  c(res.run))

res = rbind(res,res.run)

write.table(res,"sim10k.ZcurveVSbacon.dat") # write results each trial
                                   # can resume if stopped 
                                   # by loading the completed results

} # End of for loops

dim(res)  # 66 rows and 18 columns





#write the completed results / overwrites the data from the for loop
#write.table(res,"sim10k.ZcurveVSbacon.dat")


########################################################
### GET STORED RESULTS
########################################################

# load the saved data
setwd("C:/Users/ulric/Dropbox/PHPCurves/DOCS/z-curve/Tutorial")
check = read.table("sim10k.ZcurveVSbacon.dat",row.names=NULL)
check = data.frame(check)
check = check[,2:dim(check)[2]]
dim(check)
check[,1]

# careful, if you run this code any results in the res file will be replaced 
# with the saved data. 
# table(round(res[,c(1:14)],5) == round(check[,c(1:14)],5) )
# round(cbind(res[,7],check[,7]),5)

#res = check  # use only if sure to replace res; run analysis with "res" matrix 

dim(res)

summary(res)


########################################################################
### Evaluation of EDR estimates: z-curve (OF) VS z-curve (EXT free ZSD)
########################################################################

# columns
# -  5: true EDR
# -  6: res.1
# -  7: res.2 
# -  8: res.3
# -  9: res.4
# - 10: res.5 - bacon

round(res[,c(1:4,5:10)],3)

cor(res[,5:10])

rmse.edr.1 = sqrt(mean((res[,5]-res[,6])^2));rmse.edr.1
rmse.edr.2 = sqrt(mean((res[,5]-res[,7])^2));rmse.edr.2
rmse.edr.3 = sqrt(mean((res[,5]-res[,8])^2,na.rm=TRUE));rmse.edr.3
rmse.edr.4 = sqrt(mean((res[,5]-res[,9])^2,na.rm=TRUE));rmse.edr.4
rmse.edr.5 = sqrt(mean((res[,5]-res[,10])^2,na.rm=TRUE));rmse.edr.5

rmse.edr.1
rmse.edr.2
rmse.edr.3
rmse.edr.4
rmse.edr.5


res = res[order(res[,5]),]

lwd.ci = .5

plot(res[,5],res[,10],pch=15,cex=1,xlim=c(0,1),ylim=c(0,1),col="purple3",ylab="",xlab="")  # plot OF and ER estimates
par(new=TRUE)
plot(res[,5],res[,7],xlim=c(0,1),ylim=c(0,1),pch=16,cex=1,col="forestgreen",ylab="Estimated ERR",
	xlab = "Simulated True EDR")  # plot OF and ER estimates
abline(a = 0, b = 1,lty=2)
legend(.1,.9,legend=c("Default Z-Curve EM", "Full Bacon"),lwd=2,col=c("forestgreen","purple3"))


########################################################################
### Evaluation of ERR estimates: z-curve (OF) VS z-curve (EXT free ZSD)
########################################################################

# columns
# - 11: true ERR
# - 12: res.1
# - 13: res.2 
# - 14: res.3
# - 15: res.4
# - 16: res.5

rmse.err.1 = sqrt(mean((res[,11]-res[,12])^2));rmse.err.1
rmse.err.2 = sqrt(mean((res[,11]-res[,13])^2));rmse.err.2
rmse.err.3 = sqrt(mean((res[,11]-res[,14])^2));rmse.err.3
rmse.err.4 = sqrt(mean((res[,11]-res[,15])^2));rmse.err.4
rmse.err.5 = sqrt(mean((res[,11]-res[,16])^2));rmse.err.5

rmse.err.1
rmse.err.2
rmse.err.3
rmse.err.4
rmse.err.5

	
res = res[order(res[,11]),]

lwd.ci = .5

plot(res[,11],res[,16],pch=15,cex=1,xlim=c(0,1),ylim=c(0,1),col="purple3",ylab="",xlab="")  # plot OF and ER estimates
par(new=TRUE)
plot(res[,11],res[,12],xlim=c(0,1),ylim=c(0,1),pch=16,cex=1,col="forestgreen",ylab="Estimated ERR",
	xlab = "Simulated True EDR")  # plot OF and ER estimates
abline(a = 0, b = 1,lty=2)
legend(.1,.9,legend=c("Default Z-Curve", "Full Bacon"),lwd=2,col=c("forestgreen","purple3"))




