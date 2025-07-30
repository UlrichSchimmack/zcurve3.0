

###############################################################################
# Chapter 5 – Comparing z-curve and p-curve 
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
pcurve <- "Pcurve.Function.R"
source(pcurve)


# Alternatively, load z-curve function directly from GitHub
zcurve3 <- "https://raw.githubusercontent.com/UlrichSchimmack/zcurve3.0/refs/heads/main/Zing.25.07.11.test.R"
source(zcurve3)
pcurve <- "https://github.com/UlrichSchimmack/zcurve3.0/raw/refs/heads/main/Pcurve.Function.R"
source(pcurve)


### 1.0 – ILLUSTRATION WITH OSC REP PROJECT

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


### Convert two-sided p-values to z-values
z.ori <- qnorm(1 - pval$ori / 2)


source(zcurve3)
ymax = 1
ymin = -.05
hist.bar.width = .4
Est.Method = "EM"
boot.iter = 500
Zing(z.ori)


pcurve.input = paste0("z = ",z.ori)
pcurve_app(pcurve.input,SHOW.PLOT=TRUE)



### 2.0 – SIMULATION 

sim.k = 10000
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

### Illustrate mixture of (folded) normal distributions 

x1 = seq(0,6,.01)
x2 = seq(0,6,.01)*-1
d0 = dnorm(x1,0) + dnorm(x2,0)
d2 = dnorm(x1,2) + dnorm(x2,2)
d4 = dnorm(x1,4) + dnorm(x2,4)

w = weights[53,];w
d.mix = w[1]*d0 + w[2]*d2 + w[3]*d4

plot(x1,d0,xlim=c(0,6),ylim=c(0,1),type="l",lwd=2,col="purple3",
	xlab="",ylab="",lty=2)
par(new=TRUE)
plot(x1,d2,xlim=c(0,6),ylim=c(0,1),type="l",lwd=2,col="forestgreen",
	xlab="",ylab="",lty=2)
par(new=TRUE)
plot(x1,d4,xlim=c(0,6),ylim=c(0,1),type="l",lwd=2,col="royalblue",
	xlab="",ylab="",lty=2)
par(new=TRUE)
plot(x1,d.mix,xlim=c(0,6),ylim=c(0,1),type="l",lwd=4,
	col="black",lty=1,ylab="Density",
  main = paste0("Weights: w0 = ",w[1],"  w2 = ",w[2],"  w4 = ",w[3]) )





# run the simulation and collect the z-curve estimates for old-fashioned (OF)
# and slower EM methods of estimation

res = c()

run.i = 64
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

print(sim.err)
if(sim.err > 1) print("FUCK")

components = length(sim.ncz)

sim.z = c()
i = 3
for (i in 1:components) sim.z = c(sim.z,
	rnorm(sim.k*21*sim.w.all[i],sim.ncz[i],sim.zsd[i]) 
)

sim.z = abs(sim.z)
sim.z = sim.z[sim.z > z.crit]
sim.z = sim.z[order(runif(length(sim.z),0,1))]
sim.z = sim.z[1:sim.k]
table(sim.z > z.crit)
summary(sim.z)

# Run basic z-curve model with "EM" method
source(zcurve3)    # load zcurve3 every time to get default settings
Title = paste0("True Average Power of Significant Results (ERR): ",round(sim.err*100));Title
ymax = 1
ymin = -.05
#Est.Method = "EM"
boot.iter = 500
res.1 = Zing(sim.z)$res[2,]
res.1

# Run Extended z-curve model with 1 component and free SD
# A one-component model with fixed SD is equivalent to p-curve
# The free SD parameter allows z-curve to model heterogeneity, ...
# assuming a normal distribution of non-central z-values
Est.Method = "EXT"    # specify extended method with free parameters (zcurve.3.0)
ncz = 2				# specify only one component with any starting value
zsd = 5				# specify SD with maximum value to allow for free estimation
res.2 = Zing(sim.z)$res[2,]

pcurve.input = paste0("z = ",sim.z)
res.3 = pcurve_app(pcurve.input,SHOW.PLOT=TRUE)
res.3 

res.run = c(run.i,weights[run.i,],sim.edr,sim.err,res.1,res.2,res.3)
res.run

res = rbind(res,res.run)

write.table(res,"sim10k.ZcurveVSPcurve.dat") # write results each trial
                                   # can resume if stopped 
                                   # by loading the completed results

} # End of for loops

dim(res)  # 66 rows and 18 columns

#write the completed results / overwrites the data from the for loop
#write.table(res,"sim.ZcurveVSPcurve.dat")

########################################################
### GET STORED RESULTS
########################################################

# load the saved data
setwd("C:/Users/ulric/Dropbox/PHPCurves/DOCS/z-curve/Tutorial")
check = read.table("sim10k.ZcurveVSPcurve.dat",row.names=NULL)
check = data.frame(check)
check = check[,2:dim(check)[2]]
dim(check)

# careful, if you run this code any results in the res file will be replaced 
# with the saved data. 
#res = check  # use only if sure to replace res; run analysis with "res" matrix 

dim(res)


#########################################################
### Evaluation of ERR estimates: z-curve (EM) VS p-curve
#########################################################

# columns
# -  6: true ERR
# -  7: res.1
# - 10: res.2 # not used
# - 13: res.3

round(res[,c(1:4,6,7,10,13)],3)

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
legend(.1,.9,legend=c("Z-curve","P-curve"),lwd=2,col=c("forestgreen","purple3"))



# compute root mean square error for Z-Curve
rmse.err.1 = sqrt(mean((res[,6]-res[,7])^2))
rmse.err.1

# compute root mean square error for P-Curve
rmse.err.3 = sqrt(mean((res[,6]-res[,13])^2))
rmse.err.3

# compute directional bias
dir.err.1 = res[,7] - res[,6];summary(dir.err.1)
dir.err.3 = res[,13] - res[,6];summary(dir.err.3)


# confidence interval coverage 
tab.z = table(res[,8] < res[,6] & res[,9] > res[,6])
tab.z / sum(tab.z)

tab.p = table(res[,14] < res[,6] & res[,15] > res[,6])
tab.p / sum(tab.p)





### Not relevant for blog post


#####################################
### Evaluation of z-curve EXT 
#####################################

round(res[,c(1:4,6,7,10)],3)

plot(res[,6],res[,7],xlim=c(0,1),ylim=c(0,1),pch=16,chx=1.2,col="forestgreen",ylab="Estimated ERR",
	xlab = "Simulated True EDR")  # plot OF and ER estimates
par(new=TRUE)
plot(res[,6],res[,10],pch=15,xlim=c(0,1),ylim=c(0,1),col="purple3",ylab="",xlab="")  # plot OF and ER estimates
abline(a = 0, b = 1,lty=2)
legend(.1,.9,legend=c("Z-curve EM","Z-curve EXT"),lwd=2,col=c("forestgreen","purple3"))


# compute root mean square error for Z-Curve EM
rmse.err.1 = sqrt(mean((res[,6]-res[,7])^2))
rmse.err.1

# compute root mean square error for Z-Curve SD > 1
rmse.err.2 = sqrt(mean((res[,6]-res[,10])^2))
rmse.err.2

# compute directional bias
dir.err.1 = res[,7] - res[,6];summary(dir.err.1)
dir.err.2 = res[,10] - res[,6];summary(dir.err.2)


# regress directional bias on weights

reg.err.2 = summary(lm(dir.err.2 ~ res[,2] * res[,3] * res[,4] ))
reg.err.2

reg.err.2 = summary(lm(dir.err.2 ~ res[,2] * res[,3] * res[,4] + ncz.var ))
reg.err.2


err.bias = cbind(res[,1:4],ncz.var,dir.err.1,dir.err.2,dir.err.3)
round(err.bias,3)

### p-curve does better than z-curve
round(err.bias[abs(err.bias[,8]) < abs(err.bias[,6]),],3)

### z-curve does better than p-curve
round(err.bias[abs(err.bias[,6]) < abs(err.bias[,8]),],3)


### z-curve 1SD better than z-curve EM
round(err.bias[abs(err.bias[,7]) < abs(err.bias[,6]),],3)


### z-curve 1SD bad
round(err.bias[abs(err.bias[,7]) > .15,],3)




#65 65 0.9 0.1 0.0   0.316    -0.060    -0.223    -0.087

