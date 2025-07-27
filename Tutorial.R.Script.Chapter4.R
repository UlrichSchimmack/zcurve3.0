

###############################################################################
# Chapter 4 – Conduct Simulation Studies 
###############################################################################

### 0.0 – SETUP

# disable scientific notation
options(scipen = 999)

# Set work directory
setwd("C:/Users/ulric/Dropbox/PHPCurves/DOCS/z-curve/Tutorial")

# Install required packages (run once only)
install.packages(c("zcurve", "KernSmooth", "stringr", "parallel"))

# Load z-curve function from a local file
setwd("C:/Users/ulric/Dropbox/PHPCurves/DOCS/z-curve/Tutorial")
zcurve3 <- "Zing.25.07.11.test.R"
source(zcurve3)

# Alternatively, load z-curve function directly from GitHub
zcurve3 <- "https://raw.githubusercontent.com/UlrichSchimmack/zcurve3.0/refs/heads/main/Zing.25.07.11.test.R"
source(zcurve3)

### 1.0 – SIMULATION 

sim.k = 50000
sim.ncz = c(0,2,4) # change or add more components
sim.zsd = c(1,1,1) #
                   # has to stay at 1 to compute true power            
                   # other values require estimation with large sample simulations and no bias

# check that there is an equal number of means and sds of components
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

run.i = 52
weights[run.i,]

b = 1
e = 53
e = nrow(weights);e

for (run.i in b:e ) {

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
#sim.z = sim.z[sim.z > z.crit]
sim.z = sim.z[order(runif(length(sim.z),0,1))]
sim.z = sim.z[1:sim.k]
table(sim.z > z.crit)
summary(sim.z)

Title = paste0(round(sim.edr*100),"   ",round(sim.err*100));Title

# Run basic z-curve model
Est.Method = "OF"
res.1 = Zing(sim.z)

Est.Method = "EM"
res.2 = Zing(sim.z)

res.run = c(run.i,weights[run.i,],sim.edr,sim.err,
	res.1$res,
	res.2$res
	)
res.run

res = rbind(res,res.run)

write.table(res,"sim.OFvsEM.dat") # write results each trial
                                   # can resume if stopped 
                                   # by loading the completed results

} # End of for loops

dim(res)  # 66 rows and 15 columns

#write the completed results / overwrites the data from the for loop
#write.table(res,"sim.OFvsEM.dat")


########################################################
### GET STORED RESULTS
########################################################

# load the saved data
check = read.table("sim.OFvsEM.dat",row.names=NULL)
#write.table(check,"sim.edr.est.backup.dat")
check = data.frame(check)
check = check[,2:dim(check)[2]]
dim(check)

# careful, if you run this code any results in the res file will be replaced 
# with the saved data. 
#res = check  # use only if sure to replace res 

dim(res)


#####################################
### Evaluation of ERR estimates
#####################################

# columns
# -  6: true ERR
# -  8: res.1
# - 13: res.2

round(res[,c(1:4,6,8,13)],3)

plot(res[,6],res[,8],pch=15,col="purple3",ylab="",xlab="")  # plot OF and ER estimates
par(new=TRUE)
plot(res[,6],res[,13],pch=16,chx=1.2,col="forestgreen",ylab="Estimated ERR",
	xlab = "Simulated True EDR")  # plot OF and ER estimates
abline(a = 0, b = 1,lty=2)

# compute root mean square error for bw = .05
rmse.err.1 = sqrt(mean((res[,6]-res[,8])^2))
rmse.err.1

# compute root mean square error for bw = .10
rmse.err.2 = sqrt(mean((res[,6]-res[,13])^2))
rmse.err.2

# compute directional bias
dir.err.1 = res[,8] - res[,6];summary(dir.err.1)
dir.err.2 = res[,13] - res[,6];summary(dir.err.2)

# regress directional bias on weights
reg.err.1 = summary(lm(dir.err.1 ~ res[,2] * res[,3] * res[,4] ))
reg.err.1
reg.err.2 = summary(lm(dir.err.2 ~ res[,2] * res[,3] * res[,4] ))
reg.err.2


# check individual conditions
err.bias = cbind(res[,1:4],dir.err.1,dir.err.2)
round(err.bias[abs(err.bias[,5])>.02,],3)
round(err.bias[abs(err.bias[,6])>.02,],3)


### practically no bias
### no systematic effect of z-value distribution (shape of z-curve)
### very similar results for fast density and slower EM method


#####################################
### Evaluation of EDR estimates
#####################################

# repeat for EDR estimates 

round(res[,c(1:4,5,9,14)],3)


plot(res[,5],res[,9],pch=15,col="purple3",ylab="",xlab="")  # plot OF and ER estimates
par(new=TRUE)
plot(res[,5],res[,14],pch=16,chx=1.2,col="forestgreen",ylab="Estimated ERR",
	xlab = "Simulated True EDR")  # plot OF and ER estimates
abline(a = 0, b = 1,lty=2)

rmse.edr.1 = sqrt(mean((res[,5]-res[,9])^2))
rmse.edr.1

rmse.edr.2 = sqrt(mean((res[,5]-res[,14])^2))
rmse.edr.2

dir.edr.1 = res[,9] - res[,5];summary(dir.edr.1)
dir.edr.2 = res[,14] - res[,5];summary(dir.edr.2)

reg.edr.1 = summary(lm(dir.edr.1 ~ res[,2] * res[,3] * res[,4] ))
reg.edr.1

reg.edr.2 = summary(lm(dir.edr.2 ~ res[,2] * res[,3] * res[,4] ))
reg.edr.2

edr.bias = cbind(res[,1:4],dir.edr.1,dir.edr.2)
edr.bias
max.bias = .10
round(edr.bias[abs(edr.bias[,5]) > max.bias,],3)
round(edr.bias[abs(edr.bias[,6]) > max.bias,],3)

table(abs(edr.bias[,5]) + .05 < abs(edr.bias[,6]) )
table(abs(edr.bias[,6]) + .05 < abs(edr.bias[,5]) )

max.bias = .05
round(edr.bias[abs(edr.bias[,5]) > max.bias,],3)

min(edr.bias[,5])

run.i = 41 # biggest upward bias
run.i = 20 # biggest downward bias
sim.w.all = weights[run.i,];sim.w.all

set.seed(20250726)
sim.k = 200000
z.crit = qnorm(.975)
sim.ncz
sim.zsd

sim.z = c()
i = 3
for (i in 1:length(sim.ncz))  sim.z = c(sim.z,
	rnorm(sim.k*sim.w.all[i],sim.ncz[i],sim.zsd[i]) 
)

sim.z = abs(sim.z)
sim.z = sim.z[order(runif(length(sim.z),0,1))]
table(sim.z > z.crit)
summary(sim.z)

source(zcurve3)
Title = toString(sim.w.all)
hist.bar.width = .1
Zing(sim.z)

#Est.Method = "EM"; Zing(sim.z) # Warning: very slow with large k!!!
Zing(sim.z)



