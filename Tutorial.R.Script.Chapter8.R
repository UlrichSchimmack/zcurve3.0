


###############################################################################
# Chapter 8 – Transformation of Test-Statistics into z-values
###############################################################################

### 0.0 – SETUP

# disable scientific notation
options(scipen = 999)

# Set work directory
setwd("C:/Users/ulric/Dropbox/PHPCurves/DOCS/z-curve/Tutorial")

# Install required packages (run once only)
#install.packages(c("zcurve", "KernSmooth", "stringr", "parallel","pwr"))

# Load z-curve and P-curve functions from a local file
setwd("C:/Users/ulric/Dropbox/PHPCurves/DOCS/z-curve/Tutorial")
zcurve3 <- "Zing.25.07.11.test.R"
source(zcurve3)

# Alternatively, load z-curve function directly from GitHub
zcurve3b <- "https://raw.githubusercontent.com/UlrichSchimmack/zcurve3.0/refs/heads/main/Zing.25.07.11.test.R"
source(zcurve3b)

### 1.0 – ILLUSTRATION WITH EXAMPLES

sim.k = 100000
sim.p.sig = .50
sim.ncz = qnorm(sim.p.sig,qnorm(.975));sim.ncz
N = 30
df = N-2
se = 2/sqrt(N)
nct = pwr.t.test(n = N/2, 
	d=NULL, sig.level = .05, power = sim.p.sig,
	type = "two.sample", alternative = "two.sided")$d/se;nct

z.val = abs(rnorm(sim.k,sim.ncz))
t.val = abs(rt(sim.k,df,nct))

summary(z.val)
summary(t.val)

##########################

col1a = rgb(0, 100, 100, 255, max = 255)
col2a = rgb(100,0,100, max = 255)

col1b = rgb(0, 100, 100, 255, max = 255,alpha = 50)
col2b = rgb(100,0,100, max = 255, alpha = 50)


x.max = 6
graphics.off()
Title = paste0("Sample Size, N = ",N,"    P(sig) = ",round(sim.p.sig*100),"%")
hist(main=Title,z.val,freq=FALSE,ylim=c(0,.5),xlim=c(0,x.max),col=col1b,xlab="z-value")
par(new=TRUE)
hist(main="",t.val,freq=FALSE,ylim=c(0,.5),xlim=c(0,x.max),col=col2b,xlab="")
legend(0.5,.5,c("z-values","t-values"),col=c(col1a,col2a),lwd=2)
par(new=TRUE)
curve(dnorm(x,sim.ncz)+dnorm(-x,sim.ncz),0,7,ylim=c(0,.5),xlim=c(0,x.max),col=col1a,ylab="",xlab="",lwd=3)
par(new=TRUE)
curve(dt(x,df,nct)+dt(-x,df,nct),0,7,ylim=c(0,.5),xlim=c(0,x.max),col=col2a,lty=2,ylab="",xlab="",lwd=3)

graphics.off()

source(zcurve3)
Int.Beg = 0
Title = paste0("z-values, Power = ",round(sim.p.sig*100),"%")
Zing(z.val)

##########################

source(zcurve3)
Title = paste0("t-values, ",sim.df," df, Power = ",round(sim.p.sig*100),"%")
Zing(t.val)
Int.Beg = 0
Zing(t.val)

### t2p2z transform t-values to p-values and p-values to z-values 
t2p2z.val = qnorm(pt(t.val,sim.df))  # equivalent to qnorm(1-p/2) for two.tailed p-values
summary(t2p2z.val)

source(zcurve3)
Title = paste0("t->p->z values, ",sim.df," df, Power = ",round(sim.p.sig*100),"%")
Zing(t2p2z.val)
Int.Beg = 0
Zing(t2p2z.val)

##########################

### t2psig2z transform t-values to p.sig and p.sig to z-values 
median.t.sig = median(t.val[t.val > qt(.975,sim.df)]);median.t.sig
use.nct = 0
use.ncz = 0
if (median.t.sig > 2.5) {
	use.ncz = qnorm(.5,qnorm(.975));use.ncz
	use.nct = pwr.t.test(n = N/2, 
	d=NULL, sig.level = .05, power = .5,
	type = "two.sample", alternative = "two.sided")$d/se;use.nct
}
if (median.t.sig > 3) {
	use.ncz = qnorm(.75,qnorm(.975));use.ncz
	use.nct = pwr.t.test(n = N/2, 
	d=NULL, sig.level = .05, power = .75,
	type = "two.sample", alternative = "two.sided")$d/se;use.nct
}
use.nct
use.ncz
t2psig2z.val = qnorm(pt(t.val,sim.df,use.nct),use.ncz)
summary(t2psig2z.val)

source(zcurve3)
Title = paste0("t->p.sig->z values, ",sim.df," df, Power = ",round(sim.pow*100),"%")
Zing(t2psigz.val)
Int.Beg = 0
Zing(t2psig2z.val)

##########################

source(zcurve3)
Title = paste0("t-curve, ",sim.df," df, Power = ",round(sim.pow*100),"%")
Est.Method = "DF"
Zing(t.val,df=df)
Int.Beg = 0
Zing(t.val,df)



###################################################


### 2.0 – SIMULATION STUDY 

### START 

### 2.0 – SIMULATION 

sim.runs = 1       # repeat each simulation 100 times
sim.k.sig = 100000 # simulate 100 significant results per "study"
sim.ncz = c(0,2,4) # change or add more components
sim.zsd = c(1,1,1) #
                   # has to stay at 1 to compute true power            
                   # other values require estimation with large sample simulations and no bias

N = 30
sim.df = N-2
sim.se = 2/sqrt(N);se

sim.pow.2 = pnorm(2,1.96);sim.pow.2
sim.pow.4 = pnorm(4,1.96);sim.pow.4

sim.nct.2 <- pwr.t.test(n = N/2, sig.level = .025,
                         power = sim.pow.2,
                         type = "two.sample",
                         alternative = "greater")$d/sim.se;sim.nct.2
sim.nct.4 <- pwr.t.test(n = N/2, sig.level = .025,
                         power = sim.pow.4,
                         type = "two.sample",
                         alternative = "greater")$d/sim.se;sim.nct.4

sim.nct = c(0,sim.nct.2,sim.nct.4)

sim.ncz
sim.nct


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

cbind(1:66,weights)

run.i = 2
weights[run.i,]

b = 1                 # begin 
e = 5                 # end for testing
e = nrow(weights);e   # end for all runs


###

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
sim.t = c()

sim.i = 1
for (sim.i in 1:sim.components) {
	sim.z = c(sim.z,
	rnorm(sim.k.sig*21*sim.w.all[sim.i],sim.ncz[sim.i],sim.zsd[sim.i]) )
	sim.t = c(sim.t,
	rt(sim.k.sig*21*sim.w.all[sim.i],sim.df,sim.nct[sim.i]) )
}

check.z = sim.z
check.t = sim.t

table(check.z == sim.z)

length(sim.z)
length(sim.t)

sim.z = abs(sim.z)
sim.t = abs(sim.t)

sim.z = sim.z[order(runif(length(sim.z),0,1))]
sim.t = sim.t[order(runif(length(sim.t),0,1))]

tab = table(sim.z > z.crit);tab/sum(tab)
tab = table(sim.t > qt(.975,sim.df));tab/sum(tab)

sim.z = sim.z[sim.z > qnorm(.975)]
sim.t = sim.t[sim.t > qt(.975,sim.df)]

length(sim.z)
length(sim.t)

sim.z = sim.z[1:sim.k.sig]
sim.t = sim.t[1:sim.k.sig]

length(sim.z)
length(sim.t)

tab = table(sim.z > z.crit);tab / sum(tab)
tab = table(sim.t > qt(.975,sim.df)); tab/sum(tab)

summary(sim.z)
summary(sim.t)

########################

source(zcurve3)
Title = paste0("z-values, ",round(sim.edr*100),"  ",round(sim.err*100))
bw.est = .02
res.1 = Zing(sim.z);res.1
res.1 = res.1$res[2:3];res.1

###

source(zcurve3)
Title = paste0("t-values, ",round(sim.edr*100),"  ",round(sim.err*100))
bw.est = .02
res.2 = Zing(sim.t);res.2
res.2 = res.2$res[2:3];res.2

###

# t2p2z transform t-values to p-values and p-values to z-values 
t2p2z.val = qnorm(pt(sim.t,sim.df))  # equivalent to qnorm(1-p/2) for two.tailed p-values
summary(t2p2z.val)

source(zcurve3)
Title = paste0("t->p->z, ",round(sim.edr*100),"  ",round(sim.err*100))
bw.est = .02
res.3 = Zing(t2p2z.val);res.3
res.3 = res.3$res[2:3];res.3


### t2power2z transform t-values to power and power to z-values 
median.t.sig = median(sim.t[sim.t > qt(.975,sim.df)]);median.t.sig
use.nct = 0
use.ncz = 0
if (median.t.sig > 2.5) {
	use.ncz = qnorm(.5,qnorm(.975));use.ncz
	use.nct = pwr.t.test(n = N/2, 
	d=NULL, sig.level = .05, power = .5,
	type = "two.sample", alternative = "two.sided")$d/se;use.nct
}
if (median.t.sig > 3.5) {
	use.ncz = qnorm(.75,qnorm(.975));use.ncz
	use.nct = pwr.t.test(n = N/2, 
	d=NULL, sig.level = .05, power = .75,
	type = "two.sample", alternative = "two.sided")$d/se;use.nct
}
if (median.t.sig > 4.0) {
	use.ncz = qnorm(.90,qnorm(.975));use.ncz
	use.nct = pwr.t.test(n = N/2, 
	d=NULL, sig.level = .05, power = .90,
	type = "two.sample", alternative = "two.sided")$d/se;use.nct
}
use.nct
use.ncz
t2pow2z.val = qnorm(pt(sim.t,sim.df,use.nct),use.ncz)
summary(t2pow2z.val)

source(zcurve3)
Title = paste0("t->power->z, ",round(sim.edr*100),"  ",round(sim.err*100))
bw.est = .02
res.4 = Zing(t2pow2z.val);res.4
res.4 = res.4$res[2:3];res.4


source(zcurve3)
bw.est = .02
Title = paste0("t-curve, ",round(sim.edr*100),"  ",round(sim.err*100))
Est.Method = "DF"
Int.Beg = qt(.975,df);Int.Beg
res.5 = Zing(sim.t,df);res.5
res.5 = res.5$res[2:3];res.5


res.run = rbind(
	res.1,
	res.2,
	res.3,
	res.4,
	res.5
)

rownames(res.run) = c("z","t","t.p.z","t.pow.z","t-curve")

res.run

res.run = c(run.i,weights[run.i,],sim.err,res.run[,1],sim.edr,res.run[,2])
res.run

res = rbind(res,res.run)

write.table(x,"sim10k.Transformation.t30.dat") # write results each trial
                                   # can resume if stopped 
                                   # by loading the completed results

} # End of for loops

dim(res)  # 66 rows and 16 columns



#write the completed results / overwrites the data from the for loop
#write.table(res,"sim100k.Transformation.N30.dat")


########################################################
### GET STORED RESULTS
########################################################

# load the saved data
setwd("C:/Users/ulric/Dropbox/PHPCurves/DOCS/z-curve/Tutorial")
check = read.table("sim100k.Transformation.N30.dat",row.names=NULL)
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

#res = res[67:132,]

summary(res)

res[,1:4]

round(res[2,],2)

########################################################################
### Evaluation of ERR estimates: z-values versus t(28)->p->z values
########################################################################

# columns
# -  5: true ERR
# -  6: res.1
# -  7: res.2 
# -  8: res.3
# -  9: res.4
# - 10: res.5 

summary(res[,5:10])

round(res[,c(1:4,5:10)],3)

cor(res[,5:10])

rmse.err.1 = sqrt(mean((res[,5]-res[,6])^2));rmse.err.1
rmse.err.2 = sqrt(mean((res[,5]-res[,7])^2));rmse.err.2
rmse.err.3 = sqrt(mean((res[,5]-res[,8])^2,na.rm=TRUE));rmse.err.3
rmse.err.4 = sqrt(mean((res[,5]-res[,9])^2,na.rm=TRUE));rmse.err.4
rmse.err.5 = sqrt(mean((res[,5]-res[,10])^2,na.rm=TRUE));rmse.err.5

rmse.err.1
rmse.err.2
rmse.err.3
rmse.err.4
rmse.err.5

res = res[order(res[,5]),]

lwd.ci = .5

plot(res[,5],res[,7],pch=15,cex=1,xlim=c(0,1),ylim=c(0,1),col="purple3",ylab="",xlab="")  # plot OF and ER estimates
par(new=TRUE)
plot(res[,5],res[,6],xlim=c(0,1),ylim=c(0,1),pch=16,cex=1,col="forestgreen",ylab="Estimated ERR",
	xlab = "Simulated True ERR")  # plot OF and ER estimates
abline(a = 0, b = 1,lty=2)
legend(.1,.9,legend=c("z-values", "t-values"),lwd=2,col=c("forestgreen","purple3"))

###

plot(res[,5],res[,8],pch=15,cex=1,xlim=c(0,1),ylim=c(0,1),col="purple3",ylab="",xlab="")  # plot OF and ER estimates
par(new=TRUE)
plot(res[,5],res[,6],xlim=c(0,1),ylim=c(0,1),pch=16,cex=1,col="forestgreen",ylab="Estimated ERR",
	xlab = "Simulated True ERR")  # plot OF and ER estimates
abline(a = 0, b = 1,lty=2)
legend(.1,.9,legend=c("z-values", "t->p->z"),lwd=2,col=c("forestgreen","purple3"))

###

plot(res[,5],res[,9],pch=15,cex=1,xlim=c(0,1),ylim=c(0,1),col="purple3",ylab="",xlab="")  # plot OF and ER estimates
par(new=TRUE)
plot(res[,5],res[,6],xlim=c(0,1),ylim=c(0,1),pch=16,cex=1,col="forestgreen",ylab="Estimated ERR",
	xlab = "Simulated True ERR")  # plot OF and ER estimates
abline(a = 0, b = 1,lty=2)
legend(.1,.9,legend=c("z-values", "t->power->z"),lwd=2,col=c("forestgreen","purple3"))

###

plot(res[,5],res[,10],pch=15,cex=1,xlim=c(0,1),ylim=c(0,1),col="purple3",ylab="",xlab="")  # plot OF and ER estimates
par(new=TRUE)
plot(res[,5],res[,6],xlim=c(0,1),ylim=c(0,1),pch=16,cex=1,col="forestgreen",ylab="Estimated ERR",
	xlab = "Simulated True ERR")  # plot OF and ER estimates
abline(a = 0, b = 1,lty=2)
legend(.1,.9,legend=c("z-values", "t-curve"),lwd=2,col=c("forestgreen","purple3"))


########################################################################
### Evaluation of EDR estimates: z-values versus t(28)->p->z values
########################################################################

# columns
# - 11: true EDR
# - 12: res.1
# - 13: res.2 
# - 14: res.3
# - 15: res.4
# - 16: res.5

rmse.edr.1 = sqrt(mean((res[,11]-res[,12])^2));rmse.edr.1
rmse.edr.2 = sqrt(mean((res[,11]-res[,13])^2));rmse.edr.2
rmse.edr.3 = sqrt(mean((res[,11]-res[,14])^2));rmse.edr.3
rmse.edr.4 = sqrt(mean((res[,11]-res[,15])^2));rmse.edr.4
rmse.edr.5 = sqrt(mean((res[,11]-res[,16])^2));rmse.edr.5

rmse.edr.1
rmse.edr.2
rmse.edr.3
rmse.edr.4
rmse.edr.5


###

res = res[order(res[,11]),]

lwd.ci = .5

plot(res[,11],res[,13],pch=15,cex=1,xlim=c(0,1),ylim=c(0,1),col="purple3",ylab="",xlab="")  # plot OF and ER estimates
par(new=TRUE)
plot(res[,11],res[,12],xlim=c(0,1),ylim=c(0,1),pch=16,cex=1,col="forestgreen",
	ylab="Estimated EDR",xlab = "Simulated True EDR")  # plot OF and ER estimates
abline(a = 0, b = 1,lty=2)
legend(.1,.9,legend=c("z-values", "t-values"),lwd=2,col=c("forestgreen","purple3"))

###

plot(res[,11],res[,14],pch=15,cex=1,xlim=c(0,1),ylim=c(0,1),col="purple3",ylab="",xlab="")  # plot OF and ER estimates
par(new=TRUE)
plot(res[,11],res[,12],xlim=c(0,1),ylim=c(0,1),pch=16,cex=1,col="forestgreen",
	ylab="Estimated EDR",xlab = "Simulated True EDR")  # plot OF and ER estimates
abline(a = 0, b = 1,lty=2)
legend(.1,.9,legend=c("z-values", "t->p->z"),lwd=2,col=c("forestgreen","purple3"))

###

plot(res[,11],res[,15],pch=15,cex=1,xlim=c(0,1),ylim=c(0,1),col="purple3",ylab="",xlab="")  # plot OF and ER estimates
par(new=TRUE)
plot(res[,11],res[,12],xlim=c(0,1),ylim=c(0,1),pch=16,cex=1,col="forestgreen",
	ylab="Estimated EDR",xlab = "Simulated True EDR")  # plot OF and ER estimates
abline(a = 0, b = 1,lty=2)
legend(.1,.9,legend=c("z-values", "t->p.sig->z"),lwd=2,col=c("forestgreen","purple3"))
abline(v=.5,lty=2)
abline(h=.5,lty=2)

dir.bias = cbind(res[,1:4],res[,15]-res[,11])
dir.bias[abs(dir.bias[,5]) > .10,]

###

plot(res[,11],res[,16],pch=15,cex=1,xlim=c(0,1),ylim=c(0,1),col="purple3",ylab="",xlab="")  # plot OF and ER estimates
par(new=TRUE)
plot(res[,11],res[,12],xlim=c(0,1),ylim=c(0,1),pch=16,cex=1,col="forestgreen",
	ylab="Estimated EDR",xlab = "Simulated True EDR")  # plot OF and ER estimates
abline(a = 0, b = 1,lty=2)
legend(.1,.9,legend=c("z-values", "t-curve"),lwd=2,col=c("forestgreen","purple3"))

round(res[,c(1:4,11,16)],2)

dir.bias = cbind(res[,1:4],res[,16]-res[,11])
dir.bias
dir.bias[abs(dir.bias[,5]) > .20,]





