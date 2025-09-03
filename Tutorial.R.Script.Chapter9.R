

#### Co-Authord with Edward Lee 2025/09/02

###############################################################################
# Chapter 9 – P-Hacking Simulation 1: Sample Patchwork
###############################################################################

### 0.0 – SETUP

# Disable scientific notation
options(scipen = 999)

# Set work directory
setwd("C:/Users/ulric/Dropbox/PHPCurves/DOCS/z-curve/Tutorial")

# Install required packages (run once only)
#install.packages(c("zcurve", "KernSmooth", "stringr", "parallel","pwr"))

# Load z-curve and P-curve functions from a local file
setwd("C:/Users/ulric/Dropbox/PHPCurves/DOCS/z-curve/Tutorial")
zcurve3 <- "Zing.25.07.11.test.R"
source(zcurve3)
pcurve <- "Pcurve.Function.R"
source(pcurve)


# Alternatively, load z-curve and p-curve functions directly from GitHub
zcurve3 <- "https://raw.githubusercontent.com/UlrichSchimmack/zcurve3.0/refs/heads/main/Zing.25.07.11.test.R"
source(zcurve3)


### 1.0 – ILLUSTRATION WITH EXAMPLES

### get data from run.i = 1 simulation

pcurve.input = paste0("t(",df=results$df[results$sim.sel.p == 1],") = ",
  results$abs.t[results$sim.sel.p == 1])
pcurve.input = pcurve.input[results$V7[results$sim.sel.p == 1] < .05]
pcurve_app(pcurve.input,SHOW.PLOT=TRUE)


### get data from run.i = 22 simulation

### Figure 2

pcurve.input = paste0("t(",df=results$df[results$sim.sel.p == 1],") = ",
  results$abs.t[results$sim.sel.p == 1])
pcurve.input = pcurve.input[results$V7[results$sim.sel.p == 1] < .05]
pcurve_app(pcurve.input,SHOW.PLOT=TRUE)

### Figure 3

source(zcurve3)
sims[run.i,]
Title = "Mean = 0.2, SD = 0.4, n = 15"
ymax = 1.5
ymin = -.06
bw.est = .02
Est.Method = "EXT"
ncp = 3
zsds = 2
TEST4BIAS = TRUE
just = .8
Int.Beg = crit + just
Zing(results$z[results$sim.sel.p == 1]);res.8



### 2.0 – SIMULATION 

sim.runs = 1                 # repeat each simulation 100 times
sim.k.sig = 10000            # how many significant results?

sim.es.mean = seq(0,1,.2)    # mean effect size, Cohen's d
sim.es.sd   = seq(0,.8,.2)   # heterogeneity of es (tau)
sim.n.obs   = c(15,20,30)    # sample size per cell 

# create the simulation conditions

sims = c()
sim.i = 0
for (i in 1:length(sim.es.mean) ) {
	for (j in 1:length(sim.es.sd)) {
		for (k in 1:length(sim.n.obs) ) {
			sim.i = sim.i + 1
			sims = rbind(sims,c(sim.i,sim.es.mean[i],sim.es.sd[j],sim.n.obs[k]) )
}}}
sims = data.frame(sims)
colnames(sims) = c("run","es.mean","es.sd","n.obs")
dim(sims)

sims


b = 1            # form row b
e = 5            # to row e, used for testing
e = nrow(sims);e # to last row, run all conditions

res = c()        # collect the results 

###

############################
### Start Simulation Loop
############################

for (run.i in b:e ) { # for loop to run the simulations

#run.i = 22       # used for testing without loop
print(paste("Run: ",run.i))
sims[run.i,]


sample1 = 1:sims$n.obs[run.i]
sample2 = (sims$n.obs[run.i]+1):(2*sims$n.obs[run.i])
sample3 = (2*sims$n.obs[run.i]+1):(3*sims$n.obs[run.i])
samples <- list(sample1 = sample1,
                    sample2 = sample2,
                    sample3 = sample3,
                    sample12 = c(sample1, sample2),
                    sample13 = c(sample1, sample3),
                    sample23 = c(sample2, sample3),
                    sample123 = c(sample1, sample2, sample3))

samples

results = c()
count.sig = 0
count.run = 0

set.seed(run.i)

while (count.sig < sim.k.sig) {

count.run = count.run + 1

pop.es = rnorm(1,sims$es.mean[run.i],sims$es.sd[run.i]);pop.es
exp_group <- rnorm(n = 3*sims$n.obs[run.i], mean = pop.es, sd = 1)
summary(exp_group)

control <- rnorm(n = 3*sims$n.obs[run.i], mean = 0, sd = 1)
summary(control)

i.patch = 1   
res.run = c()
for (i.patch in 1:7) {
      sub_exp <- exp_group[unlist(samples[i.patch])]
      sub_control <- control[unlist(samples[i.patch])] 
      length(sub_exp)
      length(sub_control)
      tval <- t.test(sub_exp, sub_control, paird=FALSE,var.equal=TRUE)$statistic
      dfs <- t.test(sub_exp, sub_control, var.equal=TRUE)$parameter
      pval <- t.test(sub_exp, sub_control, var.equal=TRUE)$p.value
      se <- t.test(sub_exp, sub_control, var.equal=TRUE)$stderr
      es <- t.test(sub_exp, sub_control, var.equal=TRUE)$estimate
      res.run = rbind(res.run,c(count.run,pop.es,es[1]-es[2],se,tval,dfs,pval)	)
    }

res.run

results = rbind(results,res.run)

if (min(res.run[,7]) < .05) count.sig = count.sig + 1

#print(count.run)
#print(count.sig)


} # EOF while loop 

count.run
count.sig

dim(results)

results = data.frame(results)
results$N = results$df+2
results$se = 2/sqrt(results$N)
results$nct = abs(results[,2]/results$se) 
results$obs.es = results$nct*results$se 
results$abs.t = abs(results$t)
results$z = qnorm(pt(results$abs.t,results$N-2,log.p=TRUE),log.p=TRUE)
results$pow.dir = pt(qt(.975,results$df),results$df,results$nct,lower.tail=FALSE)
results$pow.sign.err = pt(-qt(.975,results$df),results$df,results$nct,lower.tail=TRUE)
results$pow = results$pow.dir + results$pow.sign.err

summary(results)

results$rn = rep(1:7,nrow(results)/7)

tapply(results$df,results$rn,mean)
tapply(results$pow,results$rn,mean)


true.edr = tapply(results$pow,results$rn,mean)[c(1,4,7)];true.edr

true.err.1 = tapply(results$pow*results$pow.dir,results$rn,sum)
true.err.2 = tapply(results$pow,results$rn,sum)
true.err = true.err.1/true.err.2
true.err = true.err[c(1,4,7)];true.err 

true.edr.2 = sum(results$pow[results$rn <= 3]/results$pow[results$rn <=3])/
  sum(1/results$pow[results$rn <= 3])
true.edr.2




true.edr
true.err

########################################
### Start P-Hacking 
########################################

phack = function(p) {
  sel = rep(0,7)
  if (p[1] < .05) sel[1] = 1	
  if (p[2] < .05) sel[2] = 1	
  if (p[3] < .05) sel[3] = 1	
  if(sum(as.numeric(p[c(1,2)] < .05)) == 0 & p[4] < .05) sel[4] = 1
  if(sum(as.numeric(p[c(1,3,4)] < .05)) == 0 & p[5] < .05) sel[5] = 1
  if(sum(as.numeric(p[c(2,3,4,5)] < .05)) == 0 & p[6] < .05) sel[6] = 1
  if(sum(as.numeric(p[1:6] < .05)) == 0 & p[7] < .05) sel[7] = 1
return(sel)
}
round(cbind(results$sim.sel.p,results$V7)[1:21,],2)


sim.sel.p = apply(matrix(results$V7,7,),2,function(x) phack(x) )
table(sim.sel.p)
results$sim.sel.p = c(sim.sel.p)
table(results$sim.sel.p,results$df)
table(results$sim.sel.p,results$se)


attempts1  = nrow(results)/7*3;attempts1
success1 = sum(as.numeric(results$V7[results$rn %in% 1:3] < .05)) / attempts1;success1

attempts2  = nrow(results)/7;attempts2
success2 = sim.k.sig / attempts2;success2

table(results$df,results$rn)
tapply(results$pow,results$rn,mean)
tapply(results$pow,list(sim.sel.p,results$df),mean)

table(results$sim.sel.p,results$df)
tapply(results$pow[results$sim.sel.p == 1],results$df[results$sim.sel.p == 1],mean)

check = cbind(results[,2],results$df,results$nct,results$pow,results$V7,results$sim.sel.p)
round(check[1:70,],2)
hist(results$pow[results$sim.sel.p == 1])
round(check[check[,4] > .99,],2)

table(results$sim.sel.p,results$df)


#PET-regression
results$v = results$se^2
w = 1/results$v
summary(lm(results$obs.es[results$sim.sel.p == 1] ~ results$se[results$sim.sel.p == 1]),
weights = w)


true.edr
success1

#res[run.i,]




#################################################
### Analyze the Simulated Data 
#################################################

source(zcurve3)
Title = paste(round(true.edr[1]*100),"  ",round(true.err[1]*100))
ymax = 1
ymin = -.06
bw.est = .02
res.1 = Zing(results$z[results$rn %in% 1:3]);res.1
res.1 = res.1$res
res.1

source(zcurve3)
Title = paste(round(true.edr[1]*100),"  ",round(true.err[1]*100))
ymax = 1
ymin = -.06
bw.est = .02
Est.Method = "EXT"
ncp = 3
zsds = 2
res.2 = Zing(results$z[results$rn %in% 1:3]);res.2
res.2 = res.2$res
res.2

source(zcurve3)
Title = paste(round(true.edr[1]*100),"  ",round(true.err[1]*100))
ymax = 1
ymin = -.06
bw.est = .02
TEST4BIAS = TRUE
just = .8
Int.Beg = crit + just
res.3 = Zing(results$z[results$rn %in% 1:3]);res.3
res.3 = c(res.3$res,res.3$bias)
res.3

source(zcurve3)
Title = paste(round(true.edr[1]*100),"  ",round(true.err[1]*100))
ymax = 1
ymin = -.06
bw.est = .02
Est.Method = "EXT"
ncp = 3
zsds = 2
TEST4BIAS = TRUE
just = .8
Int.Beg = crit + just
res.4 = Zing(results$z[results$rn %in% 1:3]);res.4
res.4 = c(res.4$res,res.4$bias)
res.4


###########################
### patchwork 
###########################


source(zcurve3)
Title = paste(round(true.edr[1]*100),"  ",round(true.err[1]*100))
ymax = 1
ymin = -.06
bw.est = .02
#res.5 = Zing(results$z[results$sim.sel.p == 1 & results$rn %in% 1:3]);res.5
res.5 = Zing(results$z[results$sim.sel.p == 1]);res.5
res.5 = res.5$res
res.5

source(zcurve3)
Title = paste(round(true.edr[1]*100),"  ",round(true.err[1]*100))
ymax = 1
ymin = -.06
bw.est = .02
Est.Method = "EXT"
ncp = 3
zsds = 2
#res.6 = Zing(results$z[results$sim.sel.p == 1 & results$rn %in% 1:3]);res.6
res.6 = Zing(results$z[results$sim.sel.p == 1]);res.6
res.6 = res.6$res
res.6

source(zcurve3)
Title = paste(round(true.edr[1]*100),"  ",round(true.err[1]*100))
ymax = 1
ymin = -.06
bw.est = .02
TEST4BIAS = TRUE
just = .8
Int.Beg = crit + just
#res.7 = Zing(results$z[results$sim.sel.p == 1 & results$rn %in% 1:3]);res.7
res.7 = Zing(results$z[results$sim.sel.p == 1]);res.7
res.7 = c(res.7$res,res.7$bias)
res.7

source(zcurve3)
sims[run.i,]
Title = "Mean = 0.2, SD = 0.4, n = 15"
Title = paste(round(true.edr[1]*100),"  ",round(true.err[1]*100))
ymax = 1.5
ymin = -.06
bw.est = .02
Est.Method = "EXT"
ncp = 3
zsds = 2
TEST4BIAS = TRUE
just = .8
Int.Beg = crit + just
#res.8 = Zing(results$z[results$sim.sel.p == 1 & results$rn %in% 1:3]);res.8
res.8 = Zing(results$z[results$sim.sel.p == 1]);res.8
res.8 = c(res.8$res,res.8$bias)
res.8

res.run = c(unlist(sims[run.i,]),
	success1,success2,true.err,true.edr,
	res.1,res.2,res.3,res.4,
	res.5,res.6,res.7,res.8
)

print(round(res.run,2))

res = rbind(res,res.run)

write.table(res,"sim.patchwork.bias.dat")    # write results each trial
                                             # can resume if stopped 
                                             # by loading the completed results

} # End of for loop

dim(res)  # 

#write the completed results / overwrites the data from the for loop
#write.table(res,"sim.patchwork.bias.dat")


########################################################
### GET STORED RESULTS
########################################################

# load the saved data
setwd("C:/Users/ulric/Dropbox/PHPCurves/DOCS/z-curve/Tutorial")
check = read.table("sim.patchwork.bias.dat",row.names=NULL)
check = data.frame(check)
check = check[,2:dim(check)[2]]
dim(check)

table(round(res[,c(1:14)],5) == round(check[,c(1:14)],5) )

check[,c(1:6,13)]

# careful, if you run this code any results in the res file will be replaced 
# with the saved data. 

#res = check  # use only if sure to replace res; run analysis with "res" matrix 

res = data.frame(res)
dim(res)

summary(res)



#############################################
### Inflation (False Positive Rate for EDR = 5%
#############################################

summary(res)

#ODR 
summary(res[,c(5,10,13,18,23,31,39)])

#EDR 
summary(res[,c(10,15,20,25,33,41,46,51,59)])

#ERR 
summary(res[,c(7,14,19,24,32,40,45,50,58)])

#BIAS 
#default no bias
summary(res[,c(28,29,30)])
table(res[,30] < .05)

#normal, no bias
summary(res[,c(36,37,38)])
table(res[,38] < .05)

#default, bias
summary(res[,c(54,55,56)])
table(res[,56] < .05)

#normal, bias
summary(res[,c(62,63,64)])
table(res[,64] < .05)

most.bias = which(res[,62] - res[,63] == max(res[,62] - res[,63]))
res[most.bias,c(1:4,62,63,64)]


###

summary(res[,c(10,5,6)]) # True EDR 

round(res[,c(10,6)],2)

plot(res[,10],res[,6],xlim=c(0,1),ylim=c(0,1),xlab="True Discovery Rate",
	ylab="Inflated Discovery Rate with P-hacking",pch=15)
abline(b = 1,a = 0,lty=2)

inflation = res[,6] - res[,10]
type1 = as.numeric(res[,1] <= 3)
type1
table(type1)
tapply(inflation,type1,mean)

summary(lm(inflation ~ type1 + res[,2] + res[,3] + res[,4]))

round(tapply(inflation,list(res[,2],res[,3]),mean),2)




########################################################################
### Evaluation of ERR estimates (no bias)
########################################################################

# columns
# -  7: true ERR
# - 14: z-curve default
# - 19: z-curve normal
# - 24: z-curve default Int.Beg = 2.8
# - 32: z-curve normal Int.Beg = 2.8

#ERR 
summary(res[,c(7,14,19,24,32)])

rmse.err.1 = sqrt(mean((res[,7]-res[,14])^2));rmse.err.1
rmse.err.2 = sqrt(mean((res[,7]-res[,19])^2));rmse.err.2
rmse.err.3 = sqrt(mean((res[,7]-res[,24])^2,na.rm=TRUE));rmse.err.3
rmse.err.4 = sqrt(mean((res[,7]-res[,32])^2,na.rm=TRUE));rmse.err.4
###
rmse.err.1
rmse.err.2
rmse.err.3
rmse.err.4

lwd.ci = .5
graphics.off()
plot(res[,7],res[,14],xlim=c(0,1),ylim=c(0,1),pch=16,cex=1,col="forestgreen",ylab="Estimated ERR",
	xlab = "Simulated True ERR")  # plot OF and ER estimates
par(new=TRUE)
plot(res[,7],res[,19],xlim=c(0,1),ylim=c(0,1),pch=16,cex=1,col="firebrick3",ylab="Estimated ERR",
	xlab = "Simulated True ERR")  # plot OF and ER estimates
abline(a = 0, b = 1,lty=2)
abline(v = .5,lty=2)
abline(h = .5,lty=2)
legend(.1,.9,legend=c("z-curve: default", "z-curve: normal"),lwd=2,col=c("forestgreen","purple3"))



########################################################################
### Evaluation of ERR estimates (patchwork p-hacked) 
########################################################################

summary(res[,c(10,15,20,25,33,41,46,51,59)])
# -  7: true ERR
# - 40: z-curve default
# - 45: z-curve normal
# - 50: z-curve default Int.Beg = 2.8
# - 58: z-curve normal Int.Beg = 2.8

rmse.err.sel.1 = sqrt(mean((res[,7]-res[,40])^2));rmse.err.sel.1
rmse.err.sel.2 = sqrt(mean((res[,7]-res[,45])^2));rmse.err.sel.2
rmse.err.sel.3 = sqrt(mean((res[,7]-res[,50])^2,na.rm=TRUE));rmse.err.sel.3
rmse.err.sel.4 = sqrt(mean((res[,7]-res[,58])^2,na.rm=TRUE));rmse.err.sel.4

rmse.err.sel.1
rmse.err.sel.2
rmse.err.sel.3
rmse.err.sel.4

### sel significant results

graphics.off()
plot(res[,7],res[,45],pch=15,cex=1,xlim=c(0,1),ylim=c(0,1),col="purple3",ylab="",xlab="")  # plot OF and ER estimates
par(new=TRUE)
plot(res[,7],res[,58],xlim=c(0,1),ylim=c(0,1),pch=16,cex=1,col="forestgreen",ylab="Estimated ERR",
	xlab = "Simulated True ERR")  # plot OF and ER estimates
abline(a = 0, b = 1,lty=2)
abline(v = .5,lty=2)
abline(h = .5,lty=2)
legend(.01,.9,legend=c("z-curve - normal, > 1.96 ", "z-curve - normal, z > 2.8"),lwd=2,col=c("purple3","forestgreen"))


dim(res)
res$dir.bias.1 = (res[,45] - res[,7])
summary(res$dir.bias.1)

res[res$dir.bias.1 < -.10,c(1:4,7,45,58)]



########################################################################
### Evaluation of EDR estimates (no bias)
########################################################################

# columns
# - 10: true EDR
# - 15: z-curve default
# - 20: z-curve normal
# - 25: z-curve default Int.Beg = 2.8
# - 33: z-curve normal Int.Beg = 2.8

#EDR 
summary(res[,c(10,15,20,25,33)])

rmse.edr.1 = sqrt(mean((res[,10]-res[,15])^2));rmse.edr.1
rmse.edr.2 = sqrt(mean((res[,10]-res[,20])^2));rmse.edr.2
rmse.edr.3 = sqrt(mean((res[,10]-res[,25])^2,na.rm=TRUE));rmse.edr.3
rmse.edr.4 = sqrt(mean((res[,10]-res[,33])^2,na.rm=TRUE));rmse.edr.4
###
rmse.edr.1
rmse.edr.2
rmse.edr.3
rmse.edr.4

lwd.ci = .5
graphics.off()
plot(res[,10],res[,15],xlim=c(0,1),ylim=c(0,1),pch=16,cex=1,col="purple3"
	,ylab="Estimated EDR",xlab = "Simulated True EDR") 
par(new=TRUE)
plot(res[,10],res[,20],xlim=c(0,1),ylim=c(0,1),pch=16,cex=1,col="forestgreen",ylab="",
	xlab = "")  # plot OF and ER estimates
abline(a = 0, b = 1,lty=2)
abline(v = .5,lty=2)
abline(h = .5,lty=2)
legend(.1,.9,legend=c("z-curve: default", "z-curve: normal"),lwd=2,col=c("forestgreen","purple3"))



########################################################################
### Evaluation of ERR estimates (patchwork p-hacked) 
########################################################################

# - 10: true EDR
# - 41: z-curve default
# - 46: z-curve normal
# - 51: z-curve default Int.Beg = 2.8
# - 59: z-curve normal Int.Beg = 2.8

rmse.edr.sel.1 = sqrt(mean((res[,10]-res[,41])^2));rmse.edr.sel.1
rmse.edr.sel.2 = sqrt(mean((res[,10]-res[,46])^2));rmse.edr.sel.2
rmse.edr.sel.3 = sqrt(mean((res[,10]-res[,51])^2,na.rm=TRUE));rmse.edr.sel.3
rmse.edr.sel.4 = sqrt(mean((res[,10]-res[,59])^2,na.rm=TRUE));rmse.edr.sel.4

rmse.edr.sel.1
rmse.edr.sel.2
rmse.edr.sel.3
rmse.edr.sel.4

### sel significant results

graphics.off()
plot(res[,10],res[,46],pch=15,cex=1,xlim=c(0,1),ylim=c(0,1),col="purple3",ylab="",xlab="")  # plot OF and ER estimates
par(new=TRUE)
plot(res[,10],res[,59],xlim=c(0,1),ylim=c(0,1),pch=16,cex=1,col="forestgreen",ylab="Estimated ERR",
	xlab = "Simulated True ERR")  # plot OF and ER estimates
abline(a = 0, b = 1,lty=2)
abline(v = .5,lty=2)
abline(h = .5,lty=2)
legend(.1,.9,legend=c("z-curve: default", "z-curve, z > 2.8"),lwd=2,col=c("forestgreen","purple3"))


dim(res)
res$dir.bias.edr.sel = (res[,46] - res[,10])
summary(res$dir.bias.edr.sel)

res[res$dir.bias.edr.sel < -.10,c(1:4,10,46,59)]



