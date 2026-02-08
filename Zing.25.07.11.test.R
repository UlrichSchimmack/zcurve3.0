

########################################################################
### SETTING PARAMETERS FOR Z-CURVE MODEL
########################################################################

version <- "Version 25.12.14"   # Version label to appear on plots

# Optional cleanup 
# rm(list = ls())
# options(scipen = 999)  # Disable scientific notation

### INSTALL PACKAGES (only once – manually run if needed)
if (1 == 2) {  # This block is ignored unless manually changed to (1 == 1)
  install.packages("pwr")
  install.packages("zcurve")
  install.packages("KernSmooth")
  install.packages("parallel")
  install.packages("stringr")
} # END install block

### LOAD LIBRARIES
library(parallel)
library(KernSmooth)
library(zcurve)
library(stringr)
library(pwr)

# Optional: Setup parallel processing (currently not functional)
# cl <- parallel::makeCluster(parallel::detectCores() - 10)


########################################################################
### GLOBAL PARAMETERS
########################################################################

TESTING <- FALSE             # Toggle for development/debugging mode

CURVE.TYPE = "z"             # Set to "t" for t-distribution
df = c()

### SPEED CONTROL

parallel <- FALSE            # Placeholder – parallel functionality not yet implemented
max_iter <- 1e6              # Max iterations for model estimation
max_iter_boot <- 1e5         # Max iterations for bootstrapped estimates

EM.criterion <- 1e-3         # Convergence threshold for EM algorithm
EM.max.iter <- 1000          # Max iterations for EM

Plot.Fitting <- FALSE        # Plot fitting curve (only for Est.Method = "OF" or "EXT")

### PLOT SETTINGS

Title <- ""                  # Optional plot title

letter.size <- 1             # Text size in plots
letter.size.1 <- letter.size # Used for version labels in plot
y.line.factor <- 3           # Controls spacing of plot text

x.lim.min <- 0               # X-axis lower bound
x.lim.max <- 6               # X-axis upper bound
ymax <- 0.6                  # Y-axis upper bound
ymin <- 0                    # OUTDATED Y-axis lower bound (for label space)

Show.Histogram <- TRUE       # Toggle histogram in plot
Show.Text <- TRUE            # Toggle model results in plot
Show.Curve.All <- TRUE       # Show predicted z-curve
Show.Curve.Sig <- FALSE      # Option: show z-curve only for significant values
Show.Significance <- TRUE    # Show z = critical value line
Show.KD <- FALSE             # Toggle kernel density overlay (density method only)

sig.levels <- c()            # Optional: mark additional p-value thresholds on plot

int.loc <- 0.5               # Plot local power intervals below x-axis (set 0 to disable)
hist.bar.width <- 0.2        # Width of histogram bars
bw.draw <- 0.10              # Smoothing for kernel density display

### CONSOLE OUTPUT

Show.Iterations <- TRUE      # Show iterations for slow procedures (e.g., EXT, TEST4HETEROGENEITY)

### MODEL PARAMETERS

alpha <- 0.05                        # Significance level
crit <- qnorm(1 - alpha / 2)       # Corresponding two-sided critical z

two.sided <- TRUE                   # Assume two-sided z-values (use abs(z)); not yet compatible with signed z-values

# Color scheme
col.curve <- "violetred3"
col.hist <- "blue3"
col.kd <- "green3"

Est.Method <- "OF"                  # Estimation method: "OF", "EM", or "EXT"
                                    # Clustered Data: "CLU-W" (weighted),"CLU-B" (bootstrap)   
Int.Beg <- 1.96                     # Default: critical value for alpha = .05
Int.End <- 6                        # End of modeling interval (z > 6 = power = 1)

ncp <- 0:6                          # Component locations (z-values at which densities are centered)
components <- length(ncp)           # Number of components
zsd <- 1                            # SD of standard normal z-distribution
zsds = rep(zsd,components)          # one SD for each component

just <- 0.8                         # Cutoff for "just significant" z-values (used in optional bias test)

ZSDS.FIXED <- FALSE                 # Fix SD values for EXT method 
NCP.FIXED <- FALSE                  # Fix non-central parameter(NCP) means values for EXT method
W.FIXED   <- FALSE                  # Fix weights for EXT method

fixed.false.positives <- 0          # If > 0, constrains proportion of false positives (e.g., weight for z = 0 component)

### DENSITY-BASED SETTINGS (Only used with Est.Method = "OF")

n.bars <- 512                       # Number of bars in histogram

Augment <- TRUE                     # Apply correction for bias at lower bound
Augment.Regression <- FALSE         # Use Slope for Augmentation
Augment.Factor <- 1                 # Amount of augmentation

bw.est <- 0.05                      # Bandwidth for kernel density (lower = less smoothing, higher = more smoothing)
bw.aug <-  .20					 # Width of Augmentation interval

### INPUT RESTRICTIONS

MAX.INP.Z <- Inf                    # Optionally restrict very large z-values (set Inf to disable)

### CONFIDENCE INTERVALS / BOOTSTRAPS

boot.iter <- 0                      # Number of bootstrap iterations (suggest 500+ for final models)
ERR.CI.adjust <- 0.03               # Conservative widening of confidence intervals for ERR
EDR.CI.adjust <- 0.05               # Conservative widening for EDR

CI.ALPHA <- 0.05                    # CI level (default = 95%)

### CI levels for Heterogeneity Test
fit.ci <- c(.01, .025, .05, .10, .17, .20, .50, .80, .83, .90, .95, .975, .99)  # CI levels for model fit test

TEST4BIAS <- FALSE                  # Enable optional bias test
TEST4HETEROGENEITY <- 0             # Optional heterogeneity test (slow) — set number of bootstrap iterations


### DISPLAY FINAL STATUS

print(version)
print("Parameter OK")

##################################################################
### END OF SETTING DEFAULT PARAMETERs
##################################################################


##################################################################
##################################################################
### AAA Start of Zing Code
##################################################################
##################################################################

Zing = function(val.input,df=c(),cluster.id=c(),p=c(), lp=c()) {

##################################################################

### Functions used in Zing

#######################################################################

get.t.density = function(tx) {     

n.bars = length(tx);n.bars
bar.width = tx[2]-tx[1]

Dens	= c()
for(i in 1:n.bars) {
	for (j in 1:length(ncp)) {
		Dens = c(Dens,
		dt(tx[i],df,ncp[j]) +
		dt(-tx[i],df,ncp[j])
		)
	}
}
Dens = matrix(Dens,length(ncp),byrow=FALSE)
row.sum.dens = rowSums(Dens)
Dens = Dens/(row.sum.dens * bar.width)

return(Dens)

} # EOF get.t.density 

####################################

get.z.density = function(zx) {     

n.bars = length(zx);n.bars
bar.width = zx[2] - zx[1]

Dens	= c()
for(i in 1:n.bars) {
	for (j in 1:length(ncp)) {
		Dens = c(Dens,
		dnorm(zx[i],ncp[j],zsds[j]) +
		dnorm(-zx[i],ncp[j],zsds[j])
		)
	}
}
Dens = matrix(Dens,length(ncp),byrow=FALSE)
row.sum.dens = rowSums(Dens)
Dens = Dens/(row.sum.dens * bar.width)
dim(Dens)

return(Dens)

} #EOF get.z.density


###########################################################


### FUN.2 - HETEROGENEITY TEST

run.heterogeneity.test = function(val.input,boot.run = 500, 
	fit.ci = c(.01,.025,.05,.10,.17,.20,.50,.80,.83,.90,.95,.975,.99)
	) {


res.boot = c()

boot.i = 1
for (boot.i in 1:boot.run ) {
	print(boot.i)

	zboot = sample(val.input,replace=TRUE)
	summary(zboot)

	ncp = 2
	zsds = 1
	components = length(ncp)
	para.est.EXT = extended.curve(zboot,ncp,zsds);para.est.EXT
	fit.hom = para.est.EXT[(3*components+1)];fit.hom

	NCP.FIXED = FALSE
	ncp = c(0,3)
	zsds = c(1,1)
	components = length(ncp)
	para.est.EXT = extended.curve(zboot,ncp,zsds);para.est.EXT
	fit.het.1 = para.est.EXT[length(para.est.EXT)];fit.het.1

	ncp = 2
	zsds = 5
	components = length(ncp)
	para.est.EXT = extended.curve(zboot,ncp,zsds);para.est.EXT
	fit.het.2 = para.est.EXT[4];fit.het.2
	sd.het.2 = para.est.EXT[3];sd.het.2

	delta.fit.1 = fit.hom - fit.het.1
	delta.fit.2 = fit.hom - fit.het.2
	delta.fit.3 = fit.het.2 - fit.het.1

	res.boot = rbind(res.boot,c(fit.hom,fit.het.1,fit.het.2,sd.het.2,
		delta.fit.1,delta.fit.2,delta.fit.3))

} # EOF boot iterations

fit.hom.ci = quantile(res.boot[,1],fit.ci);fit.hom.ci
fit.het.1.ci = quantile(res.boot[,2],fit.ci);fit.het.1.ci
fit.het.2.ci = quantile(res.boot[,3],fit.ci);fit.het.2.ci

sd.het.2.ci = quantile(res.boot[,4],fit.ci);sd.het.2.ci

delta.fit.1.ci = quantile(res.boot[,5],fit.ci);delta.fit.1.ci
delta.fit.2.ci = quantile(res.boot[,6],fit.ci);delta.fit.2.ci
delta.fit.3.ci = quantile(res.boot[,7],fit.ci);delta.fit.3.ci

return.res = cbind(
	fit.hom.ci,fit.het.1.ci,fit.het.2.ci,
	sd.het.2.ci,
	delta.fit.1.ci,delta.fit.2.ci,delta.fit.3.ci)

dim(return.res)

colnames(return.res) = c(
	"fit.hom","fit.het.1","fit.het.2",
	"sd.het.2",	
	"delta.fit.2-1","delta.fit.3-1","delta.fit.2-3"
)

return(return.res)

} # EOF run.heterogeneity.test



##############################################################
##############################################################
##############################################################


### function to detect bias

test.bias = function(w.all) {

#Int.Beg = 2 + .4

D.X = seq(x.lim.min,x.lim.max,.01)
n.bars = length(D.X)
bar.width = D.X[2]-D.X[1]

components = length(ncp)
means = ncp
sds = zsds
w = w.all


### get the densities for each interval and each non-centrality parameter
Dens	= c()
for(i in 1:n.bars) {
	for (j in 1:components) {
		if (CURVE.TYPE == "z") {
			Dens = c(Dens,
			dnorm(D.X[i],means[j],sds[j]) +
			dnorm(-D.X[i],means[j],sds[j]) 
			)
		} else {
			sds = sds - 1
			sds[sds <= 0] = 0
			sds
			if (max(sds) > 0) {
				Dens = c(Dens,
				dt(D.X[i],df,means[j]) +
				dt(-D.X[i],df,means[j])
				)
			} else { 
				Dens = c(Dens,
				dt(D.X[i],df,means[j]) +
				dt(-D.X[i],df,means[j])
				)
			}
		}
	}
}

Dens = matrix(Dens,length(ncp),byrow=FALSE)

### rescale the densities for the range of z-values to 1
row.sum.dens = rowSums(Dens)
Dens = Dens/(row.sum.dens * bar.width)
dim(Dens)

### compute the new estimated density distribution
E.D.Y = colSums(Dens*w)
length(E.D.Y)

E.D.Y = E.D.Y / sum(E.D.Y * bar.width)
sum(E.D.Y * bar.width)

prob1 = sum(E.D.Y[D.X > crit & D.X < crit + just]*bar.width);prob1
prob2 = sum(E.D.Y[D.X > crit + just & D.X < Int.End]*bar.width);prob2
prob3 = sum(E.D.Y[D.X < crit]*bar.width);prob3
sum(prob1,prob2,prob3)
prob = prob1 / (prob1 + prob2)
prob
adj = .00
prob.adj = prob + adj 

k = length(val.input[val.input < Int.End])

sig.k = length(val.input[val.input > crit & val.input < Int.End])
sig.k

#Int.Beg = 0

just.sig.k = length(val.input[val.input > crit & val.input < crit + just])
just.sig.k

just.sig.k/sig.k

### binomial
p.bias.binomial = 1 - pbinom(just.sig.k - 1, sig.k, prob = prob.adj);p.bias.binomial
p.bias.binomial
print(p.bias.binomial)


### for all
#prob = prob1 / (prob1 + prob2 + prob3)
#just.sig.k/k
#p.bias.binomial = 1 - pbinom(just.sig.k - 1, k, prob = prob.adj);p.bias.binomial


### chi2 approximation
#O = just.sig.k
#E = sig.k*prob;E
#chi2.val <- (O - E)^2 / (E * (1 - prob));chi2.val
#p.bias.chi2 <- 1 - pchisq(chi2.val, df = 1)
#print(p.bias.chi2)

bias.res = c(just.sig.k/sig.k,prob.adj,p.bias.binomial)
print(bias.res)

return(bias.res)

} #EOF bias test



### function to run the zcurve package to get EM weights and fit

run.zcurve = function(Est.Method = "OF",kd.model="KD2",K=6,
	alpha = .05,Int.Beg = 1.96,Int.End = 6, boot.iter = 0,parallel=TRUE) {

	meth = Est.Method
	if (Est.Method == "OF") meth = "density"

	biter = boot.iter
	if (boot.iter == 0) biter = FALSE

	z.res = zcurve(val.input,bootstrap=biter,method=meth,
		control=list(parallel = parallel,
		sig_level=alpha,a = Int.Beg,b = Int.End,mu=ncp))

	return(z.res)

} ### End run.zcurve

if (TESTING) {
	meth = "density"
	z.res = run.zcurve(val.input,Est.Method=meth,boot.iter=0)
	z.res
	#str(z.res)
	round(z.res$fit$weights,3)
	w.inp = z.res$fit$weights
	round(w.inp,3)
}


############################################
############################################
############################################

EXT.boot = function()	{

		#INT = z.rep

		boot = 1
		boot.res = c()
		
		#boot.iter = 20
		components = length(ncp)
		for (boot in 1:boot.iter) {

			### Get Bootstrap Sample
   		  	val.sample = sample(INT, size=length(INT), replace=TRUE)

			### Submit Bootstrap Sample to Parameter Estimation Function

			if (Show.Iterations) print(paste0("EXT.boot Iteration: ",boot))

			para.est.EXT = extended.curve(val.sample,ncp,zsds);para.est.EXT
		
			fit = para.est.EXT[(3*components+1)];fit
			para.est.EXT = para.est.EXT[1:(3*components)]

			w.inp = para.est.EXT[1:components]
			w.inp = w.inp/sum(w.inp)
			w.inp

			ncp = para.est.EXT[(components+1):(2*components)]
			ncp

			zsds = para.est.EXT[(2*components+1):(3*components)]
			zsds

			cp.input = c(w.inp,ncp,zsds)

			if (mean(zsds) > 1.1) { cp.res = 
				Compute.Power.SDG1(cp.input,BOOT=TRUE)
			} else { cp.res = Compute.Power.Z(cp.input) }

			cp.res

			EDR = cp.res[1];EDR
			ERR = cp.res[4];ERR

			w.all = cp.res[which(substring(names(cp.res),1,5) == "w.all")];w.all
			w.sig = cp.res[which(substring(names(cp.res),1,5) == "w.sig")];w.all

			w.all[is.na(w.all)] = 0
			w.sig[is.na(w.sig)] = 0

			boot.res = rbind(boot.res,c(ERR,EDR,fit,ncp,zsds,w.all,w.sig))
			round(boot.res,3)

		}

		#round(boot.res,3)

		CIs = c()
	
		CIs = rbind(CIs,quantile(boot.res[,1],c(CI.ALPHA/2,1-CI.ALPHA/2)) )
		CIs[1,1] = CIs[1,1]-ERR.CI.adjust
		CIs[1,2] = CIs[1,2]+ERR.CI.adjust
		if (CIs[1,1] < alpha/2) CIs[1,1] = alpha/2
		if (CIs[1,2] > 1) CIs[1,2] = 1

		CIs = rbind(CIs,quantile(boot.res[,2],c(CI.ALPHA/2,1-CI.ALPHA/2)) )
		CIs[2,1] = CIs[2,1]-EDR.CI.adjust
		CIs[2,2] = CIs[2,2]+EDR.CI.adjust
		if (CIs[2,1] < alpha) CIs[2,1] = alpha
		if (CIs[2,2] > 1) CIs[2,2] = 1

		FDR = (1/CIs[2,] - 1)*(alpha/(1-alpha))
		CIs = rbind(CIs,FDR[2:1])

		i = 6
		boot.res[,6]
		for (i in 3:dim(boot.res)[2]) CIs = rbind(CIs,quantile(boot.res[,i],
			c(CI.ALPHA/2,1-CI.ALPHA/2)) )

		rownames(CIs) = c("ERR","EDR","FDR","FIT",
			rep("NCP",components),
			rep("ZSD",components),
			rep("WALL",components),
			rep("WSIG",components)
		)

		if (components == 1) CIs[which(substring(rownames(CIs),1,1) == "W"),] = 1

		#dim(CIs)
		#print(CIs)

		return(CIs)


} # End function EXT.boot

#####################################
#####################################
#####################################

get.ci.info = function(Est.Method = "EM") {

#val.input = abs(c(rnorm(1000,1.1),rnorm(1000,2.8)))
#boot.iter = 0

if (Est.Method %in% c("OF","EM")) {

#Int.Beg = 2
zres = run.zcurve(val.input, Est.Method=Est.Method, alpha = alpha,boot.iter = boot.iter,
	Int.Beg = Int.Beg, Int.End = Int.End,parallel=parallel)

zres
ERR = summary(zres)$coefficients[1,];ERR
EDR = summary(zres)$coefficients[2,];EDR
FDR = (1/EDR - 1)*(alpha/(1-alpha));FDR
FDR = FDR[c(1,3,2)];FDR

w.inp = summary(zres, type="parameters")$coefficients[,2]
w.inp.low = summary(zres, type="parameters")$coefficients[,3]
w.inp.high = summary(zres, type="parameters")$coefficients[,4]

round(cbind(w.inp,w.inp.low,w.inp.high),3)

pow.sel = pnorm(ncp,Int.Beg) + pnorm(-ncp,Int.Beg);pow.sel

sum(w.all)
w.all = w.inp/(pow.sel)
w.all = w.all/sum(w.all)
round(w.all,3)

w.inp.low = w.inp.low/sum(w.inp.low)
w.all.low = w.inp.low/(pow.sel)
w.all.low = w.all.low/sum(w.all.low)
round(w.all.low,3)

w.inp.high = w.inp.high/sum(w.inp.high)
w.all.high = w.inp.high/(pow.sel)
w.all.high = w.all.high/sum(w.all.high)
w.all.high

fit.val = summary(zres)$model$fit_index
fit.val

w.all = cbind(w.all,w.all.low,w.all.high)

res.ci = rbind(EDR,ERR,FDR,fit.val,w.all)
print(round(res.ci,3))


} 

if (Est.Method == "EXT") {

res.ci = EXT.boot()

res.ci

#print("METHOD EXT")
#print("RETURN res.ci")
#print(res.ci)

}

return(res.ci)

}

#######################################################################


Write.Local.Power = function(loc.power) {
	names(loc.power) = paste0("LP",seq(1,length(loc.power)))
	int = seq(x.lim.min,x.lim.max-int.loc,int.loc)+int.loc/2

	old_mar <- par("mar")
	par(mar = old_mar + c(1, 0, 0, 0))  # increase bottom margin by 2 lines

	new.y = 0

#for (i in 1:length(int)) text(int[i],new.y,
#	paste0(format(round(loc.power[i]*100),nsmall=0),"%"),srt=90,pos=1
#                              ,cex=letter.size)
#par(xpd = FALSE)

lab = c()
for (i in 1:length(int)) lab[i] = paste0("   ",format(round(loc.power[i]*100),nsmall=0),"%")
lab

x_pos  <- seq(x.lim.min, x.lim.max, length.out = length(lab))
x_pos

#axis(1)

mtext(lab, side=1, line=.8, at=x_pos, cex=1.0, las=2)

par(mar = old_mar) 

}


###################################################
#### Begin Draw Histogram
###################################################

Draw.Histogram = function(w,cola="blue3",
	results,Write.CI = FALSE) {

#	int.start = round(Int.Beg,1)
#	if (round(Int.Beg,2) == 1.96) int.start = 2


#w = w.all;cola = "blue3"; results = res.text


make_breaks_straddle <- function(Int.Beg=1.96, bw=.2, x_min, x_max) {
  # Put Int.Beg at the center of a bin
  start0 <- Int.Beg - bw/2

  # Shift left so that start <= x_min
  n_left <- ceiling((start0 - x_min) / bw)
  start <- start0 - n_left * bw

  # Ensure we go past x_max
  end <- start + ceiling((x_max - start) / bw) * bw

  seq(start, end, by = bw)
}

#bw = .2
brks <- make_breaks_straddle(Int.Beg, bw, x.lim.min, x.lim.max)



	z.hist.val = val.input[val.input > x.lim.min & val.input < x.lim.max]
	if (round(Int.Beg,2) == 1.96) {
		z.hist = val.input[val.input > x.lim.min & val.input < x.lim.max -.04] + .04
	}
	table(z.hist > 1.96)

	n.breaks = seq(x.lim.min,x.lim.max,hist.bar.width);n.breaks

	par(cex.axis = 1)
	par(cex.lab = 1)
	par(family = "sans") 
	par(font.axis = 1) #"Microsoft Himalaya")
	par(font.lab = 1) # Microsoft Himalaya")

	col1 = adjustcolor(cola, alpha.f = 0.2)
	col2 = adjustcolor(cola, alpha.f = 0.3)

	h.res = hist(z.hist.val,breaks=brks,freq=FALSE,
		col=col1,border="white",
		xlim=c(x.lim.min,x.lim.max),ylim=c(ymin,ymax),
		ylab="Density",xlab="",axes=FALSE,main=Title,lwd=1)


	bin_left  <- h.res$breaks[-length(h.res$breaks)]
	bin_right <- h.res$breaks[-1]
	dark_bins <- bin_left >= Int.Beg

	par(new=TRUE)

op <- par(mar = c(6, 4, 2, 2)) 

plot(
  h.res,
  col = NA,
  border = NA,
  xlab = "",
  ylab = "Density",
  ylim = c(0, ymax),
  main = "",
  xaxt = "n",   # suppress x-axis
  yaxt = "n"
)

abline(h = 0, col = "black", lwd = 1)

for (i in seq_along(bin_left)) {
  rect(
    xleft   = bin_left[i],
    ybottom = 0,
    xright  = bin_right[i],
    ytop    = h.res$density[i],
    col     = if (dark_bins[i]) col2 else col1,
    border  = "black"
  )
}

axis(2, at = pretty(c(0, ymax)))
axis(1, line = 2)  # draw x-axis lower

######################################### 
######################################### 

if (Show.Text) {

#par(family = fam[2])


par(new=TRUE)
hist(c(0),main="",ylim=c(ymin,ymax),ylab="",xlab="",xlim=c(x.lim.min,x.lim.max),
	density=0,border="white",axes=FALSE)

	min.z = min(val.input)
	max.z = max(val.input)
	n.z = length(val.input)
	n.z.sig = length(val.input[val.input > crit])
	n.not.shown = length(val.input[val.input > x.lim.max])

	### set location parameters for writing text
	y.line.factor 
	y.text = ymax
	y.line = ymax/100*y.line.factor

	### write copyright 
	i = 0
	text(x.lim.min,ymax-y.line*i,"Schimmack, Bartos, Brunner",pos=4,cex =letter.size)
	i = i + 1.5
	text(x.lim.min,ymax-y.line*i,"z-curve 3.0 (2025)",pos=4,cex=letter.size)
	i = i + 1.5
	text(x.lim.min,ymax-y.line*i,version,pos=4,cex =letter.size)

	########

	#TESTING
	#y.text = ymax;Write.CI = TRUE;results=res.text;Show.Text = FALSE; Draw.Histogram(w.all,results=res.text,cola=col.hist)

	i = 0
	results.x = x.lim.max 

#	text(results.x,y.text-y.line*i,paste0("Range: ",
#		round(min.z,2)," to ",round(max.z,2)),
#		pos=2,cex=letter.size)
#	i = i + 2

	text(results.x+0.1,y.text-y.line*i,paste0(toString(format(n.z,big.mark=",")),
		" tests "),pos=2,cex = letter.size)
	i = i + 1.5

	if (substring(Est.Method,1,3) == "CLU") {
		text(results.x+0.1,y.text-y.line*i,paste0(toString(format(cluster.k,big.mark=",")),
			" articles "),pos=2,cex = letter.size)
		i = i + 1.5
	}

	text(results.x+0.1,y.text-y.line*i,paste0(toString(format(n.z.sig,big.mark=","))," significant "),
	  pos=2,cex = letter.size)
	i = i + 1.5
	text(results.x,y.text-y.line*i,paste0(n.not.shown," z > ",x.lim.max),
		pos=2,cex = letter.size)
	

#	par(family = fam[1])

	######

	#print("Check write CI")
	#print(Write.CI)

	#############################################
	#############################################


	if (Write.CI) {


	### results = res.text

	p.bias = results[5,3]
	if (is.na(p.bias)) TEST4BIAS = FALSE

	if(TEST4BIAS) { 
		if (p.bias < .00005) { bias.res = "EJS, p < .0001" } else {  
			bias.res = paste0("EJS, p = ",sub("^0","",formatC(p.bias,format="f",digits=4))) }
		bias.res
	}	

	results = round(results,2)
	results = round(results*100);results

	ODR = results[1,1]
	ODR.low = results[1,2]
	ODR.high = results[1,3]
	EDR = results[2,1]
	EDR.low = results[2,2]
	EDR.high = results[2,3]
	ERR = results[3.1]
	ERR.low = results[3,2]
	ERR.high = results[3,3]
	FDR = results[4.1]
	FDR.low = results[4,2]
	FDR.high = results[4,3]

	if (ODR.high > 100) ODR.high = 100
	#print("Check ERR.low")
	#print(ERR.low)
	if(ERR.low < 2.5) ERR.low = 2.5
	if(EDR.low < 5) EDR.low = 5
	
	i = i + 3
	text(results.x,y.text-y.line*i,
		paste0(round((1-CI.ALPHA)*100),"% CI     "),
		pos=2,cex=letter.size.1)

	x.loc = results.x + c(-1.5,-1.1,-1.0,-.6,-.5,-.1,0)

	ODR.a = str_pad(round(ODR), width = 3, pad = " ");ODR.a
	ODR.b = str_pad(round(ODR.low), width = 3, pad = " ");ODR.b
	ODR.c = str_pad(round(ODR.high), width = 3, pad = " ");ODR.c

	EDR.a = str_pad(round(EDR), width = 3, pad = " ");EDR.a
	EDR.b = str_pad(round(EDR.low), width = 3, pad = " ");EDR.b
	EDR.c = str_pad(round(EDR.high), width = 3, pad = " ");EDR.c

	ERR.a = str_pad(round(ERR), width = 3, pad = " ");ERR.a
	ERR.b = str_pad(round(ERR.low), width = 3, pad = " ");ERR.b
	ERR.c = str_pad(round(ERR.high), width = 3, pad = " ");ERR.c

	FDR.a = str_pad(round(FDR), width = 3, pad = " ");FDR.a
	FDR.b = str_pad(round(FDR.low), width = 3, pad = " ");FDR.b
	FDR.c = str_pad(round(FDR.high), width = 3, pad = " ");FDR.c

	res4text = c()
	res4text = rbind(res4text,c("ODR: ","EDR: ","ERR :","FDR :"))
	res4text = rbind(res4text,c(ODR.a,EDR.a,ERR.a,FDR.a))
	res4text = rbind(res4text,rep("[",4))
	res4text = rbind(res4text,c(ODR.b,EDR.b,ERR.b,FDR.b))
	res4text = rbind(res4text,rep(",",4))
	res4text = rbind(res4text,c(ODR.c,EDR.c,ERR.c,FDR.c))
	res4text = rbind(res4text,rep("]",4))

	dim(res4text)
	res4text

	x.i = 1
	for (x.i in 1:length(x.loc)) {
		for (ii in 1:4) {
			text(x.loc[x.i],y.text-y.line*(i+2*ii),res4text[x.i,ii]
			,pos=2,cex=letter.size.1)
		}
	}

	if(TEST4BIAS) {
		i = i + 10
		text(results.x,y.text-y.line*i,bias.res,
			pos=2,cex=letter.size.1)
		}

	##############################################
	# End of IF CI show results with CI
	} else {  
	# End of IF CI show results without CI

	p.bias = results[5]

	results = round(results*100);results

	ODR = results[1];ODR
	EDR = results[2];ERR
	ERR = results[3];EDR
	FDR = results[4];FDR

	if (is.na(p.bias)) TEST4BIAS = FALSE

	if(TEST4BIAS) { 
		if (p.bias < .00005) { bias.res = "EJS, p < .0001" } else { 
			bias.res = paste0("EJS, p = ",sub("^0","",formatC(p.bias,format="f",digits=4))) }
		bias.res
	}	

	# Witout CI
	#print("Writing without CI")

	i = i + 3
	ODR.string = str_pad(round(ODR), width = 3, pad = " ");ODR.string
	text(results.x,y.text-y.line*i,
		paste("ODR:",ODR.string,"%"),
		pos=2,cex=letter.size.1)

	i = i + 2
	EDR.string = str_pad(round(EDR), width = 3, pad = " ");EDR.string
	text(results.x,y.text-y.line*i,
		paste("EDR:",EDR.string,"%"),
		pos=2,cex=letter.size.1)

	i = i + 2
	ERR.string = str_pad(round(ERR), width = 3, pad = " ");ERR.string
	text(results.x,y.text-y.line*i,
		paste("ERR:",ERR.string,"%"),
		pos=2,cex=letter.size.1) 

	i = i + 2
	FDR.string = str_pad(round(FDR), width = 3, pad = " ");FDR.string
	text(results.x,y.text-y.line*i,
		paste("FDR:",FDR.string,"%"),
		pos=2,cex=letter.size.1) 

	if (is.na(p.bias)) TEST4BIAS = FALSE

	if(TEST4BIAS) {
		i = i + 2
		text(results.x,y.text-y.line*i,bias.res,
			pos=2,cex=letter.size.1)
	}


#	i = i + 2
#	text(results.x,y.text-y.line*i,
#		paste0("OSR:   ",OSR,"%"),
#		pos=2,cex=letter.size.1)

	i = i + 2

	if (Est.Method == "OF" & boot.iter > 0 & Write.CI == FALSE) text(results.x,ymax-y.line*i,
		"WAIT FOR BOOSTRAPPED CIs",
		pos=2,cex=letter.size,col="red")

	i = i + 2
	if (TEST4HETEROGENEITY) text(results.x,ymax-y.line*i,
		"WAIT FOR HETEROGENEITY BOOSTRAPPED CIs",
		pos=2,cex=letter.size,col="red")



	} # End of CI

} # End of Show.Text


### restore default graphic parameters
par(op)


} # End of Histogram
 
######################################### 
######################################### 
######################################### 

######################################### 
######################################### 
######################################### 


Draw.KD = function(z.draw,w,Write.CI=FALSE,cola="blue",Lwidth=5) {

	#ymin=-.015;ymax = .6
	#Draw.Histogram(w.all,results=res,cola=col.hist)

	#z.draw = val.input
	#cola = "blue"

	x.adj.max = x.lim.max + 3 * bw.draw

	### densitiy line

	summary(z.draw)

	d.all = Get.Densities(z.draw[z.draw >= x.lim.min & z.draw < x.lim.max],
		bw=bw.draw,d.x.min=x.lim.min,d.x.max=x.lim.max,Augment=Augment)
	summary(d.all)
	bar.width = d.all[2,1] - d.all[1,1];bar.width
	d.all = d.all[d.all[,1] > x.lim.min & d.all[,1] <= x.lim.max,]
	dim(d.all)
	sum(d.all[,2])*bar.width
	d.all.X = d.all[,1]
	d.all.Y = d.all[,2]
	summary(d.all.X)


	summary(z.draw)

	d.sig = Get.Densities(z.draw[z.draw >= Int.Beg & z.draw < x.lim.max],
		bw=bw.draw,d.x.min=Int.Beg,d.x.max=x.lim.max,Augment=Augment)
	
	bar.width = d.sig[2,1] - d.sig[1,1];bar.width

	dim(d.sig)
	d.sig = d.sig[d.sig[,1] > Int.Beg & d.sig[,1] <= x.lim.max,]

	#dim(d.sig)
	#sum(d.sig[,2])*bar.width

	d.sig.X = d.sig[,1]
	d.sig.Y = d.sig[,2]
	#summary(d.sig.X)

	d.fit = length(z.draw[z.draw > Int.Beg & z.draw < x.lim.max]) /
		length(z.draw[z.draw > x.lim.min & z.draw < x.lim.max])

	d.fit

	d.sig.Y = d.sig.Y * d.fit

	sum(d.sig[,2])*bar.width
	
#	### draw the density of the observed values in the selected region
#	par(new=TRUE)
#	plot(d.all.X[d.all.X > x.lim.min+.05 & d.all.X < x.lim.max - .05],
#		d.all.Y[d.all.X > x.lim.min+.05 & d.all.X < x.lim.max - .05],
#		type="l",col=adjustcolor(cola, alpha.f = 0.3),lty=1,lwd=Lwidth,
#		xlim =c(x.lim.min,x.lim.max),ylim=c(ymin,ymax),xlab="",ylab="")


	### draw the density of the observed values in the selected region
	par(new=TRUE)
	plot(d.sig.X[d.sig.X > Int.Beg +.05 & d.sig.X < x.lim.max - .05],
		d.sig.Y[d.sig.X > Int.Beg + .05 & d.sig.X < x.lim.max - .05],
		type="l",col=adjustcolor(cola, alpha.f = 0.4),lty=1,lwd=Lwidth,
		xlim =c(x.lim.min,x.lim.max),ylim=c(ymin,ymax),xlab="",ylab="",
		axes=FALSE)
	### draw vertical line for Beginning

	Int.Beg
	height = max(d.sig.Y[which(round(d.sig.X,1) == round(Int.Beg+.02,1))]);height
	segments(Int.Beg+.02,0,Int.Beg+.02,height,lty=1,lwd=3,col=cola)




} # End of Draw.KD

###################################################
#### End Draw.KD
###################################################



###################################################
#### Begin Draw Zcurve
###################################################


Draw.Curve.All = function(w,cola="black",Lwidth=4,
	Ltype=1,x.start=x.lim.min,x.end=x.lim.max) {


#Draw.Histogram()

bar.width = .01
D.X = seq(x.lim.min,x.lim.max,bar.width)

if (CURVE.TYPE == "z") {
  Dens = get.z.density(D.X)
} else {
  Dens = get.t.density(D.X)
}

D.Y = colSums(Dens*w)

d.dense.sel = sum(D.Y[D.X > x.lim.min & D.X >= Int.Beg]*bar.width)
d.dense.sel

d.hist.sel = mean(as.numeric(val.input > x.lim.min & val.input > Int.Beg & val.input < Int.End)) /
	mean(as.numeric(val.input > x.lim.min & val.input < Int.End))
d.hist.sel

check = D.Y[D.X > crit & D.X < Int.Beg]
length(check)
EPJS = sum(D.Y[D.X > crit & D.X < Int.Beg])/sum(D.Y[D.X > crit])
EPJS

OFJS = sum(as.numeric(val.input > crit & val.input < Int.Beg))
OFJS

OFS = sum(as.numeric(val.input > crit & val.input < Int.End))
OFS

OPJS = OFJS/OFS
OPJS



scale = d.hist.sel/d.dense.sel;scale

lines(D.X[which(D.X >= x.start & D.X < x.end)],
	D.Y[which(D.X >= x.start & D.X < x.end)]*scale,
	lty=Ltype,col=cola,lwd=Lwidth,xlim=c(x.lim.min,x.lim.max),ylim=c(ymin,ymax))

boundary = Int.Beg
if (round(Int.Beg,2) == 1.96) {boundary = 2}

segments(boundary,0,boundary,
	D.Y[which(round(D.X,2) == round(boundary,2))[1]]*scale,
	col=cola,lwd=3)


if (Show.Significance) {

	if (length(sig.levels) == 0) sig.levels = alpha

		sig.level = crit
		if(two.sided == FALSE) sig.level = c(-sig.level,sig.level)

		sig.lty = rep(2,length(sig.level));sig.lty

		i = 1
		for (i in 1:length(sig.level)) {
			height = max(scale*D.Y[which(round(D.X,1) == round(sig.level[i],1))]);height
			segments(sig.level[i],0,sig.level[i],height,lty=sig.lty,lwd=2,col="firebrick3")
		}

} # End of Show.Significance


} # End of Draw.Curve.All


###################################################
#### Begin Draw Zcurve SDG1
###################################################

Draw.Curve.All.SDG1 = function(w.inp,ncp=ncp,zsds=zsds,cola=col.curve,
	Ltype=1,Lwidth=4,x.start=x.lim.min,x.end=x.lim.max) {

#x.start = x.lim.min; x.end = x.lim.max;Ltype = 3;Lwidth = 4;cola = col.hist

#Draw.Histogram(w.all,results=res.text,cola=col.hist)

components = length(ncp)

if (components == 1) { 
	D.X = seq(x.lim.min,x.lim.max,.01)
	if (CURVE.TYPE == "z") {
		D.Y = dnorm(D.X,ncp,zsds)
	} else {
		D.Y = dt(D.Y,df,ncp)
	}
	d.sim = cbind(D.X,D.Y)
} else {
	k.sim = 500000
	sim.val = c()
	for (i in 1:components) {
		if (CURVE.TYPE == "z") { 
			sim.val = c(sim.val,rnorm(w.inp[i]*k.sim,ncp[i],zsds[i]))
		} else {
			sim.val = c(sim.val,rt(w.inp[i]*k.sim,df,ncp[i]))
		}
	}
	sim.val = abs(sim.val)
	d.sim = Get.Densities(sim.val[sim.val >= x.lim.min & sim.val < x.lim.max],
		bw=bw.draw,d.x.min=x.lim.min,d.x.max=x.lim.max,Augment=Augment)
}
dim(d.sim)

bar.width = d.sim[2,1] - d.sim[1,1];bar.width
d.sim = d.sim[d.sim[,1] > x.lim.min & d.sim[,1] <= x.lim.max,]
dim(d.sim)
d.sim[,2] = d.sim[,2]/(sum(d.sim[,2])*bar.width)
sum(d.sim[,2])*bar.width

d.sim.X = d.sim[,1]
d.sim.Y = d.sim[,2]

#plot(d.sim.X,d.sim.Y,ylim=c(ymin,ymax))

summary(d.sim.X)
sum(d.sim.Y*bar.width)

d.hist.sel = mean(as.numeric(val.input > x.lim.min & val.input > Int.Beg & val.input < Int.End)) /
	mean(as.numeric(val.input > x.lim.min & val.input < Int.End))
d.hist.sel

d.dense.sel = sum(d.sim.Y[d.sim.X > x.lim.min & d.sim.X >= Int.Beg]*bar.width)
d.dense.sel

scale = d.hist.sel/d.dense.sel;scale

lines(d.sim.X[which(d.sim.X >= x.start & d.sim.X < x.end)],
	d.sim.Y[which(d.sim.X >= x.start & d.sim.X < x.end)]*scale,
	lty=Ltype,col=cola,lwd=Lwidth,xlim=c(x.lim.min,x.lim.max)
	,ylim=c(ymin,ymax))


if (Show.Significance) {

	if (length(sig.levels) == 0) sig.levels = alpha

		cz.alpha = -qnorm(sig.levels/2);cz.alpha
		if(two.sided == FALSE) cz.alpha = c(-cz.alpha,cz.alpha);cz.alpha

		sig.lty = rep(2,length(cz.alpha));sig.lty

		i = 1
		for (i in 1:length(cz.alpha)) {
			cz.alpha[i]
			height = max(d.sim.Y[which(round(d.sim.X,1) == round(cz.alpha[i],1))]);height
			segments(cz.alpha[i],0,cz.alpha[i],height,lty=sig.lty,lwd=2,col="firebrick3")
		}

} # End of Show.Significance



} # End of Draw.Curve.SDG1

###################################################




##############################################
### Get Densities
##############################################

Get.Densities = function(zval,bw="nrd0",d.x.min=0,d.x.max=6,Augment=TRUE) {

#bw = bw.est; zval = val.input[val.input > Int.Beg];d.x.max = 6;d.x.min = Int.Beg
#hist(zval)
#summary(zval)

### find the maximum z-score. This is only needed if the maximum z-score is below Int.End
max.z = Int.End
if (max(zval) < d.x.max) max.z = max(zval)
max.z

#max.z = 6

z.check = length(val.input[val.input > round(Int.Beg,1)  
	& val.input < round(Int.Beg,1) + 1.5 + 4*bw.est])
z.check

if(Augment.Regression & z.check > 20) { 

Z.Density = bkde(val.input[val.input > round(Int.Beg,1)  
	& val.input < round(Int.Beg,1) + 2 + 4*bw.est],bandwidth=bw.est,
	range=c(Int.Beg,Int.Beg+2))
D = data.frame(Z.Density$x,Z.Density$y)
dim(D)
colnames(D) = c("ZX","ZY")
bw.reg = 1.5
D = D[D$ZX > Int.Beg & D$ZX < Int.Beg + bw.reg,]
dim(D)
summary(D$ZX)

#plot(D$ZX,D$ZY)

d.reg = -lm(scale(D$ZY) ~ scale(D$ZX))$coefficients[2];d.reg

if (Int.Beg > 1) Augment.Factor = 1 + d.reg;Augment.Factor

} 


### Augment z-scores on the left side of Interval to avoid downward trend 
### of kernal density function (assumes values go to 0)

if (Augment == TRUE) { 

	AUG = c()
	n.AUG = round(Augment.Factor*length(zval[zval > d.x.min & zval < d.x.min+bw.aug]));n.AUG
	if (n.AUG > 0) AUG = seq(d.x.min-bw.aug,d.x.min,bw.aug/n.AUG)

	summary(zval)
	Z.INT.USE = c(zval,AUG)
	#hist(zval[zval & 1.9 & zval < 6],breaks=seq(1.9,6.1,.1),freq=FALSE)
	#hist(Z.INT.USE[Z.INT.USE < 6],breaks=seq(1.9,6,.1),freq=FALSE)
	#summary(AUG)

} else {
	#print("Not Augmenting")
	Z.INT.USE = zval[zval > d.x.min & zval <= max.z]
}

#print(summary(Z.INT.USE))

Z.Density = bkde(Z.INT.USE,bandwidth=bw,range=c(d.x.min-bw.aug,d.x.max)) 
val.max = d.x.max
D = data.frame(Z.Density$x,Z.Density$y)
colnames(D) = c("ZX","ZY")
#plot(D$ZX,D$ZY);abline(v=d.x.min)


D = D[D$ZX > d.x.min & D$ZX < val.max,]
dim(D)
#plot(D$ZX,D$ZY)

bar.width = D$ZX[2] - D$ZX[1]
D$ZY = D$ZY/(sum(D$ZY*bar.width)) 
sum(D$ZY*bar.width)

#print(summary(D))

#print("End Augment check")

return(D)


}  ### End of Get Densities 


#######################################################
### End of Get Densities
#######################################################



#######################################################
### Old Fashioned Density Method
#######################################################

#weight = startval

old.fashioned = function(val.input,cola = "springgreen2") {

curve.fitting = function(theta,RetEst=FALSE)    {

### get the weights and rescale 
weight = theta
weight = weight/sum(weight)

if (fixed.false.positives > 0) weight = c(fixed.false.positives,weight*(1-fixed.false.positives))
sum(weight)

### compute the new estimated density distribution
E.D.Y = colSums(Dens*weight)

### compare to observed density distribution
rmse = sqrt(mean((E.D.Y-O.D.Y)^2))

### return either fit if continue or estimates if finished
value = rmse
if(RetEst) value = E.D.Y


### showing the fitting of the function in a plot
if(Plot.Fitting) {

	rval = runif(1)
	if (rval > .9) {

	tit = ""
	xL = ""
	yL = ""
	plot(D.X,O.D.Y*scale,type='l',
		xlim=c(x.lim.min,x.lim.max),ylim=c(ymin,ymax),
		main=tit,xlab=xL,ylab=yL,axes=FALSE)
	lines(D.X,E.D.Y*scale,lty=1,col="red1",
		xlim=c(x.lim.min,x.lim.max),ylim=c(ymin,ymax),axes=FALSE)
	#points(D.X,z.est,pch=20,col="red1",ylim=c(0,ymax),)

	}

} ### End of Plot Fitting

### return value to optimization function
return(value)

} ### End of curve.fitting

############################
#### End of Fitting Function
############################

components = length(ncp)

### get the densities for each interval and each non-centrality parameter

INT = val.input[val.input >= Int.Beg & val.input <= Int.End]

densy = Get.Densities(INT,bw=bw.est,d.x.min=Int.Beg,d.x.max=Int.End,Augment=Augment)

D.X = densy[,1]
O.D.Y = densy[,2]

#plot(D.X,O.D.Y,type="l")

n.bars = length(D.X)

bar.width = D.X[2] - D.X[1]
bar.width

### Finish getting observed densities 

#print("Create Dens")
Dens	= c()
for(i in 1:n.bars) {
	for (j in 1:length(ncp)) {
		if (CURVE.TYPE == "z") {
			Dens = c(Dens,
			dnorm(D.X[i],ncp[j],zsds[j]) +
			dnorm(-D.X[i],ncp[j],zsds[j])
			)
		} else {
			Dens = c(Dens,
			dt(D.X[i],df,ncp[j]) +
			dt(-D.X[i],df,ncp[j])
			)
		}
	}
}
summary(Dens)
Dens = matrix(Dens,length(ncp),byrow=FALSE)
sum.dens = rowSums(Dens)
Dens = Dens/(sum.dens * bar.width)
dim(Dens)

startval = rep(1/(components),components)
#startval[1] = 1
startval = startval/sum(startval)
startval

lowlim = c(rep(0,components),ncp,zsds);lowlim
highlim = c(rep(1,components),ncp,zsds);highlim

if (fixed.false.positives > 0 & 0 %in% ncp) {
  startval = rep(1/(components-1),components-1)
  lowlim = c(rep(0,components-1),ncp,zsds);lowlim
  highlim = c(rep(1,components-1),ncp,zsds);highlim
}


d.hist.sel = mean(as.numeric(val.input > x.lim.min & val.input > Int.Beg & val.input < Int.End)) /
	mean(as.numeric(val.input > x.lim.min & val.input < Int.End))
d.hist.sel

scale = d.hist.sel

#TESTING = TRUE
if (TESTING) Plot.Fitting = TRUE else Plot.Fitting = FALSE


auto = nlminb(startval,curve.fitting,lower=lowlim,upper=highlim,
	control=list(eval.max=1000,abs.tol = 1e-20))

fit = auto$objective;fit

### get the estimated weights 
WT = auto$par
WT = WT/sum(WT)
if (fixed.false.positives > 0) WT = c(fixed.false.positives,WT*(1-fixed.false.positives))
sum(WT)

res = c(WT,fit);res

return(res)

} ### End of old.fashioned

######################################################
### End of Old Fashioned Density Method
#######################################################







################################################################
### USE EXTENDED CURVE (Est.Method = "EXT"  (slower than OF)
################################################################

extended.curve = function(val.input,ncp=ncp,zsds=zsds) {

#w.inp = 1;ncp = 4;zsds = 2;
#weights = 1;means = 0;sds = 1;

### this is the actual fitting function
extended.curve.fitting = function(para,RetEst=FALSE,Fixed.Null=FALSE)    {

### get the weights
weights = para[1:components]
means = para[(components+1):(2*components)]
sds = para[(components*2+1):(3*components)]

i = 1
j = 1
k = 1
### get the densities for each interval and each non-centrality parameter
Dens	= c()
for(i in 1:n.bars) {
	for (j in 1:components) {
		if (CURVE.TYPE == "z") {
			Dens = c(Dens,
			dnorm(D.X[i],means[j],sds[j]) +
			dnorm(-D.X[i],means[j],sds[j]) 
			)
		} else {
			sds = sds - 1
			sds[sds <= 0] = 0
			sds
			if (max(sds) > 0) {
				Dens = c(Dens,
				dt(D.X[i],df,means[j]) +
				dt(-D.X[i],df,means[j])
				)
			} else { 
				Dens = c(Dens,
				dt(D.X[i],df,means[j]) +
				dt(-D.X[i],df,means[j])
				)
			}
		}
	}
}

Dens = matrix(Dens,length(ncp),byrow=FALSE)

### rescale the densities for the range of z-values to 1
row.sum.dens = rowSums(Dens)
Dens = Dens/(row.sum.dens * bar.width)
dim(Dens)
summary(Dens)


### compute the new estimated density distribution
E.D.Y = colSums(Dens*weights)

### compare to observed density distribution
rmse = sqrt(mean((E.D.Y-O.D.Y)^2));rmse


### return either fit if continue or estimates if finished
value = rmse
if(RetEst) value = est


### showing the fitting of the function in a plot
if(Plot.Fitting) {

	rval = runif(1)
	if (rval > .4) {

	tit = ""
	xyy = cbind(D.X,O.D.Y,E.D.Y)
	plot(xyy[,1],xyy[,2],type='l',
		xlim=c(x.lim.min,x.lim.max),ylim=c(ymin,ymax),
		main=tit,xlab='Z')
	lines(xyy[,1],xyy[,3],lty=1,col="red1",
		xlim=c(x.lim.min,x.lim.max),ylim=c(ymin,ymax))
	}

}


### return value to optimization function
return(value)

} 

####

### create set with z-scores in the interval used for model fitting
INT = val.input[val.input >= Int.Beg & val.input <= Int.End]
summary(INT)

components = length(ncp);components

#print(Augment)
densy = Get.Densities(INT,bw=bw.est,d.x.min=Int.Beg,d.x.max=Int.End,Augment=Augment)

D.X = densy[,1]
O.D.Y = densy[,2]

#plot(D.X,D.Y)

n.bars = length(D.X)
n.bars

bar.width = D.X[2] - D.X[1]
bar.width


if (W.FIXED) {
	startval = w.fix
	lowlim = w.fix
	highlim = w.fix
} else {
	startval = rep(1/components,components)
	lowlim = rep(0,components)
	highlim = rep(1,components)
}

if (NCP.FIXED) {
	startval = c(startval,ncp)
	lowlim = c(lowlim,ncp)
	highlim = c(highlim,ncp)
} else {
	startval = c(startval,ncp)
	lowlim = c(lowlim,rep(0,components))
	highlim = c(highlim,rep(Int.End,components))
}

if (ZSDS.FIXED) { 
	startval = c(startval,zsds)
	lowlim = c(lowlim,zsds)
	highlim = c(highlim,zsds)
} else {
	startval = c(startval,zsds)
	lowlim = c(lowlim,rep(1,components))
	highlim = c(highlim,zsds)
}

if(components == 1) lowlim[1] = 1

#print(startval);print(lowlim);print(highlim)

#TESTING = TRUE
if (TESTING == TRUE) Plot.Fitting = TRUE

auto = nlminb(startval,extended.curve.fitting,lower=lowlim,upper=highlim,control=list(eval.max=1000))

para = auto$par
para[1:components] = para[1:components]/sum(para[1:components])
para

fit = auto$objective;fit

res = c(para,fit)
res

return(c(res))

} 

######################################################
### End of Extended Zcurve
#######################################################


#######################################################
### Begin Old Fashioned Zcurve (Est.Method = "OF" 
#######################################################


#########################################################################
### This Function Computes Power from Weights and Non-Centrality Parameters
#########################################################################

Compute.Power.Z = function(para,Int.Beg=crit,BOOT=FALSE) {

crit = qnorm(alpha/2,lower.tail=FALSE)

#para = para.est.OF
#para = res.1$w.all

ext.all = length(val.input[val.input > Int.End]) / 
	length(val.input)
ext.all

ext.sig = length(val.input[val.input > Int.End]) / 
	length(val.input[val.input > crit])
ext.sig

ext.inp = length(val.input[val.input > Int.End]) / 
	length(val.input[val.input > Int.Beg])
ext.inp


### the input weights based on z.curve method
### these are the weights based on the Int.Beg value
### Int.Beg could be 1.96 (all significant)
### but it can also be other values
### Therefore the weights cannot be directly used
### to estimate power/replicability
w.inp = para[1:components]
ncp = para[(1+components):(2*components)]
zsds = para[(1+2*components):(3*components)]

### get power values for the components (ncp)
pow.dir = pnorm(abs(ncp),crit);pow.dir 

### get the opposite sign probability
sign.error = 1-pnorm(ncp,-crit);round(sign.error,3)

pow = pow.dir + sign.error
pow.ext = c(pow,1)

### this gives the power with the Int.Beg selection criterion as alpha 
### this power is used as a weight to get the weights for the full distribution
### using Jerry's insight that weight before selection is weight after selection divided by power
### get power values for the components (ncp)
pow.sel = pnorm(ncp,Int.Beg) + pnorm(-ncp,Int.Beg);pow.sel
pow.sel.ext = c(pow.sel,1)

w.inp.ext = c(w.inp*(1-ext.inp), ext.inp)          # estimated weights

w.all.ext   <- (w.inp.ext + w.inp.ext*(1-pow.sel.ext)
   /(pow.sel.ext)) / sum(w.inp.ext + w.inp.ext*(1-pow.sel.ext)
   /(pow.sel.ext))

EDR = sum(w.all.ext*pow.ext);EDR

w.sig.ext = w.all.ext * pow.ext
w.sig.ext = w.sig.ext / EDR
sum(w.sig.ext)

pow.dir.ext = c(pow.dir,1)

w.all = w.all.ext[1:components]
w.sig = w.sig.ext[1:components]


ERR = sum(w.sig.ext*pow.dir.ext);ERR

if (ERR > 1) ERR = 1
if (EDR > 1) EDR = 1

if (ERR < alpha/2) ERR = alpha/2
if (EDR < alpha) EDR = alpha 


ERR.pos = NA
ERR.neg = NA

EDR.pos = NA
EDR.neg = NA

res.est = c(EDR,EDR.pos,EDR.neg,ERR,ERR.pos,ERR.neg)
res.est
res = c(res.est,w.all,w.sig)
names(res) = c("EDR","EDR.pos","EDR.neg","ERR","ERR.pos","ERR.neg",
paste0("w.all.",ncp[1:components]),paste0("w.sig",ncp[1:components]) )
res

### now we compute mean power as a function of z-scores continuously
### this is only performed if local power is requested (int.loc > 0)
if (int.loc > 0) {

	bar.width = .01 # how fine should be the resolution
	X = seq(x.lim.min,Int.End,bar.width);summary(X)	 # set of z-scores 

	W.D.Sum = unlist(lapply(X,function(x) sum(dnorm(x,ncp)*w.all) ))
	#plot(Z.X,Z.W.D.Sum)

	loc.p = unlist(lapply(X,function(x) 
		sum(dnorm(x,ncp)*w.all*(pow.dir+sign.error)) /	sum(dnorm(x,ncp)*w.all)	))

	int = seq(x.lim.min,x.lim.max,int.loc)
	local.power = c()
	i = 1
	for (i in 1:(length(int)-1)) local.power = c(local.power,
		sum(loc.p[X > int[i] & X < int[i+1]]*
			W.D.Sum[X > int[i] & X < int[i+1]])/
		sum(W.D.Sum[X > int[i] & X < int[i+1]])	 )		

	local.power
	names(local.power) = paste0("lp.",1:length(local.power))
	res = c(res,local.power)
}


#print("Check res local power")
#print(length(res))
#print(res)

### to be past back to the main program from this function
return(res)

}

### Finish Compute.Power.Z

#####################################
#####################################
#####################################


Compute.Power.T = function(para,Int.Beg=crit,BOOT=FALSE) {

#para = para.est.OF

t.crit = qt(1-alpha/2,df);t.crit

t.ext.all = length(val.input[val.input > Int.End]) / 
	length(val.input)
t.ext.all

t.ext.sig = length(val.input[val.input > Int.End]) / 
	length(val.input[val.input > crit])
t.ext.sig

t.ext.inp = length(val.input[val.input > Int.End]) / 
	length(val.input[val.input > Int.Beg])
t.ext.inp

w.inp = para[1:components]
nct = para[(1+components):(2*components)]

### get power values for the components (ncp)
pow.dir = pt(t.crit,df,nct,lower.tail=FALSE);pow.dir 

### get the opposite sign probability
sign.error = pt(-t.crit,df,nct,lower.tail=TRUE);round(sign.error,3)

pow = pow.dir + sign.error
pow.ext = c(pow,1)

pow.sel = pt(Int.Beg,df,nct,lower.tail=FALSE) + pt(-Int.Beg,df,nct);pow.sel
pow.sel.ext = c(pow.sel,1)
pow.sel

w.inp.ext = c(w.inp*(1-t.ext.inp), t.ext.inp)          # estimated weights
w.inp.ext

w.all.ext   <- (w.inp.ext + w.inp.ext*(1-pow.sel.ext)
   /(pow.sel.ext)) / sum(w.inp.ext + w.inp.ext*(1-pow.sel.ext)
   /(pow.sel.ext))

EDR = sum(w.all.ext*pow.ext);EDR

w.sig.ext = w.all.ext * pow.ext
w.sig.ext = w.sig.ext / EDR
sum(w.sig.ext)

pow.dir.ext = c(pow.dir,1)

w.all = w.all.ext[1:components]
w.sig = w.sig.ext[1:components]

ERR = sum(w.sig.ext*pow.dir.ext);ERR

if (ERR > 1) ERR = 1
if (EDR > 1) EDR = 1

if (ERR < alpha/2) ERR = alpha/2
if (EDR < alpha) EDR = alpha 


ERR.pos = NA
ERR.neg = NA

EDR.pos = NA
EDR.neg = NA

res.est = c(EDR,EDR.pos,EDR.neg,ERR,ERR.pos,ERR.neg)
res.est
res = c(res.est,w.all,w.sig)
names(res) = c("EDR","EDR.pos","EDR.neg","ERR","ERR.pos","ERR.neg",
paste0("w.all.",ncp[1:components]),paste0("w.sig",ncp[1:components]) )
res


if (int.loc > 0) {

	bar.width = .01 # how fine should be the resolution
	Z.X = seq(x.lim.min,Int.End,bar.width);summary(Z.X)	 # set of z-scores 

	Z.W.D.Sum = unlist(lapply(Z.X,function(x) sum(dt(x,df,ncp)*w.all) ))
	#plot(Z.X,Z.W.D.Sum)

	loc.p = unlist(lapply(Z.X,function(x) 
		sum(dt(x,df,ncp)*w.all*(pow.dir+sign.error)) /	
			sum(dt(x,df,ncp)*w.all)	))

	loc.p
	int = seq(x.lim.min,x.lim.max,int.loc)
	local.power = c()
	i = 1
	for (i in 1:(length(int)-1)) local.power = c(local.power,
		sum(loc.p[Z.X > int[i] & Z.X < int[i+1]]*
			Z.W.D.Sum[Z.X > int[i] & Z.X < int[i+1]])/
		sum(Z.W.D.Sum[Z.X > int[i] & Z.X < int[i+1]])	 )		

	local.power
	names(local.power) = paste0("lp.",1:length(local.power))
	res = c(res,local.power)
}

print(res)

return(res)

} # EOF Compute.Power.T


#####################################
#####################################
#####################################




###############################################################

Compute.Power.SDG1 = function(para,BOOT=FALSE,Int.Beg=Int.Beg) {

ext.all = length(val.input[val.input > Int.End]) / 
	length(val.input)
ext.all

ext.sig = length(val.input[val.input > Int.End]) / 
	length(val.input[val.input > crit])
ext.sig

ext.inp = length(val.input[val.input > Int.End]) / 
	length(val.input[val.input > Int.Beg])
ext.inp


#para = para.est.EXT
#para = c(w.inp,ncp,zsds)

components = length(ncp)

w.inp = para[1:components];w.inp
if (components == 1) w.all = w.inp
if (Int.Beg == 0) w.all = w.inp;w.all
est.cw.all = w.all

est.ncp = para[(1+components):(2*components)];est.ncp
est.zsds = para[(1+2*components):(3*components)];est.zsds
est.zsds[est.zsds < 1] = 1  
est.zncp.sd = sqrt(est.zsds^2 - 1);est.zncp.sd
est.components = length(est.ncp)

#curve(dnorm(x,est.ncp,est.zncp.sd),0,6)

###

zx.bw = .01
zx = seq(0,Int.End,zx.bw)
pow.zx.dir = pnorm(abs(zx),1.96) 
pow.zx.sign.error = + pnorm(-1.96,abs(zx))
pow.zx = pow.zx.dir + pow.zx.sign.error

i = 1
est.wd.all = c()
for (i in 1:est.components) {
	if (est.zncp.sd[i] == 0) {
		wd = rep(0,length(zx))
		wd[which(round(zx,2) == round(est.ncp[i],2))] = 1
	} else {
		wd = dnorm(zx,est.ncp[i],est.zncp.sd[i])
	}
	wd = wd/sum(wd)
	wd = wd*est.cw.all[i]
	sum(wd)
	est.wd.all = rbind(est.wd.all,wd)
}

dim(est.wd.all)

est.wd.all = colSums(est.wd.all)
est.wd.all = est.wd.all/sum(est.wd.all)
sum(est.wd.all)

#plot(zx,est.wd.all)

est.wd.all.ext = c(est.wd.all*(1-ext.all),ext.all)
pow.zx.dir.ext = c(pow.zx.dir,1)

#cbind(est.wd.all.ext,pow.zx.dir.ext)

pow.zx.ext = c(pow.zx,1)
EDR = sum(pow.zx.ext*est.wd.all.ext);EDR

est.wd.sig.ext = est.wd.all.ext*pow.zx.ext
est.wd.sig.ext = est.wd.sig.ext/sum(est.wd.sig.ext)

ERR = sum(pow.zx.dir.ext*est.wd.sig.ext)/sum(est.wd.sig.ext);ERR

if (ERR > 1) ERR = 1
if (EDR > 1) EDR = 1

if (ERR < alpha/2) ERR = alpha/2
if (EDR < alpha) EDR = alpha 

ERR
EDR

EDR.pos = NA
EDR.neg = NA

ERR.pos = NA
ERR.neg = NA

res.est = c(EDR,EDR.pos,EDR.neg,ERR,ERR.pos,ERR.neg)
res.est
res = c(res.est,rep(NA,components),rep(NA,components))
names(res) = c("EDR","EDR.pos","EDR.neg","ERR","ERR.pos","ERR.neg",
paste0("w.all.",ncp[1:components]),paste0("w.sig.",ncp[1:components]))
res

if (int.loc > 0 & BOOT==FALSE) {

	zx.bw = .01
	zx = seq(x.lim.min,x.lim.max,zx.bw)
	pow.zx.dir = pnorm(abs(zx),1.96) 
	pow.zx.sign.error = + pnorm(-1.96,abs(zx))
	pow.zx = pow.zx.dir + pow.zx.sign.error
	
	i = 1
	est.wd.all = c()
	for (i in 1:est.components) {
		if (est.zncp.sd[i] == 0) {
			wd = rep(0,length(zx))
			wd[which(round(zx,2) == round(est.ncp[i],2))] = 1
		} else {
			wd = dnorm(zx,est.ncp[i],est.zncp.sd[i])
		}
		wd = wd/sum(wd)
		wd = wd*est.cw.all[i]
		sum(wd)
		est.wd.all = rbind(est.wd.all,wd)
	}

	dim(est.wd.all)

	est.wd.all = colSums(est.wd.all)
	est.wd.all = est.wd.all/sum(est.wd.all)
	sum(est.wd.all)

	#plot(zx,est.wd.all)

	sim.z = seq(x.lim.min,x.lim.max,.01)
	x = sim.z[1]
	loc.p = unlist(lapply(sim.z,function(x) 
		sum(pow.zx*dnorm(x,est.wd.all))/sum(dnorm(x,est.wd.all)) 
	))

	loc.p


	### compute local mean power for different intervals
	### each local power value is weighted by the density 
	int = seq(x.lim.min,x.lim.max,int.loc)
	local.power = c()
	i = 1
	for (i in 1:(length(int)-1)) local.power = c(local.power,
		sum(loc.p[zx > int[i] & zx < int[i+1]]*
			est.wd.all[zx > int[i] & zx < int[i+1]])/
		sum(est.wd.all[zx > int[i] & zx < int[i+1]])	 )		
	names(local.power) = paste0("lp.",1:length(local.power))
	local.power
	res = c(res,local.power)

} #EOF int.loc


res

#print("Compute Extended Power End")

### to be past back to the main program from this function
return(res)

} ### EOF Compute.Power for SDG1


#####################################
#####################################
#####################################



#####################################
### BBB ZingStart #START #Begin of Main Program 
### BBB ZingStart #START #Begin of Main Program 
### BBB ZingStart #START #Begin of Main Program 
### BBB ZingStart #START #Begin of Main Program 
### BBB ZingStart #START #Begin of Main Program 
### BBB ZingStart #START #Begin of Main Program 
#####################################


if (length(p) > 0) {
  val.input = qnorm(1-p/2)
}

if (substring(Est.Method,1,3) == "CLU") {
	if (length(cluster.id > 0)) {
	   clu.inp = data.frame(val.input,cluster.id)
	   clu.inp = clu.inp[complete.cases(clu.inp),]
	   z.clu.input = zcurve_data(paste("z = ", clu.inp[,1]), id = clu.inp[,2])
	   tab = table(cluster.id);tab
	   cluster.k = length(tab);cluster.k
	} else{
       print("Cluster IDs are missing!")
	}
} 

if (CURVE.TYPE == "z") { 
  crit = qnorm(alpha/2,lower.tail=FALSE)
 } else {
  crit = qt(alpha/2,df,lower.tail=FALSE)
}
crit

if (length(Int.Beg) == 0) Int.Beg = crit


components = length(ncp)


### remove NA
val.input = val.input[!is.na(val.input)];length(val.input)


if (two.sided) val.input = abs(val.input)

### limit range
val.input[val.input > MAX.INP.Z] = MAX.INP.Z
length(val.input)

if (length(val.input) > 10) print("!START!") else print("Insufficient Data")

print(paste0("Title: ",Title))

### create set with z-scores in the interval used for model fitting
INT = val.input[val.input >= Int.Beg & val.input <= Int.End]
summary(INT)

n.z = length(val.input);n.z
n.z.sig = length(val.input[abs(val.input) > crit])

ODR = n.z.sig/n.z;ODR

ODR.se = sqrt(ODR * (1 - ODR) / n.z)
ODR.se

ODR.low = round((ODR-1.96*ODR.se),3);ODR.low
ODR.high = round((ODR+1.96*ODR.se),3);ODR.high


ODR.res = c(ODR,ODR.low,ODR.high)
ODR.res


extreme = c(
	length(val.input[val.input < -Int.End])/length(val.input[val.input < -crit]),
	length(val.input[val.input > Int.End])/length(val.input[val.input > crit])
)
names(extreme) = c("Ext.Neg","Ext.Pos")
(extreme*100)


### End of Preparation
###############################################
###############################################
###############################################


### OLD FASHIONED Z-CURVE 

#zsds = rep(1,components)

if (Est.Method %in% c("OF","EM","CLU")) {

#summary(val.input)
para.est.OF = old.fashioned(val.input)

fit = para.est.OF[components+1];fit

para.est.OF = c(para.est.OF[1:components],ncp,zsds)

if (CURVE.TYPE == "z") {
  cp.res = Compute.Power.Z(para.est.OF,Int.Beg=Int.Beg)
} else {
  cp.res = Compute.Power.T(para.est.OF,Int.Beg=Int.Beg)
}
round(cp.res,3)

EDR = cp.res[1];EDR
ERR = cp.res[4];ERR

w.all = cp.res[which(substring(names(cp.res),1,5) == "w.all")]
w.all
sum(w.all)

w.sig = cp.res[which(substring(names(cp.res),1,5) == "w.sig")]
round(w.sig,3)

loc.power = cp.res[which(substring(names(cp.res),1,2) == "lp")]
round(loc.power,3)

#print("Finish Old Fashioned")

} # EOF if Est.Method OF
 
#########################################
### EXTENDED Z-CURVE

if (Est.Method == "EXT") {

	components = length(ncp)

	para.est.EXT = extended.curve(val.input,ncp,zsds);para.est.EXT
	fit = para.est.EXT[(3*components+1)];fit

	para.est.EXT = para.est.EXT[1:(3*components)]
	para.est.EXT

	w.inp = para.est.EXT[1:components];w.inp
	sum(w.inp)
	round(w.inp,3)

	if (W.FIXED) w.inp = w.fix

	ncp = para.est.EXT[(components+1):(2*components)]
	ncp

	zsds = para.est.EXT[(2*components+1):(3*components)]
	zsds

	#print("zsds")
	#print(zsds)

	if (max(zsds) < 1.05) {
		print("SD==1")
		if (CURVE.TYPE == "z") {
			cp.res = Compute.Power.Z(para.est.EXT,Int.Beg=Int.Beg)
		} else {
			cp.res = Compute.Power.T(para.est.EXT,Int.Beg=Int.Beg)
		}
	} else {
		print("SD > 1")
		cp.res = Compute.Power.SDG1(para.est.EXT,Int.Beg=Int.Beg)
	}			

	round(cp.res,3)	

	EDR = cp.res[1];EDR
	ERR = cp.res[4];ERR

	#print("Check")
	#print(cp.res)

	w.all = cp.res[which(substring(names(cp.res),1,5) == "w.all")]
	if (components == 1) w.all = 1

	if (Int.Beg == 0) w.all = w.inp

	w.sig = cp.res[which(substring(names(cp.res),1,5) == "w.sig")]
	if (components == 1) w.sig = 1
	round(w.sig,3)

	loc.power = cp.res[which(substring(names(cp.res),1,2) == "lp")]
	loc.power

#	print("Extended Version Completed")

} # EOF Extended Version 


##########################################################

 
### if requested, run slower EM
### Est.Method = "EM"
if(Est.Method == "EM" & boot.iter == 0) {

	components = length(ncp)	
	#summary(val.input)
	#print("Fitting EM")

	#Int.Beg = 2.4; alpha = .005

	res.em = run.zcurve(val.input,Est.Method="EM",boot.iter=boot.iter,alpha=alpha,
		Int.Beg=Int.Beg,Int.End=Int.End,parallel=parallel)
	summary(res.em)     

	fit = res.em$fit$Q
	fit

	para.est.EM = res.em$fit$weights
	para.est.EM

	para.est.EM = c(para.est.EM,ncp,zsds)

	cp.res = Compute.Power.Z(para.est.EM,Int.Beg=Int.Beg)
	cp.res

	if (1 == 2) {

		para.est.OF = old.fashioned(val.input)
		fit = para.est.OF[components+1];fit
		para.est.OF = c(para.est.OF[1:components],ncp,zsds)
		cp.res.of = Compute.Power.Z(para.est.OF,Int.Beg=Int.Beg)

		summary(res.em)$coefficients[2:1]
		cp.res[c(1,4)]
		cp.res.of[c(1,4)]

	}


	EDR = cp.res[1];EDR
	ERR = cp.res[4];ERR

	w.all = cp.res[which(substring(names(cp.res),1,5) == "w.all")]
	w.all
	sum(w.all)

	w.sig = cp.res[which(substring(names(cp.res),1,5) == "w.sig")]
	round(w.sig,3)

	loc.power = cp.res[which(substring(names(cp.res),1,2) == "lp")]
	round(loc.power,3)


	#print("Finished EM")
} 



### use only for testing, replaced by OF
if(Est.Method == "density" & boot.iter == 0) {

	print(paste("Est.Method = ",Est.Method))

	components = length(ncp)	
	#summary(val.input)
	#print("Fitting EM")
	z.res = run.zcurve(val.input,Est.Method="density",boot.iter=boot.iter,
		Int.Beg=Int.Beg,Int.End=Int.End,parallel=parallel)
	summary(z.res)     
	para.est.EM = c(summary(z.res, type="parameters")$coefficients)
	para.est.EM = para.est.EM[c((components+1):(2*components),1:components)]
	para.est.EM = c(para.est.EM,ncp)
	para.est.EM
	
	fit = z.res$fit$objective

	w.inp =	 summary(z.res, type="parameters")$coefficients[(components+1):(2*components)]
	round(w.inp,3)

	#w.all = w.sig/(pow.dir+sign.error)
	#w.all = w.all/sum(w.all)
	#round(w.all,3)

	cp.res = Compute.Power.Z(Int.Beg = Int.Beg,c(w.inp,ncp,zsds))
	round(cp.res,3)

	w.all = cp.res[which(substring(names(cp.res),1,5) == "w.all")]
	round(w.all,3)
	sum(w.all)

	loc.power = cp.res[which(substring(names(cp.res),1,2) == "lp")]
	round(loc.power,3)

	EDR = z.res$coefficients[2];EDR
	ERR = z.res$coefficients[1];ERR

	print("Finished density")

} #EOF density


### ADD RESULTS

FDR = round((1/EDR - 1)*(alpha/(1-alpha)),2);
names(FDR) = "FDR"
FDR

res = c(ODR,EDR,ERR,FDR,fit)
names(res) = c("ODR","EDR","ERR","FDR","FIT")

#print("Check res")
#print(res)

#print("Finished Computations")
#res is needed for text in plot

p.bias = c(NA,NA,NA)
if(TEST4BIAS) { 
	p.bias = test.bias(w.all) 
}
names(p.bias) = c("OBS.JS","EXP.JS","EJS.p")
#print("bias test")
#print(p.bias)

res.text = c(res[1:4],p.bias[3])

res.zsds = zsds


#print("RESULTS")
#print(res.text)

### END OF ADD RESULTS


####################################################################

if (Est.Method %in% c("CLU", "CLU-W","CLU-B") & boot.iter >= 0) {

	bb = Int.End
    if (bb > 6) bb = 6
	
	print("USING CLUSTER OPTION BOOTSTRAP")

	method = "b"

    if (Est.Method == "CLU-W") print("CLU-W NOT SUPPORTED") #method = "w"

    z.clu = zcurve_clustered(z.clu.input,method = method, bootstrap=boot.iter,
            control=list(method="EM",a = Int.Beg,b = bb,
            sig_level=alpha,max_iter = max_iter,max_iter_boot = max_iter_boot))

	summary(z.clu)     

	fit = z.clu$fit$Q

	w.inp =	 summary(z.clu, type="parameters")$coefficients[(components+1):(2*components)]
	round(w.inp,3)

	cp.res = Compute.Power.Z(Int.Beg = Int.Beg,c(w.inp,ncp,zsds))
	round(cp.res,3)

	w.all = cp.res[which(substring(names(cp.res),1,5) == "w.all")]
	round(w.all,3)
	sum(w.all)

	loc.power = cp.res[which(substring(names(cp.res),1,2) == "lp")]
	round(loc.power,3)

	#str(summary(z.clu))

	EDR = summary(z.clu)$coefficients[2,];EDR
	ERR = summary(z.clu)$coefficients[1,];ERR
	FDR = round((1/EDR - 1)*(alpha/(1-alpha)),2)[c(1,3,2)];FDR

} # EOF Cluster Method 


########################################################## 


##########################################
### This Code is Used to Create Graphic
##########################################

if (Show.Histogram & sum(extreme,na.rm=TRUE) < .95) { 

	print("SHOW RESULTS")
	print(res.text)

	Draw.Histogram(w.all,results=res.text,cola=col.hist)

	if (Show.KD) Draw.KD(val.input,w.all,cola=col.kd)

	if (Show.Curve.All & max(zsds) < 1.05) {

		Draw.Curve.All(w=w.all,cola=col.curve,
			Ltype=3,Lwidth = 4,x.start=x.lim.min,x.end=x.lim.max)
		Draw.Curve.All(w=w.all,cola=col.curve,
			Ltype=1,Lwidth = 4,x.start=Int.Beg,x.end=Int.End)

		}

	if (Show.Curve.All & max(zsds) > 1.05) {
		Draw.Curve.All.SDG1(w=w.all,ncp=ncp,zsds=zsds,cola=col.curve,
			Ltype = 3,Lwidth = 4,x.start=x.lim.min,x.end=x.lim.max)
		Draw.Curve.All.SDG1(w=w.all,ncp=ncp,zsds=zsds,cola=col.curve,
			Ltype=1,Lwidth = 4,x.start=Int.Beg,x.end=Int.End)
		}

	if (length(loc.power > 0)) Write.Local.Power(loc.power)

	} # End of Show.Histogram


########################################################################
########################################################################
########################################################################

res.het = NA
if (TEST4HETEROGENEITY > 0) {
	boot.run = TEST4HETEROGENEITY
	if (!is.numeric(TEST4HETEROGENEITY)) print("!!! Need to provide number of iterations!!!")
	else res.het = run.heterogeneity.test(val.input,boot.run=boot.run,fit.ci)
	}



print(round(res,3))

return.results = list(
    res = res,
    ncp = ncp,
    zsds = res.zsds,
    w.all = w.all,
    bias = p.bias,
    fit.comp = res.het
  )




### code for confidence intervals 

########################################################################
########################################################################
########################################################################

### If Confidence Intervals are requested, compute CI (boot.iter > 0)
if (boot.iter > 0 & Est.Method %in% c("OF","EM") )  {

	res.with.ci = get.ci.info(Est.Method = Est.Method)

	res.with.ci = rbind(ODR.res,res.with.ci,p.bias)
	round(res.with.ci,3)


	res.text = rbind(res.with.ci[1:4,],p.bias)
	rownames(res.text) = c("ODR","EDR","ERR","FDR","bias")
	round(res.text,3)

	res = c(res,w.all)
	round(res,3)	 
	length(res)

}


#boot.iter = 50	# Testing
if (boot.iter > 0 & Est.Method == "EXT") {

	res.ci = get.ci.info(Est.Method = Est.Method)
	res.ci = rbind(c(ODR.low,ODR.high),res.ci)
	rownames(res.ci)[1] = "ODR"
	round(res.ci,3)
	dim(res.ci)

	FDR = round((1/EDR - 1)*(alpha/(1-alpha)),2);
	names(FDR) = "FDR"
	FDR

	res.text = res[1:4]
	names(res.text) = c("ODR","ERR","EDR","FDR")

	res.text = cbind(res.text,res.ci[1:4,])
	res.text = round(res.text,3)

	res.text = rbind(res.text,c(p.bias))
	rownames(res.text)[5] = "BIAS"
	round(res.text,3)
	
	res = res.text
	res

	res.zsds = c(zsds,res.ci[which(rownames(res.ci) == "ZSD"),])

}


if (boot.iter > 0 & substring(Est.Method,1,3) == "CLU") {

	#print("MAKE RESTEXT 4 CLU")

	results = rbind(c(ODR,ODR.low,ODR.high),EDR,ERR,FDR)
	rownames(results) = c("ODR","EDR","ERR","FDR")
	round(results,3)

	res.text = rbind(results,c(p.bias))
	res.text = round(res.text,3)

	rownames(res.text)[5] = "BIAS"


	#print(round(res.text,3))

	res = res.text
	res

}

if (boot.iter > 0 & Show.Histogram) {

	Draw.Histogram(w.all,results=res.text,cola=col.hist,Write.CI = TRUE)

	if (Show.KD) Draw.KD(val.input,w.all,cola=col.kd)

	if (Show.Curve.All & Est.Method != "DF" & max(zsds) < 1.05) {
		Draw.Curve.All(w=w.all,cola=col.curve,
			Ltype=3,Lwidth = 4,x.start=x.lim.min,x.end=x.lim.max)
		Draw.Curve.All(w=w.all,cola=col.curve,
			Ltype=1,Lwidth = 4,x.start=Int.Beg,x.end=Int.End)
		}

	if (Show.Curve.All & Est.Method != "DF" & max(zsds) > 1.05) {
		Draw.Curve.All.SDG1(w=w.all,ncp=ncp,zsds=zsds,cola=col.curve,
			Ltype = 3,Lwidth = 4,x.start=x.lim.min,x.end=x.lim.max)
		Draw.Curve.All.SDG1(w=w.all,ncp=ncp,zsds=zsds,cola=col.curve,
			Ltype=1,Lwidth = 4,x.start=Int.Beg,x.end=Int.End)
		}

	if (Show.Curve.Sig) {
		Draw.Curve.Sig(val.input,ncp=ncp,zsds=zsds,w.sig,cola=col.curve,Ltype=3,x.start=x.lim.min)
		Draw.Curve.Sig(val.input,ncp=ncp,zsds=zsds,w.sig,cola=col.curve,Ltype=1,x.start=Int.Beg)
		}
		#	par(family = fam[1])	
		if (length(loc.power > 0)) Write.Local.Power(loc.power)

} # End of Show.Histogram	for bootstrap


if (boot.iter > 0) {

return.results = list(
    res = res.text[1:4,],
    ncp = ncp,
    zsds = res.zsds,
    w.all = w.all,
    bias = p.bias,
    fit.comp = res.het
  )

}


return(return.results)

} #End of Zing


######################################################################
######################################################################
######################################################################


zcurve.es.estimation = function(N,pow) {

ZCURVE.EM.HOM = c()

pow.i = 1
for (pow.i in 1:length(pow)) {

	
	if (pow[pow.i] <= .05) { ZCURVE.EM.HOM = c(ZCURVE.EM.HOM,0) } else {

		n = round(N/2)

		avg.n = mean(N)/2;avg.n
		d.1 = pwr.t.test(
		n = avg.n, d = NULL, sig.level = 0.05, power = pow[pow.i],
		type = "two.sample", alternative = "two.sided")$d 
		d.1

		d.2 = unlist(lapply(n,function(x) pwr.t.test(
		n = x, d = NULL, sig.level = 0.05, power = pow[pow.i],
		type = "two.sample", alternative = "two.sided")$d ))
		d.2 = mean(d.2)
		d.2

		### do not use d.1
		ZCURVE.EM.HOM = c(ZCURVE.EM.HOM,d.2)
	} # end else

} # End pow.i

#print(ZCURVE.EM.HOM)

return(ZCURVE.EM.HOM)

} # End of Function

