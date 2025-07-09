

#rm(list = ls())

#options(scipen = 999)

#install.packages("pwr")
#install.packages("zcurve")
#install.packages("KernSmooth")

#install.packages("extrafont")
#install.packages("showtext")

library(parallel)
library(KernSmooth)
library(zcurve)
library(stringr)
library(pwr)

#library(extrafont)
#library(showtext)
#cl <- parallel::makeCluster(parallel::detectCores() - 10)


#setwd("C:\\Users\\ulric\\Dropbox\\PHPCurves\\R-Code\\JetBrainsMono-2.304\\fonts\\ttf\\")
#font_add("JetBrains","JetBrainsMono-Medium.ttf")
#font_add("JetBrains","JetBrainsMono-Regular.ttf")
#par(family = "JetBrainsMonoNL")
#par(family = "JetBrainsMono")
#par(family = "mono")

#### Start Default Settings
#### Gobal Variables Can be Changed in R-Code Anywhere

TESTING = FALSE

#windowsFonts()

#fam = c("mono", "mono")
#fam = c("mono", "Monotype Corsiva")
#fam = c("Microsoft Himalaya","Monotype Corsiva")
#fam = rep("Microsoft Himalaya",2);fam
#fam = c("Cancadia Mono SemiBold","JetBrainsMono")

letter.size = 1
letter.size.1 = letter.size

Title = ""

version = "Version 25.06"

alpha = .05
z.crit = qnorm(1-alpha/2)

two.sided = TRUE
#two.sided == FALSE means one-sided tests in both directions

### settings
x.lim.min = 0
x.lim.max = 6

fixed.false.positives = 0

Int.Beg = z.crit
Int.End = x.lim.max
Int.Ext = Int.End

ncz.int = 1
ncz = seq(0,Int.End,ncz.int);ncz
components = length(ncz)
zsd = 1
SD.FIXED = TRUE

lowlim = c(rep(0,length(ncz)),ncz,rep(1,length(ncz)));lowlim
highlim = c(rep(1,length(ncz)),ncz,rep(1,length(ncz)));lowlim


### is used for the point estimate
max_iter = 1e6
max_iter_boot = 1e5

fast_iter = max_iter

### faster iteration for bootstrap
fast_iter_boot = 1e4

bw.draw = .10
bw.est = .05


y.line.factor = 3

ymax = .6
ymin = -ymax/15;ymin

tiva.x = 1.9

int.loc = .5

### resolution of density function (doesn't seem to matter much)	
n.bars = 512

fit.crit = .02
fit.crit.BKD = .01
fit.crit.EM = .006

ERR.CI.adjust = .03
EDR.CI.adjust = .05

MAX.INP.Z = 100

Plot.Fitting = FALSE
BOOT = FALSE

FAST.BOOT = TRUE

Show.Fitted = FALSE
Show.Histogram = TRUE
Show.Zcurve.All = TRUE
Show.Zcurve.Sig = FALSE
Show.KD = FALSE
Show.Text = TRUE
Show.Iterations = TRUE
Show.Significance = TRUE

WRITE.ADJ.Z = FALSE

sig.levels = c()

col.zcurve = "violetred3"
col.hist = "blue3"
col.kd = "black"


Augment = TRUE
Augment.Factor = 1.4
Augment.bw = .20

USE.SLOPE = FALSE
SLOPE.crit = 1
slope.width = .5

hist.bar.width = .2

Return.Weights = TRUE

Est.Method = "OF"
#Est.Method = "EM"

boot.iter = 0

TEST4HETEROGENEITY = FALSE

EM.criterion = 1e-3
EM.max.iter = 1000

parallel = TRUE

CI.ALPHA = .05

TEST4BIAS = FALSE
just = .6 ### ~ p = .01 two-tailed

print("Zing 25.06.09.test")
print("Parameter OK")
print(ncz)

###################################################
#### End of Default Settings
###################################################


###################################################
#### Include Functions in Zing 
###################################################

Zing = function(z.val.input,lp=c() ) {

#######################################################################
#######################################################################
### AAA Start Functions
#######################################################################
#######################################################################

### function to detect bias

test.bias = function(w.all) {

#w.all = c(1,rep(0,6))
print(w.all)

sig.k = length(z.val.input[z.val.input > z.crit & z.val.input < Int.End])
sig.k

just.sig.k = length(z.val.input[z.val.input > z.crit & z.val.input < z.crit + just])
just.sig.k

bar.width = .01 # how fine should be the resolution
Z.X = seq(x.lim.min,Int.End,bar.width);summary(Z.X)	 # set of z-scores 

w.all
Z.W.D = unlist(lapply(Z.X,function(x) sum(dnorm(x,ncz)*w.all) ))
DD = rbind(Z.X,Z.W.D)

DD = DD[,DD[1,] > z.crit]
#plot(DD[1,],DD[2,])

DD[2,] = DD[2,]/(sum(DD[2,])*bar.width)

prob = sum(DD[2,DD[1,] < z.crit + just])*bar.width;prob
prob2 = sum(DD[2,DD[1,] > z.crit + just])*bar.width;prob2
prob + prob2

### binomial
p.bias.binomial = 1 - pbinom(just.sig.k - 1, sig.k, prob = prob)
#print(p.bias.binomial)

### chi2 approximation
#O = just.sig.k
#E = sig.k*prob;E
#chi2.val <- (O - E)^2 / (E * (1 - prob));chi2.val
#p.bias.chi2 <- 1 - pchisq(chi2.val, df = 1)
#print(p.bias.chi2)

return(c(just.sig.k/sig.k,prob,p.bias.binomial))

}


### function to run the zcurve package to get EM weights and fit

run.zcurve = function(Est.Method = "OF",kd.model="KD2",K=6,
	alpha = .05,Int.Beg = 1.96,Int.End = 6, boot.iter = 0,parallel=TRUE) {

	print("RUNNING ZCURVE")
	print(Est.Method)
	print(parallel)

	if (Est.Method == "OF") {
		z.res = zcurve(z.val.input,bootstrap=boot.iter,method="density",
			control=list(parallel = parallel,
			a = Int.Beg,b = Int.End))
	} else {
		z.res = zcurve(z.val.input,bootstrap=boot.iter,method="EM",
			control=list(parallel = parallel,
			alpha = alpha, a = Int.Beg,b = Int.End))
	}

	z.res

return(z.res)
} ### End run.zcurve

if (TESTING) {
	z.res = run.zcurve(z.val.input,Est.Method="OF",
		kd.model="KD2",Int.Beg=0,parallel=parallel)
	z.res
}

### testing 
### ncz; z.res = run.zcurve(zval,Est.Method = "EM",Int.Beg = 0,Int.End = 6);z.res

### testing 
### Est.Method = "OF"; boot.iter = 50; z.val.input = rnorm(100,2.8)

############################################
############################################
############################################

EXT.boot = function()	{

		boot = 1
		boot.res = c()

		components = length(ncz)
		for (boot in 1:boot.iter) {

			### Get Bootstrap Sample
   		  	z.sample = sample(Z.INT, size=length(Z.INT), replace=TRUE)

			### Submit Bootstrap Sample to Parameter Estimation Function

			if (Show.Iterations) print(paste0("EXT.boot Iteration: ",boot))

			para.est.EXT = extended.zcurve(z.sample,ncz,zsd);para.est.EXT
		
			fit = para.est.EXT[(3*components+1)];fit
			para.est.EXT = para.est.EXT[1:(3*components)]

			WT = para.est.EXT[1:components];WT
			WT = WT/sum(WT)
			WT

			ncz = para.est.EXT[(components+1):(2*components)]
			ncz

			zsds = para.est.EXT[(2*components+1):(3*components)]
			zsds = 1 #testing
			zsds

			para.est.EXT = c(WT,ncz,zsds)

			if (mean(zsds) == 1) cp.res = Compute.Power(para.est.EXT,Int.Beg=Int.Beg,BOOT=TRUE)
			#if (mean(zsds) != 1) cp.res = Compute.Power.EXT(para.est.EXT)

			cp.res

			EDR = cp.res[1];EDR
			ERR = cp.res[4];ERR

			w.all = cp.res[(6+components*1):(6+components*2-1)];w.all
			w.sig = cp.res[(6+components*2):(6+components*3-1)];w.sig

			boot.res = rbind(boot.res,c(ERR,EDR,fit,ncz,zsds,w.all,w.sig))
			boot.res

		}

		dim(boot.res)

		CIs = c()
	
		CIs = rbind(CIs,quantile(boot.res[,1],c(CI.ALPHA/2,1-CI.ALPHA/2)) )

		CIs = rbind(CIs,quantile(boot.res[,2],c(CI.ALPHA/2,1-CI.ALPHA/2)) )

		CIs = rbind(CIs,quantile(boot.res[,3],c(CI.ALPHA/2,1-CI.ALPHA/2)) )

		for (i in 1:components) CIs = rbind(CIs, 
			quantile(boot.res[,(3+components*1):(3+(components*2-1))],c(CI.ALPHA/2,1-CI.ALPHA/2)) )

		for (i in 1:components) CIs = rbind(CIs, 
			quantile(boot.res[,(3+components*2):(3+(components*3-1))],c(CI.ALPHA/2,1-CI.ALPHA/2)) )

		for (i in 1:components) CIs = rbind(CIs,
			quantile(boot.res[,(3+components*3):(3+(components*4-1))],c(CI.ALPHA/2,1-CI.ALPHA/2)) )

		for (i in 1:components) CIs = rbind(CIs,
			quantile(boot.res[,(3+components*4):(3+(components*5-1))],c(CI.ALPHA/2,1-CI.ALPHA/2)) )


		print(CIs)

		CIs[1,1] = CIs[1,1]-ERR.CI.adjust
		CIs[1,2] = CIs[1,2]+ERR.CI.adjust

		CIs[2,1] = CIs[2,1]-EDR.CI.adjust
		CIs[2,2] = CIs[2,2]+EDR.CI.adjust

		colnames(CIs) = c("low","high")
		rownames(CIs) = c("ERR","EDR","FIT",rep("NCZ",components),
			rep("ZSD",components),rep("WALL",components),
			rep("WSIG",components) )
		print(CIs)

		return(CIs)


} # End function EXT.boot

#####################################
#####################################
#####################################


#z.val.input = z

get.ci.info = function(Est.Method = "EM") {

if (Est.Method %in% c("OF","EM")) {

zres = run.zcurve(z.val.input, Est.Method=Est.Method, boot.iter = boot.iter,
	Int.Beg = Int.Beg, Int.End = Int.End,parallel=parallel)

zres
ERR = summary(zres)$coefficients[1,];ERR
EDR = summary(zres)$coefficients[2,];EDR
FDR = round((1/EDR - 1)*(alpha/(1-alpha)),3);FDR

w.all = summary(zres, type="parameters")$coefficients[,2]
w.all.low = summary(zres, type="parameters")$coefficients[,3]
w.all.high = summary(zres, type="parameters")$coefficients[,4]

#round(cbind(w.all,w.all.low,w.all.high),3)	  

fit.val = summary(zres)$model$fit_index
fit.val

res.ci = c(ERR,EDR,FDR,fit.val,w.all,w.all.low,w.all.high)

} 

if (Est.Method == "EXT") {

res.ci = EXT.boot()

FDR = round((1/res.ci[2,2:1] - 1)*(alpha/(1-alpha)),2);FDR

res.ci = rbind(res.ci[1:2,],FDR,res.ci[3:7,])

res.ci

#print("METHOD EXT")
#print("RETURN res.ci")
#print(res.ci)

}

return(res.ci)

}

#################################

test.heterogeneity = function(z.val.input,boot.runs = 500, 
	fit.ci = c(.01,.025,.05,.10,.17,.20,.50,.80,.83,.90,.95,.975,.99)
	) {


#boot.i = 2; boot.runs = 50
res.boot = c()

for (boot.i in 1:boot.runs ) {
	print(boot.i)
	zboot = sample(z.val.input,replace=TRUE)
	ncz = 2
	zsd = 1
	components = length(ncz)
	para.est.EXT = extended.zcurve(zboot,ncz,zsd);para.est.EXT
	fit.hom = para.est.EXT[(3*components+1)];fit.hom
	ncz = c(1,3)
	zsd = 1
	components = length(ncz)
	para.est.EXT = extended.zcurve(zboot,ncz,zsd);para.est.EXT
	fit.het.1 = para.est.EXT[(3*components+1)];fit.het.1
	ncz = 2
	zsd = 5
	components = length(ncz)
	para.est.EXT = extended.zcurve(zboot,ncz,zsd);para.est.EXT
	fit.het.2 = para.est.EXT[4];fit.het.2
	sd.het.2 = para.est.EXT[3];sd.het.2

	delta.fit.1 = fit.hom - fit.het.1
	delta.fit.2 = fit.hom - fit.het.2
	delta.fit.3 = fit.het.2 - fit.het.1

	res.boot = rbind(res.boot,c(fit.hom,fit.het.1,fit.het.2,sd.het.2,
		delta.fit.1,delta.fit.2,delta.fit.3))
	print(res.boot[boot.i,])

} # end of boot

### diagnostics
if (1 == 2) {
	dim(res.boot)
	summary(res.boot)
	cor(res.boot[,1:3])
}

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
#	"hom.odr","hom.edr","hom.err","hom.fdr",
	"fit.hom","fit.het.1","fit.het.2",
	"sd.het.2",	
	"delta.fit.2-1","delta.fit.3-1","delta.fit.2-3"
)

print(return.res)
	
return(return.res)

}


#######################################################################


Write.Local.Power = function(loc.power) {
	names(loc.power) = paste0("LP",seq(1,length(loc.power)))
	int = seq(x.lim.min,x.lim.max-int.loc,int.loc)+int.loc/2
	for (i in 1:length(int)) text(int[i],ymin/10,paste0(format(round(loc.power[i]*100),nsmall=0),"%"),srt=90,pos=2
                               ,cex=letter.size)
}


###################################################
#### Begin Draw Histogram
###################################################

Draw.Histogram = function(z.draw,w,cola="blue3",
	results,Write.CI = FALSE) {


	#z.draw = sim.dat$z
	#cola="blue3"
	z.hist = z.draw[z.draw > x.lim.min & z.draw < x.lim.max]
	if (round(Int.Beg,2) == 1.96) {
		z.hist = z.draw[z.draw > x.lim.min & z.draw < x.lim.max -.04] + .04
	}
	z.hist

	scale = 
	length(z.draw[z.draw > x.lim.min & z.draw < x.lim.max]) /
	length(z.draw[z.draw > x.lim.min & z.draw < x.lim.max])
	scale

	Int.Beg
	ymin/scale
	ymax/scale
	
	n.breaks = seq(x.lim.min,x.lim.max,hist.bar.width);n.breaks

	par(cex.axis = 1)
	par(cex.lab = 1)
	par(family = "sans") 
	par(font.axis = 1) #"Microsoft Himalaya")
	par(font.lab = 1) # Microsoft Himalaya")

	###z.hist = z.val.input
	### draw a histogram of observed z-scores
	hist(z.hist[z.hist < x.lim.max],breaks=n.breaks,freq=FALSE,
		col=adjustcolor(cola, alpha.f = 0.2),border="white",
		xlim=c(x.lim.min,x.lim.max),ylim=c(ymin,ymax),
		ylab="Density",xlab="absolute z-value",main=Title,lwd=1)

	par(new=TRUE)


	scale = 
	length(z.draw[z.draw > Int.Beg & z.draw < Int.End]) /
	length(z.draw[z.draw > x.lim.min & z.draw < x.lim.max])
	scale

	#n.breaks = seq(round(Int.Beg,1),x.lim.max,hist.bar.width);n.breaks

	hist(z.hist[z.hist > round(Int.Beg,1) & z.hist < Int.End],
		breaks=n.breaks,freq=FALSE,
		col=adjustcolor(cola, alpha.f = 0.3),border="white",
		xlim=c(x.lim.min,x.lim.max),ylim=c(ymin/scale,ymax/scale),
		ylab=,xlab="absolute z-value",main=Title,lwd=1,axes=FALSE)




######################################### 
######################################### 
######################################### 

if (Show.Text) {

#par(family = fam[2])

par(new=TRUE)
hist(c(0),main="",ylim=c(ymin,ymax),ylab="",xlab="",xlim=c(x.lim.min,x.lim.max),
	density=0,border="white",axes=FALSE)

	#results = res
	#results = res.with.ci


	if(TEST4BIAS) p.bias = results[length(results)];p.bias

	results = round(results*100);results

	ODR = results[1];ODR
	EDR = results[2];EDR
	ERR = results[3];ERR
	FDR = results[4];FDR
	
	min.z = min(z.draw)
	max.z = max(z.draw)
	n.z = length(z.draw)
	n.z.sig = length(z.draw[z.draw > z.crit])
	n.not.shown = length(z.draw[z.draw > x.lim.max])

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

	i = 0

	results.x = x.lim.max 



#	text(results.x,y.text-y.line*i,paste0("Range: ",
#		round(min.z,2)," to ",round(max.z,2)),
#		pos=2,cex=letter.size)
#	i = i + 2

	text(results.x+0.1,y.text-y.line*i,paste0(toString(format(n.z,big.mark=",")),
		" tests "),pos=2,cex = letter.size)
	i = i + 1.5
	text(results.x+0.1,y.text-y.line*i,paste0(toString(format(n.z.sig,big.mark=","))," significant "),
	  pos=2,cex = letter.size)
	i = i + 1.5
	text(results.x,y.text-y.line*i,paste0(n.not.shown," z > ",x.lim.max),
		pos=2,cex = letter.size)
	

#	par(family = fam[1])

	######

	#print("Check write CI")
	#print(Write.CI)

	if (Write.CI) {

	ODR = results[1]
	ODR.low = results[2]
	ODR.high = results[3]
	ERR = results[4]
	ERR.low = results[5]
	ERR.high = results[6]
	EDR = results[7]
	EDR.low = results[8]
	EDR.high = results[9]
	FDR = results[10]
	FDR.low = results[12]
	FDR.high = results[11]

	if (ODR.high > 100) ODR.high = 100
	#print("Check ERR.low")
	#print(ERR.low)
	if(ERR.low < 2.5) ERR.low = 2.5
	if(EDR.low < 5) EDR.low = 5
	
	i = i + 3
	text(results.x,y.text-y.line*i,
		paste0(round((1-CI.ALPHA)*100),"% CI     "),
		pos=2,cex=letter.size.1)

	i = i + 2

	ODR.a = str_pad(round(ODR), width = 3, pad = " ");ODR.a
	ODR.b = str_pad(round(ODR.low), width = 3, pad = " ");ODR.b
	ODR.c = str_pad(round(ODR.high), width = 3, pad = " ");ODR.c
	text(results.x,y.text-y.line*i,
		paste0("ODR: ",ODR.a,"% [",	ODR.b,"%;",ODR.c,"%]")
		,pos=2,cex=letter.size.1)

	i = i + 2

	EDR.a = str_pad(round(EDR), width = 3, pad = " ");EDR.a
	EDR.b = str_pad(round(EDR.low), width = 3, pad = " ");EDR.b
	EDR.c = str_pad(round(EDR.high), width = 3, pad = " ");EDR.c
	text(results.x,y.text-y.line*i,
		paste0("EDR: ",EDR.a,"% [",EDR.b,"%;",EDR.c,"%]")
		,pos=2,cex=letter.size.1)

	i = i + 2

	ERR.a = str_pad(round(ERR), width = 3, pad = " ");ERR.a
	ERR.b = str_pad(round(ERR.low), width = 3, pad = " ");ERR.b
	ERR.c = str_pad(round(ERR.high), width = 3, pad = " ");ERR.c
	text(results.x,y.text-y.line*i,
		paste0("ERR: ",ERR.a,"% [",ERR.b,"%;",ERR.c,"%]")
		,pos=2,cex=letter.size.1)


	i = i + 2

	FDR.a = str_pad(round(FDR), width = 3, pad = " ");FDR.a
	FDR.b = str_pad(round(FDR.low), width = 3, pad = " ");FDR.b
	FDR.c = str_pad(round(FDR.high), width = 3, pad = " ");FDR.c
	text(results.x,y.text-y.line*i,
		paste0("FDR: ",FDR.a,"% [",FDR.b,"%;",FDR.c,"%]")
		,pos=2,cex=letter.size.1)

	if(TEST4BIAS) {
		i = i + 2
		text(results.x,y.text-y.line*i,
			paste0("EJS p = ",sub("^0","",format(round(p.bias,4),nsmall=4))),
			pos=2,cex=letter.size.1)
	}


#	i = i + 2
#	text(results.x,y.text-y.line*i,
#		paste0("OSR:   ",	OSR,"%   [ ",OSR.low,"% ; ",OSR.high,"% ]")
#		,pos=2,cex=letter.size)

	# End of IF CI show.ci.results
	
	} else {  

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

	if(TEST4BIAS) {
		i = i + 2
		text(results.x,y.text-y.line*i,
			paste0("EJS, p = ",sub("^0","",format(round(p.bias,4),nsmall=4))),
			pos=2,cex=letter.size.1)
	}


#	i = i + 2
#	text(results.x,y.text-y.line*i,
#		paste0("OSR:   ",OSR,"%"),
#		pos=2,cex=letter.size.1)

	i = i + 2
	#print("Check Waiting for Bootstrap")
	#print(boot.iter)
	#print(i)
	#print(results.x)
	#print(ymax-y.line*i)

	if (boot.iter > 0 & Write.CI == FALSE) text(results.x,ymax-y.line*i,
		"WAIT FOR BOOSTRAPPED CIs",
		pos=2,cex=letter.size,col="red")


	} # End of CI

} # End of Show.Text


abline(h=0)


} # End of Histogram
 
######################################### 
######################################### 
######################################### 

######################################### 
######################################### 
######################################### 


Draw.KD = function(z.draw,w,Write.CI=FALSE,cola="blue",Lwidth=5) {

	#ymin=-.015;ymax = .6
	#Draw.Histogram(z.val.input,w.all,results=res,cola=col.hist)

	#z.draw = z.val.input
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

Draw.Zcurve.All.SDG1 = function(w.inp,ncz=ncz,zsds=zsds,cola="black",
	Ltype=1,Lwidth=3,x.start=x.lim.min,x.end=x.lim.max) {

components = length(ncz)

nczs = seq(x.lim.min,x.lim.max,.5)
pows = pnorm(nczs,z.crit) + pnorm(-nczs,z.crit)
cbind(nczs,pows)

comp.i = 1
set = c()

for (comp.i in 1:components) {
	w.nczs = dnorm(nczs,ncz[comp.i],zsds[comp.i]-1) 
	w.nczs = w.nczs/sum(w.nczs)
	#plot(nczs,w.nczs)
	subset = data.frame(nczs,pows,w.nczs)
	subset = subset[round(w.nczs,4) > 0,]
	subset$w.nczs = subset$w.nczs / sum(subset$w.nczs)
	subset$comp.i = comp.i
	subset$w.inp = w.inp[comp.i] 
	set = rbind(set,subset)
}
dim(set)

set$w.all = set$w.inp/set$pows
set$w.all = set$w.all/sum(set$w.all)
set$w.all

set

z.x.i = 1
z.y = c()
for (z.x.i in 1:length(z.x)) {
	w = set$w.nczs * set$w.all
	w = w / sum(w)
	z.y = c(z.y,sum(dnorm(z.x[z.x.i],set$nczs)*w)/sum(w) )
} 
z.y

#plot(z.x,z.y)

z.y = z.y/sum(z.y*bar.width)

Z.Density.X = z.x
z.est = z.y

sum(z.est*bar.width)

#plot(Z.Density.X,z.est)

#############################

d.hist.sel = mean(as.numeric(z.val.input > x.lim.min & z.val.input > Int.Beg & z.val.input < Int.End)) /
	mean(as.numeric(z.val.input > x.lim.min & z.val.input < Int.End))
d.hist.sel

table(Z.Density.X >= Int.Beg)
d.dense.sel = sum(z.est[Z.Density.X > x.lim.min & Z.Density.X >= Int.Beg]*bar.width)
d.dense.sel

cbind(Z.Density.X,z.est)

scale = d.hist.sel/d.dense.sel;scale

par(new=TRUE)

lines(Z.Density.X[which(Z.Density.X >= x.start & Z.Density.X < x.end)],
	z.est[which(Z.Density.X >= x.start & Z.Density.X < x.end)]*scale,
	lty=Ltype,col=cola,lwd=4,xlim=c(x.lim.min,x.lim.max),ylim=c(ymin,ymax))

z.boundary = Int.Beg
if (round(Int.Beg,2) == 1.96) {z.boundary = 2}

Z.Density.X	
segments(z.boundary,0,z.boundary,
	z.est[which(round(Z.Density.X,2) == round(z.boundary,2))[1]]*scale,
	col=cola,lwd=3)


if (Show.Significance) {

	if (length(sig.levels) == 0) sig.levels = alpha

		cz.alpha = -qnorm(sig.levels/2);cz.alpha
		if(two.sided == FALSE) cz.alpha = c(-cz.alpha,cz.alpha);cz.alpha

		sig.lty = rep(2,length(cz.alpha));sig.lty

		i = 1
		for (i in 1:length(cz.alpha)) {
			cz.alpha[i]
			height = max(scale*z.est[which(round(Z.Density.X,1) == round(cz.alpha[i],1))]);height
			segments(cz.alpha[i],0,cz.alpha[i],height,lty=sig.lty,lwd=2,col="firebrick3")
		}

} # End of Show.Significance


} # End of Draw.Zcurve.SDG1

################################################
################################################
################################################

###################################################
#### Begin Draw Zcurve
###################################################

Draw.Zcurve.All = function(z.val.input,w,ncz,zsds,cola="black",Lwidth=3,
	Ltype=1,x.start=x.lim.min,x.end=x.lim.max) {

print("Draw.Zcurve.All")

if (1 == 2) {
	x.start = x.lim.min
	x.start = Int.Beg
	Lwidth = 2
	cola = col.zcurve
	z.draw = z.val.input
	w = WT;w
	w = w.all
	Ltype = 1
}

bar.width = .01
Z.Density.X = seq(0,x.lim.max,bar.width)

n.bars = length(Z.Density.X);n.bars

Dens	= c()
for(i in 1:n.bars) {
	for (j in 1:length(ncz)) {
		Dens = c(Dens,dnorm(Z.Density.X[i],ncz[j],zsds[j]))
	}
}
Dens = matrix(Dens,length(ncz),byrow=FALSE)
sum.dens = rowSums(Dens)
Dens = Dens/(sum.dens * bar.width)

z.est = c()
for (i in 1:n.bars) z.est[i] = sum(Dens[,i]*w)
summary(z.est)
sum(z.est*bar.width)

d.hist.sel = mean(as.numeric(z.val.input > x.lim.min & z.val.input > Int.Beg & z.val.input < Int.End)) /
	mean(as.numeric(z.val.input > x.lim.min & z.val.input < Int.End))
d.hist.sel

table(Z.Density.X >= Int.Beg)
d.dense.sel = sum(z.est[Z.Density.X > x.lim.min & Z.Density.X >= Int.Beg]*bar.width)
d.dense.sel

cbind(Z.Density.X,z.est)

scale = d.hist.sel/d.dense.sel;scale

par(new=TRUE)

lines(Z.Density.X[which(Z.Density.X >= x.start & Z.Density.X < x.end)],
	z.est[which(Z.Density.X >= x.start & Z.Density.X < x.end)]*scale,
	lty=Ltype,col=cola,lwd=4,xlim=c(x.lim.min,x.lim.max),ylim=c(ymin,ymax))

z.boundary = Int.Beg
if (round(Int.Beg,2) == 1.96) {z.boundary = 2}

Z.Density.X	
segments(z.boundary,0,z.boundary,
	z.est[which(round(Z.Density.X,2) == round(z.boundary,2))[1]]*scale,
	col=cola,lwd=3)


if (Show.Significance) {

	if (length(sig.levels) == 0) sig.levels = alpha

		cz.alpha = -qnorm(sig.levels/2);cz.alpha
		if(two.sided == FALSE) cz.alpha = c(-cz.alpha,cz.alpha);cz.alpha

		sig.lty = rep(2,length(cz.alpha));sig.lty

		i = 1
		for (i in 1:length(cz.alpha)) {
			cz.alpha[i]
			height = max(scale*z.est[which(round(Z.Density.X,1) == round(cz.alpha[i],1))]);height
			segments(cz.alpha[i],0,cz.alpha[i],height,lty=sig.lty,lwd=2,col="firebrick3")
		}

} # End of Show.Significance


} # End of Draw.Zcurve.All

################################################
################################################
################################################



Draw.Zcurve.Sig = function(z.draw,w,ncz,zsds,cola="black",Lwidth=3,Ltype=1,x.start=x.lim.min) {

#print("Draw.Zcurve.Sig")

if (1 == 2) {
	x.start = x.lim.min
	x.start = Int.Beg
	Lwidth = 2
	cola = col.zcurve
	z.draw = z.val.input
	w = WT;w
	w = w.sig
	Ltype = 1
}

bar.width = .01
Z.Density.X = seq(Int.Beg,x.lim.max,bar.width)
n.bars = length(Z.Density.X);n.bars


Dens	= c()
for(i in 1:n.bars) {
	for (j in 1:length(ncz)) {
		Dens = c(Dens,dnorm(Z.Density.X[i],ncz[j],zsds[j]))
	}
}
Dens = matrix(Dens,length(ncz),byrow=FALSE)
sum.dens = rowSums(Dens)
Dens = Dens/(sum.dens * bar.width)

z.est = c()
for (i in 1:n.bars) z.est[i] = sum(Dens[,i]*w)
summary(z.est)
sum(z.est*bar.width)

#plot(Z.Density.X,z.est)

par(new=TRUE)
lines(Z.Density.X[which(Z.Density.X >= x.start)],z.est[which(Z.Density.X >= x.start)]*scale,lty=Ltype,col=cola,lwd=4,
 xlim=c(x.lim.min,x.lim.max),ylim=c(ymin,ymax))

z.boundary = Int.Beg
if (round(Int.Beg,2) == 1.96) {z.boundary = 2}

#Z.Density.X	
#segments(z.boundary,0,z.boundary,
#	z.est[which(round(Z.Density.X,2) == round(z.boundary,2))[1]]*scale,
#	col=cola,lwd=3)


if (Show.Significance) {

	if (length(sig.levels) == 0) sig.levels = alpha

		cz.alpha = -qnorm(sig.levels/2);cz.alpha
		if(two.sided == FALSE) cz.alpha = c(-cz.alpha,cz.alpha);cz.alpha

		sig.lty = rep(2,length(cz.alpha));sig.lty

		i = 1
		for (i in 1:length(cz.alpha)) {
			cz.alpha[i]
			height = max(scale*z.est[which(round(Z.Density.X,1) == round(cz.alpha[i],1))]);height
			segments(cz.alpha[i],0,cz.alpha[i],height,lty=sig.lty,lwd=2,col="firebrick3")
		}

} # End of Show.Significance


} # End of Draw.Zcurve Sig

################################################
################################################
################################################


################################################
################################################
################################################

##############################################
### Get Densities
##############################################

Get.Densities = function(zval,bw="nrd0",d.x.min=0,d.x.max=6,Augment=TRUE) {


#zval = z.val.input;d.x.min = Int.Beg; d.x.max = 6;bw=bw.draw
#zval = z.val.input[z.val.input > Int.Beg];d.x.max = 6;d.x.min = z.crit;bw = .2

### find the maximum z-score. This is only needed if the maximum z-score is below Int.End
max.z = d.x.max
if (max(zval) < d.x.max) max.z = max(zval)
max.z

### Augment z-scores on the left side of Interval to avoid downward trend 
### of kernal density function (assumes values go to 0)

if (Augment) { 

	AUG = c()
	n.AUG = round(Augment.Factor*length(zval[zval > d.x.min & zval < d.x.min+Augment.bw]));n.AUG
	if (n.AUG > 0) AUG = seq(d.x.min-Augment.bw,d.x.min-.01,Augment.bw/n.AUG)

	Z.INT.USE = c(zval,AUG)
	#hist(Z.INT.USE[Z.INT.USE < 6],breaks=20)
	#summary(AUG)

} else {
	#print("Not Augmenting")
	Z.INT.USE = zval[zval > d.x.min & zval <= max.z]
}

#print(summary(Z.INT.USE))

Z.Density = bkde(Z.INT.USE,bandwidth=bw,range=c(d.x.min,d.x.max)) 
val.max = d.x.max
D = data.frame(Z.Density$x,Z.Density$y)
colnames(D) = c("ZX","ZY")
#plot(D$ZX,D$ZY)
D = D[D$ZX > d.x.min & D$ZX < val.max,]
dim(D)

bar.width = D$ZX[2] - D$ZX[1]
D$ZY = D$ZY/(sum(D$ZY*bar.width)) 
sum(D$ZY*bar.width)

#print(summary(D))

#print("End Augment check")

return(D)

}  ### End of Get Densities 


### Test Get.Densities 
if(1==2){
Int.End = 6
z.val.input = z
Z.INT = z.val.input[z.val.input >= Int.Beg & z.val.input <= Int.End]
Augment = FALSE
densy = Get.Densities(zval=Z.INT,bw=bw.est,d.x.min=Int.Beg,d.x.max=Int.End,Augment=Augment);dim(densy);plot(densy[,1],densy[,2])
dim(densy)
densy[densy[,1] > 5.8,]
}

#######################################################
### End of Get Densities
#######################################################




#######################################################
### Begin Old Fashioned Fast Zcurve 
#######################################################

old.fashioned.zcurve = function(z.val.input,cola = "springgreen2") {

### create set with z-scores in the interval used for model fitting
### THIS IS THE FUNCTION THAT COMPARES OBSERVED TO PREDICTED Z-VALUE DISTRIBUTIONS

zcurve.fitting = function(theta,RetEst=FALSE)    {

### get the weights and rescale 
weight = theta
weight = weight/sum(weight)

if (fixed.false.positives > 0) weight = c(fixed.false.positives,weight*(1-fixed.false.positives))
sum(weight)

### compute the new estimated density distribution
z.est = c()
#i = 1
#weight = WT
#cbind(Dens[,i],weight)
#dim(Dens)
for (i in 1:n.bars) z.est[i] = sum(Dens[,i]*weight)

### compare to observed density distribution
#misfit = mean(abs(z.est-Z.Density.Y))
rmse = sqrt(mean((z.est-Z.Density.Y)^2))

### return either fit if continue or estimates if finished
value = rmse
if(RetEst) value = z.est

### showing the fitting of the function in a plot
### showing the fitting of the function in a plot
if(Plot.Fitting) {

	rval = runif(1)
	if (rval > .9) {

	tit = ""
	xL = ""
	yL = ""
	plot(Z.Density.X,Z.Density.Y*scale,type='l',
		xlim=c(x.lim.min,x.lim.max),ylim=c(ymin,ymax),
		main=tit,xlab=xL,ylab=yL,axes=FALSE)
	lines(Z.Density.X,z.est*scale,lty=1,col="red1",
		xlim=c(x.lim.min,x.lim.max),ylim=c(ymin,ymax),axes=FALSE)
	#points(Z.Density.X,z.est,pch=20,col="red1",ylim=c(0,ymax),)

	}

} ### End of Plot Fitting

### return value to optimization function
return(value)

} ### End of zcurve.fitting

############################
#### End of Fitting Function
############################

### TEST VALUES
#cola = "orange"
#Augment.Factor = 1
#Int.Beg = 0.6
#zsd = 1.7
#x.lim.min = 0.6
#z.val.input = abs(c(rnorm(10000,0,zsd),rnorm(5000,2.8,zsd),rnorm(0,5)))
#ncz = c(0,1,2,3,4,5,6);zsd = 1;fixed.false.positives = 0;components = length(ncz)
#zsd = 1

components = length(ncz)


### get the densities for each interval and each non-centrality parameter

Z.INT = z.val.input[z.val.input >= Int.Beg & z.val.input <= Int.End]
summary(Z.INT)
#hist(Z.INT)

densy = Get.Densities(Z.INT,bw=bw.est,d.x.min=Int.Beg,d.x.max=Int.End,Augment=Augment)

Z.Density.X = densy[,1]
Z.Density.Y = densy[,2]

#plot(Z.Density.X,Z.Density.Y,type="l")

n.bars.int.beg = length(Z.Density.X)
n.bars = n.bars.int.beg;n.bars

bar.width = Z.Density.X[2] - Z.Density.X[1]
bar.width

### Finish getting observed densities 

zsds = rep(zsd,components)

#print("Create Dens")
Dens	= c()
for(i in 1:n.bars) {
	for (j in 1:length(ncz)) {
		Dens = c(Dens,dnorm(Z.Density.X[i],ncz[j],zsds[j]))
	}
}
summary(Dens)
Dens = matrix(Dens,length(ncz),byrow=FALSE)
sum.dens = rowSums(Dens)
Dens = Dens/(sum.dens * bar.width)
dim(Dens)

startval = rep(1/(components),components)
startval[1] = 1
startval = startval/sum(startval)
startval

lowlim = c(rep(0,components),ncz,zsds);lowlim
highlim = c(rep(1,components),ncz,zsds);highlim

if (fixed.false.positives > 0 & 0 %in% ncz) {
  startval = rep(1/(components-1),components-1)
  lowlim = c(rep(0,components-1),ncz,zsds);lowlim
  highlim = c(rep(1,components-1),ncz,zsds);highlim
}

startval
lowlim
highlim
ncz
Int.Beg
x.lim.min

#TESTING = TRUE
if (TESTING) Plot.Fitting = TRUE else Plot.Fitting = FALSE

d.hist.sel = mean(as.numeric(z.val.input > x.lim.min & z.val.input > Int.Beg & z.val.input < Int.End)) /
	mean(as.numeric(z.val.input > x.lim.min & z.val.input < Int.End))
d.hist.sel

scale = d.hist.sel

auto = nlminb(startval,zcurve.fitting,lower=lowlim,upper=highlim,
	control=list(eval.max=1000,abs.tol = 1e-20))

fit = auto$objective;fit

### get the estimated weights 
WT = auto$par
WT = WT/sum(WT)
if (fixed.false.positives > 0) WT = c(fixed.false.positives,WT*(1-fixed.false.positives))
sum(WT)

res = c(WT,fit);res


if (TESTING) {

cola = "forestgreen"

summary(Z.Density.X)
n.bars = length(Z.Density.X);n.bars
bar.width = Z.Density.X[2]-Z.Density.X[1]

Dens	= c()
for(i in 1:n.bars) {
	for (j in 1:length(ncz)) {
		Dens = c(Dens,dnorm(Z.Density.X[i],ncz[j],zsds[j]))
	}
}
Dens = matrix(Dens,length(ncz),byrow=FALSE)
sum.dens = rowSums(Dens)
Dens = Dens/(sum.dens * bar.width)

z.est = c()
for (i in 1:n.bars) z.est[i] = sum(Dens[,i]*WT)
z.est = z.est*scale

if (Show.Fitted) Draw.Histogram(z.val.input,WT,
	results=res,col=col.hist,Write.CI = FALSE)

par(new=TRUE)

lines(Z.Density.X,z.est,lty=1,col=cola,lwd=4,
 xlim=c(x.lim.min,x.lim.max),ylim=c(ymin,ymax))

}


###

if (TESTING | Show.Fitted) {

cp.res = Compute.Power(WT);cp.res;length(cp.res)
w.all = cp.res[8:(7+components)]
w.sig = cp.res[15:(14+components)]
WT
w.sig
w.all

Z.Density.X = seq(0,x.lim.max,bar.width)

n.bars = length(Z.Density.X);n.bars


Dens	= c()
for(i in 1:n.bars) {
	for (j in 1:length(ncz)) {
		Dens = c(Dens,dnorm(Z.Density.X[i],ncz[j],zsds[j]))
	}
}
Dens = matrix(Dens,length(ncz),byrow=FALSE)
sum.dens = rowSums(Dens)
Dens = Dens/(sum.dens * bar.width)

z.est = c()
for (i in 1:n.bars) z.est[i] = sum(Dens[,i]*w.all)
summary(z.est)
sum(z.est*bar.width)

d.hist.sel = mean(as.numeric(z.val.input > x.lim.min & z.val.input > Int.Beg & z.val.input < Int.End)) /
	mean(as.numeric(z.val.input > x.lim.min & z.val.input < Int.End))
d.hist.sel

table(Z.Density.X >= Int.Beg)
d.dense.sel = sum(z.est[Z.Density.X > x.lim.min & Z.Density.X >= Int.Beg]*bar.width)
d.dense.sel

scale = d.hist.sel/d.dense.sel;scale

par(new=TRUE)
lines(Z.Density.X,z.est*scale,lty=2,col="red3",lwd=4,
 xlim=c(x.lim.min,x.lim.max),ylim=c(ymin,ymax))

#print("Finished Fitting with Line")


} # End of Testing


res

return(res)

} ### End of old.fashioned.zcurve


######################################################
### End of Old Fashioned Zcurve
#######################################################



#######################################################
### extended.zcurve allows for free M and SD parameters
### slower!!!
#######################################################

extended.zcurve = function(z.val.input,ncz,zsd) {

##################
### THIS IS THE FUNCTION THAT COMPARES OBSERVED TO PREDICTED Z-VALUE DISTRIBUTIONS
##################

ext.zcurve.fitting = function(para,RetEst=FALSE,Fixed.Null=FALSE)    {


### get the weights
weights = para[1:length(ncz)];weights
means = para[(length(ncz)+1):(2*length(ncz))];means
sds = para[(length(ncz)*2+1):(length(ncz)*3)];sds

### Finish getting observed densities for Int.Beg to Maximum Z-Score in Interval

### get the densities for each interval and each non-centrality parameter
Dens	= c()
for(i in 1:n.bars) {
	for (j in 1:length(ncz)) {
		Dens = c(Dens,dnorm(Z.Density.X[i],means[j],sds[j]))
	}
}
Dens = matrix(Dens,length(ncz),byrow=FALSE)
dim(Dens)

### rescale the densities for the range of z-values to 1
sum.dens = rowSums(Dens)
Dens = Dens/(sum.dens * bar.width)

### compute the new estimated density distribution
z.est = c()
for (i in 1:n.bars) z.est[i] = sum(Dens[,i]*weights)

### compare to observed density distribution
#misfit = mean(abs(z.est-Z.Density.Y))
rmse = sqrt(mean((z.est-Z.Density.Y)^2))


### return either fit if continue or estimates if finished
value = rmse
if(RetEst) value = z.est

### showing the fitting of the function in a plot
### showing the fitting of the function in a plot
if(Plot.Fitting) {

	rval = runif(1)
	if (rval > .4) {

	tit = ""
	plot(Z.Density.X,Z.Density.Y,type='l',
		xlim=c(x.lim.min,x.lim.max),ylim=c(ymin,ymax),
		main=tit,xlab='Z')
	lines(Z.Density.X,z.est,lty=1,col="red1",
		xlim=c(x.lim.min,x.lim.max),ylim=c(ymin,ymax))
	#points(Z.Density.X,z.est,pch=20,col="red1",ylim=c(0,ymax),)

	}
}


### return value to optimization function
return(value)

} 

####

###ncz = c(0,2);zsd=1 

#z.val.input = rnorm(500,-.5,1)

### create set with z-scores in the interval used for model fitting
Z.INT = z.val.input[z.val.input >= Int.Beg & z.val.input <= Int.End]
summary(Z.INT)

#ncz = 2;zsd = 1
components = length(ncz);components

zsds = rep(zsd,components);zsds

#print(Augment)
densy = Get.Densities(Z.INT,bw=bw.est,d.x.min=Int.Beg,d.x.max=Int.End,Augment=Augment)

Z.Density.X = densy[,1]
Z.Density.Y = densy[,2]

n.bars = length(Z.Density.X)
n.bars

bar.width = Z.Density.X[2] - Z.Density.X[1]
bar.width


wstart = rep(1/components,components)
wstart[1] = 1
wstart = wstart/sum(wstart)
wstart

if (SD.FIXED) { 
	startval = c(wstart,ncz,rep(1,components));startval
	lowlim = c(rep(0,components),rep(0,components),rep(1,components))
	highlim = c(rep(1,components),rep(Int.End,components),rep(1,components))
} else {
	startval = c(wstart,ncz,zsds);startval
	lowlim = c(rep(0,components),rep(0,components),rep(1,components))
	highlim = c(rep(1,components),rep(Int.End,components),zsds)
}

#print(startval)

if(components == 1) { 
	startval = c(1,2,1)
	lowlim = c(1,0,0.5)
	highlim = c(1,6,zsds)
}

startval
lowlim
highlim

#ymax = 1.5
#TESTING = TRUE
if (TESTING == TRUE) Plot.Fitting = TRUE

auto = nlminb(startval,ext.zcurve.fitting,lower=lowlim,upper=highlim,control=list(eval.max=1000))
para = auto$par
para

fit = auto$objective;fit

res = c(para,fit)
res

return(c(res))

} 

### Testing 
### test = extended.zcurve(z.val.input,ncz,zsd);test
######################################################
### End of Extended Zcurve
#######################################################


#########################################################################
### This Function Computes Power from Weights and Non-Centrality Parameters
#########################################################################

Compute.Power = function(para,Int.Beg=z.crit,BOOT=FALSE) {

#para = para.est.OF

z.ext = z.extreme
if (two.sided) z.ext = z.ext[2] else z.ext = sum(z.ext);z.ext
#print("Check Extreme")
#print(z.ext)

### the input weights based on z.curve method
### these are the weights based on the Int.Beg value
### Int.Beg could be 1.96 (all significant)
### but it can also be other values
### Therefore the weights cannot be directly used
### to estimate power/replicability
w.inp = para[1:components]
ncz = para[(1+components):(2*components)]
zsds = para[(1+2*components):(3*components)]

#print(w.inp);print(ncz);print(zsds);print(z.ext)

### get power values for the components (ncz)
pow.dir = pnorm(abs(ncz),z.crit);pow.dir 

### get the opposite sign probability
sign.error = 1-pnorm(ncz,-z.crit);round(sign.error,3)

### this gives the power with the Int.Beg selection criterion as alpha 
### this power is used as a weight to get the weights for the full distribution
### using Jerry's insight that weight before selection is weight after selection divided by power
### get power values for the components (ncz)
pow.sel = pnorm(ncz,Int.Beg) + pnorm(-ncz,Int.Beg);pow.sel

### now we apply Jerry's third theorem in reverse to compute 
### the weights before selection (w.all)
### once we have the weights, we devided by sum of all weights 
### so that they add up to 1
w.all = w.inp / pow.sel
w.all = w.all / sum(w.all)
round(cbind(pow.sel,w.inp,w.all),3)
	
### now we are ready to compute the weights after selection for significance
### using Jerry's formula going from before selection to after selection
### by multiplying by power (w)
### again all the weights are standardized by dividing by the sum of all weights
w.sig = w.all * (pow.dir + sign.error)
w.sig = w.sig / sum(w.sig)
round(cbind(pow.sel,w.inp,w.all,w.sig),3)

#print("Check Two-Sided")
#print(two.sided)

if (two.sided) {

	### compute expected replication rate
	### this is easy, replicabilty is simply the weighted sum of power
	### using the weights after selection for significance, w.sig
	### [but limited to values below Int.End
	w.sig.ext = c(w.sig * (1-z.ext),z.ext);round(w.sig.ext,3)
	w.sig.ext
	sum(w.sig.ext)

	pow.sig.ext = c(pow.dir,1);pow.sig.ext

	ERR = sum(pow.sig.ext*w.sig.ext);ERR
	
	ERR.pos = NA
	ERR.neg = NA

	### compute expected discovery rate 
	### weights are reweighted so that extreme values are included with power = 1
	### as the weighted average of power using the weights before selection, w.all

	pow.ext = c(pow.dir+sign.error,1)
	round(pow.ext,3)

	w.all.ext = c(w.all*(1-z.ext),z.ext)
	cbind(round(pow.ext,3),round(w.all.ext,3))

	EDR = sum(w.all.ext*pow.ext);EDR
	EDR.pos = NA
	EDR.neg = NA

	#print("Check Compute Power")
	#print(c(EDR,ERR))
	#print("Check Compute Power")

} else {
	w.pos = c(w.all[ncz >= 0]*(1-sum(z.ext[2])),z.ext[2])
	w.pos = w.pos/sum(w.pos)
	pow.pos = c(pow[ncz >= 0],1)
	cbind(c(ncz[ncz >= 0],10),pow.pos,w.pos)
	EDR.pos = sum(w.pos*pow.pos)/sum(w.pos);EDR.pos

	EDR.neg = NA
	w.neg = c(z.ext[1],w.all[ncz <= 0]*(1-sum(z.ext[1])))	
	w.neg = w.neg/sum(w.neg)	
	pow.neg = c(1,pow[ncz <= 0])
	cbind(c(-10,ncz[ncz <= 0]),pow.neg,w.neg)
	EDR.neg = sum(w.neg*pow.neg);EDR.neg
}


### res stores the results to be past back from the function

res.est = c(EDR,EDR.pos,EDR.neg,ERR,ERR.pos,ERR.neg)
res.est
res = c(res.est,w.all,w.inp)
names(res) = c("EDR","EDR.pos","EDR.neg","ERR","ERR.pos","ERR.neg",
paste0("w.all.",ncz[1:components]),paste0("w.sig",ncz[1:components]) )
res

### now we compute mean power as a function of z-scores continuously
### this is only performed if local power is requested (int.loc > 0)
if (int.loc > 0) {

	bar.width = .01 # how fine should be the resolution
	Z.X = seq(x.lim.min,Int.End,bar.width);summary(Z.X)	 # set of z-scores 

	Z.W.D.Sum = unlist(lapply(Z.X,function(x) sum(dnorm(x,ncz)*w.all) ))
	#plot(Z.X,Z.W.D.Sum)

	loc.p = unlist(lapply(Z.X,function(x) 
		sum(dnorm(x,ncz)*w.all*(pow.dir+sign.error)) /	sum(dnorm(x,ncz)*w.all)	))

	#loc.p = .5
	#adj.z = qnorm(loc.p,z.crit)
	#plot(Z.X,adj.z,ylim=c(0,6),xlim=c(0,6))

	### compute local mean power for different intervals
	### each local power value is weighted by the density 
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


#print("Check res local power")
#print(length(res))
#print(res)

### to be past back to the main program from this function
return(res)

}

### Finish Compute.Power

#####################################
#####################################
#####################################


### Testing
if (1 == 2) {
w.inp = c(.5,.5)
ncz = c(1.5,2.5)
zsds = c(1.1,1.1)
}


Compute.Power.SDG1 = function(para,BOOT=FALSE) {

#para = para.est.EXT;para;BOOT=TRUE
#TEST
#ncz = c(1.1,2.8);zsds=c(1.05,1.05);w.inp=c(.5,.5);Int.Beg = 1.96;BOOT=FALSE


#print("Compute Power Start")


#print("Compute Power")
z.crit = qnorm(alpha/2,lower.tail=FALSE); z.crit
#print(z.crit)


z.ext = z.extreme
if (two.sided) z.ext = z.ext[2] else z.ext = sum(z.ext)
print("Check Extreme")
print(z.ext)


### the input weights based on z.curve method
### these are the weights based on the Int.Beg value
### Int.Beg could be 1.96 (all significant)
### but it can also be other values
### Therefore the weights cannot be directly used
### to estimate power/replicability
w.inp = para[1:components]
ncz = para[(1+components):(2*components)]
zsds = para[(1+2*components):(3*components)]

bar.width = .01 # how fine should be the resolution
zx = seq(x.lim.min,x.lim.max,bar.width)

###

components = length(ncz)

nczs = seq(x.lim.min,x.lim.max,.5)
pows = pnorm(nczs,z.crit) + pnorm(-nczs,z.crit)
cbind(nczs,pows)

comp.i = 1
set = c()

for (comp.i in 1:components) {
	w.nczs = dnorm(nczs,ncz[comp.i],zsds[comp.i]-1) 
	w.nczs = w.nczs/sum(w.nczs)
	#plot(nczs,w.nczs)
	subset = data.frame(nczs,pows,w.nczs)
	subset = subset[round(w.nczs,4) > 0,]
	subset$w.nczs = subset$w.nczs / sum(subset$w.nczs)
	subset$comp.i = comp.i
	subset$w.inp = w.inp[comp.i] 
	set = rbind(set,subset)
}
dim(set)
set

set$w.all = set$w.inp/set$pows
set$w.all = set$w.all/sum(set$w.all)
set$w.all


#local power is the true average power for an observed z-score 
#we need to get the power values for all components
#we need to get the weights of components 
#we also need to get the densities of the z-value for all components
#we create a double weighted weight using the weights of components and densities
#if there is more than one component, we also need to take the model weights into account
#the product of all 3 weights has to be scaled to 1
#we then compute the weighted average power 
#this is the local power for each observed z-value
#this can be used to compute edr and err

loc.pow = c()
for (zx.i in 1:length(zx)) {
	w.z = dnorm(zx[zx.i],set$nczs)
	w.z = w.z/sum(w.z)
	w = w.z * set$w.nczs * set$w.all
	w = w / sum(w)
	loc.pow = c(loc.pow,sum(w*set$pows)/sum(w))
} 

#plot(zx,loc.pow,ylim=c(0,1))

edr = mean(loc.pow);edr
err = mean(loc.pow[zx > z.crit]);err

ERR = err*(1-z.ext)+z.ext
EDR = edr*(1-z.ext)+z.ext

ERR
EDR

EDR.pos = NA
EDR.neg = NA

ERR.pos = NA
EDR.neg = NA

res.est = c(EDR,EDR.pos,EDR.neg,ERR,ERR.pos,ERR.neg)
res.est
res = c(res.est,rep(NA,components),rep(NA,components))
names(res) = c("EDR","EDR.pos","EDR.neg","ERR","ERR.pos","ERR.neg",
paste0("w.all.",ncz[1:components]),paste0("w.sig.",ncz[1:components]))
res

if (int.loc > 0 & BOOT==FALSE) {

	local.power = tapply(loc.pow,cut(zx,seq(x.lim.min-.0001,x.lim.max,int.loc),
		include.lowest=TRUE),mean)
	local.power
	length(local.power)
	names(local.power) = paste0("lp.",1:length(local.power))
	res = c(res,local.power)

} #EOF int.loc

res

#print("Compute Extended Power End")

### to be past back to the main program from this function
return(res)

} ### EOF Compute.Power


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

#z.val.input = fig1.zval

#Zing = function(z.val.input,lp=c() ) {

#lp = c()
if (length(lp) > 0) z.val.input = -qnorm(lp -log(2),log.p=TRUE)

z.crit = qnorm(alpha/2,lower.tail=FALSE); z.crit

components = length(ncz)

zsds = rep(zsd,length(ncz))[1:components];zsds


### remove NA
z.val.input = z.val.input[!is.na(z.val.input)];length(z.val.input)
z.val.input = z.val.input[z.val.input >= x.lim.min]
### limit range
z.val.input[z.val.input > MAX.INP.Z] = MAX.INP.Z
length(z.val.input)

if (length(z.val.input) > 10) print("!START!") else print("Insufficient Data")

print(paste0("Title: ",Title))

### create set with z-scores in the interval used for model fitting
Z.INT = z.val.input[z.val.input >= Int.Beg & z.val.input <= Int.End]
summary(Z.INT)

n.z = length(z.val.input);n.z
n.z.sig = length(z.val.input[abs(z.val.input) > z.crit])

ODR = n.z.sig/n.z;ODR

ODR.se = sqrt(ODR * (1 - ODR) / n.z)
ODR.se

ODR.low = round((ODR-1.96*ODR.se),2);ODR.low
ODR.high = round((ODR+1.96*ODR.se),2);ODR.high

z.extreme = c(
	length(z.val.input[z.val.input < -Int.Ext])/length(z.val.input),
	length(z.val.input[z.val.input > Int.Ext])/length(z.val.input)
)
names(z.extreme) = c("Ext.Neg","Ext.Pos")
(z.extreme*100)


### End of Preparation
###############################################
###############################################
###############################################


### OLD FASHIONED Z-CURVE 

if (Est.Method == "OF") {

#print("Fitting Old Fashioned")

para.est.OF = old.fashioned.zcurve(z.val.input)
fit = para.est.OF[components+1];fit

para.est.OF = c(para.est.OF[1:components],ncz,zsds)
para.est.OF
#weights = para.est.OF[1:components]

cp.res = Compute.Power(para.est.OF,Int.Beg=Int.Beg)
round(cp.res,3)

#print("Check cp.res")
#print(cp.res)

EDR = cp.res[1];EDR
ERR = cp.res[4];ERR

w.all = cp.res[which(substring(names(cp.res),1,5) == "w.all")]
w.all
sum(w.all)

w.sig = cp.res[which(substring(names(cp.res),1,5) == "w.sig")]
round(w.sig,3)

loc.power = cp.res[which(substring(names(cp.res),1,2) == "lp")]
loc.power

#print("Finish Old Fashioned")

}

#########################################
### EXTENDED Z-CURVE

if (Est.Method == "EXT") {
	#print("TESTING EXTENDED")
	print(ncz)
	components = length(ncz)
	zsds
	#zsd = 1
	para.est.EXT = extended.zcurve(z.val.input,ncz,zsd);para.est.EXT
	fit = para.est.EXT[(3*components+1)];fit
	para.est.EXT = para.est.EXT[1:(3*components)]

	para.est.EXT

	WT = para.est.EXT[1:components];WT
	WT = WT/sum(WT)
	WT

	ncz = para.est.EXT[(components+1):(2*components)]
	ncz

	zsds = para.est.EXT[(2*components+1):(3*components)]
	zsds

	if (max(zsds) < 1.1) {
		cp.res = Compute.Power(para.est.EXT)
	} else {
		cp.res = Compute.Power.SDG1(para.est.EXT)
	}			

	round(cp.res,3)	

	EDR = cp.res[1];EDR
	ERR = cp.res[4];ERR
	#OSR = cp.res[7];OSR

	#print("Check")
	#print(cp.res)

	w.all = cp.res[which(substring(names(cp.res),1,5) == "w.all")]
	w.all
	sum(w.all)

	w.sig = cp.res[which(substring(names(cp.res),1,5) == "w.sig")]
	round(w.sig,3)

	loc.power = cp.res[which(substring(names(cp.res),1,2) == "lp")]
	loc.power

#	print("Extended Version Completed")

}

### end of Extended Version 


### if requested, run slower EM
### Est.Method = "EM"
if(Est.Method == "EM" & boot.iter == 0) {
	#summary(z.val.input)
	#print("Fitting EM")
	z.res = run.zcurve(z.val.input,Est.Method="EM",boot.iter,
		Int.Beg=Int.Beg,Int.End=Int.End,parallel=parallel)
	summary(z.res)     
	para.est.EM = c(summary(z.res, type="parameters")$coefficients)
	para.est.EM = para.est.EM[c((components+1):(2*components),1:components)]
	para.est.EM = c(para.est.EM,rep(1,7))
	para.est.EM
	
	fit = z.res$fit$Q

	w.all =	 summary(z.res, type="parameters")$coefficients[(components+1):(2*components)]
	w.all

	cp.res = Compute.Power(para.est.EM,Int.Beg=0)

	w.sig = cp.res[which(substring(names(cp.res),1,5) == "w.sig")]
	round(w.sig,3)

	loc.power = cp.res[which(substring(names(cp.res),1,2) == "lp")]
	loc.power

	str(z.res)
	EDR = z.res$coefficients[2];EDR
	ERR = z.res$coefficients[1];ERR

	#print("Finished EM")
} 


FDR = round((1/EDR - 1)*(alpha/(1-alpha)),2);FDR

res = c(ODR,EDR,ERR,FDR,fit)
names(res) = c("ODR","EDR","ERR","FDR","FIT")

#print("Check res")
#print(res)
#print(WT)

#print("Finished Computations")
#res is needed for text in plot

p.bias = c(NA,NA,NA)
if(TEST4BIAS) { 
	p.bias = test.bias(w.all) 
	res = c(res,p.bias)
}
#print("bias test")
#print(p.bias)

##########################################
### This Code is Used to Create Graphic
##########################################

if (Show.Histogram) { 

	Draw.Histogram(z.val.input,w.all,results=res,cola=col.hist)

	if (Show.KD) Draw.KD(z.val.input,w.all,cola=col.kd)

	if (Show.Zcurve.All) {
		Draw.Zcurve.All(z.val.input,w=w.all,ncz=ncz,zsds=zsds,cola=col.zcurve,
			Ltype=3,Lwidth = 1,x.start=x.lim.min,x.end=x.lim.max)
		Draw.Zcurve.All(z.val.input,w=w.all,ncz=ncz,zsds=zsds,cola=col.zcurve,
			Ltype=1,Lwidth = 3,x.start=Int.Beg,x.end=Int.End)
		}

	if (Show.Zcurve.Sig) {
		Draw.Zcurve.Sig(z.val.input,ncz=ncz,zsds=zsds,w.sig,cola=col.zcurve,
			Ltype=3,x.start=x.lim.min,x.end=x.lim.max)
		Draw.Zcurve.Sig(z.val.input,ncz=ncz,zsds=zsds,w.sig,cola=col.zcurve,
			Ltype=1,x.start=Int.Beg,x.end=Int.End)
		}

	if (length(loc.power > 0)) Write.Local.Power(loc.power)

	} # End of Show.Histogram


########################################################################
########################################################################
########################################################################

### code for confidence intervals 

########################################################################
########################################################################
########################################################################







### If Confidence Intervals are requested, compute CI (boot.iter > 0)
if (boot.iter > 0 & Est.Method != "EXT") {

	res = get.ci.info(Est.Method = Est.Method)
	res = c(ODR,ODR.low,ODR.high,res[1:9],p.bias)
	names(res)[13:15] = c("pbias obs","pbias exp","pbias p")
	round(res,3)
}


#Testing: boot.iter = 20	
if (boot.iter > 0 & Est.Method == "EXT") {

	res.ci = get.ci.info(Est.Method = Est.Method)
	res.ci = rbind(c(ODR.low,ODR.high),res.ci)
	rownames(res.ci)[1] = "ODR"
	res.ci
	dim(res.ci)

	res = c(ODR,EDR,ERR,FDR,fit,ncz,zsds,w.all,w.sig)
	names(res) = rownames(res.ci)
	length(res)

	res = cbind(res,res.ci)

	res = rbind(res,c(p.bias))
	rownames(res)[10] = "BIAS"
	round(res,3)

}


if (boot.iter > 0) {

	if (Show.Histogram) { 
		Draw.Histogram(z.val.input,w.all,
		results=res,cola=col.hist,Write.CI = TRUE)
		if (Show.KD) Draw.KD(z.val.input,w.all,cola=col.kd)
		if (Show.Zcurve.All) {
			Draw.Zcurve.All(z.val.input,w=w.all,ncz=ncz,zsds=zsds,cola=col.zcurve,
				Ltype=3,x.start=x.lim.min,x.end=x.lim.max)
			Draw.Zcurve.All(z.val.input,w=w.all,ncz=ncz,zsds=zsds,cola=col.zcurve,
				Ltype=1,x.start=Int.Beg,x.end=Int.End)
			}
		if (Show.Zcurve.Sig) {
			Draw.Zcurve.Sig(z.val.input,ncz=ncz,zsds=zsds,w.sig,cola=col.zcurve,Ltype=3,x.start=x.lim.min)
			Draw.Zcurve.Sig(z.val.input,ncz=ncz,zsds=zsds,w.sig,cola=col.zcurve,Ltype=1,x.start=Int.Beg)
			}
		#	par(family = fam[1])	
		if (length(loc.power > 0)) Write.Local.Power(loc.power)
	} # End of Show.Histogram	

} # End of boot.iter condition

### return results for further analysis 
if (Return.Weights) {
	start = length(res)	
	res = c(res,ncz,zsds,w.all,w.sig);length(res)
	names(res)[(start+1):length(res)] = c(paste0("M",1:length(ncz)),paste0("SD",1:length(ncz)),
		paste0("w.all.",ncz[1:length(ncz)]),paste0("w.sig.",ncz[1:length(ncz)]) )
	}

#print(round(res,3))

if (TEST4HETEROGENEITY) {
	res.het = test.heterogeneity(z.val.input)
	res = c(res,c(res.het))
}

return(res)


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

