

#rm(list = ls())


options(scipen=999)

home_dir = "C:/Users/ulric/Dropbox/PHPCurves/DOCS/z-curve/Simulations"
setwd(home_dir)
getwd()

library(ashr)
library(weightr)
library(metafor)


zcurve = "zcurve.V3.90.R"
source(zcurve)

source("SimpleMix.R")


### helper function

trunc_moments <- function(mu, tau, min_es = 0) {

    lower = min_es
    a   <- (lower - mu) / tau
    Z   <- 1 - pnorm(a)                 # P(X > lower)
    lam <- dnorm(a) / Z                 # inverse Mills ratio
    m   <- mu + tau * lam               # truncated mean
    med <- mu + tau * qnorm(pnorm(a) + 0.5 * Z)   # truncated median
    v   <- tau^2 * (1 + a * lam - lam^2) # truncated variance
    ret = c(mean = m[1], med = med[1], sd = sqrt(v)[1]) #, p_above = Z)
    ret
}




##########################################################

#Load Simulation design (based on Z-curve.2.0 simulations U) 

##########################################################

block <- readRDS(file = "settings_U.RDS")
dim(block)

sel = (block$N == 50) & (block$k == 1000)
block = block[sel,]

level.k = length(names(table(block$k)))
level.d = length(names(table(block$d)))
level.SD = length(names(table(block$SD)))
level.N = length(names(table(block$N)))
level.FDR = length(names(table(block$FDR)))
cells = level.k * level.d * level.SD * level.N * level.FDR
cells

dim(block)

settings = block

settings$seed = 1:nrow(settings)

settings$p_h0 = settings$FDR / 3 * 4
table(settings$p_h0)

settings$N_low   = 99
settings$N_high  = 101
settings$N       = 100
settings$N_shape = .8

settings$bias_shape = "no_bias"
settings$shape = 30

settings$pop_dist = "normal"
settings$dist_shape = 5        # df for t.dist, shape for gamma
settings$beta1 = 1
settings$beta2 = 10

settings$pub_bias = .7

level.k = length(names(table(settings$k)))
level.d = length(names(table(settings$d)))
level.SD = length(names(table(settings$SD)))
level.N = length(names(table(settings$N)))
level.p_h0 = length(names(table(settings$p_h0)))
level.bias_shape = length(names(table(settings$bias)))
level.pop_dist = length(names(table(settings$pop_dist)))

cells = level.k * level.d * level.SD * level.N * level.p_h0 * level.bias_shape * level.pop_dist
cells

dim(settings)

z_bins = seq(0,6,1)

only.pos = FALSE
sel_crit = 0
min_z = -8
max_z = 8
min_d = -2.5
max_d = 2.5

which(settings$d[1:64] == .4)
which(settings$SD[1:64] == .6)


fn_dat = "SimpleMix_all_normal_ZC_ASHR_VW_RMA.dat"


###################################################################
### Run Simulation 
###################################################################

source(zcurve)

Est.Method = "EM"
boot.iter = 0
Directional = TRUE
Folded = FALSE

Int.Beg = -8
Int.End = 8
x.lim.min = Int.Beg
x.lim.max = Int.End
ncp = seq(Int.Beg,Int.End,1);length(ncp)
components = length(ncp)
zsds = rep(1,components)

###

res = c()

### you can break and look at results and then continue with the next run
dim(res)

i = b
i

b = nrow(res) + 1
b

e = cells
e = 700

b = 1
e = nrow(settings)

b
e

i = 45
i = b

### 
if (1 == 2) {
  pick = which(settings$k == 1000 & settings$d == .4 & settings$SD == .4 & settings$N == 100 
    & settings$p_h0 == 0 & settings$dist == "t.dist" & settings$bias_shape == "no_bias")
  i = pick[1]
  i 
}
###

### Start of simulation run 
for (i in b:e) {
   
  print(paste0("Working on ",i))
  print(settings[i,])

  set.seed(settings$seed[i])   # reproducible per replication

  k        = settings$k[i]
  d        = settings$d[i]
  SD       = settings$SD[i]
  p_h0     = settings$p_h0[i]

  N        = settings$N[i]
  N_low    = settings$N_low[i]
  N_high   = settings$N_high[i]
  N_shape  = settings$N_shape[i]

  pop_dist             = settings$pop_dist[i]
  pop_dist_beta_shape1 = settings$beta1[i]
  pop_dist_beta_shape2 = settings$beta2[i]
  bias                 = settings$bias[i]
  bias_shape           = settings$bias_shape[i]
  shape                = settings$shape[i]
  pub_bias             = settings$pub_bias[i]
  dist_shape           = settings$dist_shape[i]

  #source("SimpleMix.R")

  sim.res <- SimpleMix(k=k,d=d,SD=SD,N=N,N_low=N_low,
       N_high=N_high,p_h0=p_h0,only.pos=only.pos,
       min_z = min_z,max_z = max_z, min_d = min_d, max_d = max_d,
       pop_dist = pop_dist,dist_shape = dist_shape,
       pop_dist_beta_shape1 = pop_dist_beta_shape1,
       pop_dist_beta_shape2 = pop_dist_beta_shape2,
       )

  sim.dat = data.frame(sim.res$dat)

  sel1 = sim.dat$z > min_z & sim.dat$z < max_z
  sel2 = sim.dat$pop_d >= min_d & sim.dat$pop_d < max_d
  table(sel1,sel2)
  sel = sel1 & sel2
  table(sel)

  sim.dat = sim.dat[sel,]

  #hist(sim.dat$pop_d)
  summary(sim.dat$obs_d)
  summary(sim.dat$pop_d)
  sd(sim.dat$pop_d)
  summary(sim.dat$z)

  ### No abs if only.pos FALSE and Directional TRUE
  z.crit = qnorm(1-alpha/2)
  true_mean_all   = mean(sim.dat$pop_d)
  true_median_all = median(sim.dat$pop_d)
  true_mean_sig   = mean(sim.dat$pop_d[sim.dat$z > z.crit])
  true_tau        = sd(sim.dat$pop_d)
  table(cut(sim.dat$z,seq(Int.Beg,Int.End,1)))
  true_local_es = tapply(sim.dat$pop_d,cut(sim.dat$z,seq(Int.Beg,Int.End,1)),mean,na.rm=TRUE)

  sim_run = c(
         run = i,
         true_mean_all   = true_mean_all,
         true_median_all = true_median_all,   
         true_mean_sig   = true_mean_sig,
         true_tau        = true_tau,
         true_local_es   = true_local_es
        ) 

  #sim_run
  #########

  Title = paste0("Run: ",i)

  val.input   = sim.dat$z
  yi          = sim.dat$obs_d
  sei         = sim.dat$se 

  #hist(val.input)
  #summary(val.input) #ncp
  #sd(sim.dat$pop_d)
  #sqrt(var(val.input)-1)*mean(sim.dat$se)
  #var(val.input)-1

  #source(zcurve)

  do.zcurve = TRUE
  if(do.zcurve) {
  z.res  <- Zing(
    val.input = val.input,
    yi        = yi, 
    sei       = sei
   )

  } else {

    z.res = NULL

  }


  #z.res$w.all
  #z.res$es_tau
  #sum(z.res$w.all[1,])

  if(!is.null(z.res)) {
    z_run <- c(
      z.res$es_mean_all, #3
      z.res$es_mean_sig, #6
      z.res$local_es[1,],#22
      z.res$local_es[2,],#38
      z.res$local_es[3,],#54
      z.res$w.all[1,], # 71
      z.res$es_tau     # 74  
    )#
  } else { 
      z_run = rep(NA,74)
  }

  #round(z_run,2);length(z_run)

  ###
  betahat = sim.dat$obs_d 
  sebetahat = sim.dat$se 

  #z-curve version of ash
  ash.zcurve = TRUE
 
  if(ash.zcurve) {

    K = (length(ncp)-1)
    ncp2 = seq(-1,1,1/K)
    ncp2
    K2 = length(ncp2)

    ash_res = ash(   
      betahat = betahat,
      sebetahat = sebetahat,
      g = normalmix(
          pi          = rep(1/K2, K2),
          mean        = ncp2,
          sd          = rep(.01, K2)),  # ~point masses
          fixg        = FALSE,            # estimate the weights
          prior       = "uniform",        # <- penalty OFF, so it matches unpenalized ML
          mixcompdist = "normal",
          optmethod   = "mixEM",    
          control = list(tol = 1e-8,maxiter = 5000, trace = FALSE)
       )

     #summary(ash_res)
   
    } else { 

      ash_res = ash(   
        betahat = betahat,
        sebetahat = sebetahat,
        optmethod   = "mixEM",    
        control = list(tol = 1e-8,maxiter = 1000, trace = FALSE)
      )

     #summary(ash_res)

   }


   g <- ash_res$fitted_g          # the fitted mixture prior
   # g has components: g$pi (weights), g$mean (component means), g$sd (component SDs)
   # mean of g:
   g_mean <- sum(g$pi * g$mean)
   # variance of g = within-component variance + between-component variance
   g_var  <- sum(g$pi * (g$sd^2 + g$mean^2)) - g_mean^2
   ash_tau <- sqrt(g_var)

   pm <- get_pm(ash_res)
   length(pm)

   ash_local_es = tapply(pm,cut(sim.dat$z,seq(Int.Beg,Int.End,1)),mean)
   ash_mean_all = mean(pm)
   ash_mean_sig = mean(pm[abs(sim.dat$z) > z.crit])

    ash_run = c(
      ash_mean_all = ash_mean_all,
      ash_mean_sig = ash_mean_sig,
      ash_local_es,
      ash_tau = ash_tau
    )
    #ash_run["ash_tau"]
    #true_tau


  ### vevea * woods weightr

  yi = sim.dat$obs_d
  vi = sim.dat$se^2

  #hist(yi)

  res.weightr <- tryCatch(
      {
        suppressWarnings(
          weightr::weightfunct(yi, vi, steps = c(.5,.025))
        )
      },
       error = function(e) {
        message("weightr failed: ", conditionMessage(e))
        return(NULL)
      }
    )

  #res.weightr

  if (!is.null(res.weightr)) {

    tau2_hat <- res.weightr$adj_est[1, 1]
    tau2_se  <- res.weightr$adj_se[1, 1]
    tau_hat <- sqrt(tau2_hat)
    tau2_ci <- tau2_hat + c(-1, 1) * 1.96 * tau2_se
    tau_ci <- sqrt(pmax(0, tau2_ci))
    tau <- c(tau_hat,tau_ci)

    vw_run = t(rbind(res.weightr$adj_est[2:4],
        res.weightr$ci.lb_adj[2:4],
        res.weightr$ci.ub_adj[2:4]))
    vw_run = c(vw_run[1,],tau,vw_run[2,],vw_run[3,])

    names(vw_run) = c(
    "vw_mean_pe","vw_mean_lb","vw_mean_ub",
    "vw_tau_pe","vw_tau_lb","vw_tau_ub",
    "vw_weight1_pe","vw_weight1_lb","vw_weight1_ub",
    "vw_weight2_pe","vw_weight2_lb","vw_weight2_ub")

    #vw_run
    min_es = 0
    moments.pe = trunc_moments(vw_run[1],vw_run[4],min_es = min_es)
    moments.lb = trunc_moments(vw_run[2],vw_run[5],min_es = min_es)
    moments.ub = trunc_moments(vw_run[3],vw_run[6],min_es = min_es)
    moments = c(rbind(moments.pe,moments.lb,moments.ub)) 
    names(moments) = c(
        "vw_mean_pos_pe","vw_mean_pos_lb","vw_mean__pos_ub",
        "vw_median_pos_pe","vw_median_pos_lb","vw_median_pos_ub",
        "vw_tau_pos_pe","vw_tau_pos_lb","vw_tau_pos_ub")
    #moments  
    
    vw_run = c(vw_run,moments)
    vw_run[is.nan(vw_run)] = NA

    if(length(vw_run) < 21) vw_run = rep(NA,21)


  } else {
  
    vw_run = rep(NA,21)

  }

  #vw_run["vw_tau_pe"]

  #### 

  # plain random-effects MA on the same signed effect sizes
  # yi = observed effects, vi = sampling variances (se^2)
  rma_fit <- rma(yi = sim.dat$obs_d, vi = sim.dat$se^2, method = "REML")
  rma_mean_pe <- as.numeric(rma_fit$b)      # estimated mean effect
  rma_tau    <- sqrt(rma_fit$tau2)               # estimated heterogeneity
  rma_run = c(rma_mean_pe = rma_mean_pe,rma_tau = rma_tau)
  rma_run


  ####

  res_run = c(
          sim_run,
          z_run,
          ash_run,
          vw_run,
          rma_run
         )   

  print(paste0("Run: ",i))

  print(res_run[c("true_mean_all","es_mean_all_pe",
         "ash_mean_all","vw_mean_pe","rma_mean_pe")])

  print(res_run[c("true_tau","es_tau_pe",
         "ash_tau","vw_tau_pe","rma_tau")])

  #res_run;length(res_run)

  if(length(res_run) != 137) stop("length wrong")
  
  ### dim(res)
  ### res = c()
  res = rbind(res,res_run)
  res = data.frame(res)

  colnames(res) = names(res_run)


  if(sum(!is.na(res$es_mean_all_pe)) > 2) print(
      mean(abs(res$true_mean_all - res$es_mean_all_pe),na.rm=TRUE)
     )
  if(sum(!is.na(res$ash_mean_all)) > 2) print(
      mean(abs(res$true_mean_all - res$ash_mean_all),na.rm=TRUE)
     )
  if(sum(!is.na(res$vw_mean_pe)) > 2) print(
      mean(abs(res$true_mean_all - res$vw_mean_pe),na.rm=TRUE)
     )
  if(sum(!is.na(res$rma_mean_pe)) > 2) print(
      mean(abs(res$true_mean_all - res$rma_mean_pe),na.rm=TRUE)
     )

  print("tau, tau, tau, tau")

  if(sum(!is.na(res$es_mean_all_pe)) > 2) print(
      mean(abs(res$true_tau - res$es_tau_pe),na.rm=TRUE)
     )
  if(sum(!is.na(res$ash_mean_all)) > 2) print(
      mean(abs(res$true_tau - res$ash_tau),na.rm=TRUE)
     )
  if(sum(!is.na(res$vw_mean_pe)) > 2) print(
      mean(abs(res$true_tau - res$vw_tau_pe),na.rm=TRUE)
     )
  if(sum(!is.na(res$rma_mean_pe)) > 2) print(
      mean(abs(res$true_tau - res$rma_tau),na.rm=TRUE)
     )

 
  if (i %% 64 == 0) {
    write.table(
      res,
      fn_dat,
      row.names = FALSE
    )
  } 

} # End of sim loop 


res = data.frame(res)
dim(res)



##########################################

##########################################

##########################################

##########################################

##########################################

home_dir = "C:/Users/ulric/Dropbox/PHPCurves/DOCS/z-curve/Simulations"
setwd(home_dir)
getwd()

#rres = res
colnames(res)

fn_dat
xres = data.frame(read.table(fn_dat,,header=TRUE))
dim(xres)

#res = xres

cells

finished = trunc(nrow(res) / cells)
finished # finished = 8

keep = cells*finished
keep

res = res[1:keep,]
dim(res)

sel.set = settings[1:nrow(res),]
dim(sel.set)

comp = data.frame(sel.set,res)
dim(comp)


### point estimtes for all

sel = rep(TRUE,nrow(comp))

colnames(comp)
DV = comp$true_mean_all
tab.truth = tapply(DV[sel],list(comp$d[sel],comp$SD[sel],comp$FDR[sel]),mean,na.rm=TRUE)
tab.truth = as.data.frame.table(tab.truth)
tab.truth = as.data.frame(tab.truth)
tab.truth

colnames(comp)
DV = comp$es_mean_all_pe
tab.z.pe = tapply(DV[sel],list(comp$d[sel],comp$SD[sel],comp$FDR[sel]),mean,na.rm=TRUE)
tab.z.pe = as.data.frame.table(tab.z.pe)
tab.z.pe = as.data.frame(tab.z.pe)
tab.z.pe

DV = comp$ash_mean_all
tab.ash.all.pe = tapply(DV[sel],list(comp$d[sel],comp$SD[sel],comp$FDR[sel]),mean,na.rm=TRUE)
tab.ash.all.pe = as.data.frame.table(tab.ash.all.pe)
tab.ash.all.pe = as.data.frame(tab.ash.all.pe)
tab.ash.all.pe

DV = comp$vw_mean_pe
tab.vw.mean.pe = tapply(DV[sel],list(comp$d[sel],comp$SD[sel],comp$FDR[sel]),mean,na.rm=TRUE)
tab.vw.mean.pe = as.data.frame.table(tab.vw.mean.pe)
tab.vw.mean.pe = as.data.frame(tab.vw.mean.pe)
tab.vw.mean.pe

DV = comp$rma_mean_pe
tab.rma.mean.pe = tapply(DV[sel],list(comp$d[sel],comp$SD[sel],comp$FDR[sel]),mean,na.rm=TRUE)
tab.rma.mean.pe = as.data.frame.table(tab.rma.mean.pe)
tab.rma.mean.pe = as.data.frame(tab.rma.mean.pe)
tab.rma.mean.pe


tab.comp = data.frame(
   tab.truth[,1:4],
   mean.truth   = round(tab.truth[,4],2),
   mean.z.est   = round(tab.z.pe[,4],2),
   mean.ash.est = round(tab.ash.all.pe[,4],2),
   mean.vw.est  = round(tab.vw.mean.pe[,4],2),
   mean.rma.est = round(tab.rma.mean.pe[,4],2)
  )

summary(tab.comp[,5:9])

summary(tab.z.pe[,4] - tab.truth[,4])
summary(tab.ash.all.pe[,4] - tab.truth[,4])
summary(tab.vw.mean.pe[,4] - tab.truth[,4])
summary(tab.rma.mean.pe[,4] - tab.truth[,4])

summary(sqrt(mean((tab.z.pe[,4] - tab.truth[,4])^2)))
summary(sqrt(mean((tab.ash.all.pe[,4] - tab.truth[,4])^2)))
summary(sqrt(mean((tab.vw.mean.pe[,4] - tab.truth[,4])^2)))
summary(sqrt(mean((tab.rma.mean.pe[,4] - tab.truth[,4])^2)))

dim(tab.comp)



###

colnames(comp)
DV = comp$true_tau
tau.truth = tapply(DV[sel],list(comp$d[sel],comp$SD[sel],comp$FDR[sel]),mean,na.rm=TRUE)
tau.truth = as.data.frame.table(tau.truth)
tau.truth = as.data.frame(tau.truth)
tau.truth

colnames(comp)
DV = comp$es_tau_pe
tau.z.pe = tapply(DV[sel],list(comp$d[sel],comp$SD[sel],comp$FDR[sel]),mean,na.rm=TRUE)
tau.z.pe = as.data.frame.table(tau.z.pe)
tau.z.pe = as.data.frame(tau.z.pe)
tau.z.pe

DV = comp$ash_tau
tau.ash = tapply(DV[sel],list(comp$d[sel],comp$SD[sel],comp$FDR[sel]),mean,na.rm=TRUE)
tau.ash = as.data.frame.table(tau.ash)
tau.ash = as.data.frame(tau.ash)
tau.ash

DV = comp$vw_tau_pe
tau.vw = tapply(DV[sel],list(comp$d[sel],comp$SD[sel],comp$FDR[sel]),mean,na.rm=TRUE)
tau.vw = as.data.frame.table(tau.vw)
tau.vw = as.data.frame(tau.vw)
tau.vw

DV = comp$rma_tau
tau.rma = tapply(DV[sel],list(comp$d[sel],comp$SD[sel],comp$FDR[sel]),mean,na.rm=TRUE)
tau.rma = as.data.frame.table(tau.rma)
tau.rma = as.data.frame(tau.rma)
tau.rma

summary(sqrt(mean((tau.z.pe[,4] - tau.truth[,4])^2)))
summary(sqrt(mean((tau.ash.all.pe[,4] - tau.truth[,4])^2)))
summary(sqrt(mean((tau.vw[,4] - tau.truth[,4])^2)))
summary(sqrt(mean((tau.rma[,4] - tau.truth[,4])^2)))


tab.comp = data.frame(tab.comp,
   tau.truth   = round(tau.truth[,4],2),
   tau.z.est   = round(tau.z.pe[,4],2),
   tau.ash.est = round(tau.ash[,4],2),
   tau.vw.est  = round(tau.vw[,4],2),
   tau.rma.est = round(tau.rma[,4],2)
  )

dim(tab.comp)

colnames(tab.comp)


###

#png("sim.signed.tdist5.mean.png", width = 8, height = 8, units = "in", res = 300)

par(mfrow = c(2, 2), mar = c(4, 4, 3, 1))

plot_panel <- function(x, y, main, col = "black") {
  plot(
    x, y,
    pch = 21,
    bg = adjustcolor(col, alpha.f = .35),
    col = adjustcolor("black", alpha.f = .35),
    xlim = c(0, 1),
    ylim = c(-.1, .8),
    xlab = "True target mean",
    ylab = "Estimated mean",
    main = main
  )
  abline(0, 1, lty = 2)
  abline(h = 0, lty = 1)
}


colnames(tab.comp)

plot_panel(tab.comp$mean.truth, tab.comp$mean.z.est,
           "Z-curve (signed-EM)", "black")

plot_panel(tab.comp$mean.truth, tab.comp$mean.ash.est,
           "ash-package (EM)", "darkgreen")

plot_panel(tab.comp$mean.truth, tab.comp$mean.vw.est,
           "weightr-package", "purple")

plot_panel(tab.comp$mean.truth, tab.comp$mean.rma.est,
           "rma (metafor-package", "darkred")


#dev.off()


###

#png("sim.signed.tdist5.tau.png", width = 8, height = 8, units = "in", res = 300)

par(mfrow = c(2, 2), mar = c(4, 4, 3, 1))

plot_panel <- function(x, y, main, col = "black") {
  plot(
    x, y,
    pch = 21,
    bg = adjustcolor(col, alpha.f = .35),
    col = adjustcolor("black", alpha.f = .35),
    xlim = c(0, 1),
    ylim = c(-.1, .8),
    xlab = "True target tau",
    ylab = "Estimated tau",
    main = main
  )
  abline(0, 1, lty = 2)
  abline(h = 0, lty = 1)
}


plot_panel(tab.comp$tau.truth, tab.comp$tau.z.est,
           "Z-curve (signed-EM)", "black")

plot_panel(tab.comp$tau.truth, tab.comp$tau.ash.est,
           "ash-package (EM)", "darkgreen")

plot_panel(tab.comp$tau.truth, tab.comp$tau.vw.est,
           "weightr-package", "purple")

plot_panel(tab.comp$tau.truth, tab.comp$tau.rma.est,
           "rma (metafor-package", "darkred")


#dev.off()



###




