

#rm(list = ls())

options(scipen = 999) 

library(weightr)
source("https://raw.githubusercontent.com/UlrichSchimmack/zcurve3.0/refs/heads/main/run_boot_cluster_weightr.R")

zcurve3 <- "https://raw.githubusercontent.com/UlrichSchimmack/zcurve3.0/refs/heads/main/zcurve.V3.6.R"
source(zcurve3)

dat_file <- "https://raw.githubusercontent.com/UlrichSchimmack/zcurve3.0/refs/heads/main/data.csv"

#Read in data, making "" = NA
dat_all = data.frame(read.csv(dat_file,header=T, na.strings=c("",".","NA"))) 
table(dat_all$d > 2.5)
dat <- dat_all[which(dat_all$d < 2.5),]

dat$cluster.id = as.character(dat$uniquestudy)
tab = table(dat$cluster.id)
length(tab)

table(is.na(dat$cluster.id))
dat$se = sqrt(dat$mn)
dat$z = dat$d / dat$se
table(is.na(dat$z))
table(cut(dat$z,c(-10,-1.96,0,1.65,1.96,6,100)))

zval = data.frame(dat$z,dat$cluster.id)
table(complete.cases(zval))

### ALL

dat$sel = rep(TRUE,nrow(dat))
table(dat$sel)
res.es.all = run_boot_cluster_weightr(yi = dat$d[dat$sel], vi = dat$mn[dat$sel],
  cluster = dat$cluster.id[dat$sel],
  steps = c(.025,.05,.5,.975) )





### Subliminal
dat$sel = dat$Liminality.Binary == 0
table(dat$sel)
res.es.sub = run_weightr_cluster_boot(yi = dat$d[dat$sel], vi = dat$mn[dat$sel],
  cluster = dat$cluster.id[dat$sel],
  steps = c(.05, .025) )


source(zcurve3)
Title = "Subliminal Priming (Dai et al., 2023)"
set.seed = 1865
boot.iter = 1000
Est.Method = "EM"
int.loc = 2
ymax = 1
dat$sel = dat$Liminality.Binary == 0
res.sub = Zing(val.input = dat$z[dat$sel],
   cluster.id = dat$cluster.id[dat$sel])
res.sub$BIAS
res.sub$local.power

dat$sel = dat$Liminality.Binary == 0 & dat$z > 4
dat[dat$sel,c("Short.Title","n","d","se","z"),]





#install.packages("remotes")
#remotes::install_github("jepusto/metaselection")
#Sys.which("make")
#remotes::install_github("jepusto/metaselection", force = TRUE, verbose = TRUE)
#remotes::install_github("jepusto/metaselection", build_vignettes = FALSE, INSTALL_opts = "--no-multiarch")



res.es.all

dat$sel = dat$se > .4
mean(dat$d[dat$sel])
sd(dat$d[dat$sel])

dat$sel = dat$se > .3 & dat$se <= .4
mean(dat$d[dat$sel])
sd(dat$d[dat$sel])

dat$sel = dat$se > .2 & dat$se <= .3
mean(dat$d[dat$sel])
sd(dat$d[dat$sel])

dat$sel = dat$se > .1 & dat$se <= .2
mean(dat$d[dat$sel])
sd(dat$d[dat$sel])

dat$sel = dat$se > 0 & dat$se <= .1
mean(dat$d[dat$sel])
sd(dat$d[dat$sel])

table(cut(dat$se,c(0,.1,.2,.3,.4,1)))


ymax = 1.2
title = "Effect Size Distribution" 
curve(main=title,dnorm(x,.2044,.3993),xlim=c(-1,2),ylim=c(0,ymax),col="grey0"
 ,ylab="Density",lwd=4,lty=2)
par(new = TRUE)
curve(dnorm(x,.4044,.5723),xlim=c(-1,2),ylim=c(0,ymax),col="grey20",ylab="",lty=1)
par(new = TRUE)
curve(dnorm(x,.4499,.6063),xlim=c(-1,2),ylim=c(0,ymax),col="grey40",ylab="",lty=1)
par(new = TRUE)
curve(dnorm(x,.3393,.4258),xlim=c(-1,2),ylim=c(0,ymax),col="grey40",ylab="",lty=1)
par(new = TRUE)
curve(dnorm(x,.2856,.3478),xlim=c(-1,2),ylim=c(0,ymax),col="grey40",ylab="",lty=1)
par(new = TRUE)
curve(dnorm(x,.2856,.3478),xlim=c(-1,2),ylim=c(0,ymax),col="grey40",ylab="",lty=1)
#par(new = TRUE)
#curve(dnorm(x,.048,.077),xlim=c(-1,2),ylim=c(0,ymax),col="grey40",ylab="",lty=2)
legend(1,1,legend=c("Latent True Effects","Observed Effects"),lty=c(2,1),lwd=c(4,1))



## =========================================================
## Effect-size distribution plot with colored SE-band fills
## =========================================================

#pop.mu = .2; pop.tau = .3

run_prediction_interval_plot = function{
    pop.mu, 
    pop.tau, 
    yi,
    show_band <- c(FALSE, TRUE, TRUE, TRUE, TRUE) ) {


## SE bands
## first band (0 to .1) is included here but not plotted by default
band_labels <- c("SE 0-.1", "SE .1-.2", "SE .2-.3", "SE .3-.4", "SE > .4")

## Representative SE for each band
## For the open-ended band (> .4), choose a reasonable representative value
se_rep <- c(0.05, 0.15, 0.25, 0.35, 0.50)

## Colors (dark = more precise, light = less precise)
band_cols <- c("#66A61E", "#E6AB02", "#1B9E77",  "#E7298A", "#7570B3")

## Transparency
fill_alpha <- rep(.1,5)

## ---- x grid ----
x <- seq(-1, 2, length.out = 1200)

## ---- Population distribution ----
pop_y <- dnorm(x, mean = pop.mu, sd = pop.tau)

dat$se_band <- cut(
  dat$se,
  breaks = c(0, .1, .2, .3, .4, Inf),
  labels = c("SE 0-.1", "SE .1-.2", "SE .2-.3", "SE .3-.4", "SE > .4"),
  right = TRUE,
  include.lowest = TRUE
)


## ---- Observed mean and SD by SE band ----

band_stats <- aggregate(
  yi ~ se_band,
  data = dat,
  FUN = function(x) c(
    n = sum(!is.na(x)),
    mean = mean(x, na.rm = TRUE),
    sd = sd(x, na.rm = TRUE)
  )
)

## ---- Clean up aggregate output ----

band_stats <- do.call(
  data.frame,
  band_stats
)

names(band_stats) <- c("se_band", "n", "mean", "sd")

band_stats

## ---- Observed density curves from empirical band means/SDs ----

sample_y <- lapply(seq_len(nrow(band_stats)), function(i) {
  dnorm(
    x,
    mean = band_stats$mean[i],
    sd   = band_stats$sd[i]
  )
})

names(sample_y) <- band_stats$se_band

## ---- y-axis limits ----
y_max <- max(
  pop_y,
  unlist(sample_y[show_band])
)

## ---- Plot setup ----
par(mar = c(5, 5, 4, 2) + 0.1)

plot(
  x, pop_y,
  type = "n",
  xlab = "Effect Size (d)",
  ylab = "Density",
  main = "Effect Size Distribution",
  ylim = c(0, y_max * 1.10),
  xaxs = "i",
  yaxs = "i"
)

## ---- Draw filled sample distributions ----
## Draw widest curves first so smaller ones remain visible
draw_order <- rev(which(show_band))

for (i in draw_order) {
  polygon(
    x = c(x, rev(x)),
    y = c(sample_y[[i]], rep(0, length(x))),
    col = adjustcolor(band_cols[i], alpha.f = fill_alpha[i]),
    border = NA
  )
}

## ---- Add thin outlines for sample distributions ----
for (i in which(show_band)) {
  lines(x, sample_y[[i]], col = band_cols[i], lwd = 4)
}

## ---- Add population curve ----
lines(x, pop_y, col = "black", lwd = 3, lty = 2)

## ---- Legend ----
legend(
  "topright",
  inset = 0.03,
  legend = c("Population", band_labels[show_band]),
  col    = c("black", band_cols[show_band]),
  lty    = c(2, rep(1, sum(show_band))),
  lwd    = c(4, rep(1.8, sum(show_band))),
  pch    = c(NA, rep(21, sum(show_band))),
  pt.cex = 2,
  pt.bg  = c(NA, adjustcolor(band_cols[show_band], alpha.f = 1)),
  bty    = "n"
)





table(dat$d > .5 & dat$d < 1.5 & dat$se < .2, dat$z > 4)

table(!is.na(dat$d))
mean(dat$d > .5 & dat$se < .2)


sel = !(dat$d > .5 & dat$se < .2) & dat$z > 4
sel = (dat$d > .5 & dat$se < .2) & dat$z <4
table(sel)
tab = dat[sel,c("Short.Title","d","se","z")]
tab = tab[order(tab$z),]
tab

table(dat$d < 1.5, dat$se > .2)

sel = dat$d < 1.5 & dat$d > .5 & dat$se < .2 & dat$z > 4
table(sel)

colnames(dat)
tab = dat[sel,c("Short.Title","n","d","se","z")]
tab = tab[order(tab$z),]
tab

dat[dat$Short.Title == "Effect of watermarks",c("n","d","se","z","Title")]


#9 756 Effect of watermarks 0.5569614 0.1330380  4.186483

#8 170 Prime and Performanc 0.7767042 0.1627577  4.772151

#7 767 Hedonic Eating Goals 0.7938469 0.1873172  4.237983
#7 765 Hedonic Eating Goals 0.9281469 0.1873172  4.954948

#6 745 An (un)healthy poste 0.9029184 0.1781742  5.067617

#5 749 Call to (In)Action:  0.7467889 0.1290994  5.784602

#4 415 Why Does the ????Sin 0.6297085 0.1087857  5.788525

#3 140 Money Cues Increase  1.2214930 0.1760902  6.936747

#2 389 To lead or to be lik 0.7286567 0.1276867  5.706596
#2 818 To lead or to be lik 1.0951340 0.1276867  8.576724

#1 759 Exposure to national 1.3647889 0.1280369 10.659342



## ---- Optional direct labels near the right tail ----
## Uncomment if you want labels on the curves instead of relying only on legend

# label_x <- c(0.95, 1.05, 1.15, 1.25, 1.35)
# shown <- which(show_band)
# k <- 1
# for (i in shown) {
#   lx <- label_x[k]
#   ly <- approx(x, sample_y[[i]], xout = lx)$y
#   text(lx, ly, labels = band_labels[i], pos = 4, cex = 0.85, col = band_cols[i], xpd = NA)
#   k <- k + 1
# }



dat$sel = dat$se > .4
mean(dat$d[dat$sel])
sd(dat$d[dat$sel])

### Z-CURVE 

val.input = dat$z
cluster.id = dat$cluster.id

dat$sel = z > 0 & sqrt(dat$mn) < .3
table(dat$sel)

source(zcurve3)
Title = "Priming Behavior (Dai et al., 2023)"
int.loc = 1
boot.iter = 500
Est.Method = "EM"
res.all = Zing(val.input = dat$z,
   cluster.id = dat$cluster.id)


source(zcurve3)
Title = "Subliminal Priming (Dai et al., 2023)"
boot.iter = 500
Est.Method = "EM"
int.loc = 1
ymax = 1
dat$sel = dat$Liminality.Binary == 0
res.sub = Zing(val.input = dat$z[dat$sel],
   cluster.id = dat$cluster.id[dat$sel])


res.sub$BIAS

table(dat$z > 4,dat$sel)

dat$Title <- iconv(dat$Title, from = "", to = "UTF-8", sub = "?")
dat$Short.Title = substring(dat$Title,1,20)

dat$sel = dat$PrimeContent == 4
table(dat$sel)
source(zcurve3)
boot.iter = 500
Est.Method = "EM"
Title = "Morality Priming (Dai et al., 2023)"
ymax = 1
res.moral = Zing(val.input = dat$z[dat$sel],
   cluster.id = dat$cluster.id[dat$sel])





772  4.115077       2852 Men Percieve Anti-ma
691  4.129307        882 When Conflicts Are G
756  4.186483        621 Effect of watermarks

767  4.237983       2704 Hedonic Eating Goals
765  4.954948       2704 Hedonic Eating Goals

851  4.154398        693 The Effects of Chron
636  4.475843        694 The Effects of Chron

277  4.743429       3213 The Effect of a Cont
170  4.772151       3152 Prime and Performanc

647  4.796065       1481 The Influence of Alc

574  4.860575        702 Motivating Goal Purs

###

490  4.027501       1015 Changing, priming, a
495  5.096095       1015 Changing, priming, a

745  5.067617       2652 An (un)healthy poste


dat[355,]

colnames(dat)
tab = dat[dat$sel & dat$z > 4,c("z","cluster.id","Title")]
tab = tab[order(tab$z),]
tab

colnames(dat)

dat$PrimeContent = as.factor(dat$PrimeContent)
table(dat$PrimeContent)

summary(lm(z ~ 
  + se 
  + Year
  + Funneled.Debriefing.Present
  + Liminality.Binary
  + PrimeContent
  + Filler.Task
, data = dat))


#1. Chan (2019) (1 z > 4)
dat[grepl("Exposure to national flags reduces tax evasion",dat$Title),
  c("Short.Title","z","effectsize","d","se")]

#2. Case (2018) (2 z > 4)
seg = "To lead or to be liked"
dat[grepl(seg,dat$Title),
  c("Short.Title","z","effectsize","d","se")]

#3. 
seg = "The Impact of Priming on Speed Reduction on a Ski Slope"
dat[grepl(seg,dat$Title),
  c("Short.Title","z","effectsize","d","se")]

#4. Gasiorowska...Vohs (2016) (4 z > 4)
seg = "Money Cues Increase Agency and Decrease Prosociality Among Children"
dat[grepl(seg,dat$Title),
  c("Short.Title","z","effectsize","d","se")]


#5. Versluis and Papies (2016) (1 z > 4) error 
seg = "Eating less from big" 
dat[grepl(seg,dat$Title),
  c("Short.Title","z","effectsize","d","se")]


#6. Ding et al. (2016)
seg = "Why Does" 
dat[grepl(seg,dat$Title),
  c("Short.Title","z","effectsize","d","se")]

#7. Hassell (2015)
seg = "Call to" 
dat[grepl(seg,dat$Title),
  c("Short.Title","z","effectsize","d","se")]

#8. Ginsberg (2012)
seg = "Priming of disabilit"
dat[grepl(seg,dat$Title),
  c("Short.Title","z","effectsize","d","se")]

#9. Lin et al. (2016)
seg = "The Prosocial Impact"
dat[grepl(seg,dat$Title),
  c("Short.Title","z","effectsize","d","se")]

library(readxl)
lin = read_xlsx("Lin_focal_z_values.xlsx")
table(lin$Z > 4)
source(zcurve3)
Int.Beg = 0
ymax = 1
hist.bar.width = .5
NCP.FIXED = FALSE
ncp = 4
Zing(lin$Z)

