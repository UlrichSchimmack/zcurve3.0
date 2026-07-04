
rm(list = ls())

library(readxl)
library(weightr)
library(metafor)


#home_dir = "C:/Users/ulric/Dropbox/PHPCurves/DOCS/z-curve/Simulations"
#setwd(home_dir)
zcurve = "zcurve.V3.86.R"

# load directly from guithub
zcurve <- "https://raw.githubusercontent.com/UlrichSchimmack/zcurve3.0/main/zcurve.V3.86.R"

setwd("C:\\Users\\ulric\\Dropbox\\PHPCurves\\DOCS\\z-curve\\TerrorManagement")

### download data from OSF or my github

### Meta-analysis data file with wrong effect sizes
tm <- read_xlsx("data_meta.xlsx")
tm = data.frame(tm)
dim(tm)

### p-curve dataset with actual test results
tmz <- read_xlsx("data_pcurve.xlsx")
tmz = data.frame(tmz)
dim(tmz)


### Diagnose the coding mistake 

tm.merged = merge(tm,tmz,by="ID",all=TRUE)
dim(tm.merged)

sel = !is.na(tm.merged$yi.x) & !is.na(tm.merged$yi.y)
table(sel)

sel = tm.merged$Statistics == "F"
plot(tm.merged$yi.y[sel],tm.merged$yi.x[sel])
abline(v=c(0,1))
abline(h=c(0,1))

sel = tm.merged$Statistics == "t"
plot(tm.merged$yi.y[sel],tm.merged$yi.x[sel])
abline(v=c(0,1))
abline(h=c(0,1))

tapply(tm.merged$yi.x,tm.merged$Statistics,mean,na.rm=TRUE)
tapply(tm.merged$yi.y,tm.merged$Statistics,mean,na.rm=TRUE)


### analysis of meta data

hist(tm$yi,main="Histogram of Effect Size Estimates",xlab="Standardized Mean Difference")

mean(tm$yi > 0)

PET_RE <- rma(
  yi = yi,
  vi = vi,
  mods = ~ sei,
  method = "REML",
  data = tm
)

summary(PET_RE)

pet_tau = sqrt(PET_RE$tau2)
pet_tau


plot(
  tm$sei, tm$yi,
  xlim = c(0, max(tm$sei, na.rm = TRUE)),
  ylim = range(c(tm$yi, -0.11), na.rm = TRUE),
  xlab = "Standard error",
  ylab = "Effect size"
)

abline(PET)
abline(h = 0, lty = 2)
abline(v = 0, lty = 2)


tm$z = tm$yi / tm$sei
hist(tm$z)
abline(v = 0,lwd=3,col="red")

yi = tm$yi
vi = tm$vi
vw.rep = weightfunct(yi, vi, steps=c(.025,.05,.975) )
vw.rep

.36 + 1.95 * .705
.36 - 1.95 * .705


yi = tm$yi
vi = tm$vi
vw.fix = weightfunct(yi, vi, steps=c(.025,.05,.5) )
vw.fix

-.336 - 1.96 + .956
-.336 + 1.96 + .956


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

mom = trunc_moments(-.336,.956, 0)

mom[1] + 1.96* mom[3]


source(zcurve)
Directional = TRUE
boot.iter = 500
ymax = 1
sel = tm$z > 0
table(sel)
zres = Zing(
  val.input = tm$z[sel],
  yi        = tm$yi[sel],
  sei       = tm$sei[sel]
)


### get z-values from tmz file
stat_string_to_z <- function(x,
                             alternative = "two.sided",
                             signed = TRUE,
                             sign = 1,
                             compute_es = TRUE) {
  
  alternative <- match.arg(alternative, c("two.sided", "greater", "less"))
  
  x0 <- x
  x <- tolower(x)
  x <- gsub("\\s+", "", x)
  
  # normalize formats
  x <- gsub("^z\\([^)]*\\)=", "z=", x)
  x <- gsub("χ2", "chi2", x)
  x <- gsub("chisq", "chi2", x)
  x <- gsub("χ\\^2", "chi2", x)
  
  # defaults
  z <- p <- df_error <- yi <- sei <- NA_real_
  test <- NA_character_
  
  get_es_from_t <- function(tval, df) {
    d <- tval * sqrt(4 / (df + 2))
    vi <- 4 / (df + 2) + d^2 / (2 * df)
    c(yi = d, sei = sqrt(vi))
  }
  
  # z = value
  if (grepl("^z=", x)) {
    
    m <- regexec("^z=([-+]?[0-9.]+)$", x)
    out <- regmatches(x, m)[[1]]
    if (length(out) == 0) stop("Could not parse z string: ", x0)
    
    test <- "z"
    val <- as.numeric(out[2])
    
    if (alternative == "two.sided") {
      p <- 2 * pnorm(abs(val), lower.tail = FALSE)
      z <- abs(val)
      if (signed) z <- sign(val) * z
    } else if (alternative == "greater") {
      p <- pnorm(val, lower.tail = FALSE)
      z <- qnorm(1 - p)
    } else {
      p <- pnorm(val, lower.tail = TRUE)
      z <- qnorm(1 - p)
    }
  }
  
  # t(df) = value
  else if (grepl("^t\\(", x)) {
    
    m <- regexec("^t\\(([0-9.]+)\\)=([-+]?[0-9.]+)$", x)
    out <- regmatches(x, m)[[1]]
    if (length(out) == 0) stop("Could not parse t-test string: ", x0)
    
    test <- "t"
    df_error <- as.numeric(out[2])
    tval <- as.numeric(out[3])
    
    if (alternative == "two.sided") {
      p <- 2 * pt(abs(tval), df = df_error, lower.tail = FALSE)
      z <- qnorm(1 - p / 2)
      if (signed) z <- sign(tval) * z
    } else if (alternative == "greater") {
      p <- pt(tval, df = df_error, lower.tail = FALSE)
      z <- qnorm(1 - p)
    } else {
      p <- pt(tval, df = df_error, lower.tail = TRUE)
      z <- qnorm(1 - p)
    }
    
    if (compute_es) {
      es <- get_es_from_t(tval, df_error)
      yi <- es["yi"]
      sei <- es["sei"]
    }
  }
  
  # F(df1, df2) = value
  else if (grepl("^f\\(", x)) {
    
    m <- regexec("^f\\(([0-9.]+),([0-9.]+)\\)=([0-9.]+)$", x)
    out <- regmatches(x, m)[[1]]
    if (length(out) == 0) stop("Could not parse F-test string: ", x0)
    
    test <- "F"
    df1 <- as.numeric(out[2])
    df_error <- as.numeric(out[3])
    fval <- as.numeric(out[4])
    
    p <- pf(fval, df1 = df1, df2 = df_error, lower.tail = FALSE)
    z <- qnorm(1 - p / 2)
    if (signed) z <- sign * z
    
    # Only F(1, df2) is equivalent to t^2
    if (compute_es && df1 == 1) {
      tval <- sqrt(fval) * sign
      es <- get_es_from_t(tval, df_error)
      yi <- es["yi"]
      sei <- es["sei"]
    }
  }
  
  # chi2(df) = value
  else if (grepl("^chi2\\(", x)) {
    
    m <- regexec("^chi2\\(([0-9.]+)\\)=([0-9.]+)$", x)
    out <- regmatches(x, m)[[1]]
    if (length(out) == 0) stop("Could not parse chi-square string: ", x0)
    
    test <- "chi2"
    df_error <- as.numeric(out[2])
    val <- as.numeric(out[3])
    
    p <- pchisq(val, df = df_error, lower.tail = FALSE)
    z <- qnorm(1 - p / 2)
    if (signed) z <- sign * z
  }
  
  # r(df) = value, where df = n - 2
  else if (grepl("^r\\(", x)) {
    
    m <- regexec("^r\\(([0-9.]+)\\)=([-+]?[0-9.]+)$", x)
    out <- regmatches(x, m)[[1]]
    if (length(out) == 0) stop("Could not parse correlation string: ", x0)
    
    test <- "r"
    df_error <- as.numeric(out[2])
    r <- as.numeric(out[3])
    
    tval <- r * sqrt(df_error / (1 - r^2))
    p <- 2 * pt(abs(tval), df = df_error, lower.tail = FALSE)
    z <- qnorm(1 - p / 2)
    if (signed) z <- sign(r) * z
  }
  
  # p = value
  else if (grepl("^p[<=>]", x)) {
    
    m <- regexec("^p[<=>]([0-9.]+)$", x)
    out <- regmatches(x, m)[[1]]
    if (length(out) == 0) stop("Could not parse p-value string: ", x0)
    
    test <- "p"
    p <- as.numeric(out[2])
    z <- qnorm(1 - p / 2)
    if (signed) z <- sign * z
  }
  
  else {
    stop("Unknown or unsupported test-statistic string: ", x0)
  }
  
  c(
    z = z,
    p = p,
    df = df_error,
    yi = yi,
    sei = sei
  )
}


tmz$Statistics.1[810:820]	

which(grepl("103",tmz$Statistics.1))

i = 816
x = tmz$Statistics.1[i]
x
stat_string_to_z(x)

tmz[i,]

dim(tmz)
tmz.add = data.frame(t(sapply(tmz$Statistics.1,function(x) stat_string_to_z(x))))
colnames(tmz.add)

tmz$z = tmz.add$z
tmz$yi = tmz.add$yi.yi
tmz$sei = tmz.add$sei.sei

source(zcurve)
Directional = TRUE
boot.iter = 500
ymax = 7
sel = tmz$z > 0
table(sel)
zres = Zing(
  val.input = tmz$z[sel],
  yi        = tmz$yi[sel],
  sei       = tmz$sei[sel]
)

zres$es_estimate


yi = tmz$yi[sel]
vi = tmz$sei[sel]^2
vw.pcurve = weightfunct(yi, vi, steps=c(.025) )
vw.pcurve



sel = tmz$sei < .3 & tmz$yi < .8 & !is.na(tmz$yi) & !is.na(tmz$sei)
table(sel)

tmz$yi
tmz$sei

dat = tmz[sel,]

PET_RE <- rma(
  yi = yi,
  sei = sei,
  mods = ~ sei,
  method = "REML",
  data = dat
)

summary(PET_RE)

sqrt(PET_RE$tau2)


