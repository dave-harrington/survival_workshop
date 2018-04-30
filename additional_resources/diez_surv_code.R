rm(list=ls())

#===> loading packages and such <===#
#install.packages("OIsurv")
library(OIsurv)
data(aids)
aids
attach(aids)
infect
detach(aids)
aids$infect

#===> survival object <===#
data(tongue); attach(tongue)   # the following will not affect computations
# create a subset for just the first group by using [type==1]
my.surv.object <- Surv(time[type==1], delta[type==1])
my.surv.object
detach(tongue)

data(psych); attach(psych)
my.surv.object <- Surv(age, age+time, death)
my.surv.object
detach(psych)

#===> K-M Estimate <===#
data(tongue); attach(tongue)
my.surv <- Surv(time[type==1], delta[type==1])
survfit(my.surv ~ 1)
my.fit <- survfit(my.surv ~ 1)
summary(my.fit)$surv     # returns the Kaplan-Meier estimate at each t_i
summary(my.fit)$time     # {t_i}
summary(my.fit)$n.risk   # {Y_i}
summary(my.fit)$n.event  # {d_i}
summary(my.fit)$std.err  # standard error of the K-M estimate at {t_i}
summary(my.fit)$lower    # lower pointwise estimates (alternatively, $upper)
str(my.fit)              # full summary of the my.fit object
str(summary(my.fit))     # full summary of the my.fit object
pdf("../figures/kmPlot.pdf", 7, 4.5)
plot(my.fit, main="Kaplan-Meier estimate with 95% confidence bounds",
   xlab="time", ylab="survival function")
dev.off()
my.fit1 <- survfit( Surv(time, delta) ~ type )   # here the key is "type"
detach(tongue)

#===> confidence bands <===#
data(bmt); attach(bmt)
my.surv <- Surv(t2[group==1], d3[group==1])
my.cb <- confBands(my.surv, confLevel=0.95, type="hall")
pdf("../figures/confBand.pdf", 8, 5)
plot(survfit(my.surv ~ 1), xlim=c(100, 600), xlab="time",
  ylab="Estimated Survival Function",
  main="Reproducing Confidence Bands for Example 4.2 in Klein/Moeschberger")

lines(my.cb$time, my.cb$lower, lty=3, type="s")
lines(my.cb$time, my.cb$upper, lty=3, type="s")
legend(100, 0.3, legend=c("K-M survival estimate",
  "pointwise intervals","confidence bands"), lty=1:3)
dev.off()
detach(bmt)


#===> cumulative hazard <===#
data(tongue); attach(tongue)
my.surv <- Surv(time[type==1], delta[type==1])
my.fit  <- summary(survfit(my.surv ~ 1))
H.hat   <- -log(my.fit$surv)
H.hat   <- c(H.hat, tail(H.hat, 1))
h.sort.of <- my.fit$n.event / my.fit$n.risk
H.tilde   <- cumsum(h.sort.of)
H.tilde   <- c(H.tilde, tail(H.tilde, 1))
pdf("../figures/cumHazard.pdf", 6, 4)
plot(c(my.fit$time, 250), H.hat, xlab="time", ylab="cumulative hazard",
  main="comparing cumulative hazards", ylim=range(c(H.hat, H.tilde)), type="s")
points(c(my.fit$time, 250), H.tilde, lty=2, type="s")
legend("topleft", legend=c("H.hat","H.tilde"), lty=1:2)
dev.off()
detach(tongue)


#===> mean/median <===#
data(drug6mp); attach(drug6mp)
my.surv <- Surv(t1, rep(1, 21))   # all placebo patients observed
survfit(my.surv ~ 1)
print(survfit(my.surv ~ 1), print.rmean=TRUE)
detach(drug6mp)

#===> test for 2+ samples <===#
data(btrial); attach(btrial)
survdiff(Surv(time, death) ~ im)   # output omitted
survdiff(Surv(time, death) ~ im, rho=1)   # output omitted
detach(btrial)

#===> coxph, time-independent <===#
data(burn); attach(burn)
my.surv   <- Surv(T1, D1)
coxph.fit <- coxph(my.surv ~ Z1 + as.factor(Z11), method="breslow")
coxph.fit
co <- coxph.fit$coefficients  # may use coxph.fit$coeff instead
va <- coxph.fit$var           # I^(-1), estimated cov matrix of the estimates
ll <- coxph.fit$loglik        # log-likelihood for alt and null MLEs, resp.
my.survfit.object <- survfit(coxph.fit)
hold <- survfit(my.surv ~ 1)
#source("http://www.stat.ucla.edu/~david/teac/surv/local-coxph-test.R")
coxph.fit
C   <- matrix(c(0, 1, 0, 0,
                0, 0, 1, 0,
                0, 0, 0, 1), nrow=3, byrow=TRUE)
d   <- rep(0, 3)
t1  <- C %*% co - d
t2  <- C %*% va %*% t(C)
XW2 <- c(t(t1) %*% solve(t2) %*% t1)
pchisq(XW2, 3, lower.tail=FALSE)
#local.coxph.test(coxph.fit, 2:4)
my.survfit.object <- survfit(coxph.fit)
detach(burn)


#===> coxph, time-dependent <===#
data(relapse)
relapse

N  <- dim(relapse)[1]
t1 <- rep(0, N+sum(!is.na(relapse$int)))  # initialize start time at 0
t2 <- rep(-1, length(t1))                 # build vector for end times
d  <- rep(-1, length(t1))                 # whether event was censored
g  <- rep(-1, length(t1))                 # gender covariate
i  <- rep(FALSE, length(t1))              # initialize intervention at FALSE

j  <- 1
for(ii in 1:dim(relapse)[1]){
  if(is.na(relapse$int[ii])){      # no intervention, copy survival record
    t2[j] <- relapse$event[ii]
    d[j]  <- relapse$delta[ii]
    g[j]  <- relapse$gender[ii]
    j <- j+1
  } else {                         # intervention, split records
    g[j+0:1] <- relapse$gender[ii] # gender is same for each time
    d[j]     <- 0                  # no relapse observed pre-intervention
    d[j+1]   <- relapse$delta[ii]  # relapse occur post-intervention?
    i[j+1]   <- TRUE               # intervention covariate, post-intervention
    t2[j]    <- relapse$int[ii]-1  # end of pre-intervention
    t1[j+1]  <- relapse$int[ii]-1  # start of post-intervention
    t2[j+1]  <- relapse$event[ii]  # end of post-intervention
    j <- j+2                       # two records added
  }
}

mySurv <- Surv(t1, t2, d)        # pg 3 discusses left-trunc. right-cens. data
myCPH  <- coxph(mySurv ~ g + i)

#data(burn); attach(burn)
##source("http://www.stat.ucla.edu/~david/teac/surv/time-dep-coxph.R")
#td.coxph <- timeDepCoxph(burn, "T1", "D1", 2:4, "Z1", verbose=FALSE)
#td.coxph   # some model output is omitted for brevity
#detach(burn)

#===> AFT models <===#
data(larynx)
attach(larynx)
srFit <- survreg(Surv(time, delta) ~ as.factor(stage) + age, dist="weibull")
summary(srFit)
srFitExp <- survreg(Surv(time, delta) ~ as.factor(stage) + age, dist="exponential")
summary(srFitExp)
srFitExp$coeff    # covariate coefficients
srFitExp$icoef    # intercept and scale coefficients
srFitExp$var      # variance-covariance matrix
srFitExp$loglik   # log-likelihood
srFit$scale       # not using srFitExp (defaulted to 1)
detach(larynx)
