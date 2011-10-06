##### Analysis of host preference experiment #####
## put 2 branches from a high density plant 
## and 2 branches from a low density plant
## into a cage with 4 females
## recorded locations and ovipositions at intervals
# 5 October 2011

library(bbmle)

d = read.csv('~/Documents/DATA/2011 DATA/choice experiments/pref_expt_2011.csv')

d$diff = d$galls10 - d$galls4

dr = data.frame(p10 = rep(NA,15), p4 = NA, diff=NA)

dr$p10 = as.vector(by(d$p10, d$rep, sum))
dr$p4 = as.vector(by(d$p4, d$rep, sum))
dr$diff = as.vector(by(d$diff, d$rep, mean))
dr= dr[-1,]
dr$p = dr$p10 / (dr$p10 + dr$p4)



m1 = mle2(dr$p10 ~ dbinom(size=dr$p10 + dr$p4, prob=1/(1+exp(a))), 
	start=list(a=1), data=dr)
m2 = mle2(dr$p10 ~ dbinom(size=dr$p10 + dr$p4, prob=1/(1+exp(a + b * dr$diff))), 
	start=list(a=1, b=0), data=dr)
AICtab(m1,m2, weights=TRUE)
plot(p ~ diff, data=dr, ylim=c(0,1), xlim=c(0,30))
abline(h=0.5)

mx1 = mle2(dr$p ~ dnorm(mean=prob, sd=sd), start=list(prob=0.5, sd=1),
	data=dr)
mx2 = mle2(dr$p ~ dnorm(mean=a + b * dr$diff, sd=sd), start=list(a=0.5, b=0, sd=1),
	data=dr)
abline(mx2)
AICtab(mx1, mx2)
mone = mle2(87 ~ dbinom(size=137, prob=prob), start=list(prob=0.5), data=dr)
summary(mone)