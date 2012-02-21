##### Analysis of host preference experiment #####
## put 2 branches from a high density plant 
## and 2 branches from a low density plant
## into a cage with 4 females
## recorded locations and ovipositions at intervals
# 5 October 2011

library(bbmle)
library(ggplot2)

d = read.csv('~/Documents/DATA/2011 DATA/choice experiments/pref_expt_2011.csv')

## data business: ##
# add a variable for the difference in gall numbers on branches
d$diff = d$galls10 - d$galls4 

# make a dataframe to hold data summarized by replicate
dr = data.frame(p10 = rep(NA,15), p4 = NA, diff=NA)

dr$p10 = as.vector(by(d$p10, d$rep, sum))
dr$p4 = as.vector(by(d$p4, d$rep, sum))
dr$diff = as.vector(by(d$diff, d$rep, mean))
dr= dr[-1,] # didn't keep track of number of galls on 10+ and 4- plants
# PROPORTION of observations of flies on branches that they were on 10+ branches:
dr$p = dr$p10 / (dr$p10 + dr$p4)  

### I SHould do a hierarchical model where observations within a rep are nested
## models: ##
# binomial with number of observations of flies on 10+ branches as outcome
m1 = mle2(dr$p10 ~ dbinom(size=dr$p10 + dr$p4, prob=1/(1+exp(a))), 
	start=list(a=1), data=dr)
m2 = mle2(dr$p10 ~ dbinom(size=dr$p10 + dr$p4, prob=1/(1+exp(a + b * dr$diff))), 
	start=list(a=1, b=0), data=dr)
AICtab(m1,m2, weights=TRUE)
plot(p ~ diff, data=dr, ylim=c(0,1), xlim=c(0,30))
abline(h=0.5)

# normal models with p per replicate
mx1 = mle2(dr$p ~ dnorm(mean=prob, sd=sd), start=list(prob=0.5, sd=1),
	data=dr)
mx2 = mle2(dr$p ~ dnorm(mean=a + b * dr$diff, sd=sd), start=list(a=0.5, b=0, sd=1),
	data=dr)
mx3 = mle2(dr$p ~ dnorm(mean=a * dr$diff / (b + dr$diff ), sd=sd), start=list(a=0.5, b=1, sd=1),
	data=dr)

abline(mx2)
AICtab(mx1, mx2, mx3)
anova(mx3, mx1)

pp = predict(mx3)

p1 = ggplot(data=dr, aes(x = diff, y = p)) + geom_point(alpha=1, position='jitter') +
	#stat_smooth(method='loess', fill='NA') +
	#stat_smooth(method="loess",fill=NA,colour="black",linetype=2,geom="ribbon") +
	scale_x_continuous('Difference in gall abundance') + 
	scale_y_continuous('Proportion on high density branches') +
	theme_bw() +
	opts( panel.grid.minor=theme_blank(), panel.grid.major=theme_blank(),
	axis.title.x = theme_text(vjust = 0)) +
	geom_hline(aes(yintercept=c(0.5)), linetype=2, lwd=0.5, alpha=0.5) +
	geom_smooth(aes(x=dr$diff, y=pp), data=dr, stat='identity')
	#geom_abline(intercept=0.1959, slope=0.0287139)
print(p1)
ggsave('~/Documents/Analysis-repos/pref-expt/figs/p-diff.pdf', width=3.5, height=3)


hist(dr$p)
curve(dbeta(x, shape1=86, shape2=49), 0, 1, add=TRUE)
mbeta1 = mle2(I(dr$p-0.01) ~ dbeta(shape1 = a, shape2 = b), start=list(a=2, b=2),
	data=dr)
curve(dbeta(x, shape1=coef(mbeta1)[1], shape2=coef(mbeta1)[2]), 0, 1, add=TRUE)

logistic = function(b1, b2){
	1 / (1 + exp(b1 + b2 * dr$diff))
}

mbeta2 = mle2(I(dr$p-0.01) ~ dbeta(shape1 = a, 
	shape2 = (a * logistic(b1, b2)) / logistic(b1, b2)), 
	start=list(a=0.6, b1=-0.8, b2=3.6), data=dr, method='SANN')
summary(mbeta2)


## proper beta regression with parm for mean and precision ##
## 26 Jan 2012 ##

hist(dr$p, xlim=c(0,1))
lines(density(dr$p))

logistic = function(x){ 1 / (1 + exp(-x))}

dr$p.corr = dr$p
dr$p.corr[dr$p.corr==1] = 0.999
#dr$p.corr = dr$p - 0.001

mu = 0.9
theta = 1
curve(dbeta(x, shape1=mu * theta, shape2=(1 - logistic(mu)) * theta, log=TRUE), 0,1)
dbeta(c(0,0.0001,0.01, 0.9999, 0.99, 1), shape1=mu * theta, shape2=(1 - logistic(mu)) * theta)
dbeta(dr$p, shape1=mu * theta, shape2=(1 - logistic(mu)) * theta, log=TRUE)
dbeta(dr$p.corr, shape1=mu * theta, shape2=(1 - logistic(mu)) * theta, log=TRUE)

# models

m0 = mle2(p.corr ~ dbeta(shape1 = logistic(mu) * theta, 
	shape2 = (1 - logistic(mu)) * theta), 
	start=list(mu=0, theta=0.001),
	data=dr, method='Nelder-Mead')

m1 = mle2(p.corr ~ dbeta(shape1 = logistic(b0 + b1 * diff) * theta, 
	shape2 = (1 - logistic(b0 + b1 * diff)) * theta), 
	start=list(b0 = 0.854, b1=0.1, theta=1.468),
	data=dr, method='Nelder-Mead')

m1.notlogistic = mle2(p.corr ~ dbeta(shape1 = (b0 + b1 * diff) * theta, 
	shape2 = (1 - (b0 + b1 * diff)) * theta), 
	start=list(b0 = 0.5, b1=0.01, theta=1),
	data=dr, method='Nelder-Mead')

m2 = mle2(p.corr ~ dbeta(shape1 = logistic(b0 + b1 * log(diff)) * theta, 
	shape2 = (1 - logistic(b0 + b1 * log(diff))) * theta), 
	start=list(b0 = 0, b1=0.01, theta=10),
	data=dr, method='Nelder-Mead')

AICtab(m0,m1,m2,mx3)
anova(m0,m1)

plot(p.corr ~ diff, data=dr)
newdiff = min(dr$diff):max(dr$diff)
m1predicted = logistic(coef(m1)['b0'] + coef(m1)['b1'] * newdiff)
points(m1predicted ~ newdiff, type='l')

m1.notlogistic.predicted = coef(m1.notlogistic)['b0'] + coef(m1.notlogistic)['b1'] * newdiff
points(m1.notlogistic.predicted ~ newdiff, type='l', lty=2)

plot(p.corr ~ log(diff), data=dr)
m2predicted = logistic(coef(m2)['b0'] + coef(m2)['b1'] * log(newdiff))
points(m2predicted ~ log(newdiff), type='l')



m1.profile = profile(m1)
m2.profile = profile(m2)

m1.ci = confint(m1.profile)
m2.ci = confint(m2.profile)

logistic(m1.ci)
logistic(m2.ci)