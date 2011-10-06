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


## models: ##
# binomial with number of observations of flies on 10+ branches as outcome
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
mx2log = mle2(dr$p ~ dnorm(mean=a + b * log(dr$diff), sd=sd), start=list(a=0.5, b=0, sd=1),
	data=dr)
abline(mx2)
AICtab(mx1, mx2, mx2log)
mone = mle2(87 ~ dbinom(size=137, prob=prob), start=list(prob=0.5), data=dr)
summary(mone)


p1 = ggplot(data=dr, aes(x = diff, y = p)) + geom_point(alpha=1, position='jitter') +
	stat_smooth(method='loess') +
	scale_x_continuous('Difference in gall abundance') + 
	scale_y_continuous('Proportion on high density branches') +
	theme_bw() +
	opts( panel.grid.minor=theme_blank(), panel.grid.major=theme_blank(),
	axis.title.x = theme_text(vjust = 0)) +
	geom_hline(aes(yintercept=c(0.5)), linetype=2, lwd=0.5, alpha=0.5)
print(p1)
ggsave('~/Documents/Analysis repos/pref-expt/figs/p~diff.pdf', width=3.5, height=3)


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

	
	