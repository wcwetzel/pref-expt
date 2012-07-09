##### Analysis of host preference experiment #####
## put 2 branches from a high density plant 
## and 2 branches from a low density plant
## into a cage with 4 females
## recorded locations and ovipositions at intervals
# 5 October 2011

library(bbmle)
library(ggplot2)
library(lme4)
library(rethinking)

d = read.csv('~/Documents/DATA/2011 DATA/choice experiments/pref_expt_2011.csv')

logistic = function(x) 1/(1 + exp(-x))


## data business: ##
# add a variable for the difference in gall numbers on branches
d$diff = d$galls10 - d$galls4 
d$prop = d$p10 / (d$p10 + d$p4)
with(d, mean(prop[!is.nan(prop)]))
d$ttime = strptime(d$time, format = '%H:%M')

d = d[!is.nan(d$prop),]

with(d, pairs(cbind(rep, days, time, prop, fliesfrom)))

plot(prop ~ g10, data=d)

plot(d$ttime, d$diff, xaxt='n', cex=2*d$prop)
axis.POSIXct(side=1, x=d$ttime, format="%H:%M")

ts(prop ~ ttime, data=d)


ggplot(data=d, aes(x=days, y=prop)) +
	geom_point()+
	stat_smooth()

ggplot(data=d, aes(x=days, y=prop)) +
	geom_point()+
	stat_smooth()



ggplot(data=d, aes(x=diff, y=prop)) +
	geom_point()+
	stat_smooth()


obs = with(d, cbind(p10, p4))

with(d, p10 / (p10 + p4))

ggplot(data=d[d$rep>1,], aes(x=diff, y=prop, color=factor(rep))) +
	geom_point(alpha=1/2, position = position_jitter(w = 0.2, h = 0.02))+
	#stat_smooth() +
	geom_smooth(aes(x=newdiff, y=ppred), 
	data=d, stat='identity', colour='royalblue')





#-----------binomial mixed models of decisions as a function of diff----------#

m0 = lmer(obs[d$rep>1,] ~ (1|d$rep[d$rep>1]) + 1, data=d[d$rep>1,], family=binomial)
m1 = lmer(obs[d$rep>1,] ~ (1|d$rep[d$rep>1]) + diff, data=d[d$rep>1,], family=binomial)

AICtab(m0, m1, weights=TRUE)
AICctab(m0, m1, nobs=14, weights=TRUE)

anova(m0, m1)

### plotting

post.m1 = sample.naive.posterior(m1)

newdiff = 9:30

m1.mu = sapply(newdiff,
	function(z) mean(logistic(post.m1[,1] + post.m1[,2] * z)))

m1.ci = sapply(newdiff,
	function(z) HPDI(logistic(post.m1[,1] + post.m1[,2] * z)))

# I want to plot by replicate
# make a dataframe to hold data summarized by replicate
dr = data.frame(p10 = rep(NA,15), p4 = NA, diff=NA)
dr$p10 = as.vector(by(d$p10, d$rep, sum))
dr$p4 = as.vector(by(d$p4, d$rep, sum))
dr$diff = as.vector(by(d$diff, d$rep, mean))
dr= dr[-1,] # didn't keep track of number of galls on 10+ and 4- plants
# PROPORTION of observations of flies on branches that they were on 10+ branches:
dr$p = dr$p10 / (dr$p10 + dr$p4)  


ggplot(data=dr, aes(x=diff, y=p)) +
	scale_x_continuous('Difference in gall abundance') + 
	scale_y_continuous('Proportion on high abundance branches') +
	geom_point() +
	geom_line(aes(x=newdiff, y=m1.mu)) +
	geom_line(aes(x=newdiff, y=m1.ci[1,]), lty=2) +
	geom_line(aes(x=newdiff, y=m1.ci[2,]), lty=2) +
	theme_bw() +
	opts( panel.grid.minor=theme_blank(), panel.grid.major=theme_blank(),
	axis.title.x = theme_text(vjust = 0))





