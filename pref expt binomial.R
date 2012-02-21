##### Analysis of host preference experiment #####
## put 2 branches from a high density plant 
## and 2 branches from a low density plant
## into a cage with 4 females
## recorded locations and ovipositions at intervals
# 5 October 2011

library(bbmle)
library(ggplot2)
library(lme4)

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

m0 = lmer(obs ~ (1|rep) + 1, data=d, family=binomial)
m1 = lmer(obs ~ (1|rep) + diff, data=d, family=binomial)

AIC(m0, m1)
anova(m0, m1)

coef(m1)
mean(coef(m1)$rep[,1])


newdiff = 9:30

ppred = logistic(mean(coef(m1)$rep[,1]) + coef(m1)$rep[1,2] * newdiff)

ggplot(data=d, aes(x=diff, y=prop, color=factor(rep))) +
	geom_point(alpha=1/2, position = position_jitter(w = 0.2, h = 0.02))+
	#stat_smooth() +
	geom_smooth(aes(x=newdiff, y=ppred), 
	data=d, stat='identity', colour='royalblue')




