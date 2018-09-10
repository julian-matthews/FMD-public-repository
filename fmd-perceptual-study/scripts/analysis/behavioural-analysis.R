## Behavioural analysis rm-anova

## Mixed ANOVA

alldata <- JASP.export

objper <- alldata[ , c(1,2,9:14)]

mahal <- mahalanobis(objper[ , 3:8],
                     colMeans(objper[ , 3:8]),
                     cov(objper[ , 3:8]))
cutoff <- qchisq(1-.001, ncol(objper[ , 3:8]))

mahal < cutoff

noout <- subset(objper, mahal < cutoff)

## assumptions
random <- rchisq(nrow(noout), 7)
fake <- lm(random ~ ., data = noout[ , -1])
standardized <- rstudent(fake)
fitted <- scale(fake$fitted.values)

##normality
hist(standardized)

##linearity
qqnorm(standardized)
abline(0,1)

##homog
plot(fitted, standardized)
abline(0,0)
abline(v=0)

library(reshape)
longdata = melt(noout,
                id = c('Subject','Group'),
                variable_name = 'Attention')
longdata$Expectations <- longdata$Attention
longdata$Subject <- as.factor(longdata$Subject)

levels(longdata$Attention) <- list(Full=c("X1D.E1AF","X1D.E2AF","X1D.E3AF"), 
                                   Diverted=c("X1D.E1AD","X1D.E2AD","X1D.E3AD"))
levels(longdata$Expectations) <- list(Low=c("X1D.E1AF","X1D.E1AD"), 
                                   Medium=c("X1D.E2AF","X1D.E2AD"),
                                   High=c("X1D.E3AF","X1D.E3AD"))

longdata <- rename(longdata, c(value = 'Dprime'))

library(ez)
fullnov <- ezANOVA(data = longdata,
        dv = Dprime,
        within = .(Attention,Expectations),
        between = Group,
        wid = Subject,
        type = 3)

# Expectations follow-up
exp_des <- ezStats(
  data = longdata,
  dv = Dprime,
  wid = Subject,
  within = Expectations,
  within_full = .(Attention,Expectations)
)
print(exp_des)

ezStats(
  data = longdata,
  dv = Dprime,
  wid = Subject,
  within = Attention,
  within_full = .(Attention,Expectations),
  between = Group
)

ezPlot(
  data = longdata,
  dv = .(Dprime),
  wid = .(Subject),
  within = .(Expectations),
  within_full = .(Attention,Expectations),
  between_full = .(Group),
  x = .(Expectations),
  do_lines = TRUE,
  x_lab = 'Expectations',
  y_lab = 'Type I d'
)

with(longdata, pairwise.t.test(Dprime, Expectations, p.adj = 'holm',paired = T))
with(longdata, pairwise.t.test(Dprime, Attention, p.adj = 'holm',paired = T))

##levenes
ezANOVA(data = longdata,
        dv = Dprime,
        between = Group,
        wid = Subject,
        type = 3)

options(scipen = 999)

# Only significant is Expectations F(2,114)=11.30, p<.001, geta2=.02
# 

with(longdata, tapply(Dprime, list(Expectations,Group), mean))
with(longdata, tapply(pCON, list(Group), mean))
with(longdata, tapply(pCON, list(attention,group), mean))
with(longdata, tapply(pCON, list(attention,group), sd))

cont <- subset(longdata, Group=='Control')
fmd <- subset(longdata, Group=='FMD')
omd <- subset(longdata, Group=='OMD')

fullatt <- subset(longdata, Attention=='Full')
divatt <- subset(longdata, Attention=='Diverted')

with(cont, pairwise.t.test(Dprime, Attention,
                               paired = T, p.adj = 'holm'))
with(fmd, pairwise.t.test(Dprime, Attention,
                           paired = T, p.adj = 'holm'))
with(omd, pairwise.t.test(Dprime, Attention,
                           paired = T, p.adj = 'holm'))

with(fullatt, pairwise.t.test(Dprime, Group,
                                   paired = F,
                                   p.adjust.method = 'holm'))
with(divatt, pairwise.t.test(Dprime, Group,
                                  paired = F,
                                  p.adjust.method = 'holm'))

## Posthoc comparisons with multcomp

library(nlme)
library(multcomp)
lme_dprime <- lme(Dprime ~ Attention*Expectations, data=longdata, random = ~1 | Group)
anova(lme_dprime)
summary(glht(lme_dprime, linfct=mcp(Expectations = "Tukey")), test = adjusted(type = "bonferroni"))


