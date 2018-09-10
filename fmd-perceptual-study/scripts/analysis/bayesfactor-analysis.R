# Bayesian analysis of contrast data using BayesFactor
# Uses expsum data <- c( subjID, group, attention, pCON)

library(BayesFactor)

# Try traditional NHST mixed ANOVA

##
summary(aov(pCON ~ attention*group + Error(subjID / attention), data=expsum))

summary(aov(pCON ~ attention * group, data = expsum))
bf <- anovaBF(pCON ~ attention * group, data = expsum, whichModels = 'withmain')

plot(bf[1:3] / bf[4])

bfMainEffects = lmBF(pCON ~ attention + group, data = expsum)
bfInteraction = lmBF(pCON ~ attention + group + attention*group, data = expsum)
# Comparison to examine interaction only
bf = bfInteraction / bfMainEffects

# Recompute to reduce error (here from 1.9% to ~0.9%), default iterations = 10,000
newbf <- recompute(bf, iterations = 50000)

chains = posterior(bfInteraction, iterations = 10000)
summary(chains[ 1:2])
plot(chains [4:6])

##

summary(aov(Dprime ~ Attention*Expectations + Error(Group+Subject / (Attention*Expectations)), 
                 data=longdata))
anovaBF(Dprime ~ Attention*Expectations + Group + Subject, data = longdata, whichRandom = 'Subject')

bfBest <- lmBF(Dprime ~ Attention + Expectations + Group + Attention*Group + Subject
               , whichRandom = 'Subject'
               , data = longdata)

anovaBF(pCON ~ attention * group + subjID, data = expsum, whichRandom = c('group','subjID'), progress = T)

bfBest <- lmBF(pCON ~ attention + group + attention*group + subjID
               , whichRandom = c('subjID','group')
               , data = expsum
               , iterations = 10000)
bfSubj <- lmBF(pCON ~ subjID
               , whichRandom = 'subjID'
               , data = expsum)

bf = bfBest / bfSubj

## 
library(ez)

fullcon <- ezANOVA(expsum
                   , dv = pCON
                   , wid = subjID
                   , within = attention
                   , between = group
                   , type = 3)

# follow up
con_des <- ezStats(expsum
                   , dv = pCON
                   , wid = subjID
                   , within = attention
                   , between = group
)
 
print(con_des)

ezPlot(expsum
       , dv = .(pCON)
       , wid = .(subjID)
       , within = .(attention)
       , between = .(group)
       , x = .(attention)
       , do_lines = TRUE
       , levels = list(
         attention = list(
           new_order = c('full','diverted')
         )
       )
       , col = .(group)
       , x_lab = 'Attention'
       , y_lab = 'Gabor Contrast Threshold'
)

cont <- subset(expsum, group=='control')
fmd <- subset(expsum, group=='fmd')
omd <- subset(expsum, group=='organic')

with(cont, pairwise.t.test(pCON, attention,
                           paired = T, p.adj = 'holm'))
with(fmd, pairwise.t.test(pCON, attention,
                          paired = T, p.adj = 'holm'))
with(omd, pairwise.t.test(pCON, attention,
                          paired = T, p.adj = 'holm'))

fullatt <- subset(expsum, attention=='full')
divatt <- subset(expsum, attention=='diverted')

with(fullatt, pairwise.t.test(pCON, group,
                              paired = F,
                              p.adjust.method = 'holm'))
with(divatt, pairwise.t.test(pCON, group,
                             paired = F,
                             p.adjust.method = 'holm'))

diffScores <- fmd[fmd$attention=='full','pCON'] - fmd[fmd$attention=='diverted','pCON']
t.test(diffScores)
ttestBF(x = diffScores)

ttestBF(formula = pCON ~ group, data = expsum, paired = F)

condiff <- expsum[expsum$attention=='full','pCON'] - expsum[expsum$attention=='diverted','pCON']

condiff <- read.csv('summarise-expdat.csv')

with(condiff,pairwise.t.test(diffCON, group,
                paired = F,
                p.adj = 'holm'))

anovaBF(diffCON ~ group, data = condiff)

bfMainEffects <- lmBF(pCON ~ attention + group
                      , data = expsum)
bfAttention <- lmBF(pCON ~ attention, data = expsum)

bf = bfMainEffects / bfAttention

chains <- posterior(bfMainEffects, iterations = 15000)
summary(chains)
