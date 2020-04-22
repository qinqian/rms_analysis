library(lmer)
data(Oats)
str(Oats)

model=lm(yield~Variety*nitro, data=Oats)
summary(model)

model2 = lme(yield~Variety*nitro, data=Oats, random=~1|Block/Variety)
summary(model2)
