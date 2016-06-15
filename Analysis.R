#########################################################################
### Factors determining distribution and abundance of Striped Lychnis ###
#########################################################################

dframe1 <- read.csv("Data/Sites.csv")


library(lme4)
library(arm)
library(MuMIn)

dframe1$Site <- factor(dframe1$Site)
dframe1$WeevilDensity <- ordered(dframe1$WeevilDensity, levels = c("N","R","O","F","A","D"))
dframe1$PlantDensity <- ordered(dframe1$PlantDensity, levels = c("R","O","F","A","D"))
dframe1$WildflowerVerge <- factor(dframe1$WildflowerVerge)

summary(dframe1)

### Larva presence/absence

plot(LarvaePresent ~ PlantCount, data = dframe1)
plot(LarvaePresent ~ NNSmoothed, data = dframe1)
plot(LarvaePresent ~ WeevilDensity, data = dframe1)

model1 <- glm(LarvaePresent ~ 
                PlantCount + MeanSward + 
                WeevilDensity + NNSmoothed,
                family = binomial (link = "logit"),
                data = dframe1)

summary(model1)
drop1(model1, test="Chisq")


### Try an information theoretic approach ###

global.model <- glm(LarvaePresent ~
                    PlantCount + MeanSward +
                    WildflowerVerge + WindShelter +
                    WeevilDensity + NNSmoothed,
                    family = binomial(link = "logit"),
                    na.action = na.fail,
                    data = dframe1)

summary(global.model)
drop1(global.model, test = "Chi")

stdz.model <- standardize(global.model,
                          standardize.y = FALSE)
summary(stdz.model)

model.set <- dredge(stdz.model)

top.models <- get.models(model.set, 
                         subset = delta <2)
top.models

averaged.model <- model.avg(top.models)
summary(averaged.model)



### Larval abundance

plot(MaxLarvae ~ PlantCount, data = dframe1)
plot(MaxLarvae ~ NNSmoothed, data = dframe1)
plot(MaxLarvae ~ WeevilDensity, data = dframe1)

mean(dframe1$MaxLarvae)
var(dframe1$MaxLarvae)


model2 <- glm(MaxLarvae ~ 
              PlantCount + MeanSward + 
              WeevilDensity + NNSmoothed,
              family = poisson (link = "log"),
              data = dframe1)

summary(model2)
drop1(model2, test="Chisq")


### Larval abundance - present only

dframe1P <- subset(dframe1,LarvaePresent==1)

plot(MaxLarvae ~ PlantCount, data = dframe1P)
plot(MaxLarvae ~ NNSmoothed, data = dframe1P)
plot(MaxLarvae ~ WeevilDensity, data = dframe1P)

mean(dframe1P$MaxLarvae)
var(dframe1P$MaxLarvae)


model2P <- glm(MaxLarvae ~ 
               PlantCount + MeanSward + 
               WeevilDensity + NNSmoothed,
               family = poisson (link = "log"),
               data = dframe1P)

summary(model2P)
drop1(model2P, test="Chisq")


