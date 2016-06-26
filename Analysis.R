#########################################################################
### Factors determining distribution and abundance of Striped Lychnis ###
#########################################################################

# dframe2 <- read.csv(file.choose())

dframe1 <- read.csv("Data/Sites.csv")


library(lme4)
library(arm)
library(MuMIn)
library(bbmle)

source("CheckResidsFunction.R")

# source(file.choose())

dframe1$Site <- factor(dframe1$Site)
dframe1$WeevilDensity <- ordered(dframe1$WeevilDensity, levels = c("N","R","O","F","A","D"))
dframe1$PlantDensity <- ordered(dframe1$PlantDensity, levels = c("R","O","F","A","D"))
dframe1$WildflowerVerge <- factor(dframe1$WildflowerVerge)
dframe1$TetradIsolation <- ordered(dframe1$TetradIsolation, levels = c(0,1,2,3,4,5,6,7,8))

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

### test the different variables for...
# patch quality

model1a <- glm(LarvaePresent ~ 
                 PlantDensity + MeanSward + 
                 WeevilDensity + NNSmoothed,
               family = binomial (link = "logit"),
               data = dframe1)

summary(model1a)
drop1(model1a, test = "Chi")

AIC(model1,model1a) # model1 (PlantCount) is better

# patch isolation

model1b <- glm(LarvaePresent ~ 
                PlantCount + MeanSward + 
                WeevilDensity + NN,
              family = binomial (link = "logit"),
              data = dframe1)

summary(model1b)
drop1(model1b, test = "Chi")

model1c <- glm(LarvaePresent ~ 
                 PlantCount + MeanSward + 
                 WeevilDensity + TetradIsolation,
               family = binomial (link = "logit"),
               data = dframe1)

summary(model1c)
drop1(model1c, test = "Chi")

AIC(model1,model1b,model1c) # again model1 is best (NNSmoothed)
AICtab(model1,model1b,model1c)

summary(model1)
drop1(model1, test = "Chi")


### Try an information theoretic approach ###

global.model1 <- glm(LarvaePresent ~
                    PlantCount + MeanSward +
                    WildflowerVerge + WindShelter +
                    WeevilDensity + NNSmoothed,
                    family = binomial(link = "logit"),
                    na.action = na.fail,
                    data = dframe1)

summary(global.model1)
drop1(global.model1, test = "Chi")

stdz.model1 <- standardize(global.model1,
                          standardize.y = FALSE)
summary(stdz.model1)

model.set1 <- dredge(stdz.model1)

top.models1 <- get.models(model.set1, 
                         subset = delta <2)
top.models1

averaged.model1 <- model.avg(top.models1,rank="AICc")
summary(averaged.model1)

confint(averaged.model1,
        level=0.95)

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

chkres(model2) # pretty bad - but there are loads of zeroes skewing these data

### Larval abundance - present only

dframe1P <- subset(dframe1,LarvaePresent==1)

plot(MaxLarvae ~ PlantCount, data = dframe1P)
plot(MaxLarvae ~ PlantDensity, data = dframe1P)
plot(MaxLarvae ~ NNSmoothed, data = dframe1P)
plot(MaxLarvae ~ WeevilDensity, data = dframe1P)
plot(MaxLarvae ~ MeanSward, data = dframe1P)

mean(dframe1P$MaxLarvae)
var(dframe1P$MaxLarvae)

model2P <- glm(MaxLarvae ~
                 PlantCount + MeanSward +
                 WeevilDensity + NNSmoothed,
               family = poisson (link = "log"),
               data = dframe1P)

summary(model2P)
drop1(model2P, test = "Chi")

### test the different variables for...
# patch quality

model2a <- glm(MaxLarvae ~ 
                 PlantDensity + MeanSward + 
                 WeevilDensity + NNSmoothed,
               family = poisson (link = "log"),
               data = dframe1P)

summary(model2a)
drop1(model2a, test = "Chi")

AIC(model2P,model2a) # model2a (PlantDensity) is better

# patch isolation

model2b <- glm(MaxLarvae ~ 
                 PlantDensity + MeanSward + 
                 WeevilDensity + NN,
               family = poisson (link = "log"),
               data = dframe1P)

summary(model2b)
drop1(model2b, test = "Chi")

model2c <- glm(MaxLarvae ~ 
                 PlantDensity + MeanSward + 
                 WeevilDensity + TetradIsolation,
               family = poisson (link = "log"),
               data = dframe1P)

summary(model2c)
drop1(model2c, test = "Chi")

AIC(model2a,model2b,model2c) # again model2a is best (NNSmoothed)


summary(model2a)
drop1(model2a, test = "Chi")

# information theoretic approach

global.model2 <- glm(MaxLarvae ~ 
                     PlantCount + MeanSward + 
                     WeevilDensity + NNSmoothed,
                     family = poisson (link = "log"),
                     na.action = na.fail,
                     data = dframe1P)

summary(global.model2)
drop1(global.model2, test="Chisq")

chkres(global.model2)


stdz.model2 <- standardize(global.model2,
                          standardize.y = FALSE)
summary(stdz.model2)

model.set2 <- dredge(stdz.model2)

top.models2 <- get.models(model.set2, 
                         subset = delta <2)
top.models2

averaged.model2 <- model.avg(top.models2,rank="AICc")
summary(averaged.model2)

confint(averaged.model2,
        level=0.95)




#################################################################################
### Factors determining plant-level variation in occupancy by Striped Lychnis ###
#################################################################################

dframe2 <- read.csv("Data/Plants.csv")

dframe2$Site <- factor(dframe2$Site)

summary(dframe2)

### Larva presence/absence

plot(LarvaePresent ~ PlantHeight, data = dframe2)
plot(LarvaePresent ~ Flowerspikes, data = dframe2)

model3 <- glmer(LarvaePresent ~ 
                PlantHeight + Flowerspikes
                + (1|Site),
              family = binomial (link = "logit"),
              na.action=na.fail,
              data = dframe2)

summary(model3)
drop1(model3, test="Chisq")


### Try an information theoretic approach ###

global.model3 <- model3
  
stdz.model3 <- standardize(global.model3,
                           standardize.y = FALSE)
summary(stdz.model3)

model.set3 <- dredge(stdz.model3)

top.models3 <- get.models(model.set3, 
                          subset = delta <2)
top.models3


# only global.model selected as a top model so just use that

confint(global.model3,
        level=0.95)




#######################################################

### figures

library(coefplot)
library(arm)
library(ggplot2)

coefplot(model1)
coefplot(model2a)


names(dframe1)

### occupancy

### some predicted datasets

# probability of patch occupancy
newdata<-expand.grid(NNSmoothed=(seq(0.5,12,0.1)),PlantCount=(seq(1,350,1)),
                     MeanSward=23.69,WeevilDensity=c("N","R","O","F","A","D"),LarvaePresent=1)
newdata$LarvaePresent <- predict(model1,newdata=newdata,type="response")
preddat <- predict(model1,newdata=newdata,type="response",se.fit=TRUE)
preddat

newdata <- cbind(newdata,preddat)
newdata <- data.frame(
  newdata
  , plo = newdata$fit-1.96*newdata$se.fit
  , phi = newdata$fit+1.96*newdata$se.fit
)


# site size at 50% occupancy
newdata0.5 <- subset(newdata, LarvaePresent < 0.51,
                   select=c(NNSmoothed,PlantCount,LarvaePresent))
newdata0.5 <- subset(newdata0.5, LarvaePresent > 0.49,
                   select=c(NNSmoothed,PlantCount,LarvaePresent))


# site size at 90% occupancy
newdata0.9 <- subset(newdata, LarvaePresent < 0.91,
                   select=c(NNSmoothed,PlantCount,LarvaePresent))
newdata0.9 <- subset(newdata0.9, LarvaePresent > 0.89,
                   select=c(NNSmoothed,PlantCount,LarvaePresent))


# site size at 10% occupancy
newdata0.1 <- subset(newdata, LarvaePresent < 0.11,
                     select=c(NNSmoothed,PlantCount,LarvaePresent))
newdata0.1 <- subset(newdata0.1, LarvaePresent > 0.09,
                     select=c(NNSmoothed,PlantCount,LarvaePresent))

# site size at 97.5% occupancy
newdata0.975 <- subset(newdata, LarvaePresent < 0.98,
                     select=c(NNSmoothed,PlantCount,LarvaePresent))
newdata0.975 <- subset(newdata0.975, LarvaePresent > 0.97,
                     select=c(NNSmoothed,PlantCount,LarvaePresent))

# site size at 2.5% occupancy
newdata0.025 <- subset(newdata, LarvaePresent < 0.03,
                     select=c(NNSmoothed,PlantCount,LarvaePresent))
newdata0.025 <- subset(newdata0.025, LarvaePresent > 0.02,
                     select=c(NNSmoothed,PlantCount,LarvaePresent))


g1 <- ggplot(dframe1,
             aes(x=NNSmoothed, y=PlantCount))+
             scale_x_log10(limits=c(0.5,12),
                           breaks=c(1, 2, 5, 10))+
             scale_y_log10(limits=c(0.5,350),
                           breaks=c(1,10,50,100,200,400))+
             geom_point(size=4,colour="black",fill="white",stroke=2,
                        aes(shape=factor(LarvaePresent)))+
                            scale_shape_manual(values=c(21,16),
                                               labels=c("Absent","Present"),
                                               name=expression(paste(italic("S. lychnitis "),"larvae")))+
             geom_smooth(colour="black",
                         data = newdata0.5, aes(x=NNSmoothed, y=PlantCount),
                         method=lm,se=FALSE,
                         fullrange=TRUE)+
             geom_smooth(linetype="dashed",
                         colour="black",
                         data = newdata0.975, aes(x=NNSmoothed, y=PlantCount),
                         method=lm,se=FALSE,
                         fullrange=TRUE)+
             geom_smooth(linetype="dashed",
                         colour="black",
                         data = newdata0.025, aes(x=NNSmoothed, y=PlantCount),
                         method=lm,se=FALSE,
                         fullrange=TRUE)+
             xlab("Isolation (smoothed nearest neighbour distance, km)")+
             ylab(expression(paste("No. ",
                                   italic("V. nigrum "), "flowerspikes")))+
             theme(panel.background=element_rect(fill="white"),
                   panel.grid.major.x=element_line(colour="gray90"),
                   panel.grid.major.y=element_line(colour="gray90"),
                   panel.grid.minor=element_blank(),
                   panel.border=element_rect(color="black",fill=F,size=1),
                   text=element_text(size=15),
                   axis.text=element_text(color="black"),
                   legend.title=element_text(),
                   legend.background=element_rect(fill="white", colour="black"),
                   legend.key = element_blank(),
                   legend.justification=c(1,0), legend.position=c(0.95,0.15))

g1


### abundance


### plant-level occupancy
names(dframe2)

summary(model3)

plot(PlantHeight ~ Flowerspikes, data = dframe2)


# probability of plant occupancy
newdata2<-expand.grid(Flowerspikes=(seq(1,12,1)),PlantHeight=(seq(32,198,1)),
                     Site=c(7,16,27,32,38,40,47), LarvaePresent=1)
newdata2$LarvaePresent <- predict(model3,newdata=newdata2,type="response")


# site size at 50% occupancy
newdata2_0.5 <- subset(newdata2, LarvaePresent < 0.51,
                       select=c(Flowerspikes,PlantHeight,LarvaePresent))
newdata2_0.5 <- subset(newdata2_0.5, LarvaePresent > 0.49,
                       select=c(Flowerspikes,PlantHeight,LarvaePresent))


# site size at 90% occupancy
newdata2_0.9 <- subset(newdata2, LarvaePresent < 0.91,
                       select=c(Flowerspikes,PlantHeight,LarvaePresent))
newdata2_0.9 <- subset(newdata2_0.9, LarvaePresent > 0.89,
                       select=c(Flowerspikes,PlantHeight,LarvaePresent))


# site size at 10% occupancy
newdata2_0.1 <- subset(newdata2, LarvaePresent < 0.11,
                       select=c(Flowerspikes,PlantHeight,LarvaePresent))
newdata2_0.1 <- subset(newdata2_0.1, LarvaePresent > 0.09,
                       select=c(Flowerspikes,PlantHeight,LarvaePresent))

# site size at 97.5% occupancy
newdata2_0.975 <- subset(newdata2, LarvaePresent < 0.98,
                         select=c(Flowerspikes,PlantHeight,LarvaePresent))
newdata2_0.975 <- subset(newdata2_0.975, LarvaePresent > 0.97,
                         select=c(Flowerspikes,PlantHeight,LarvaePresent))

# site size at 2.5% occupancy
newdata2_0.025 <- subset(newdata2, LarvaePresent < 0.03,
                         select=c(Flowerspikes,PlantHeight,LarvaePresent))
newdata2_0.025 <- subset(newdata2_0.025, LarvaePresent > 0.02,
                         select=c(Flowerspikes,PlantHeight,LarvaePresent))


g3 <- ggplot(dframe2,
             aes(x=Flowerspikes, y=PlantHeight))+
  geom_point(size=2,colour="black",fill="white",stroke=1, position=position_jitter(w=0.5),
             aes(shape=factor(LarvaePresent)))+
  scale_shape_manual(values=c(21,16),
                     labels=c("Absent","Present"),
                     name=expression(paste(italic("S. lychnitis "),"larvae")))+
  scale_y_continuous(limits = c(0, 220))+
  scale_x_continuous(breaks = seq(1, 12, 1))+
  geom_smooth(colour="black",
              data = newdata2_0.5, aes(x=Flowerspikes, y=PlantHeight),
              method=lm,se=FALSE,
              fullrange=TRUE)+
  geom_smooth(linetype="dashed",
              colour="black",
              data = newdata2_0.975, aes(x=Flowerspikes, y=PlantHeight),
              method=lm,se=FALSE,
              fullrange=TRUE)+
  geom_smooth(linetype="dashed",
              colour="black",
              data = newdata2_0.025, aes(x=Flowerspikes, y=PlantHeight),
              method=lm,se=FALSE,
              fullrange=TRUE)+
  xlab(expression(paste("No. ",
                        italic("V. nigrum "), "flowerspikes per plant")))+
  ylab(expression(paste(italic("V. nigrum "), "plant height (cm)")))+
  theme(panel.background=element_rect(fill="white"),
        panel.grid.major.x=element_line(colour="gray90"),
        panel.grid.major.y=element_line(colour="gray90"),
        panel.grid.minor=element_blank(),
        panel.border=element_rect(color="black",fill=F,size=1),
        text=element_text(size=15),
        axis.text=element_text(color="black"),
        legend.title=element_text(),
        legend.background=element_rect(fill="white", colour="black"),
        legend.key = element_blank(),
        legend.justification=c(1,1), legend.position=c(0.95,0.95))

g3








# sample code

newdata3<-expand.grid(Light=(c("HPS","LED")),Regime=(c("AllNight","Midnight")),Distance=0,SeedSetYN=1)
mm3<-model.matrix(terms(model3YN),newdata3)
newdata3$SeedSetYN = exp(mm3 %*% fixef(model3YN))
pvar3 <- diag(mm3 %*% tcrossprod(vcov(model3YN),mm3))
newdata3 <- data.frame(
  newdata3
  , plo = newdata3$SeedSetYN-1.96*sqrt(pvar3)
  , phi = newdata3$SeedSetYN+1.96*sqrt(pvar3)
)
newdata3 



#Plot


g3 <- ggplot(newdata3,
             aes(x=Light, y=SeedSetYN, fill=Regime))+
  geom_bar(colour="black",stat="identity",position=position_dodge())+
  scale_fill_manual(values=c("goldenrod","gray30"))+
  guides(fill=FALSE)+
  xlab("Pollinator treatment")+ ylab("Average seed count per seed capsule")+ 
  geom_errorbar(aes(ymin = plo, ymax = phi),
                colour="black",size=0.8, width=0.5,
                position=position_dodge(width=0.9))+
  theme(panel.background=element_rect(fill="white"),
        panel.grid.major=element_blank(),
        panel.grid.minor=element_blank(),
        panel.border=element_rect(color="black",fill=F,size=1))

g3

#alternate method

model3YNs <- glm(SeedSetYN ~ Light + Regime + Distance,
                 family = binomial (link = logit),
                 data = dframeYNa)

summary(model3YNs)
drop1(model3YNs, test = "Chi")   # very similar outputs so fine to proceed

newdata3a<-expand.grid(Light=(c("HPS","LED")),Regime=(c("AllNight","Midnight")),Distance=0,SeedSetYN=1)
newdata3a$SeedSetYN <- predict(model3YNs,newdata=newdata3a,type="response")
preddat <- predict(model3YNs,newdata=newdata3a,type="response",se.fit=TRUE)
preddat

newdata3a <- cbind(newdata3a,preddat)
newdata3a

newdata3a <- data.frame(
  newdata3a
  , plo = newdata3a$fit-1.96*newdata3a$se.fit
  , phi = newdata3a$fit+1.96*newdata3a$se.fit
)
newdata3a

newdata3a$Regime <- revalue(newdata3a$Regime, c("AllNight"="Full night","Midnight"="Part night"))
newdata3a

#Plot
