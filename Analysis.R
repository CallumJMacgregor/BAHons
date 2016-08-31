#########################################################################
### Factors determining distribution and abundance of Striped Lychnis ###
#########################################################################

### Clear the workspace
rm(list=ls())


### install if necessary and then load the libraries you need

j <- c("lme4","bbmle","ggplot2","RVAideMemoire","arm","MuMIn","plyr","coefplot","scales","multcomp","sp","maps","ggmap","mapproj","rworldmap","ordinal")

new.packages <- j[!(j %in% installed.packages()[,"Package"])]
if(length(new.packages)) install.packages(new.packages)

lapply(j, require, character.only = TRUE)  # loads up any libraries that aren't already loaded


### load up Callum's custom set of functions
f <- c("CheckResidsFunction.R")
lapply(f, source)

dframe1 <- read.csv("Data/Sites.csv")

dframe1$Site <- factor(dframe1$Site)
dframe1$fWeevilDensity <- factor(dframe1$WeevilDensity)
dframe1$WeevilDensity <- ordered(dframe1$WeevilDensity, levels = c("N","R","O","F","A","D"))
dframe1$PlantDensity <- ordered(dframe1$PlantDensity, levels = c("R","O","F","A","D"))
dframe1$WildflowerVerge <- factor(dframe1$WildflowerVerge)
dframe1$TetradIsolation <- ordered(dframe1$TetradIsolation, levels = c(0,1,2,3,4,5,6,7,8))
dframe1$WindShelter <- ordered(dframe1$WindShelter, levels = c("P","G","E","C"))

summary(dframe1)

### Larva presence/absence

plot(LarvaePresent ~ PlantCount, data = dframe1)
plot(LarvaePresent ~ NNSmoothed, data = dframe1)
plot(LarvaePresent ~ WeevilDensity, data = dframe1)
plot(LarvaePresent ~ MeanSward, data = dframe1)
plot(LarvaePresent ~ WindShelter, data = dframe1)

model1 <- glm(LarvaePresent ~ 
                PlantCount + MeanSward + 
                WeevilDensity + NNSmoothed +
                WindShelter,
                family = binomial (link = "logit"),
                data = dframe1)

summary(model1)
drop1(model1, test="Chisq")

### test the different variables for...
# patch quality

model1a <- glm(LarvaePresent ~ 
                 PlantDensity + MeanSward + 
                 WeevilDensity + NNSmoothed +
                 WindShelter,
               family = binomial (link = "logit"),
               data = dframe1)

summary(model1a)
drop1(model1a, test = "Chi")

AIC(model1,model1a) # model1 (PlantCount) is better

# patch isolation

model1b <- glm(LarvaePresent ~ 
                PlantCount + MeanSward + 
                WeevilDensity + NN +
                 WindShelter,
              family = binomial (link = "logit"),
              data = dframe1)

summary(model1b)
drop1(model1b, test = "Chi")

model1c <- glm(LarvaePresent ~ 
                 PlantCount + MeanSward + 
                 WeevilDensity + TetradIsolation +
                 WindShelter,
               family = binomial (link = "logit"),
               data = dframe1)

summary(model1c)
drop1(model1c, test = "Chi")

AIC(model1,model1b,model1c) # again model1 is best (NNSmoothed), as model1c did not converge
AICtab(model1,model1b,model1c)

summary(model1)
drop1(model1, test = "Chi")

cor.test(dframe1$PlantCount,dframe1$NNSmoothed) # no correlation between variables

# try sequentially dropping variables

model1.1 <- glm(LarvaePresent ~ 
                PlantCount + MeanSward + 
                NNSmoothed +
                WindShelter,
              family = binomial (link = "logit"),
              data = dframe1)

summary(model1.1)
drop1(model1.1, test="Chisq")


model1.2 <- glm(LarvaePresent ~ 
                  PlantCount + MeanSward + 
                  NNSmoothed,
                family = binomial (link = "logit"),
                data = dframe1)

summary(model1.2)
drop1(model1.2, test="Chisq")


model1.3 <- glm(LarvaePresent ~ 
                  PlantCount + NNSmoothed,
                family = binomial (link = "logit"),
                data = dframe1)

summary(model1.3)
drop1(model1.3, test="Chisq")



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
plot(MaxLarvae ~ MeanSward, data = dframe1)
plot(MaxLarvae ~ WindShelter, data = dframe1)

mean(dframe1$MaxLarvae)
var(dframe1$MaxLarvae)


model2 <- glm(MaxLarvae ~ 
              PlantCount + MeanSward + 
              WeevilDensity + NNSmoothed +
                WindShelter,
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
plot(MaxLarvae ~ WindShelter, data = dframe1P)
plot(WindShelter ~ WeevilDensity, data = dframe1P)

mean(dframe1P$MaxLarvae)
var(dframe1P$MaxLarvae)

model2P <- glm(MaxLarvae ~
                 PlantCount + MeanSward +
                 WeevilDensity + NNSmoothed +
                 WindShelter,
               family = poisson (link = "log"),
               data = dframe1P)

summary(model2P)
drop1(model2P, test = "Chi")

### test the different variables for...
# patch quality

model2a <- glm(MaxLarvae ~ 
                 PlantDensity + MeanSward + 
                 WeevilDensity + NNSmoothed +
                 WindShelter,
               family = poisson (link = "log"),
               data = dframe1P)

summary(model2a)
drop1(model2a, test = "Chi")

AIC(model2P,model2a) # model2a (PlantDensity) has a lower AIC, but is based on only 2 categories
# therefore I'll choose to go with the continuous PlantCount anyway

# patch isolation

model2b <- glm(MaxLarvae ~ 
                 PlantCount + MeanSward + 
                 WeevilDensity + NN +
                 WindShelter,
               family = poisson (link = "log"),
               data = dframe1P)

summary(model2b)
drop1(model2b, test = "Chi")

model2c <- glm(MaxLarvae ~ 
                 PlantCount + MeanSward + 
                 WeevilDensity + TetradIsolation +
                 WindShelter,
               family = poisson (link = "log"),
               data = dframe1P)

summary(model2c)
drop1(model2c, test = "Chi")

AIC(model2P,model2b,model2c) # again model2P is best (NNSmoothed) (marginally)



summary(model2P)
drop1(model2P, test = "Chi")

# now start dropping variables

model2P1 <- glm(MaxLarvae ~
                 PlantCount + 
                 WeevilDensity + NNSmoothed +
                 WindShelter,
               family = poisson (link = "log"),
               data = dframe1P)

summary(model2P1)
drop1(model2P1, test = "Chi")


model2P2 <- glm(MaxLarvae ~
                  WeevilDensity + NNSmoothed +
                  WindShelter,
                family = poisson (link = "log"),
                data = dframe1P)

summary(model2P2)
drop1(model2P2, test = "Chi")


model2P3 <- glm(MaxLarvae ~
                  WeevilDensity + NNSmoothed,
                family = poisson (link = "log"),
                data = dframe1P)

summary(model2P3)
drop1(model2P3, test = "Chi")

plot(MaxLarvae ~ NNSmoothed, data = dframe1P)

# on visual inspection of the data, one point appears to be skewing the results here, so let's drop it and repeat the whole process 
# site 32 is the one causing the problem
dframe1Px <- dframe1P[ -which(dframe1P$Site=='32'), ]

# now:
plot(MaxLarvae ~ NNSmoothed, data = dframe1Px)
plot(MaxLarvae ~ WeevilDensity, data = dframe1Px)


model2Px <- glm(MaxLarvae ~
                 PlantCount + MeanSward +
                 WeevilDensity + NNSmoothed +
                 WindShelter,
               family = poisson (link = "log"),
               data = dframe1Px)

summary(model2Px)
drop1(model2Px, test = "Chi")

### test the different variables for...
# patch quality

model2ax <- glm(MaxLarvae ~ 
                 PlantDensity + MeanSward + 
                 WeevilDensity + NNSmoothed +
                 WindShelter,
               family = poisson (link = "log"),
               data = dframe1Px)

summary(model2ax)
drop1(model2ax, test = "Chi")

AIC(model2Px,model2ax) # models have the same AIC so let's go with continuous PlantCount variable


# patch isolation

model2bx <- glm(MaxLarvae ~ 
                 PlantCount + MeanSward + 
                 WeevilDensity + NN +
                 WindShelter,
               family = poisson (link = "log"),
               data = dframe1Px)

summary(model2bx)
drop1(model2bx, test = "Chi")

model2cx <- glm(MaxLarvae ~ 
                 PlantCount + MeanSward + 
                 WeevilDensity + TetradIsolation +
                 WindShelter,
               family = poisson (link = "log"),
               data = dframe1Px)

summary(model2cx)
drop1(model2cx, test = "Chi")

AIC(model2Px,model2bx,model2cx) # again all models the same (not enough data, probably), so let's go with NNSmoothed for consistency



summary(model2Px)
drop1(model2Px, test = "Chi")

# now start dropping variables

# mean sward
model2P1x <- glm(MaxLarvae ~
                  PlantCount + 
                  WeevilDensity + NNSmoothed +
                  WindShelter,
                family = poisson (link = "log"),
                data = dframe1Px)

summary(model2P1x)
drop1(model2P1x, test = "Chi")

# next plant count
model2P2x <- glm(MaxLarvae ~
                  WeevilDensity + NNSmoothed +
                  WindShelter,
                family = poisson (link = "log"),
                data = dframe1Px)

summary(model2P2x)
drop1(model2P2x, test = "Chi")

# now wind shelter
model2P3x <- glm(MaxLarvae ~
                  WeevilDensity + NNSmoothed,
                family = poisson (link = "log"),
                data = dframe1Px)

summary(model2P3x)
drop1(model2P3x, test = "Chi")

plot(MaxLarvae ~ NNSmoothed, data = dframe1Px)

# surprisingly, the result is the SAME as with the outlier included (including direction of effects);
# so let's go back to the model with all in

summary(model2P3)
drop1(model2P3, test = "Chi")


# posthoc comparisons between categories


summary(glht(model2P, linfct = mcp(WeevilDensity="Tukey")))
summary(glht(model2a, linfct = mcp(WeevilDensity="Tukey")))




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



#############################################
### What predicts C. nigritarsis density? ###
#############################################

plot(WeevilDensity ~ NNSmoothed, data = dframe1)
plot(WeevilDensity ~ PlantCount, data = dframe1)


model4 <- clm(WeevilDensity ~ NNSmoothed + PlantCount
                              + MeanSward + WindShelter,
              data=dframe1, 
              link = "logit")

summary(model4)
drop1(model4, test = "Chi")
confint(model4)


# drop variables: wind shelter
model4.1 <- clm(WeevilDensity ~ NNSmoothed + PlantCount
                + MeanSward,
                data=dframe1, 
                link = "logit")

summary(model4.1)
drop1(model4.1, test = "Chi")


# drop variables: sward height
model4.2 <- clm(WeevilDensity ~ NNSmoothed + PlantCount,
                data=dframe1, 
                link = "logit")

summary(model4.2)
drop1(model4.2, test = "Chi")


#######################################################

### figures

coefplot(model1.3)
coefplot(model2P3)

names(dframe1)

### occupancy

### some predicted datasets

summary(model1.3)

# probability of patch occupancy
newdata<-expand.grid(NNSmoothed=(seq(0.5,12,0.1)),PlantCount=(seq(1,350,1)),LarvaePresent=1)
newdata$LarvaePresent <- predict(model1.3,newdata=newdata,type="response")
preddat <- predict(model1.3,newdata=newdata,type="response",se.fit=TRUE)
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
                   legend.justification=c(1,0), legend.position=c(0.98,0.05))

g1


### abundance
names(dframe1)
summary(model2P3)
drop1(model2P3, test="Chi")

plot(MaxLarvae ~ WeevilDensity, data = dframe1P)
plot(MaxLarvae ~ NNSmoothed, data = dframe1P)

model2af <- glm(MaxLarvae ~ 
                 fWeevilDensity + NNSmoothed,
               family = poisson (link = "log"),
               data = dframe1P)

summary(model2af)
drop1(model2af, test = "Chi")


# larval abundance - predicted
newdata1<-expand.grid(NNSmoothed=2.038,WeevilDensity=c("N","R","O","F","A","D"),MaxLarvae=1)
preddat1 <- predict(model2P3,newdata=newdata1,type="response",se.fit=TRUE)
preddat1

newdata1 <- cbind(newdata1,preddat1)
newdata1 <- data.frame(
  newdata1
  , plo = newdata1$fit-1.96*newdata1$se.fit
  , phi = newdata1$fit+1.96*newdata1$se.fit
  , selo = newdata1$fit-newdata1$se.fit
  , sehi = newdata1$fit+newdata1$se.fit
)

newdata1$WeevilDensity <- ordered(newdata1$WeevilDensity,levels=c("N","R","O","F","A","D"))
newdata1$order <- ifelse(newdata1$WeevilDensity=="N",1,
                         ifelse(newdata1$WeevilDensity=="R",2,
                                ifelse(newdata1$WeevilDensity=="O",3,
                                       ifelse(newdata1$WeevilDensity=="F",4,
                                              ifelse(newdata1$WeevilDensity=="A",5,
                                                     ifelse(newdata1$WeevilDensity=="D",6,
                                                            "error"))))))
newdata1$order <- as.numeric(as.character(newdata1$order))
summary(newdata1$order)

# larvae abundance - model S.E. on 
newdata1a <- data.frame(
  newdata1
  , plo = newdata1$fit-1.96*newdata1$se.fit
  , phi = newdata1$fit+1.96*newdata1$se.fit
  , selo = newdata1$fit-newdata1$se.fit
  , sehi = newdata1$fit+newdata1$se.fit
)


g2 <- ggplot(newdata1,
              aes(x=WeevilDensity, y=fit))+
  scale_y_continuous(limits = c(0,25), oob=squish)+
  scale_x_discrete(breaks=c("N", "R", "O","F","A","D"),
                   labels=c("None", "Rare", "Occasional","Frequent","Abundant","Dominant"))+
  guides(fill=FALSE)+
  xlab(expression(paste("Density of ", italic("C. nigritarsis"), " larvae")))+
  ylab(expression(paste("Maximum no. ", italic("S. lychnitis"), " larvae")))+ 
  theme(panel.background=element_rect(fill="white"),
        panel.grid.major.x=element_blank(),
        panel.grid.minor=element_blank(),
        panel.grid.major.y=element_line(colour="gray60"),
        panel.border=element_rect(color="black",fill=F,size=1),
        axis.text.x=element_text(angle=45, hjust=1),
        axis.text=element_text(color="black"),
        text=element_text(size=15))+
  geom_point(data = dframe1P, aes(x=WeevilDensity, y=MaxLarvae),
             size=3,shape=21,colour="black",fill=NA,stroke=1.5, position=position_jitter(w=0.2,h=0.2))+
  geom_smooth(data = newdata1, aes(x=order, y=fit),
              colour="black",
              method = "lm", se = FALSE)+
  geom_smooth(data = newdata1, aes(x=order, y=selo),
              colour="black",
              linetype="dashed",
              method = "lm", se = FALSE)+
  geom_smooth(data = newdata1, aes(x=order, y=sehi),
              colour="black",
              linetype="dashed",
              method = "lm", se = FALSE)

g2



# larval abundance - predicted by NNSmoothed
newdata1i<-expand.grid(NNSmoothed=(seq(0.5,7.5,0.1)),WeevilDensity=c("N","R","O","F","A","D"),MaxLarvae=1)
preddat1i <- predict(model2P3,newdata=newdata1i,type="response",se.fit=TRUE)
preddat1i

newdata1i <- cbind(newdata1i,preddat1i)
newdata1i <- data.frame(
  newdata1i
  , plo = newdata1i$fit-1.96*newdata1i$se.fit
  , phi = newdata1i$fit+1.96*newdata1i$se.fit
  , selo = newdata1i$fit-newdata1i$se.fit
  , sehi = newdata1i$fit+newdata1i$se.fit
)




g2i<- ggplot(newdata1i,
             aes(x=NNSmoothed, y=fit))+
  scale_y_continuous(limits = c(0,25), oob=squish)+
  scale_x_continuous(limits = c(0,8), oob=squish)+
  guides(fill=FALSE)+
  xlab("Isolation (smoothed nearest neighbour distance, km)")+
  ylab(expression(paste("Maximum no. ", italic("S. lychnitis"), " larvae")))+ 
  theme(panel.background=element_rect(fill="white"),
        panel.grid.major.x=element_blank(),
        panel.grid.minor=element_blank(),
        panel.grid.major.y=element_line(colour="gray60"),
        panel.border=element_rect(color="black",fill=F,size=1),
        axis.text=element_text(color="black"),
        text=element_text(size=15))+
  geom_point(data = dframe1P, aes(x=NNSmoothed, y=MaxLarvae),
             size=3,shape=21,colour="black",fill=NA,stroke=1.5, position=position_jitter(w=0.2,h=0.2))+
  geom_smooth(data = newdata1i, aes(x=NNSmoothed, y=fit),
              colour="black",
              method = "lm", se = FALSE)+
  geom_smooth(data = newdata1i, aes(x=NNSmoothed, y=selo),
              colour="black",
              linetype="dashed",
              method = "lm", se = FALSE)+
  geom_smooth(data = newdata1i, aes(x=NNSmoothed, y=sehi),
              colour="black",
              linetype="dashed",
              method = "lm", se = FALSE)

g2i




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


g3a <- ggplot(dframe2,
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

g3a



########### map plot ######################

summary(dframe1$GridRef)
summary(dframe1$lat)
summary(dframe1$lon)


gadm <- readRDS("Data/GBR_adm2.rds")
plot(gadm, xlim = c(-1.8, -0.7), ylim = c(50.6, 51.4))
points(dframe1$lon,dframe1$lat,col=as.numeric(dframe1$LarvaePresent),cex=2,pch=c(21,16))




ggmap(gadm)

map9 <- get_map(location = 'Peninsula Road,Winchester,England', maptype = "roadmap",zoom = 9, color='bw')
ggmap(map9)




mapPoints9 <- ggmap(map9, extent = "device", crop = T) +
               geom_point(data = dframe1, size=5,colour="black",fill="white",stroke=1,alpha=0.6,
                          aes(x = lon, y = lat,shape=factor(LarvaePresent)))+
               scale_shape_manual(values=c(21,16),
                                  labels=c("Absent","Present"),
                                  name=expression(paste(italic("S. lychnitis "),"larvae")))+
               theme(legend.justification=c(0,0), legend.position=c(0.07,0.77))+ 
               theme(legend.text = element_text(size=15))+
               theme(legend.title = element_text(size=15))+
               theme(legend.background = element_rect(fill="white",size=0.5, linetype="solid",colour ="black"))
                

mapPoints9


map9c <- get_map(location = 'Peninsula Road,Winchester,England', maptype = "roadmap",zoom = 9)
ggmap(map9c)




mapPoints9c <- ggmap(map9c, extent = "device", crop = T) +
  geom_point(data = dframe1, size=3,colour="black",fill="white",stroke=1,alpha=0.8,
             aes(x = lon, y = lat,shape=factor(LarvaePresent)))+
  scale_shape_manual(values=c(21,16),
                     labels=c("Absent","Present"),
                     name=expression(paste(italic("S. lychnitis "),"larvae")))+
  theme(legend.justification=c(0,0), legend.position=c(0.07,0.77))+ 
  theme(legend.background = element_rect(fill="white",size=0.5, linetype="solid",colour ="black"))


mapPoints9c














map10 <- get_map(location = 'Peninsula Road,Winchester,England', zoom = 10, maptype="roadmap", color='bw')
ggmap(map10)




mapPoints10 <- ggmap(map10, extent = "device") +
  geom_point(data = dframe1, size=3,colour="black",fill="white",stroke=1,alpha=0.8,
             aes(x = lon, y = lat,shape=factor(LarvaePresent)))+
  scale_shape_manual(values=c(21,16),
                     labels=c("Absent","Present"),
                     name=expression(paste(italic("S. lychnitis "),"larvae")))+
  theme(legend.justification=c(0,0), legend.position=c(0.07,0.77))

mapPoints10















