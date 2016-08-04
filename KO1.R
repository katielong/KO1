rm(list=ls())
library(caret)
library(ROCR)
library(randomForest)
library(pROC)
library(car)
library(DMwR)
library(dplyr)
library(reshape2)
library(ggplot2)

KO1 <- memisc::spss.system.file("...") %>%
       memisc::as.data.set() %>%  
       as.data.frame()

##DATA PRE-PROCESSING
#Significant variables under significance level of 0.1
#illcit drug use/marijuana use/heroin use/cocaine use/years in jail
KO1$yearsjail <- KO1$daysjail19/365 #years in jail

# Significant variables under significance level of 0.1-0.2
# age/marital status/months since release/months served
KO1$monthsrelease <- KO1$daysrelease9/30 #months since released from prison
KO1$monthserve <- KO1$daystimeserved11/30 #months served before release 

# recode & clean
#check if the drugafterrel20 variable includes only marijuana, heroin and cocain
#if so, no need for drugafterrel20 variable
test <- subset(KO1,drugafterrel20=="No")
test[,c("drugafterrel20", "marijuanause", "heroinuse", "cocaineuse")]

sub <- select(KO1, 
              drugafterrel20, # because no marijauana use, heoin use, cocaineuse implies no drug
              marijuanause, heroinuse, cocaineuse, 
              yearsjail, 
              daysjail19, #vars with p-value below 0.1
              age2, maritalstatus5,
              #monthsrelease, monthserve, 
              daysrelease9, daystimeserved11, #vars with p-value above 0.1 but less than 0.2
              sexunprafterrel22) #outcome

# sub$drugafterrel20 <- as.factor(car::recode(sub$drugafterrel20, "1='Yes';0='No';9999=NA"))
sub$cocaineuse <- car::recode(sub$cocaineuse,"'No'='No';else='Yes'")
sub$class <- car::recode(sub$sexunprafterrel22,"NA=0") %>% #recode NA as 0 by study design
                          car::recode("1='Yes';0='No'")
sub$maritalstatus5 <- recode(sub$maritalstatus5,"'Never been married'='Never been married'; else='Others'")
sub <- select(sub,-sexunprafterrel22)
summary(sub)

#significant under .1 --> risk ratio using log-binomial regression
drug <- glm(class ~ drugafterrel20, data=sub, family=binomial(link="log"))
mari <- glm(class ~ marijuanause, data=sub, family=binomial(link="log"))
hero <- glm(class ~ heroinuse, data=sub, family=binomial(link="log"))
coca <- glm(class ~ cocaineuse, data=sub, family=binomial(link="log"))
jail <- glm(class ~ yearsjail, data=sub, family=binomial(link="log"), start=c(-3.416677e-01,-7.959915e-05))

rr <- list(drug, mari, hero, coca, jail)
lapply(rr, coef) %>% lapply(exp)
lapply(rr, confint) %>% lapply(exp)


#EDA
colnames(sub)
featurePlot(x=sub[,c("daysjail19", "age2", "daysrelease9", "daystimeserved11")],
            y=sub$class, plot="pairs", auto.key=list(columns=2),main="Unprotected Sex")
vcd::mosaic(class ~ marijuanause,
            data=sub)
vcd::mosaic(class ~ heroinuse,
            data=sub)
vcd::mosaic(class ~ cocaineuse,
            data=sub)
vcd::mosaic(class ~ maritalstatus5,
            data=sub)


###RANDOM FOREST
sub$class <- factor(sub$class, levels=rev(c("No", "Yes"))) #reverse for better interpretation of sensitivity/specificity
sub <- select(sub,-drugafterrel20, -yearsjail) #exclude variables not used in random forest

#calculate statistics on resampling subsets
SevenStat <- function(data, lev = c("No", "Yes"), model = NULL) {
  accKapp <- postResample(data[,"pred"], data[,"obs"]) #accuracy/kappa
  ROC <- twoClassSummary(data,lev=lev,model=NULL)[1] #roc-auc
  Sens <- sensitivity(data[, "pred"], data[, "obs"], "Yes") #sensitivity
  Spec <- specificity(data[, "pred"], data[, "obs"], "No") #specificity
  Prec <- sum(data[,"pred"]=="Yes" & data[,"obs"]=="Yes")/sum(data[,"pred"]=="Yes") #precision
  F1 <- 2*Prec*Sens/(Prec+Sens) #f1 score
  out <- c(accKapp,
           ROC,
           Sens, Spec,
           Prec,F1)
  names(out)[c(3:7)] <- c("ROC","Sens", "Spec","Prec","F1")
  out}

##repeated cv
control_cv <- trainControl(method="repeatedcv",
                           repeats = 5,
                           classProbs=T,
                           savePredictions = T,
                           summaryFunction = SevenStat)
##tuning
grid <- expand.grid(mtry=seq(1:6))

##no subsampling
set.seed(3239)
control_cv$sampling <- NULL
rf <- train(class ~ .,
               data = sub,
               trControl=control_cv,##
               tuneGrid = grid,
               method = "rf",
               metric = "ROC",
               ntree = 1500)

##down subsampling inside of resampling
set.seed(3239)
control_cv$sampling <- "down"
rf_d <- train(class ~ .,
            data = sub,
            trControl=control_cv,##
            tuneGrid = grid,
            method = "rf",
            metric="ROC",
            ntree=1500)

##BRF
set.seed(3239)
control_cv$sampling <- NULL
rf_id <- train(class ~ .,
              data = sub,
              trControl=control_cv,##
              tuneGrid = grid,
              method = "rf",
              metric = "ROC",
              sampsize=rep(nrow(subset(sub,class=="No")),2),
              strata = sub$class,
              ntree = 1500)

#compare accuracy across resampling and Type methods
p <- rbind(RF=cbind(rf$results[,c("mtry","ROC")], Type = rep("RF",6)),
           BRF=cbind(rf_id$results[,c("mtry","ROC")],Type = rep("BRF",6)),
           Down=cbind(rf_d$results[,c("mtry","ROC")], Type = rep("DOWN",6))) %>% 
      melt(id=c("mtry","Type"),value.name="ROC") 

#plot
ggplot(p, aes(x=mtry, y=ROC)) + geom_line(aes(color=Type)) + geom_point(aes(color=Type,shape=Type),size=2.5)

#resampling sets
models <- list(RF=rf,BRF=rf_id,Down=rf_d)
resamps <- resamples(models)
resamp <- summary(resamps)
resamp

#importance plot
lapply(models,varImp,scale=F)
lapply(models,varImp,scale=T)
