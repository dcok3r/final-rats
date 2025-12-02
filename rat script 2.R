
library(tidyverse)
library(ggplot2)
library(brms)
library(sjPlot)

options(mc.cores = parallel::detectCores())

## Load and subset data
df<-read.csv("~/NRES 746/746 Project/2015_05_06_ExternalMorphology_Middens_Distance_Adults_3species.csv",
             header=TRUE,na.strings=c("",NA))

rats<-df%>%filter(Species!='hybrid')

## Processing ====
for (x in 1:nrow(rats)) {
  if(rats$Zone.Number[x]==3){
    rats$dist[x]<-0
  } else if(rats$Zone.Number[x]==2 | rats$Zone.Number[x]==4 ){
    rats$dist[x]<-1
  } else {
    rats$dist[x]<-2
  }
}

rats$male<-ifelse(rats$Sex=="Male",1,0)

rats$macrotis<-ifelse(rats$Species=="macrotis",1,0)

rats<-rats%>%mutate(across(c(5,8:10),scale))


## PCA ====

#see if pregnant females should be excluded
df$Pregnant<-ifelse(df$Sex=="Female" & is.na(df$Pregnant),"no",df$Pregnant)
with(df%>%filter(Sex=="Female"),t.test(Weight..gram.~factor(Pregnant))) #no difference in weight

rats2<-rats%>%subset(select=c("ID","Weight..gram.","Ear..mm.","Rostrum..mm.","Hind.Foot..mm."))

cov(subset(rats2,select=-1), use = "complete.obs")
cor(subset(rats2,select=-1), use = "complete.obs")

rats2<-na.omit(rats2)

PC<-princomp(subset(rats2,select=-1),cor=T)
summary(PC)

fviz_eig(PC,addlabels=T)
fviz_pca_var(PC, col.var="cos2",gradient.cols=c("blue","orange","green"),repel=TRUE)

rats2<-cbind(rats2,PC$scores)

rats3<-rats%>%left_join(rats2)

## Morphology predictive model ====

# look at some possible interactions- do these make sense to include?
ggplot(rats3,aes(x=male,y=Comp.1,color=factor(macrotis))) + geom_jitter() + 
  geom_smooth(method="lm")
ggplot(rats3,aes(x=dist,y=Comp.1,color=factor(macrotis))) + geom_jitter() + 
  geom_smooth(method="lm")
ggplot(rats3,aes(x=dist,y=Comp.1,color=factor(male))) + geom_jitter() + geom_smooth(method="lm")


PCmodel <- brm( Comp.1 ~ male + macrotis + dist + male:macrotis + dist:macrotis + dist:male,
  data = rats3,
  family = gaussian(),
  chains = 4,
  cores = 6,
  iter = 4000
)

summary(PCmodel)

stancode(PCmodel)

male=1
macrotis=1
dist=2
0.56+(1.47*male)-(1.35*macrotis)-(0.16*dist)-(0.69*male*macrotis)-(0.14*dist*macrotis)-
  (0.23*male*dist)

plot_model(model, type="pred",ppd=T,ci.lvl=0.90,terms = c("dist","male","macrotis")) + 
  ylab("Principal component")


library(tidybayes)
model %>% spread_draws(r_condition[condition,term]) %>%
  head(10)

rats3 %>%
  add_predicted_draws(model) %>%
  ggplot(aes(x = dist, y = Comp.1, color = factor(male), fill = factor(male))) +
  stat_lineribbon(aes(y = .prediction), .width = c(.95), alpha = 0.1) +
  geom_jitter(data=rats3) +
  facet_grid(. ~ macrotis)


rats3 %>%
  add_epred_draws(model, ndraws = 200) %>%
  ggplot(aes(x = dist, y = Comp.1, color = factor(male))) +
  geom_line(aes(y = .epred, group = paste(male, .draw)), alpha = 0.25) +
  geom_jitter(data = rats3) +
  facet_grid(. ~ macrotis)




## Relative morphology -----------------------------------------------
  # get mean weight and weight SD for each species in allopatric zones
df$dist<-NA
for (x in 1:nrow(df)) {
  if(df$Zone.Number[x]==3){
    df$dist[x]<-0
  } else if(df$Zone.Number[x]==2 | df$Zone.Number[x]==4 ){
    df$dist[x]<-1
  } else {
    df$dist[x]<-2
  }
}

rats4<-df%>%filter(Species!='hybrid' & dist==2)%>%group_by(Species)%>%summarize(
  meanwgt=mean(Weight..gram.,na.rm=T), wgt.sd=sd(Weight..gram.,na.rm=T),
  meanear=mean(Ear..mm.,na.rm=T), ear.sd=sd(Ear..mm.,na.rm=T),
  meanfoot=mean(Hind.Foot..mm.,na.rm=T), foot.sd=sd(Hind.Foot..mm.,na.rm=T),
  meanrost=mean(Rostrum..mm.,na.rm=T), rost.sd=sd(Rostrum..mm.,na.rm=T)  ) %>% 
  right_join(.,df%>%filter(Species!='hybrid'),by="Species")

rats4$male<-ifelse(rats4$Sex=="Male",1,0)

rats4$macrotis<-ifelse(rats4$Species=="macrotis",1,0)

  # add morphology z-scaled within species in allopatric zone
rats4$z_weight<-(rats4$Weight..gram.-rats4$meanwgt)/rats4$wgt.sd
rats4$z_ear<-(rats4$Ear..mm.-rats4$meanear)/rats4$ear.sd
rats4$z_foot<-(rats4$Hind.Foot..mm.-rats4$meanfoot)/rats4$foot.sd
rats4$z_rost<-(rats4$Rostrum..mm.-rats4$meanrost)/rats4$rost.sd


zwgtmodel <- brm(z_weight ~ male + macrotis + dist + male:macrotis + dist:macrotis + dist:male +
  (1|Year), data = rats4, family = gaussian(), chains = 4, cores = 6, iter = 5000 )
summary(zwgtmodel)
plot_model(zwgtmodel, type="pred",ppd=T,ci.lvl=0.95,terms = c("dist","male","macrotis")) + 
  ylab("z-transformed weight")


zearmodel <- brm(z_ear ~ male + macrotis + dist + male:macrotis + dist:macrotis + dist:male,
  data = rats4, family = gaussian(), chains = 4, cores = 6, iter = 4000 )
summary(zearmodel)
plot_model(zearmodel, type="pred",ppd=T,ci.lvl=0.95,terms = c("dist","male","macrotis")) + 
  ylab("z-transformed ear length")


zfootmodel <- brm(z_foot ~ male + macrotis + dist + male:macrotis + dist:macrotis + dist:male,
  data = rats4, family = gaussian(), chains = 4, cores = 6, iter = 4000 )
summary(zfootmodel)
plot_model(zfootmodel, type="pred",ppd=T,ci.lvl=0.95,terms = c("dist","male","macrotis")) + 
  ylab("z-transformed hindfoot length")


zrostmodel <- brm(z_rost ~ male + macrotis + dist + male:macrotis + dist:macrotis + dist:male,
  data = rats4, family = gaussian(), chains = 4, cores = 6, iter = 4000 )
summary(zrostmodel)
plot_model(zrostmodel, type="pred",ppd=T,ci.lvl=0.95,terms = c("dist","male","macrotis")) + 
  ylab("z-transformed rostrum length")


## Interpretable measure means ====
morphosumm<-rats4%>%group_by(Species,dist,Sex)%>%summarize(meanwgt=mean(Weight..gram./10,na.rm=T),
          wgt.sd=sd(Weight..gram.,na.rm=T), meanear=mean(Ear..mm.,na.rm=T),
          ear.sd=sd(Ear..mm.,na.rm=T), meanfoot=mean(Hind.Foot..mm.,na.rm=T), 
          foot.sd=sd(Hind.Foot..mm.,na.rm=T), meanrost=mean(Rostrum..mm.,na.rm=T),
          rost.sd=sd(Rostrum..mm.,na.rm=T))

  # Weight
ggplot(rats4,aes(x=factor(dist),y=Weight..gram.,fill=Sex)) + geom_boxplot() +
  facet_wrap(~ Species)

  # Ear length
ggplot(rats4,aes(x=factor(dist),y=Ear..mm.,fill=Sex)) + geom_boxplot() +
  facet_wrap(~ Species)
  
  # Hind foot length
ggplot(rats4,aes(x=factor(dist),y=Hind.Foot..mm.,fill=Sex)) + geom_boxplot() +
  facet_wrap(~ Species)
  
  # Rostrum length
ggplot(rats4,aes(x=factor(dist),y=Rostrum..mm.,fill=Sex)) + geom_boxplot() +
  facet_wrap(~ Species)
  

## Ectoparasites ====

  # visualize
  parasite.summ<-rats4%>%filter(!is.na(Ectoparasite.Load))%>%
    group_by(Species,Ectoparasite.Load,dist)%>%summarize(count=n())
    
  parasite.summ<-rats4%>%filter(!is.na(Ectoparasite.Load))%>%
    group_by(Species,dist)%>%summarize(count=n())%>%
    right_join(.,parasite.summ,by=c("Species","dist"))
    
  parasite.summ<-parasite.summ%>%mutate(prop=count.y/count.x)
    
ggplot(parasite.summ,aes(x=Ectoparasite.Load,y=prop,fill=factor(Species))) + 
  geom_histogram(stat="identity",position="dodge") + facet_grid(~dist)

  

library(MASS)
rats4$Ectoparasite.Load<-factor(rats4$Ectoparasite.Load, levels=c("None","Low",
                                      "Medium","High"),ordered=T)

  model<-polr(Ectoparasite.Load ~ macrotis + male + dist + Weight..gram. +
                macrotis:male + macrotis:dist + male:dist, data=rats4, Hess=T)
  summary(model)
  confint(model)

exp(cbind(OR = coef(model), confint(model)))
  

newdat <- data.frame(
  male = rep(0:1, 300),
  macrotis = rep(0:1, each = 300),
  dist = rep(0:2,each=200),
  Weight..gram. = rep(seq(from=min(rats4$Weight..gram.,na.rm=T), 
                          to=max(rats4$Weight..gram.,na.rm=T), length.out = 100), 6))

newdat <- cbind(newdat, predict(model, newdat, type = "probs"))
  
newdat <- newdat%>%pivot_longer(cols=c(5:8),names_to='Ectoparasite',values_to="probability")

ggplot(newdat, aes(x = dist, y = probability, colour = Ectoparasite)) +
  geom_line() + facet_grid(macrotis ~ male, labeller="label_both")










  