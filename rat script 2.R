
library(tidyverse)
library(ggplot2)
library(ggbeeswarm)
library(cmdstanr)
library(bayesplot)
library(factoextra)

options(mc.cores = parallel::detectCores())

## Load data
df<-read.csv("~/NRES 746/746 Project/2015_05_06_ExternalMorphology_Middens_Distance_Adults_3species.csv",
             header=TRUE,na.strings=c("",NA))
df<-df%>%filter(Species!='hybrid')
df$Sex<-factor(df$Sex,levels=c("Female","Male"))
df$Species<-factor(df$Species,levels=c("fuscipes","macrotis"))

for (x in 1:nrow(df)) {
  if(df$Zone.Number[x]==3){
    df$dist[x]<-0
  } else if(df$Zone.Number[x]==2 | df$Zone.Number[x]==4 ){
    df$dist[x]<-1
  } else {
    df$dist[x]<-2   }
}

df$male<-as.integer(df$Sex)-1

df$macrotis<-as.integer(df$Species)-1


## PCA ====
#see if pregnant females should be excluded
df$Pregnant<-ifelse(df$Sex=="Female" & is.na(df$Pregnant),"no",df$Pregnant)
with(df%>%filter(Sex=="Female"),t.test(Weight..gram.~factor(Pregnant))) #no difference in weight
#-- -- --

rats<-df%>%mutate(across(c(5,8:10),scale))

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

## Morphology (PC) model ====

# look at some possible interactions- do these make sense to include?
ggplot(rats3,aes(x=male,y=Comp.1,color=factor(macrotis))) + geom_jitter() + 
  geom_smooth(method="lm")
ggplot(rats3,aes(x=dist,y=Comp.1,color=factor(macrotis))) + geom_jitter() + 
  geom_smooth(method="lm")
ggplot(rats3,aes(x=dist,y=Comp.1,color=factor(male))) + geom_jitter() + geom_smooth(method="lm")

## create model matrix for PC model
PCdf<-rats3%>%reframe(Comp.1,species=as.integer(Species),macrotis,male,dist,
              macrotis_distance=macrotis*dist, macrotis_male=macrotis*male,
              male_distance=male*dist)

PCdf<-na.omit(PCdf)
PCdf$macrotis<-NULL

#list of stan inputs
PC_df <- list(
  N = nrow(PCdf), K = 5, X = PCdf[,c(3:7)], y = PCdf$Comp.1, n_sp = nlevels(factor(PCdf$species)),
  species = PCdf$species  )

stanmodel<-cmdstan_model("rats_stan.stan")

fit_PC <- stanmodel$sample(
  data=PC_df, chains = 6, parallel_chains = 6,  
  iter_warmup = 2000, iter_sampling = 4000
)

loo(fit_PC$draws("log_lik", format="matrix"))

PCsummary<-fit_PC$summary(); PCsummary

# ANOVA ---
aov_PC_fus<-rats3%>%filter(Species=="fuscipes" & !is.na(Comp.1))%>%
  mutate(dist=factor(dist,levels=c("0","1","2")),sex=factor(Sex))
summary(aov(Comp.1 ~ dist*sex, data=aov_PC_fus))
TukeyHSD(aov(Comp.1 ~ dist*sex, data=aov_PC_fus))

    ###
aov_PC_mac<-rats3%>%filter(Species=="macrotis" & !is.na(Comp.1))%>%
  mutate(dist=factor(dist, levels=c("0","1","2")),sex=factor(Sex))
summary(aov(Comp.1 ~ dist*sex, data=aov_PC_mac))
TukeyHSD(aov(Comp.1 ~ dist*sex, data=aov_PC_mac))

# ---

PCsummary2<-as.data.frame(PCsummary)
PCsummary2$variable<-c("lp","alpha.fus","alpha.mac","male","dist",
                       "mac_dist","mac_male","male_dist", "sigma")

draws_PC<-fit_PC$draws(format="df") # "lp","alpha.fus","alpha.mac","male","dist",
                                    # "mac_dist","mac_male","male_dist", "sigma"

# slope credible interval ---
draws_PC %>% reframe(fus_male_slope = `beta[2]` + `beta[5]` ,
         fus_fem_slope = `beta[2]`,
         mac_male_slope = `beta[2]` + `beta[3]` + `beta[5]`,
         mac_fem_slope = `beta[2]`+ `beta[3]`) %>% 
  summarize(fus_male=mean(fus_male_slope),fus_male_upr=quantile(fus_male_slope,probs=0.975), 
            fus_male_lwr=quantile(fus_male_slope,probs=0.025), ##
                        fus_fem=mean(fus_fem_slope),fus_fem_upr=quantile(fus_fem_slope,probs=0.975),
            fus_fem_lwr=quantile(fus_fem_slope,probs=0.025), ##
                        mac_male=mean(mac_male_slope), mac_male_upr=quantile(mac_male_slope,probs=0.975), 
            mac_male_lwr=quantile(mac_male_slope,probs=0.025), ##
                        mac_fem=mean(mac_fem_slope), mac_fem_upr=quantile(mac_fem_slope,probs=0.975), 
            mac_fem_lwr=quantile(mac_fem_slope,probs=0.025)) %>% as.data.frame()
# ---

mcmc_trace(draws_PC, pars = c("alpha[1]","alpha[2]","beta[1]","beta[2]",
                              "beta[3]","beta[4]","beta[5]","sigma"))

# --- FUSCIPES  PC PLOT ----
draws_PC2<-draws_PC[1:3000,] #subset draws for shorter plotting times

# Warning: plots take a bit of time to run and appear (< ~1min)
ggplot(rats3%>%filter(Species=="fuscipes"),aes(x=dist,y=Comp.1,color=factor(Sex))) + 
  geom_quasirandom(width=0.2) +
  geom_abline(intercept=draws_PC2$`alpha[1]`, #fus female intercept
               slope=draws_PC2$`beta[2]`, # distance beta
               alpha=0.05,color="pink") +
  geom_abline(intercept=draws_PC2$`alpha[1]`+ draws_PC2$`beta[1]`, #fus male intercept
               slope=draws_PC2$`beta[2]`+ draws_PC2$`beta[5]`, # distance + male:dist
               alpha=0.05,color="lightblue") +
  geom_abline(intercept=PCsummary2$mean[PCsummary2$variable=="alpha.fus"], #fem fus meanline
              slope=PCsummary2$mean[PCsummary2$variable=="dist"],
              linewidth=1,color="darkred") +
  scale_color_manual(values=c("red","blue"))+
  geom_abline(intercept=PCsummary2$mean[PCsummary2$variable=="alpha.fus"]+ #male fus meanline
                PCsummary2$mean[PCsummary2$variable=="male"], 
              slope=PCsummary2$mean[PCsummary2$variable=="dist"]+
                PCsummary2$mean[PCsummary2$variable=="male_dist"],
              linewidth=1,color="darkblue") +
  scale_x_continuous(breaks = c(0,1,2), 
                      labels = c("sympatry","near\nsympatry","allopatry") ) +
  scale_y_continuous(limits=c(-4,4.5),breaks=seq(-4,4.5,by=1)) +
  coord_cartesian(xlim = c(-0.25, 2.25)) +
  ylab ("Principal component 1") + 
  guides(color = guide_legend(title = "Sex")) +
  theme_minimal() + 
  theme(  panel.grid.minor = element_line(size=0.5),
          plot.title = element_text(hjust = 0.5),
          panel.border = element_rect(color = "grey", fill = NA),
          axis.title.x = element_blank(),
          legend.title = element_text(size=12,hjust=0.5),
          legend.text = element_text(size=12),
          title = element_text(size=15),
          axis.text.x = element_text(size=13),
          axis.text.y = element_text(size=13),
          strip.background = element_rect(fill = "lightgrey", linetype = "solid"),
          strip.text = element_text(size=14,face="italic") ) + 
    facet_wrap(~Species) #put box around title
    # save as image... dims 650 x 600 

# --- MACROTIS PC PLOT ----

ggplot(rats3%>%filter(Species=="macrotis"),aes(x=dist,y=Comp.1,color=factor(Sex))) + 
  geom_quasirandom(width=0.3) +
  geom_abline (intercept=draws_PC2$`alpha[2]`, #mac female intercept
                slope=draws_PC2$`beta[2]`+
                draws_PC2$`beta[3]`,  # slope=distance + macrotis:dist
               alpha=0.05,color="pink") +
  geom_abline (intercept=draws_PC2$`alpha[2]`+ draws_PC2$`beta[1]`+
                draws_PC2$`beta[4]` , #mac male intercept
                slope=draws_PC2$`beta[2]`+ draws_PC2$`beta[3]`+
                draws_PC2$`beta[5]`, # distance + male:dist
               alpha=0.05,color="lightblue") +
  geom_abline(intercept=PCsummary2$mean[PCsummary2$variable=="alpha.mac"], #fem mac meanline
              slope=PCsummary2$mean[PCsummary2$variable=="dist"]+
                    PCsummary2$mean[PCsummary2$variable=="mac_dist"],
              linewidth=1,color="darkred") +
  scale_color_manual(values=c("red","blue"))+
  geom_abline(intercept=PCsummary2$mean[PCsummary2$variable=="alpha.mac"]+ #male mac meanline
                    PCsummary2$mean[PCsummary2$variable=="mac_male"]+
                    PCsummary2$mean[PCsummary2$variable=="male"], 
              slope=PCsummary2$mean[PCsummary2$variable=="dist"]+
                    PCsummary2$mean[PCsummary2$variable=="mac_dist"]+
                    PCsummary2$mean[PCsummary2$variable=="male_dist"],
              linewidth=1,color="darkblue") +
  scale_x_continuous(breaks = c(0,1,2), 
                     labels = c("sympatry","near\nsympatry","allopatry") ) +
  scale_y_continuous(limits=c(-4,4.5),breaks=seq(-4,4.5,by=1)) +
  coord_cartesian(xlim = c(-0.25, 2.25)) +
  ylab ("Principal component 1") + 
  guides(color = guide_legend(title = "Sex")) +
  theme_minimal() + 
  theme(  panel.grid.minor = element_line(size=0.5),
          plot.title = element_text(hjust = 0.5),
          panel.border = element_rect(color = "grey", fill = NA),
          axis.title.x = element_blank(),
          legend.title = element_text(size=12,hjust=0.5),
          legend.text = element_text(size=12),
          title = element_text(size=15),
          axis.text.x = element_text(size=13),
          axis.text.y = element_text(size=13),
          strip.background = element_rect(fill = "lightgrey", linetype = "solid"),
          strip.text = element_text(size=14,face="italic") ) + 
  facet_wrap(~Species)

  

# look at the difference between bayesian and frequentist estimates
ggplot(rats3%>%filter(Species=="fuscipes"),aes(x=dist,y=Comp.1,color=factor(Sex))) + 
  geom_quasirandom(width=0.2) +
  geom_abline(intercept=PCsummary2$mean[PCsummary2$variable=="alpha.fus"], #fem fus meanline
              slope=PCsummary2$mean[PCsummary2$variable=="dist"],
              linewidth=1,color="darkred") +
  geom_abline(intercept=PCsummary2$mean[PCsummary2$variable=="alpha.fus"]+ #male fus meanline
                PCsummary2$mean[PCsummary2$variable=="male"], 
              slope=PCsummary2$mean[PCsummary2$variable=="dist"]+
                PCsummary2$mean[PCsummary2$variable=="male_dist"],
              linewidth=1,color="darkblue") +
        geom_smooth(method='lm')


ggplot(rats3%>%filter(Species=="macrotis"),aes(x=dist,y=Comp.1,color=factor(Sex))) + 
  geom_quasirandom(width=0.3) +
  geom_abline(intercept=PCsummary2$mean[PCsummary2$variable=="alpha.mac"], #fem mac meanline
              slope=PCsummary2$mean[PCsummary2$variable=="dist"]+
                PCsummary2$mean[PCsummary2$variable=="mac_dist"],
              linewidth=1,color="darkred") +
  geom_abline(intercept=PCsummary2$mean[PCsummary2$variable=="alpha.mac"]+ #male mac meanline
                PCsummary2$mean[PCsummary2$variable=="male"]+
                PCsummary2$mean[PCsummary2$variable=="mac_male"], 
              slope=PCsummary2$mean[PCsummary2$variable=="dist"]+
                PCsummary2$mean[PCsummary2$variable=="mac_dist"]+
                PCsummary2$mean[PCsummary2$variable=="male_dist"],
              linewidth=1,color="darkblue") +
        geom_smooth(method='lm')


## Relative morphology -----------------------------------------------
  # get median weight and weight MAD for each species in allopatric zones

rats4<-df%>%filter(Species!='hybrid' & dist==2)%>%group_by(Species)%>%summarize(
  medwgt=median(Weight..gram.,na.rm=T), wgt.mad=mad(Weight..gram.,na.rm=T),
  medear=median(Ear..mm.,na.rm=T), ear.mad=mad(Ear..mm.,na.rm=T),
  medfoot=median(Hind.Foot..mm.,na.rm=T), foot.mad=mad(Hind.Foot..mm.,na.rm=T),
  medrost=median(Rostrum..mm.,na.rm=T), rost.mad=mad(Rostrum..mm.,na.rm=T)  ) %>% 
  right_join(.,df%>%filter(Species!='hybrid'),by=c("Species"))

# add morphology z-scaled within species
rats4$z_weight<-(rats4$Weight..gram.-rats4$medwgt)/rats4$wgt.mad
rats4$z_ear<-(rats4$Ear..mm.-rats4$medear)/rats4$ear.mad
rats4$z_foot<-(rats4$Hind.Foot..mm.-rats4$medfoot)/rats4$foot.mad
rats4$z_rost<-(rats4$Rostrum..mm.-rats4$medrost)/rats4$rost.mad

# --- Z_Weight model ----
ZWgtMM<-rats4%>%reframe(z_weight,species=as.integer(Species),macrotis,male,dist,
                macrotis_distance=macrotis*dist, macrotis_male=macrotis*male,
                male_distance=male*dist)
ZWgtMM<-na.omit(ZWgtMM)
ZWgtMM$macrotis<-NULL

#list of stan inputs
ZWgt_df <- list(
  N = nrow(ZWgtMM), K = 5, X = ZWgtMM[,c(3:7)], y = ZWgtMM$z_weight, 
  n_sp = nlevels(factor(ZWgtMM$species)), species = ZWgtMM$species)

#stanmodel<-cmdstan_model("rats_stan.stan")

fit_Zwgt <- stanmodel$sample(
  data=ZWgt_df, chains = 6, parallel_chains = 6,  
  iter_warmup = 2000, iter_sampling = 4000 )

Zwgtsummary<-fit_Zwgt$summary(); Zwgtsummary

# ANOVA ---
aov_Zwgt_fus<-rats4%>%filter(Species=="fuscipes" & !is.na(z_weight))%>%
  mutate(dist=factor(dist,levels=c("0","1","2")),sex=factor(Sex))
summary(aov(z_weight ~ dist*sex, data=aov_Zwgt_fus))
TukeyHSD(aov(z_weight ~ dist*sex, data=aov_Zwgt_fus))

###
aov_Zwgt_mac<-rats4%>%filter(Species=="macrotis" & !is.na(z_weight))%>%
  mutate(dist=factor(dist, levels=c("0","1","2")),sex=factor(Sex))
summary(aov(z_weight ~ dist*sex, data=aov_Zwgt_mac))
TukeyHSD(aov(z_weight ~ dist*sex, data=aov_Zwgt_mac))
# ---

Zwgtsummary2<-as.data.frame(Zwgtsummary)
Zwgtsummary2$variable<-c("lp","alpha.fus","alpha.mac","male","dist",
                       "mac_dist","mac_male","male_dist", "sigma")

draws_Zwgt<-fit_Zwgt$draws(format="df")

# slope credible interval ---
draws_Zwgt %>% reframe(fus_male_slope = `beta[2]` + `beta[5]` ,
                     fus_fem_slope = `beta[2]`,
                     mac_male_slope = `beta[2]` + `beta[3]` + `beta[5]`,
                     mac_fem_slope = `beta[2]`+ `beta[3]`) %>% 
  summarize(fus_male=mean(fus_male_slope),fus_male_upr=quantile(fus_male_slope,probs=0.975), 
            fus_male_lwr=quantile(fus_male_slope,probs=0.025), ##
            fus_fem=mean(fus_fem_slope),fus_fem_upr=quantile(fus_fem_slope,probs=0.975),
            fus_fem_lwr=quantile(fus_fem_slope,probs=0.025), ##
            mac_male=mean(mac_male_slope), mac_male_upr=quantile(mac_male_slope,probs=0.975), 
            mac_male_lwr=quantile(mac_male_slope,probs=0.025), ##
            mac_fem=mean(mac_fem_slope), mac_fem_upr=quantile(mac_fem_slope,probs=0.975), 
            mac_fem_lwr=quantile(mac_fem_slope,probs=0.025)) %>% as.data.frame()
# ---

draws_Zwgt2<-draws_Zwgt[1:3000,] #subset draws

  # Fuscipes
ggplot(rats4%>%filter(Species=="fuscipes"),aes(x=dist,y=z_weight,color=factor(Sex))) + 
  geom_quasirandom(width=0.3) +
  geom_abline(intercept=draws_Zwgt2$`alpha[1]`, #fus female intercept
              slope=draws_Zwgt2$`beta[2]`, # distance beta
              alpha=0.05,color="pink") +
  geom_abline(intercept=draws_Zwgt2$`alpha[1]`+ draws_Zwgt2$`beta[1]`, #fus male intercept
              slope=draws_Zwgt2$`beta[2]`+ draws_Zwgt2$`beta[5]`, # distance + male:dist
              alpha=0.05,color="lightblue") +
  scale_color_manual(values=c("red","blue"))+
  geom_abline(intercept=Zwgtsummary2$mean[Zwgtsummary2$variable=="alpha.fus"], #fem fus meanline
              slope=Zwgtsummary2$mean[Zwgtsummary2$variable=="dist"],
              linewidth=1,color="darkred") +
  geom_abline(intercept=Zwgtsummary2$mean[Zwgtsummary2$variable=="alpha.fus"]+ #male fus meanline
                Zwgtsummary2$mean[Zwgtsummary2$variable=="male"], 
              slope=Zwgtsummary2$mean[Zwgtsummary2$variable=="dist"]+
                Zwgtsummary2$mean[Zwgtsummary2$variable=="male_dist"],
              linewidth=1,color="darkblue") +
  scale_x_continuous(breaks = c(0,1,2), 
                     labels = c("sympatry","near\nsympatry","allopatry") ) +
  scale_y_continuous(limits=c(-2.5,3),breaks=seq(-2,3,by=1)) +
  coord_cartesian(xlim = c(-0.25, 2.25)) +
  ylab ("Z-scaled weight") + 
  geom_hline(yintercept=0, linetype="dashed",linewidth=1.5) +
  guides(color = guide_legend(title = "Sex")) +
  theme_minimal() + 
  theme(  panel.grid.minor = element_line(size=0.5),
          plot.title = element_text(hjust = 0.5),
          panel.border = element_rect(color = "grey", fill = NA),
          axis.title.x = element_blank(),
          legend.title = element_text(size=12,hjust=0.5),
          legend.text = element_text(size=12),
          title = element_text(size=15),
          axis.text.x = element_text(size=13),
          axis.text.y = element_text(size=13),
          strip.background = element_rect(fill = "lightgrey", linetype = "solid"),
          strip.text = element_text(size=14,face="italic") ) + 
  facet_wrap(~Species)


  #Macrotis
ggplot(rats4%>%filter(Species=="macrotis"),aes(x=dist,y=z_weight,color=factor(Sex))) + 
  geom_quasirandom(width=0.3) +
  geom_abline (intercept=draws_Zwgt2$`alpha[2]`, #mac female intercept
                slope=draws_Zwgt2$`beta[2]`+
                draws_Zwgt2$`beta[3]`,  # slope=distance + macrotis:dist
               alpha=0.05,color="pink") +
  geom_abline (intercept=draws_Zwgt2$`alpha[2]`+ draws_Zwgt2$`beta[1]`+
                draws_Zwgt2$`beta[4]` , #mac male intercept
                slope=draws_Zwgt2$`beta[2]`+ draws_Zwgt2$`beta[3]`+
                draws_Zwgt2$`beta[5]`, # distance + male:dist
               alpha=0.05,color="lightblue") +
  scale_color_manual(values=c("red","blue"))+
  geom_abline(intercept=Zwgtsummary2$mean[Zwgtsummary2$variable=="alpha.mac"], #fem mac meanline
              slope=Zwgtsummary2$mean[Zwgtsummary2$variable=="dist"]+
                    Zwgtsummary2$mean[Zwgtsummary2$variable=="mac_dist"],
              linewidth=1,color="darkred") +
  geom_abline(intercept=Zwgtsummary2$mean[Zwgtsummary2$variable=="alpha.mac"]+ #male mac meanline
                    Zwgtsummary2$mean[Zwgtsummary2$variable=="mac_male"]+
                    Zwgtsummary2$mean[Zwgtsummary2$variable=="male"], 
              slope=Zwgtsummary2$mean[Zwgtsummary2$variable=="dist"]+
                    Zwgtsummary2$mean[Zwgtsummary2$variable=="mac_dist"]+
                    Zwgtsummary2$mean[Zwgtsummary2$variable=="male_dist"],
              linewidth=1,color="darkblue") +
  scale_x_continuous(breaks = c(0,1,2), 
                     labels = c("sympatry","near\nsympatry","allopatry") ) +
  scale_y_continuous(limits=c(-2.5,3),breaks=seq(-2,3,by=1)) +
  coord_cartesian(xlim = c(-0.25, 2.25)) +
  ylab ("Z-scaled weight") + 
  geom_hline(yintercept=0, linetype="dashed",linewidth=1.5) +
  guides(color = guide_legend(title = "Sex")) +
  theme_minimal() + 
  theme(  panel.grid.minor = element_line(size=0.5),
          plot.title = element_text(hjust = 0.5),
          panel.border = element_rect(color = "grey", fill = NA),
          axis.title.x = element_blank(),
          legend.title = element_text(size=12,hjust=0.5),
          legend.text = element_text(size=12),
          title = element_text(size=15),
          axis.text.x = element_text(size=13),
          axis.text.y = element_text(size=13),
          strip.background = element_rect(fill = "lightgrey", linetype = "solid"),
          strip.text = element_text(size=14,face="italic") ) + 
  facet_wrap(~Species)

# --- Z_Ear model ----
ZearMM<-rats4%>%reframe(z_ear,species=as.integer(Species),macrotis,male,dist,
                        macrotis_distance=macrotis*dist, macrotis_male=macrotis*male,
                        male_distance=male*dist)
ZearMM<-na.omit(ZearMM)
ZearMM$macrotis<-NULL

#list of stan inputs
Zear_df <- list(
  N = nrow(ZearMM), K = 5, X = ZearMM[,c(3:7)], y = ZearMM$z_ear, 
  n_sp = nlevels(factor(ZearMM$species)), species = ZearMM$species)

#stanmodel<-cmdstan_model("rats_stan.stan")

fit_Zear <- stanmodel$sample(
  data=Zear_df, chains = 6, parallel_chains = 6,  
  iter_warmup = 2000, iter_sampling = 4000 )

Zearsummary<-fit_Zear$summary(); Zearsummary

# ANOVA ---
aov_Zear_fus<-rats4%>%filter(Species=="fuscipes" & !is.na(z_ear))%>%
  mutate(dist=factor(dist,levels=c("0","1","2")),sex=factor(Sex))
summary(aov(z_ear ~ dist*Sex, data=aov_Zear_fus))
TukeyHSD(aov(z_ear ~ dist*sex, data=aov_Zear_fus))
###
aov_Zear_mac<-rats4%>%filter(Species=="macrotis" & !is.na(z_ear))%>%
  mutate(dist=factor(dist, levels=c("0","1","2")),sex=factor(Sex))
summary(aov(z_ear ~ dist*sex, data=aov_Zear_mac))
TukeyHSD(aov(z_ear ~ dist*sex, data=aov_Zear_mac))
# ---

Zearsummary2<-as.data.frame(Zearsummary)
Zearsummary2$variable<-c("lp","alpha.fus","alpha.mac","male","dist",
                         "mac_dist","mac_male","male_dist", "sigma")

draws_Zear<-fit_Zear$draws(format="df")

# slope credible interval ---
draws_Zear %>% reframe(fus_male_slope = `beta[2]` + `beta[5]` ,
                       fus_fem_slope = `beta[2]`,
                       mac_male_slope = `beta[2]` + `beta[3]` + `beta[5]`,
                       mac_fem_slope = `beta[2]`+ `beta[3]`) %>% 
  summarize(fus_male=mean(fus_male_slope),fus_male_upr=quantile(fus_male_slope,probs=0.975), 
            fus_male_lwr=quantile(fus_male_slope,probs=0.025), ##
            fus_fem=mean(fus_fem_slope),fus_fem_upr=quantile(fus_fem_slope,probs=0.975),
            fus_fem_lwr=quantile(fus_fem_slope,probs=0.025), ##
            mac_male=mean(mac_male_slope), mac_male_upr=quantile(mac_male_slope,probs=0.975), 
            mac_male_lwr=quantile(mac_male_slope,probs=0.025), ##
            mac_fem=mean(mac_fem_slope), mac_fem_upr=quantile(mac_fem_slope,probs=0.975), 
            mac_fem_lwr=quantile(mac_fem_slope,probs=0.025)) %>% as.data.frame()
# ---

draws_Zear2<-draws_Zear[1:3000,] #subset draws

# Fuscipes
ggplot(rats4%>%filter(Species=="fuscipes"),aes(x=dist,y=z_ear,color=factor(Sex))) + 
  geom_quasirandom(width=0.3) +
  geom_abline(intercept=draws_Zear2$`alpha[1]`, #fus female intercept
              slope=draws_Zear2$`beta[2]`, # distance beta
              alpha=0.05,color="pink") +
  geom_abline(intercept=draws_Zear2$`alpha[1]`+ draws_Zear2$`beta[1]`, #fus male intercept
              slope=draws_Zear2$`beta[2]`+ draws_Zear2$`beta[5]`, # distance + male:dist
              alpha=0.05,color="lightblue") +
  scale_color_manual(values=c("red","blue"))+
  geom_abline(intercept=Zearsummary2$mean[Zearsummary2$variable=="alpha.fus"], #fem fus meanline
              slope=Zearsummary2$mean[Zearsummary2$variable=="dist"],
              linewidth=1,color="darkred") +
  geom_abline(intercept=Zearsummary2$mean[Zearsummary2$variable=="alpha.fus"]+ #male fus meanline
                Zearsummary2$mean[Zearsummary2$variable=="male"], 
              slope=Zearsummary2$mean[Zearsummary2$variable=="dist"]+
                Zearsummary2$mean[Zearsummary2$variable=="male_dist"],
              linewidth=1,color="darkblue") +
  scale_x_continuous(breaks = c(0,1,2), 
                     labels = c("sympatry","near\nsympatry","allopatry") ) +
  scale_y_continuous(limits=c(-5.5,3),breaks=seq(-5,3,by=1)) +
  coord_cartesian(xlim = c(-0.25, 2.25)) +
  ylab ("Z-scaled ear length") + 
  geom_hline(yintercept=0, linetype="dashed",linewidth=1.5) +
  guides(color = guide_legend(title = "Sex")) +
  theme_minimal() + 
  theme(  panel.grid.minor = element_line(size=0.5),
          plot.title = element_text(hjust = 0.5),
          panel.border = element_rect(color = "grey", fill = NA),
          axis.title.x = element_blank(),
          legend.title = element_text(size=12,hjust=0.5),
          legend.text = element_text(size=12),
          title = element_text(size=15),
          axis.text.x = element_text(size=13),
          axis.text.y = element_text(size=13),
          strip.background = element_rect(fill = "lightgrey", linetype = "solid"),
          strip.text = element_text(size=14,face="italic") ) + 
  facet_wrap(~Species) 


#Macrotis
ggplot(rats4%>%filter(Species=="macrotis"),aes(x=dist,y=z_ear,color=factor(Sex))) + 
  geom_quasirandom(width=0.3) +
  geom_abline (intercept=draws_Zear2$`alpha[2]`, #mac female intercept
               slope=draws_Zear2$`beta[2]`+
                 draws_Zear2$`beta[3]`,  # slope=distance + macrotis:dist
               alpha=0.05,color="pink") +
  geom_abline (intercept=draws_Zear2$`alpha[2]`+ draws_Zear2$`beta[1]`+
                 draws_Zear2$`beta[4]` , #mac male intercept
               slope=draws_Zear2$`beta[2]`+ draws_Zear2$`beta[3]`+
                 draws_Zear2$`beta[5]`, # distance + male:dist
               alpha=0.05,color="lightblue") +
  scale_color_manual(values=c("red","blue"))+
  geom_abline(intercept=Zearsummary2$mean[Zearsummary2$variable=="alpha.mac"], #fem mac meanline
              slope=Zearsummary2$mean[Zearsummary2$variable=="dist"]+
                Zearsummary2$mean[Zearsummary2$variable=="mac_dist"],
              linewidth=1,color="darkred") +
  geom_abline(intercept=Zearsummary2$mean[Zearsummary2$variable=="alpha.mac"]+ #male mac meanline
                Zearsummary2$mean[Zearsummary2$variable=="mac_male"]+
                Zearsummary2$mean[Zearsummary2$variable=="male"], 
              slope=Zearsummary2$mean[Zearsummary2$variable=="dist"]+
                Zearsummary2$mean[Zearsummary2$variable=="mac_dist"]+
                Zearsummary2$mean[Zearsummary2$variable=="male_dist"],
              linewidth=1,color="darkblue") +
  scale_x_continuous(breaks = c(0,1,2), 
                     labels = c("sympatry","near\nsympatry","allopatry") ) +
  scale_y_continuous(limits=c(-5.5,3),breaks=seq(-5,3,by=1)) +
  coord_cartesian(xlim = c(-0.25, 2.25)) +
  ylab ("Z-scaled ear length") + 
  geom_hline(yintercept=0, linetype="dashed",linewidth=1.5) +
  guides(color = guide_legend(title = "Sex")) +
  theme_minimal() + 
  theme(  panel.grid.minor = element_line(size=0.5),
          plot.title = element_text(hjust = 0.5),
          panel.border = element_rect(color = "grey", fill = NA),
          axis.title.x = element_blank(),
          legend.title = element_text(size=12,hjust=0.5),
          legend.text = element_text(size=12),
          title = element_text(size=15),
          axis.text.x = element_text(size=13),
          axis.text.y = element_text(size=13),
          strip.background = element_rect(fill = "lightgrey", linetype = "solid"),
          strip.text = element_text(size=14,face="italic") ) + 
  facet_wrap(~Species) 

# --- Z_Foot model ----
ZfootMM<-rats4%>%reframe(z_foot,species=as.integer(Species),macrotis,male,dist,
                         macrotis_distance=macrotis*dist, macrotis_male=macrotis*male,
                         male_distance=male*dist)
ZfootMM<-na.omit(ZfootMM)
ZfootMM$macrotis<-NULL

#list of stan inputs
Zfoot_df <- list(
  N = nrow(ZfootMM), K = 5, X = ZfootMM[,c(3:7)], y = ZfootMM$z_foot, 
  n_sp = nlevels(factor(ZfootMM$species)), species = ZfootMM$species)

#stanmodel<-cmdstan_model("rats_stan.stan")

fit_Zfoot <- stanmodel$sample(
  data=Zfoot_df, chains = 6, parallel_chains = 6,  
  iter_warmup = 2000, iter_sampling = 4000 )

Zfootsummary<-fit_Zfoot$summary(); Zfootsummary

# ANOVA ---
aov_Zfoot_fus<-rats4%>%filter(Species=="fuscipes" & !is.na(z_foot))%>%
  mutate(dist=factor(dist,levels=c("0","1","2")),sex=factor(Sex))
summary(aov(z_foot ~ dist*sex, data=aov_Zfoot_fus))
TukeyHSD(aov(z_foot ~ dist*sex, data=aov_Zfoot_fus))
###
aov_Zfoot_mac<-rats4%>%filter(Species=="macrotis" & !is.na(z_foot))%>%
  mutate(dist=factor(dist, levels=c("0","1","2")),sex=factor(Sex))
summary(aov(z_foot ~ dist*sex, data=aov_Zfoot_mac))
TukeyHSD(aov(z_foot ~ dist*sex, data=aov_Zfoot_mac))
# ---

Zfootsummary2<-as.data.frame(Zfootsummary)
Zfootsummary2$variable<-c("lp","alpha.fus","alpha.mac","male","dist",
                          "mac_dist","mac_male","male_dist", "sigma")

draws_Zfoot<-fit_Zfoot$draws(format="df")

# slope credible interval ---
draws_Zfoot %>% reframe(fus_male_slope = `beta[2]` + `beta[5]` ,
                       fus_fem_slope = `beta[2]`,
                       mac_male_slope = `beta[2]` + `beta[3]` + `beta[5]`,
                       mac_fem_slope = `beta[2]`+ `beta[3]`) %>% 
  summarize(fus_male=mean(fus_male_slope),fus_male_upr=quantile(fus_male_slope,probs=0.975), 
            fus_male_lwr=quantile(fus_male_slope,probs=0.025), ##
            fus_fem=mean(fus_fem_slope),fus_fem_upr=quantile(fus_fem_slope,probs=0.975),
            fus_fem_lwr=quantile(fus_fem_slope,probs=0.025), ##
            mac_male=mean(mac_male_slope), mac_male_upr=quantile(mac_male_slope,probs=0.975), 
            mac_male_lwr=quantile(mac_male_slope,probs=0.025), ##
            mac_fem=mean(mac_fem_slope), mac_fem_upr=quantile(mac_fem_slope,probs=0.975), 
            mac_fem_lwr=quantile(mac_fem_slope,probs=0.025)) %>% as.data.frame()
# ---

draws_Zfoot2<-draws_Zfoot[1:3000,] #subset draws

  # fuscipes
ggplot(rats4%>%filter(Species=="fuscipes"),aes(x=dist,y=z_foot,color=factor(Sex))) + 
  geom_quasirandom(width=0.3) +
  geom_abline(intercept=draws_Zfoot2$`alpha[1]`, #fus female intercept
              slope=draws_Zfoot2$`beta[2]`, # distance beta
              alpha=0.05,color="pink") +
  geom_abline(intercept=draws_Zfoot2$`alpha[1]`+ draws_Zfoot2$`beta[1]`, #fus male intercept
              slope=draws_Zfoot2$`beta[2]`+ draws_Zfoot2$`beta[5]`, # distance + male:dist
              alpha=0.05,color="lightblue") +
  scale_color_manual(values=c("red","blue"))+
  geom_abline(intercept=Zfootsummary2$mean[Zfootsummary2$variable=="alpha.fus"], #fem fus meanline
              slope=Zfootsummary2$mean[Zfootsummary2$variable=="dist"],
              linewidth=1,color="darkred") +
  geom_abline(intercept=Zfootsummary2$mean[Zfootsummary2$variable=="alpha.fus"]+ #male fus meanline
                Zfootsummary2$mean[Zfootsummary2$variable=="male"], 
              slope=Zfootsummary2$mean[Zfootsummary2$variable=="dist"]+
                Zfootsummary2$mean[Zfootsummary2$variable=="male_dist"],
              linewidth=1,color="darkblue") +
  scale_x_continuous(breaks = c(0,1,2), 
                     labels = c("sympatry","near\nsympatry","allopatry") ) +
  scale_y_continuous(limits=c(-3.5,6),breaks=seq(-3,6,by=1)) +
  coord_cartesian(xlim = c(-0.25, 2.25)) +
  ylab ("Z-scaled hind foot length") + 
  geom_hline(yintercept=0, linetype="dashed",linewidth=1.5) +
  guides(color = guide_legend(title = "Sex")) +
  theme_minimal() + 
  theme(  panel.grid.minor = element_line(size=0.5),
          plot.title = element_text(hjust = 0.5),
          panel.border = element_rect(color = "grey", fill = NA),
          axis.title.x = element_blank(),
          legend.title = element_text(size=12,hjust=0.5),
          legend.text = element_text(size=12),
          title = element_text(size=15),
          axis.text.x = element_text(size=13),
          axis.text.y = element_text(size=13),
          strip.background = element_rect(fill = "lightgrey", linetype = "solid"),
          strip.text = element_text(size=14,face="italic") ) + 
  facet_wrap(~Species) 

  # macrotis
ggplot(rats4%>%filter(Species=="macrotis"),aes(x=dist,y=z_foot,color=factor(Sex))) + 
  geom_quasirandom(width=0.3) +
  geom_abline (intercept=draws_Zfoot2$`alpha[2]`, #mac female intercept
               slope=draws_Zfoot2$`beta[2]`+
                 draws_Zfoot2$`beta[3]`,  # slope=distance + macrotis:dist
               alpha=0.05,color="pink") +
  geom_abline (intercept=draws_Zfoot2$`alpha[2]`+ draws_Zfoot2$`beta[1]`+
                 draws_Zfoot2$`beta[4]` , #mac male intercept
               slope=draws_Zfoot2$`beta[2]`+ draws_Zfoot2$`beta[3]`+
                 draws_Zfoot2$`beta[5]`, # distance + male:dist
               alpha=0.05,color="lightblue") +
  scale_color_manual(values=c("red","blue"))+
  geom_abline(intercept=Zfootsummary2$mean[Zfootsummary2$variable=="alpha.mac"], #fem mac meanline
              slope=Zfootsummary2$mean[Zfootsummary2$variable=="dist"]+
                Zfootsummary2$mean[Zfootsummary2$variable=="mac_dist"],
              linewidth=1,color="darkred") +
  geom_abline(intercept=Zfootsummary2$mean[Zfootsummary2$variable=="alpha.mac"]+ #male mac meanline
                Zfootsummary2$mean[Zfootsummary2$variable=="mac_male"]+
                Zfootsummary2$mean[Zfootsummary2$variable=="male"], 
              slope=Zfootsummary2$mean[Zfootsummary2$variable=="dist"]+
                Zfootsummary2$mean[Zfootsummary2$variable=="mac_dist"]+
                Zfootsummary2$mean[Zfootsummary2$variable=="male_dist"],
              linewidth=1,color="darkblue") +
  scale_x_continuous(breaks = c(0,1,2), 
                     labels = c("sympatry","near\nsympatry","allopatry") ) +
  scale_y_continuous(limits=c(-3.5,6),breaks=seq(-3,6,by=1)) +
  coord_cartesian(xlim = c(-0.25, 2.25)) +
  ylab ("Z-scaled hind foot length") + 
  geom_hline(yintercept=0, linetype="dashed",linewidth=1.5) +
  guides(color = guide_legend(title = "Sex")) +
  theme_minimal() + 
  theme(  panel.grid.minor = element_line(size=0.5),
          plot.title = element_text(hjust = 0.5),
          panel.border = element_rect(color = "grey", fill = NA),
          axis.title.x = element_blank(),
          legend.title = element_text(size=12,hjust=0.5),
          legend.text = element_text(size=12),
          title = element_text(size=15),
          axis.text.x = element_text(size=13),
          axis.text.y = element_text(size=13),
          strip.background = element_rect(fill = "lightgrey", linetype = "solid"),
          strip.text = element_text(size=14,face="italic") ) + 
  facet_wrap(~Species)


# --- Z_Rostrum model ----
ZrostMM<-rats4%>%reframe(z_rost,species=as.integer(Species),macrotis,male,dist,
                         macrotis_distance=macrotis*dist, macrotis_male=macrotis*male,
                         male_distance=male*dist)
ZrostMM<-na.omit(ZrostMM)
ZrostMM$macrotis<-NULL

#list of stan inputs
Zrost_df <- list(
  N = nrow(ZrostMM), K = 5, X = ZrostMM[,c(3:7)], y = ZrostMM$z_rost, 
  n_sp = nlevels(factor(ZrostMM$species)), species = ZrostMM$species)

#stanmodel<-cmdstan_model("rats_stan.stan")

fit_Zrost <- stanmodel$sample(
  data=Zrost_df, chains = 6, parallel_chains = 6,  
  iter_warmup = 2000, iter_sampling = 4000 )

Zrostsummary<-fit_Zrost$summary(); Zrostsummary

# ANOVA ---
aov_Zrost_fus<-rats4%>%filter(Species=="fuscipes" & !is.na(z_rost))%>%
  mutate(dist=factor(dist,levels=c("0","1","2")),sex=factor(Sex))
summary(aov(z_rost ~ dist*sex, data=aov_Zrost_fus))
TukeyHSD(aov(z_rost ~ dist*sex, data=aov_Zrost_fus))
###
aov_Zrost_mac<-rats4%>%filter(Species=="macrotis" & !is.na(z_rost))%>%
  mutate(dist=factor(dist, levels=c("0","1","2")),sex=factor(Sex))
summary(aov(z_rost ~ dist*sex, data=aov_Zrost_mac))
TukeyHSD(aov(z_rost ~ dist*sex, data=aov_Zrost_mac))
# ---

Zrostsummary2<-as.data.frame(Zrostsummary)
Zrostsummary2$variable<-c("lp","alpha.fus","alpha.mac","male","dist",
                          "mac_dist","mac_male","male_dist", "sigma")

draws_Zrost<-fit_Zrost$draws(format="df")

# slope credible interval ---
draws_Zrost %>% reframe(fus_male_slope = `beta[2]` + `beta[5]` ,
                        fus_fem_slope = `beta[2]`,
                        mac_male_slope = `beta[2]` + `beta[3]` + `beta[5]`,
                        mac_fem_slope = `beta[2]`+ `beta[3]`) %>% 
  summarize(fus_male=mean(fus_male_slope),fus_male_upr=quantile(fus_male_slope,probs=0.975), 
            fus_male_lwr=quantile(fus_male_slope,probs=0.025), ##
            fus_fem=mean(fus_fem_slope),fus_fem_upr=quantile(fus_fem_slope,probs=0.975),
            fus_fem_lwr=quantile(fus_fem_slope,probs=0.025), ##
            mac_male=mean(mac_male_slope), mac_male_upr=quantile(mac_male_slope,probs=0.975), 
            mac_male_lwr=quantile(mac_male_slope,probs=0.025), ##
            mac_fem=mean(mac_fem_slope), mac_fem_upr=quantile(mac_fem_slope,probs=0.975), 
            mac_fem_lwr=quantile(mac_fem_slope,probs=0.025)) %>% as.data.frame()
# ---

draws_Zrost2<-draws_Zrost[1:3000,] #subset draws

# fuscipes
ggplot(rats4%>%filter(Species=="fuscipes"),aes(x=dist,y=z_rost,color=factor(Sex))) + 
  geom_quasirandom(width=0.3) +w
  geom_abline(intercept=draws_Zrost2$`alpha[1]`, #fus female intercept
              slope=draws_Zrost2$`beta[2]`, # distance beta
              alpha=0.05,color="pink") +
  geom_abline(intercept=draws_Zrost2$`alpha[1]`+ draws_Zrost2$`beta[1]`, #fus male intercept
              slope=draws_Zrost2$`beta[2]`+ draws_Zrost2$`beta[5]`, # distance + male:dist
              alpha=0.05,color="lightblue") +
  scale_color_manual(values=c("red","blue"))+
  geom_abline(intercept=Zrostsummary2$mean[Zrostsummary2$variable=="alpha.fus"], #fem fus meanline
              slope=Zrostsummary2$mean[Zrostsummary2$variable=="dist"],
              linewidth=1,color="darkred") +
  geom_abline(intercept=Zrostsummary2$mean[Zrostsummary2$variable=="alpha.fus"]+ #male fus meanline
                Zrostsummary2$mean[Zrostsummary2$variable=="male"], 
              slope=Zrostsummary2$mean[Zrostsummary2$variable=="dist"]+
                Zrostsummary2$mean[Zrostsummary2$variable=="male_dist"],
              linewidth=1,color="darkblue") +
  scale_x_continuous(breaks = c(0,1,2), 
                     labels = c("sympatry","near\nsympatry","allopatry") ) +
  scale_y_continuous(limits=c(-2,4),breaks=seq(-2,4,by=1)) +
  coord_cartesian(xlim = c(-0.25, 2.25)) +
  ylab ("Z-scaled rostrum width") + 
  geom_hline(yintercept=0, linetype="dashed",linewidth=1.5) +
  guides(color = guide_legend(title = "Sex")) +
  theme_minimal() + 
  theme(  panel.grid.minor = element_line(size=0.5),
          plot.title = element_text(hjust = 0.5),
          panel.border = element_rect(color = "grey", fill = NA),
          axis.title.x = element_blank(),
          legend.title = element_text(size=12,hjust=0.5),
          legend.text = element_text(size=12),
          title = element_text(size=15),
          axis.text.x = element_text(size=13),
          axis.text.y = element_text(size=13),
          strip.background = element_rect(fill = "lightgrey", linetype = "solid"),
          strip.text = element_text(size=14,face="italic") ) + 
  facet_wrap(~Species)

# macrotis
ggplot(rats4%>%filter(Species=="macrotis"),aes(x=dist,y=z_rost,color=factor(Sex))) + 
  geom_quasirandom(width=0.3) +
  geom_abline (intercept=draws_Zrost2$`alpha[2]`, #mac female intercept
               slope=draws_Zrost2$`beta[2]`+
                 draws_Zrost2$`beta[3]`,  # slope=distance + macrotis:dist
               alpha=0.05,color="pink") +
  geom_abline (intercept=draws_Zrost2$`alpha[2]`+ draws_Zrost2$`beta[1]`+
                 draws_Zrost2$`beta[4]` , #mac male intercept
               slope=draws_Zrost2$`beta[2]`+ draws_Zrost2$`beta[3]`+
                 draws_Zrost2$`beta[5]`, # distance + male:dist
               alpha=0.05,color="lightblue") +
  scale_color_manual(values=c("red","blue"))+
  geom_abline(intercept=Zrostsummary2$mean[Zrostsummary2$variable=="alpha.mac"], #fem mac meanline
              slope=Zrostsummary2$mean[Zrostsummary2$variable=="dist"]+
                Zrostsummary2$mean[Zrostsummary2$variable=="mac_dist"],
              linewidth=1,color="darkred") +
  geom_abline(intercept=Zrostsummary2$mean[Zrostsummary2$variable=="alpha.mac"]+ #male mac meanline
                Zrostsummary2$mean[Zrostsummary2$variable=="mac_male"]+
                Zrostsummary2$mean[Zrostsummary2$variable=="male"], 
              slope=Zrostsummary2$mean[Zrostsummary2$variable=="dist"]+
                Zrostsummary2$mean[Zrostsummary2$variable=="mac_dist"]+
                Zrostsummary2$mean[Zrostsummary2$variable=="male_dist"],
              linewidth=1,color="darkblue") +
  scale_x_continuous(breaks = c(0,1,2), 
                     labels = c("sympatry","near\nsympatry","allopatry") ) +
  scale_y_continuous(limits=c(-2,4),breaks=seq(-2,4,by=1)) +
  coord_cartesian(xlim = c(-0.25, 2.25)) +
  ylab ("Z-scaled rostrum width") + 
  geom_hline(yintercept=0, linetype="dashed",linewidth=1.5) +
  guides(color = guide_legend(title = "Sex")) +
  theme_minimal() + 
  theme(  panel.grid.minor = element_line(size=0.5),
          plot.title = element_text(hjust = 0.5),
          panel.border = element_rect(color = "grey", fill = NA),
          axis.title.x = element_blank(),
          legend.title = element_text(size=12,hjust=0.5),
          legend.text = element_text(size=12),
          title = element_text(size=15),
          axis.text.x = element_text(size=13),
          axis.text.y = element_text(size=13),
          strip.background = element_rect(fill = "lightgrey", linetype = "solid"),
          strip.text = element_text(size=14,face="italic") ) + 
  facet_wrap(~Species) 

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
  rats4$Ectoparasite.Load<-factor(rats4$Ectoparasite.Load,levels=c("None","Low","Medium","High"))
  rats4$dist<-factor(rats4$dist,levels=c("0","1","2"))

  parasite.summ<-rats4%>%filter(!is.na(Ectoparasite.Load))%>%
    group_by(Species,Ectoparasite.Load,dist,.drop=FALSE)%>%dplyr::summarize(count=n())
    
  parasite.summ<-rats4%>%filter(!is.na(Ectoparasite.Load))%>%
    group_by(Species,dist,.drop=FALSE)%>%dplyr::summarize(count=n())%>%
    right_join(.,parasite.summ,by=c("Species","dist"))%>%mutate(prop=count.y/count.x)
  
ggplot(parasite.summ,aes(x=dist,y=prop,fill=factor(Species))) + 
  geom_histogram(stat="identity",position="dodge") + facet_grid(~Ectoparasite.Load) +
  scale_fill_manual(values=c("darkorange1","steelblue2"))+
  scale_x_discrete(labels = c("sympatry","near\nsympatry","allopatry") ) +
  ylab ("Proportion") + ggtitle("Ectoparasite Load") +
  guides(fill = guide_legend(title = "Species")) +
  theme_minimal() + 
  theme(  plot.title = element_text(size=16,hjust = 0.5),
          panel.border = element_rect(color = "grey", fill = NA),
          axis.title.x = element_blank(),
          legend.title = element_text(size=12,hjust=0.5),
          legend.text = element_text(size=12),
          title = element_text(size=15),
          axis.text.x = element_text(size=13),
          axis.text.y = element_text(size=13),
          strip.background = element_rect(fill = "lightgrey", linetype = "solid"),
          strip.text = element_text(size=14) )

library(MASS)
rats4$Ectoparasite.Load<-factor(rats4$Ectoparasite.Load, levels=c("None","Low",
                                      "Medium","High"),ordered=T)

  model<-polr(Ectoparasite.Load ~ macrotis + male + dist + Weight..gram. +
                macrotis:male + macrotis:dist + male:dist, data=rats4, Hess=T)
  summary(model)
  confint(model)

exp(cbind(OR = coef(model), confint(model)))
  

## Comparison of counts ---

ecto_df<-rats4%>%filter(!is.na(Ectoparasite.Load))%>%
  mutate(Species=factor(Species),dist=factor(dist,levels=c("0","1","2")),
         ectoparasite=factor(Ectoparasite.Load,levels=c("None","Low","Medium","High")))%>%
          dplyr::group_by(Species,Ectoparasite.Load,dist,.drop=FALSE)%>%dplyr::summarize(count=n())
  
ecto_df<-ecto_df%>%group_by(Species,dist, .drop=FALSE)%>%dplyr::summarize(total=sum(count))  %>%
          right_join(.,ecto_df)%>%mutate(diff=total-count)

## DISTANCE 0
 # "NONE"
with(ecto_df %>% filter(dist==0 & Ectoparasite.Load=="None"), 
     pairwise_fisher_test(as.table(cbind(count,diff)) ) )     ## ns

 # "LOW"
with(ecto_df %>% filter(dist==0 & Ectoparasite.Load=="Low"),
     pairwise_fisher_test(as.table(cbind(count,diff)) ) )     ## **

 # "MEDIUM"
with(ecto_df %>% filter(dist==0 & Ectoparasite.Load=="Medium"),
     pairwise_fisher_test(as.table(cbind(count,diff)) ) )     ## **

 # "HIGH"
with(ecto_df %>% filter(dist==0 & Ectoparasite.Load=="High"),
     pairwise_fisher_test(as.table(cbind(count,diff)) ) )     ## ns

## DISTANCE 1
  # "NONE"
with(ecto_df %>% filter(dist==1 & Ectoparasite.Load=="None"), 
     pairwise_fisher_test(as.table(cbind(count,diff)) ) )     ## ns

 # "LOW"
with(ecto_df %>% filter(dist==1 & Ectoparasite.Load=="Low"),
     pairwise_fisher_test(as.table(cbind(count,diff)) ) )     ## ****

 # "MEDIUM"
with(ecto_df %>% filter(dist==1 & Ectoparasite.Load=="Medium"),
     pairwise_fisher_test(as.table(cbind(count,diff)) ) )     ## **

 # "HIGH"
with(ecto_df %>% filter(dist==1 & Ectoparasite.Load=="High"),
     pairwise_fisher_test(as.table(cbind(count,diff)) ) )     ## **


## DISTANCE 2

   # "NONE"
with(ecto_df %>% filter(dist==2 & Ectoparasite.Load=="None"), 
     pairwise_fisher_test(as.table(cbind(count,diff)) ) )     ## ns

 # "LOW"
with(ecto_df %>% filter(dist==2 & Ectoparasite.Load=="Low"),
     pairwise_fisher_test(as.table(cbind(count,diff)) ) )     ## ns

 # "MEDIUM"
with(ecto_df %>% filter(dist==2 & Ectoparasite.Load=="Medium"),
     pairwise_fisher_test(as.table(cbind(count,diff)) ) )     ## ns

 # "HIGH"
with(ecto_df %>% filter(dist==2 & Ectoparasite.Load=="High"),
     pairwise_fisher_test(as.table(cbind(count,diff)) ) )     ## ns


## Sex-species
parasite.summ_sex<-rats4%>%filter(!is.na(Ectoparasite.Load))%>%
  group_by(Species,Ectoparasite.Load,dist,Sex,.drop=FALSE)%>%dplyr::summarize(count=n())

parasite.summ_sex<-rats4%>%filter(!is.na(Ectoparasite.Load))%>%
  group_by(Species,dist,Sex,.drop=FALSE)%>%dplyr::summarize(count=n())%>%
  right_join(.,parasite.summ_sex,by=c("Species","dist","Sex"))%>%mutate(prop=count.y/count.x)


ggplot(parasite.summ_sex%>%filter(Species=="fuscipes"),aes(x=dist,y=prop,fill=factor(Sex))) + 
  geom_histogram(stat="identity",position="dodge") + facet_grid(~Ectoparasite.Load) +
  scale_fill_manual(values=c("red","blue"))+
  scale_x_discrete(labels = c("sympatry","near\nsympatry","allopatry") ) +
  ylab ("Proportion") + ggtitle("Ectoparasite Load - N.fuscipes") +
  guides(fill = guide_legend(title = "Sex")) +
  theme_minimal() + 
  theme(  plot.title = element_text(size=16,hjust = 0.5),
          panel.border = element_rect(color = "grey", fill = NA),
          axis.title.x = element_blank(),
          legend.title = element_text(size=12,hjust=0.5),
          legend.text = element_text(size=12),
          title = element_text(size=15),
          axis.text.x = element_text(size=13),
          axis.text.y = element_text(size=13),
          strip.background = element_rect(fill = "lightgrey", linetype = "solid"),
          strip.text = element_text(size=14) )

ggplot(parasite.summ_sex%>%filter(Species=="macrotis"),aes(x=dist,y=prop,fill=factor(Sex))) + 
  geom_histogram(stat="identity",position="dodge") + facet_grid(~Ectoparasite.Load) +
  scale_fill_manual(values=c("red","blue"))+
  scale_x_discrete(labels = c("sympatry","near\nsympatry","allopatry") ) +
  ylab ("Proportion") + ggtitle("Ectoparasite Load - N.macrotis") +
  guides(fill = guide_legend(title = "Sex")) +
  theme_minimal() + 
  theme(  plot.title = element_text(size=16,hjust = 0.5),
          panel.border = element_rect(color = "grey", fill = NA),
          axis.title.x = element_blank(),
          legend.title = element_text(size=12,hjust=0.5),
          legend.text = element_text(size=12),
          title = element_text(size=15),
          axis.text.x = element_text(size=13),
          axis.text.y = element_text(size=13),
          strip.background = element_rect(fill = "lightgrey", linetype = "solid"),
          strip.text = element_text(size=14) )


## DF
ecto_df_sex<-rats4%>%filter(!is.na(Ectoparasite.Load))%>%
  mutate(Species=factor(Species),dist=factor(dist,levels=c("0","1","2")), Sex=factor(Sex),
         ectoparasite=factor(Ectoparasite.Load,levels=c("None","Low","Medium","High")))%>%
  dplyr::group_by(Species,Ectoparasite.Load,Sex,dist,.drop=FALSE)%>%dplyr::summarize(count=n())

ecto_df_sex<-ecto_df_sex%>%group_by(Species,dist,Sex,.drop=FALSE)%>%dplyr::summarize(total=sum(count))%>%
  right_join(.,ecto_df_sex)%>%mutate(diff=total-count)


## fuscipes

## DISTANCE 0
# "NONE"
with(ecto_df_sex %>% filter(Species=="fuscipes" & dist==0 & Ectoparasite.Load=="None"), 
     pairwise_fisher_test(as.table(cbind(count,diff)) ) )     ## ns

# "LOW"
with(ecto_df_sex %>% filter(Species=="fuscipes" & dist==0 & Ectoparasite.Load=="Low"),
     pairwise_fisher_test(as.table(cbind(count,diff)) ) )     ## ns

# "MEDIUM"
with(ecto_df_sex %>% filter(Species=="fuscipes" & dist==0 & Ectoparasite.Load=="Medium"),
     pairwise_fisher_test(as.table(cbind(count,diff)) ) )     ## ns

# "HIGH"
with(ecto_df_sex %>% filter(Species=="fuscipes" & dist==0 & Ectoparasite.Load=="High"),
     pairwise_fisher_test(as.table(cbind(count,diff)) ) )     ## ns

## DISTANCE 1
# "NONE"
with(ecto_df_sex %>% filter(Species=="fuscipes" & dist==1 & Ectoparasite.Load=="None"), 
     pairwise_fisher_test(as.table(cbind(count,diff)) ) )     ## ns

# "LOW"
with(ecto_df_sex %>% filter(Species=="fuscipes" & dist==1 & Ectoparasite.Load=="Low"),
     pairwise_fisher_test(as.table(cbind(count,diff)) ) )     ## ns

# "MEDIUM"
with(ecto_df_sex %>% filter(Species=="fuscipes" & dist==1 & Ectoparasite.Load=="Medium"),
     pairwise_fisher_test(as.table(cbind(count,diff)) ) )     ## ns

# "HIGH"
with(ecto_df_sex %>% filter(Species=="fuscipes" & dist==1 & Ectoparasite.Load=="High"),
     pairwise_fisher_test(as.table(cbind(count,diff)) ) )     ## ns


## DISTANCE 2

# "NONE"
with(ecto_df_sex %>% filter(Species=="fuscipes" & dist==2 & Ectoparasite.Load=="None"), 
     pairwise_fisher_test(as.table(cbind(count,diff)) ) )     ## ns

# "LOW"
with(ecto_df_sex %>% filter(Species=="fuscipes" & dist==2 & Ectoparasite.Load=="Low"),
     pairwise_fisher_test(as.table(cbind(count,diff)) ) )     ## ns

# "MEDIUM"
with(ecto_df_sex %>% filter(Species=="fuscipes" & dist==2 & Ectoparasite.Load=="Medium"),
     pairwise_fisher_test(as.table(cbind(count,diff)) ) )     ## ns

# "HIGH"
with(ecto_df_sex %>% filter(Species=="fuscipes" & dist==2 & Ectoparasite.Load=="High"),
     pairwise_fisher_test(as.table(cbind(count,diff)) ) )     ## ns

## macrotis

## DISTANCE 0
# "NONE"
with(ecto_df_sex %>% filter(Species=="macrotis" & dist==0 & Ectoparasite.Load=="None"), 
     pairwise_fisher_test(as.table(cbind(count,diff)) ) )     ## ns

# "LOW"
with(ecto_df_sex %>% filter(Species=="macrotis" & dist==0 & Ectoparasite.Load=="Low"),
     pairwise_fisher_test(as.table(cbind(count,diff)) ) )     ## ns

# "MEDIUM"
with(ecto_df_sex %>% filter(Species=="macrotis" & dist==0 & Ectoparasite.Load=="Medium"),
     pairwise_fisher_test(as.table(cbind(count,diff)) ) )     ## ns

# "HIGH"
with(ecto_df_sex %>% filter(Species=="macrotis" & dist==0 & Ectoparasite.Load=="High"),
     pairwise_fisher_test(as.table(cbind(count,diff)) ) )     ## ns

## DISTANCE 1
# "NONE"
with(ecto_df_sex %>% filter(Species=="macrotis" & dist==1 & Ectoparasite.Load=="None"), 
     pairwise_fisher_test(as.table(cbind(count,diff)) ) )     ## ns

# "LOW"
with(ecto_df_sex %>% filter(Species=="macrotis" & dist==1 & Ectoparasite.Load=="Low"),
     pairwise_fisher_test(as.table(cbind(count,diff)) ) )     ## ns

# "MEDIUM"
with(ecto_df_sex %>% filter(Species=="macrotis" & dist==1 & Ectoparasite.Load=="Medium"),
     pairwise_fisher_test(as.table(cbind(count,diff)) ) )     ## ns

# "HIGH"
with(ecto_df_sex %>% filter(Species=="macrotis" & dist==1 & Ectoparasite.Load=="High"),
     pairwise_fisher_test(as.table(cbind(count,diff)) ) )     ## ns


## DISTANCE 2

# "NONE"
with(ecto_df_sex %>% filter(Species=="macrotis" & dist==2 & Ectoparasite.Load=="None"), 
     pairwise_fisher_test(as.table(cbind(count,diff)) ) )     ## ns

# "LOW"
with(ecto_df_sex %>% filter(Species=="macrotis" & dist==2 & Ectoparasite.Load=="Low"),
     pairwise_fisher_test(as.table(cbind(count,diff)) ) )     ## ns

# "MEDIUM"
with(ecto_df_sex %>% filter(Species=="macrotis" & dist==2 & Ectoparasite.Load=="Medium"),
     pairwise_fisher_test(as.table(cbind(count,diff)) ) )     ## ns

# "HIGH"
with(ecto_df_sex %>% filter(Species=="macrotis" & dist==2 & Ectoparasite.Load=="High"),
     pairwise_fisher_test(as.table(cbind(count,diff)) ) )     ## ns




  
