## Figures for main immuno paper
library(dplyr)
library(ggplot2)
library(XLConnect)
library(tidyr)
require(gridExtra)
library(scales) # Need the scales package for transformed axes

## make sure you are working from the root directory as paths are relative


## read in data and make sure columns are of correct type
dat <- read.csv("data/juba_poc_ocv_immuno_dat.csv",as.is=T)
dat <- dat %>% mutate(date_of_v1=as.Date(date_of_v1),
                      date_of_v2=as.Date(date_of_v2),
                      date_of_v3=as.Date(date_of_v3)
                      )
dat$titer <- as.numeric(dat$titer)
dat <- dat %>% filter(!is.na(titer))

## add age groups
dat <- dat %>% mutate(age_gr=cut(age_yrs,breaks = c(0,6,18,100),right=FALSE,labels = c("young children","older children","adults")))

## check flow chart stats
## particpants with baseline blood sample by age group
dat %>% group_by(serum_sample_id_v1) %>% filter(day=="v1") %>% summarize(ag=age_gr[1]) %>% select(ag) %>% table
dat %>% group_by(serum_sample_id_v1) %>% filter(day=="v2") %>% summarize(ag=age_gr[1]) %>% select(ag) %>% table
dat %>% group_by(serum_sample_id_v1) %>% filter(day=="v3") %>% summarize(ag=age_gr[1]) %>% select(ag) %>% table

test = dat %>% group_by(serum_sample_id_v1) %>% summarize(age=age_yrs[1],age_gr=age_gr[1]) %>% arrange(age)
## some basic stats
dat %>% group_by(serum_sample_id_v1) %>% summarise(adult=ifelse(age_yrs[1]>=18,1,0),sex=sex[1]) %>% summarize(sum(adult==1),sum(adult==1 & sex==1),sum(adult==1 & sex==2))
dat %>% group_by(serum_sample_id_v1) %>% summarise(adult=ifelse(age_yrs[1]<18,1,0),sex=sex[1]) %>% summarize(sum(adult==1),sum(adult==1 & sex==1),sum(adult==1 & sex==2))


## inter visit timing
dat <- dat %>%
  mutate(
      time_d0d1=date_of_v2-date_of_v1,
      time_d0d3=date_of_v3-date_of_v1)


#GMTs by visit, assay and serotype
gmts_no_age <- dat %>% group_by(serotype,day,assay) %>% summarize(gmt=exp(mean(log(titer),na.rm=T)))

gmts <- dat %>% group_by(serotype,day,assay,age_gr) %>% summarize(gmt=exp(mean(log(titer),na.rm=T)))
gmts %>% data.frame


## max fold rises
dat_wide <- dat %>% filter(assay=="vibriocidals") %>% spread(day,titer)
fold_rises <- dat_wide %>% group_by(serum_sample_id_v1,serotype)%>%
  summarize(age=age_yrs[1],age_gr=age_gr[1],baseline=v1,fr=max(c(v2/v1,v3/v1,0),na.rm=T)) %>% filter(fr !=0)

fold_rises %>% group_by(serotype,age_gr) %>% summarise(minfr=min(fr),maxfr=max(fr), mfr=mean(fr))
fold_rises %>% filter(baseline<80) %>% group_by(serotype) %>% summarise(mfr=mean(fr))

## Figure 1 - Baseline Titer By Age
cols = brewer.pal(8,"Set1")

dat %>% filter(serotype=="inaba" & day=="v1" & assay=="vibriocidals")%>%
  ggplot(aes(x=age_yrs,y=titer)) + geom_smooth(col=cols[2]) + geom_point(color=AddAlpha(cols[2],.3)) +
  scale_y_continuous(trans=log2_trans()) +
  geom_hline(aes(yintercept=80),lty=2,col=AddAlpha('black',.5)) +
  xlab("age (years)") + ylab("vibriocidal titer (inaba)")

dat %>% filter(serotype=="ogawa" & day=="v1" & assay=="vibriocidals") %>%
  ggplot(aes(x=age_yrs,y=titer)) + geom_smooth(col=cols[2]) + geom_point(color=AddAlpha(cols[2],.3)) +
  scale_y_continuous(trans=log2_trans()) +
  geom_hline(aes(yintercept=80),lty=2,col=AddAlpha('black',.5)) +
  xlab("age (years)") + ylab("vibriocidal titer (ogawa)")

dat %>% filter(day=="v1" & assay=="vibriocidals") %>%
  ggplot(aes(x=age_yrs,y=titer)) + geom_smooth(col=cols[2]) + geom_point(color=AddAlpha(cols[2],.3)) +
  scale_y_continuous(trans=log2_trans()) +
  geom_hline(aes(yintercept=80),lty=2,col=AddAlpha('black',.5)) + #xlim(c(0,50)) +
  xlab("age (years)") + ylab("vibriocidal titer")  + facet_grid(serotype~.)

dat %>% filter(day=="v1" & assay=="IgM") %>%
  ggplot(aes(x=age_yrs,y=titer)) + geom_smooth(col=cols[2]) + geom_point(color=AddAlpha(cols[2],.3)) +
  scale_y_continuous(trans=log2_trans()) +
  xlab("age (years)") + ylab("OSP IgM titer")  + facet_grid(.~serotype)

dat %>% filter(day=="v1" & assay=="IgA") %>%
  ggplot(aes(x=age_yrs,y=titer)) + geom_smooth(col=cols[2]) + geom_point(color=AddAlpha(cols[2],.3)) +
  scale_y_continuous(trans=log2_trans()) +
  xlab("age (years)") + ylab("OSP IgA titer")  + facet_grid(.~serotype)

dat %>% filter(day=="v1" & assay=="IgG") %>%
  ggplot(aes(x=age_yrs,y=titer)) + geom_smooth(col=cols[2]) + geom_point(color=AddAlpha(cols[2],.3)) +
  scale_y_continuous(trans=log2_trans()) +
  xlab("age (years)") + ylab("OSP IgA titer")  + facet_grid(.~serotype)

## Figure 3
dat %>% filter(assay=='vibriocidals') %>%
  ggplot(aes(x=day,y=titer,fill=factor(vaccine_in_2014),color=factor(vaccine_in_2014),group=vaccine_in_2014)) +
  geom_jitter(position=position_jitterdodge(dodge.width=0.9),alpha=.5) +
  geom_boxplot(fill=NA,outlier.colour = NA, colour=NA,
               position = position_dodge(width=0.9)) +
                 stat_summary(fun.y = function(x) exp(mean(log(x))), aes(color=factor(vaccine_in_2014)), size = 5,  geom = "point",pch=3,
                              position = position_dodge(width=0.9)) +
                                facet_grid(serotype~age_gr) +
                                scale_y_continuous(trans=log2_trans())

## Sup Figure for
dat %>% filter(assay=='IgA') %>%
  ggplot(aes(x=day,y=titer,fill=factor(vaccine_in_2014),color=factor(vaccine_in_2014),group=vaccine_in_2014)) +
  geom_jitter(position=position_jitterdodge(dodge.width=0.9),alpha=.5) +
  geom_boxplot(fill=NA,outlier.colour = NA, colour=NA,
               position = position_dodge(width=0.9)) +
                 stat_summary(fun.y = function(x) exp(mean(log(x))), aes(color=factor(vaccine_in_2014)), size = 5,  geom = "point",pch=3,
                              position = position_dodge(width=0.9)) +
                                facet_grid(serotype~age_gr) +
                                scale_y_continuous(trans=log2_trans())

dat %>% filter(assay=='IgM') %>%
  ggplot(aes(x=day,y=titer,fill=factor(vaccine_in_2014),color=factor(vaccine_in_2014),group=vaccine_in_2014)) +
  geom_jitter(position=position_jitterdodge(dodge.width=0.9),alpha=.5) +
  geom_boxplot(fill=NA,outlier.colour = NA, colour=NA,
               position = position_dodge(width=0.9)) +
                 stat_summary(fun.y = function(x) exp(mean(log(x))), aes(color=factor(vaccine_in_2014)), size = 5,  geom = "point",pch=3,
                              position = position_dodge(width=0.9)) +
                                facet_grid(serotype~age_gr) +
                                scale_y_continuous(trans=log2_trans())

dat %>% filter(assay=='IgG') %>%
  ggplot(aes(x=day,y=titer,fill=factor(vaccine_in_2014),color=factor(vaccine_in_2014),group=vaccine_in_2014)) +
  geom_jitter(position=position_jitterdodge(dodge.width=0.9),alpha=.5) +
  geom_boxplot(fill=NA,outlier.colour = NA, colour=NA,
               position = position_dodge(width=0.9)) +
                 stat_summary(fun.y = function(x) exp(mean(log(x))), aes(color=factor(vaccine_in_2014)), size = 5,  geom = "point",pch=3,
                              position = position_dodge(width=0.9)) +
                                facet_grid(serotype~age_gr) +
                                scale_y_continuous(trans=log2_trans())

## by prior OCV
dat %>% filter(assay=='vibriocidals') %>%
  ggplot(aes(x=day,y=titer)) +
  geom_jitter(aes(col=day),height=0.02,alpha=.9) +
  stat_summary(fun.y = function(x) exp(mean(log(x))), colour = "black", size = 4,  geom = "point",pch=3) +
  facet_grid(serotype~age_gr) +
  scale_y_continuous(trans=log2_trans())


## now for some age standardization
standardize_age <- function(tmp_dat,my_breaks=c(0,10,15,20,25,30,35,101)){
  ## bring in age dist data from IOM
  require(survey)
  age_poc <- read.csv("data/agedist.csv",as.is=T)

  ## now let's get the totals for 5 year age catagories
  age_poc$age_class <- cut(age_poc$Age,breaks=my_breaks)
  age_poc <- age_poc %>% filter(Age!=0) %>% group_by(Age) %>% mutate(total=sum(F + M,na.rm=T))
  std_age_dist <- age_poc %>%
    group_by(age_class) %>%
    summarize(pct=sum(total)/sum(age_poc$total))

  tmp_dat$age_class_std <- cut(tmp_dat$age_yrs,breaks=my_breaks)

  ## now let's look at the age-distribution of our sample
  samp_age_dist <- tmp_dat %>%
    group_by(age_class_std) %>%
    summarize(pct_samp=n()/nrow(tmp_dat))

  comp_age_dists <- cbind(samp_age_dist,pct_pop=std_age_dist$pct)
  comp_age_dists <- comp_age_dists %>% mutate(w=pct_pop/pct_samp)

  mean_dat <- tmp_dat[,c('serum_sample_id_v1','age_class_std','titer')]
  mean_dat <- merge(mean_dat,comp_age_dists,by.x='age_class_std',by.y='age_class_std',all.x=T)

  ## now get the post-stratified seroprevalence
  mean_dat$titer <- log(mean_dat$titer)
  my_svy_design <-
    svydesign(
      id = ~1 ,
      data = subset(mean_dat,!is.na(titer)),
      weights=~w
    )

  rc <- c(svymean(~titer, my_svy_design,na.rm=T),confint(svymean(~titer, my_svy_design,na.rm=T), df = degf(my_svy_design)))
  return(rc)
}

exp(standardize_age(subset(dat,serotype=='inaba' & assay == 'vibriocidals' & day == 'v1')))
exp(standardize_age(subset(dat,serotype=='ogawa' & assay == 'vibriocidals' & day == 'v1')))

## baseline titer by fold-rise of vibriocidals
dat_wide <- dat %>% filter(assay=="vibriocidals") %>% spread(day,titer)
fold_rises <- dat_wide %>% group_by(serum_sample_id_v1,serotype)%>%
  summarize(age=age_yrs[1],baseline=v1,type="from 1 dose",fr=v2/v1)
fold_rises2 <- dat_wide %>% group_by(serum_sample_id_v1,serotype)%>%
  summarize(age=age_yrs[1],baseline=v1,type="from 2 doses",fr=v3/v1)
fold_rises <- rbind(fold_rises,fold_rises2)

fold_rises %>%
  ggplot() +
  geom_point(aes(x=baseline,y=fr,color=age),position=position_jitter(w = 0.1, h = 0.1)) +
  scale_x_log10(breaks=c(5,10,20,40,80,160,320,640,1280)) + scale_y_log10(breaks=c(2^(0:9)),limits=c(.2,512)) +
  geom_hline(yintercept=4,col="orange",lty=2,alpha=.5) + geom_vline(xintercept=80,col="orange",lty=2,alpha=.5) +
  xlab("baseline titer") + ylab("fold-rise in vibriocidal titer") + facet_grid(serotype~ type)

fold_rises %>%
  ggplot() +
  geom_point(aes(x=baseline,y=fr,color=age),position=position_jitter(w = 0.1, h = 0.1)) +
  scale_x_log10(breaks=c(5,10,20,40,80,160,320,640,1280)) + scale_y_log10(breaks=c(2^(0:9)),limits=c(.2,512)) +
  geom_hline(yintercept=4,col="orange",lty=2,alpha=.5) + geom_vline(xintercept=80,col="orange",lty=2,alpha=.5) +
  xlab("baseline titer") + ylab("fold-rise in vibriocidal titer") +
  facet_grid(serotype~ type) +
  geom_smooth(aes(x=baseline,y=fr),method = "lm", formula = y ~ splines::ns(x,2), se = FALSE)

fr_pred <- mutate(fold_rises,sc=(fr>=4))
glm(sc ~ baseline + I(age<25),data=fr_pred,family="binomial") %>% summary

## now the same for the OSP responses
dat_wide <- dat %>% filter(assay=="IgA") %>% spread(day,titer)
fold_rises <- dat_wide %>% group_by(serum_sample_id_v1,serotype)%>%
  summarize(age=age_yrs[1],baseline=v1,type="from 1 dose",fr=v2/v1)
fold_rises2 <- dat_wide %>% group_by(serum_sample_id_v1,serotype)%>%
  summarize(age=age_yrs[1],baseline=v1,type="from 2 doses",fr=v3/v1)
fold_rises <- rbind(fold_rises,fold_rises2)

fold_rises %>%
  ggplot() +
  geom_point(aes(x=baseline,y=fr,color=age),position=position_jitter(w = 0.05, h = 0.05)) +
  scale_x_log10() + scale_y_log10(breaks=c(2^(0:9)),limits=c(.2,32)) +
  geom_hline(yintercept=2,col="orange",lty=2,alpha=.5) +
  xlab("baseline titer") + ylab("fold-rise in OSP IgA titer") + facet_grid(serotype~ type)

## now for IgM
dat_wide <- dat %>% filter(assay=="IgM") %>% spread(day,titer)
fold_rises <- dat_wide %>% group_by(serum_sample_id_v1,serotype)%>%
  summarize(age=age_yrs[1],baseline=v1,type="from 1 dose",fr=v2/v1)
fold_rises2 <- dat_wide %>% group_by(serum_sample_id_v1,serotype)%>%
  summarize(age=age_yrs[1],baseline=v1,type="from 2 doses",fr=v3/v1)
fold_rises <- rbind(fold_rises,fold_rises2)

fold_rises %>%
  ggplot() +
  geom_point(aes(x=baseline,y=fr,color=age),position=position_jitter(w = 0.05, h = 0.05)) +
  scale_x_log10() + scale_y_log10(breaks=c(2^(0:9)),limits=c(.2,32)) +
  geom_hline(yintercept=2,col="orange",lty=2,alpha=.5) +
  xlab("baseline titer") + ylab("fold-rise in OSP IgM titer") + facet_grid(serotype~ type)

## and now for IgG
dat_wide <- dat %>% filter(assay=="IgG") %>% spread(day,titer)
fold_rises <- dat_wide %>% group_by(serum_sample_id_v1,serotype)%>%
  summarize(age=age_yrs[1],baseline=v1,type="from 1 dose",fr=v2/v1)
fold_rises2 <- dat_wide %>% group_by(serum_sample_id_v1,serotype)%>%
  summarize(age=age_yrs[1],baseline=v1,type="from 2 doses",fr=v3/v1)
fold_rises <- rbind(fold_rises,fold_rises2)

fold_rises %>%
  ggplot() +
  geom_point(aes(x=baseline,y=fr,color=age),position=position_jitter(w = 0.05, h = 0.05)) +
  scale_x_log10() + scale_y_log10(breaks=c(2^(0:9)),limits=c(.2,32)) +
  geom_hline(yintercept=2,col="orange",lty=2,alpha=.5) +
  xlab("baseline titer") + ylab("fold-rise in OSP IgG titer") + facet_grid(serotype~ type)

### Looking at ratio of IgG to IgM and IgG to IgA at baseline by age
dat %>% filter(day=="v1") %>% group_by(serum_sample_id_v1) %>% summarize()

## titer decay over time
dat <- dat %>% mutate(t_since=ifelse(day=="v1",0,
                          ifelse(day=="v2",date_of_v2-date_of_v1 %>% as.numeric,
                                 date_of_v3-date_of_v2 %>% as.numeric)))

## now plot vibriocidals by days
dat %>% ggplot(aes(x=t_since,y=log2(titer))) + geom_jitter(col='darkgrey') + geom_smooth() +  facet_wrap(~serotype) -> gg
to.pdf(print(gg),filename = "vibrio_decay.pdf",width=8,height=6)
