## ---------------------------  ##
## 
## Name: 03_CS-MARINEFF_comparison_photogrametry.R
##
## Purpose:
##
## Author: Bruno Serranito (bruno.serranito@mnhn.fr)
##
## Date: 2025-09-18
##
## ---------------------------  ##
rm(list=ls())



# libraries ---------------------------------------------------------------

library(dplyr)
library(tidyr)
library(readxl)
library(here)
library(vegan)
library(ggplot2)
library(ggrepel)
library(MetBrewer)
library(Hmisc)
library(glmmTMB)

# load files --------------------------------------------------------------

# source(here::here('code','03_0_photogrammetry_monitoring_preparation.R'))
# source(here::here('code','01_CS-MARINEFF_data_prep.R'))


CS_data=read_excel(here('data','modif','CS_MARINEFF_data.xlsx'), sheet=1)
PG_data=read_excel(here('data','modif','PG_MARINEFF_data.xlsx'), sheet=1)

# Formatting CS data ------------------------------------------------------



CS_data_format=CS_data %>%
  pivot_longer(cols=-c(Site, Date, Diver),
               names_to='Taxa',
               values_to='P') %>%
  mutate(Monit='CS')%>%
  group_by(Date,Site,Monit,Diver,Taxa)%>%
  drop_na() %>%
  reframe(P=sum(as.numeric(P)))  %>%
  mutate(P=ifelse(P>0,1,0))%>%
  # mutate(across(everything(.), ~replace_na(., 0)))%>%
  # mutate(Date=paste0(substring(Date,1,2),'/',substring(Date,3,4), '/',substring(Date,5,8)),
  #        Diver=substring(Diver,1,2),
  #        Site=substring(Site,1,2))%>%
  # mutate(Date=format(as.POSIXlt(Date,  format = "%d/%m/%Y"),   format="%Y-%m-%d", tz="GMT"))%>%
  ungroup() %>%
  rename('M_ensemble2'='Diver') %>% 
  mutate(Date2=lubridate::ymd(Date),
         Quat=lubridate::quarter(Date2,with_year = T))%>% 
  mutate(Site_quat=paste0(Site,'_',Quat)) %>% 
  group_by(Site,Quat,Monit,Site_quat,Taxa) %>% 
  reframe(P=mean(P)) %>% 
  # mutate(P=ifelse(P>0,1,0)) %>% 
  select(Site,Quat,Monit,Site_quat,Taxa,P) %>% 
  pivot_wider(names_from = Taxa,
              values_from=P) %>% 
  # pivot_longer(Taxa,P) %>% 
  mutate(
    across(everything(), ~replace_na(.x, 0))) %>% 
  
  # Filter only Bi, VB, Bu
  filter(Site %in% c('BI','BU','VB'))




CS_data_format_mat=CS_data_format %>% ungroup() %>%  select(-c(Site,Monit,Quat,Site_quat))



# Check CS - PG format and dates ------------------------------------------


Common_inter=intersect(PG_data$Site_quat,CS_data_format$Site_quat)



# PG data selection
PG_data_sel=PG_data %>% filter(Site_quat %in% Common_inter)

PG_data_sel_mat=PG_data_sel %>% ungroup() %>%  select(-c(Site,Monit,Quat,Site_quat))


# CS data selection
CS_data_format_sel=CS_data_format %>% filter(Site_quat %in% Common_inter)

CS_data_format_sel_mat=CS_data_format_sel %>% ungroup() %>%  select(-c(Site,Monit,Quat,Site_quat))





# Ordination CS data ------------------------------------------------------



CS_data_BC=vegdist(CS_data_format_sel_mat, method='bray')

pcoa_tot_CS_comp<- cmdscale(d = as.matrix(CS_data_BC), eig = T) # vegan package
# Eigen Values
Eig=pcoa_tot_CS_comp$eig[1:2]/sum(pcoa_tot_CS_comp$eig)
Eig=Eig*100




PCOA_SP_suivi_temps=data.frame(x=pcoa_tot_CS_comp$points[,1],y=pcoa_tot_CS_comp$points[,2],
                               Site=CS_data_format_sel$Site,
                               Quart=CS_data_format_sel$Quat)%>%
  group_by( Site,Quart) %>%
  mutate(Mx=mean(x), My=mean(y)) 

# Enfit
vf_SP <- envfit(pcoa_tot_CS_comp, CS_data_format_sel_mat, perm = 999)
spp.scrs_SP <- as.data.frame(scores(vf_SP, display = "vectors"))
spp.scrs2_SP <-(spp.scrs_SP[which(vf_SP$vectors$pvals<0.01),])

spp.scrs3_SP <- cbind(spp.scrs2_SP, Abb = rownames(spp.scrs2_SP))


ggplot(data=PCOA_SP_suivi_temps,aes(x,y))+
  geom_point(aes(col=Site), size=1, alpha=.6)+
  geom_point(data=PCOA_SP_suivi_temps,aes(x=Mx, y=My, col=Site), size=3, alpha=.6)+
  geom_segment(data=PCOA_SP_suivi_temps,aes(x=Mx, y=My, xend=x,yend=y,col=Site), alpha=.3)+
  geom_label(data=PCOA_SP_suivi_temps,aes(x=Mx, y=My,col=Site,
                                          label=Quart), alpha=.3)+
  
  #  geom_point(aes(Bar_x,Bar_y, fill=Strcut_groups),shape=22, size=5)+
  geom_point(data = spp.scrs3_SP, aes(Dim1/2,Dim2/2), col='black', size=3)+
  geom_label_repel(data = spp.scrs3_SP, aes(Dim1/2,Dim2/2,label=Abb))+
  geom_hline(yintercept=0)+
  geom_vline(xintercept=0)+
  xlab(sprintf("PCo1 (%.02f %%)",Eig[1]))+
  ylab(sprintf("PCo2 (%.02f %%)",Eig[2]))+
  theme_bw()+
  labs(title='SP PCoA')

#  Ordination PG data -----------------------------------------------------

PG_data_BC=vegdist(PG_data_sel_mat, method='bray')

pcoa_tot_PG<- cmdscale(d = as.matrix(PG_data_BC), eig = T) # vegan package
# Eigen Values
Eig=pcoa_tot_PG$eig[1:2]/sum(pcoa_tot_PG$eig)
Eig=Eig*100





# Procruste analysis ------------------------------------------------------

vare.proc <- procrustes(PG_data_BC, CS_data_BC)
vare.proc



summary(vare.proc)


g=protest(CS_data_BC, PG_data_BC, scores = "species", 
          permutations = how(nperm = 9999))





## Plot procruste ----------------------------------------------------------

Pred_X=predict(object=vare.proc, newdata=vare.proc$X, truemean = TRUE)
Pred_Y=predict(object=vare.proc,newdata=vare.proc$Yrot, truemean = TRUE)



# Predict trajectory
CS_pred_df=data.frame(X=Pred_X[,1],Y=Pred_X[,2],
                      PG_data_sel %>% ungroup %>%  select(Site,Monit,Quat))

PG_pred_df=data.frame(X=Pred_Y[,1],Y=Pred_Y[,2],
                      CS_data_format_sel %>% ungroup %>% 
                        select(Site,Monit,Quat))


CS_PG_df=rbind(CS_pred_df,PG_pred_df)

# 

ad_text=data.frame(X=c(-0.5,-0.5),Y=c(0.6,0.57), Site=c('BI','BI'),
                   text= c(paste0('Rp = ', round(g$scale,2),' '),
                           paste0('P = ', g$signif,' ***')))



CS_PG_df=CS_PG_df %>% 
  mutate(Monit2=ifelse(Monit== 'SS','PG',Monit))




Proc_p=ggplot(data=CS_PG_df)+
  geom_path(data=CS_PG_df,aes(X,Y, group=Monit2,col=Monit2),
            arrow = arrow(length = unit(0.50, "cm")), size=.8)+
  scale_shape(solid=T, 'Year-Quarter')+
  scale_size(range=c(.7,1.5))+
  facet_grid(.~Site)+
  geom_point(data=CS_PG_df %>% dplyr::select(-Site),
             aes(X,Y, fill=Monit2),shape=21, alpha=.4)+
  geom_point(data=CS_PG_df, aes(X,Y, fill=Monit2,
                                shape=as.factor(Quat)), size=3)+
  guides(size='none')+
  geom_vline(xintercept=0, linetype=2)+
  geom_hline(yintercept=0, linetype=2)+
  scale_color_manual(values=(met.brewer("Egypt", 2)))+
  scale_fill_manual(values=(met.brewer("Egypt", 2)))+
  labs(subtitle=c(paste0('Rp = ', round(g$scale,2),'\n',
                         'P = ', g$signif,' ***')))+

  theme_bw()+
  theme(plot.subtitle=element_text(size=11, face="bold", color="black"),
        legend.title=element_text(size=14))+
  xlab("Dimension 1")+
  ylab("Dimension 2")




## Export figure -----------------------------------------------------------

ggsave(Proc_p, file=here('outputs','Fig7_procruste_analysis.png'), dpi=350, width=10, height=6)


# Bootstrap -----------------------------------------------------


# Fish_SP= c("Sediment")

# 
# SP123_ok_temp=SP123 %>%
#   filter(!Taxa %in% Fish_SP) %>%
#   mutate(Monit='CS')%>%
#   group_by(Date,Site,Monit,Diver,Quadrat,Taxa)%>%
#   drop_na() %>%
#   summarize(P=sum(as.numeric(as.character(P))))%>%
#   mutate(P=ifelse(P>0,1,0))%>%
#   # mutate(across(everything(.), ~replace_na(., 0)))%>%
#   mutate(Date=paste0(substring(Date,1,2),'/',substring(Date,3,4), '/',substring(Date,5,8)),
#          Diver=substring(Diver,1,2),
#          Site=substring(Site,1,2))%>%
#   mutate(Date=format(as.POSIXlt(Date,  format = "%d/%m/%Y"),   format="%Y-%m-%d", tz="GMT"))%>%
#   ungroup() %>%
#   rename('M_ensemble2'='Diver') %>%
#   mutate(Date2=lubridate::ymd(Date),
#          Quat=lubridate::quarter(Date2,with_year = T))%>%
#   mutate(Site_quat=paste0(Site,'_',Quat)) %>%
#   ungroup() %>%
#   spread(Taxa,P) %>%
#   filter(Site %in% c('BI','BU','VB'))



CS_data_boot=CS_data %>% 
  mutate(Date2=lubridate::ymd(Date),
         Quat=lubridate::quarter(Date2,with_year = T),
         Site_quat=paste0(Site,'_',Quat)) %>% 
  select(-c(Date2,Date))


set.seed(123)
L=list()

for (i in 1: 4){
  LL=list()
  for (j in 1:50){
    # SP123_ok_rand=SP123_ok_temp %>% ungroup() %>% 
    #   mutate(Site_quat=paste0(Site,'_',Quat)) %>% 
    #   select(-c(Monit,M_ensemble2,Quadrat,Date,Date2)) %>% 
    
    CS_data_rand=CS_data_boot %>% 
      group_by(Site_quat, Site, Quat ) %>% sample_n(i) %>% 
      gather(Taxa,P, -c(Diver,Site_quat,Site,Quat)) %>% 
      
      group_by(Site_quat, Site ,Taxa,Quat) %>% 
      summarize(Mp=mean(P)) %>% 
      spread(Taxa,Mp) 
    
    
    # Common_inter=intersect(Data_PQ_mod$Site_quat,CS_data_format$Site_quat)
    
    
    SP123_ok_rand_sel=CS_data_rand %>% filter(Site_quat %in% Common_inter)
    
    
    SP123_ok_rand_sel_mat=SP123_ok_rand_sel %>% ungroup() %>%  select(-c(Site,Quat,Site_quat))
    
    
    
    SP_rand_BC=vegdist(SP123_ok_rand_sel_mat, method='bray')
    
    # 
    # vare.proc <- procrustes(SS_dist_BC, SP_rand_BC)
    # 
    # 
    # summary(vare.proc)
    
    g=protest(PG_data_BC, SP_rand_BC, scores = "species", 
              permutations = how(nperm = 9999))
    
    LL[[j]]=data.frame(R=round(g$scale,2), P=g$signif, rep=j, n=i)
    
  }
  print(i)
  LL_temp=do.call(rbind,LL)
  L[[i]] <-LL_temp
}



DF=do.call(rbind,L) %>% 
  mutate(pval=ifelse(P<0.05,"Signif.","N.S"))

DF_mean = DF %>%  group_by(n) %>% 
  summarize(M=mean(R),
            sd=sd(R))

DF_pval = DF %>%  group_by(n) %>% 
  mutate(tot=n()) %>% 
  ungroup() %>% 
  group_by(n,pval) %>% 
  summarize(per=n()/tot,
            nb=n()) %>% distinct()



Boot_proc_p=ggplot(DF, aes(n,R, group=pval,col=pval))+
  geom_jitter(width = 0.2, size=2, alpha=.5)+
  # geom_point(position=position_jitterdodge(dodge.width = 0.8),
  #             size=2, alpha=.5)+
  theme_bw()+
  theme(axis.text=element_text(size=12),
        text = element_text(size = 14),
        strip.text.x = element_text(size = 20))+
  scale_color_viridis_d(direction=-1)+
  guides(color="none",shape="none")+
  scale_color_manual(values=met.brewer("Egypt", 2))+
  stat_summary(fun='mean',geom='line',aes(group=1), linetype=2)+
  stat_summary(fun.data=mean_sdl,geom='pointrange', aes(group=n))+
  stat_summary(fun='mean',geom='point', aes(group=n), size=3, shape=21, fill='black')+
  xlab('n surveys')


Boot_prop_p=ggplot(DF_pval, aes(n,per, fill=pval))+
  geom_bar(stat='identity', col='black')+
  theme_bw()+
  xlab('n surveys')+
  ylab('Proportion of Bootstrapped samples')+
  theme(axis.text=element_text(size=12),
        text = element_text(size = 12),
        strip.text.x = element_text(size = 20))+
  guides(fill='none')+
  scale_fill_manual(values=met.brewer("Egypt", 2))



Boot_prop_double_p=ggarrange(Boot_proc_p,Boot_prop_p,ncol=2,
                             labels=c('A','B'))


DF_mean=DF %>% 
  group_by(n) %>% 
  dplyr::summarize(R_mean=mean(R))



setwd('~/R/MARINEFF/Science-Participative/Fig_papier')

ggsave(Boot_prop_double_p, file='Boot_prop_double_p.png', dpi=350)



## Model increasing--------------------------


mod1=glmmTMB(R  ~ n, 
             family = gaussian(), data = DF)


summary(mod1)




Resids <- residuals(mod1, type = "pearson")

p1<-ggplot(data.frame(resid = Resids), aes(sample = resid)) +
  stat_qq(alpha = 0.6) +
  stat_qq_line(color = "red") +
  xlab('Theorical quantiles')+
  ylab('Sample quantiles')+
  # labs(title = "QQ-plot") +
  theme_bw()+
  theme(axis.text=element_text(size=14),
        axis.title=element_text(size=14),
        legend.text = element_text(size=13),
        legend.title= element_text(size=14),
        text = element_text(size = 13),
        strip.text.x = element_text(size = 20))


#fitted versus residuals
df_res<- data.frame(log_fitted = fitted(mod1),
                             residuals  = resid(mod1, type = "deviance"))


p2 <- ggplot(df_res, aes(x = log_fitted, y = residuals))+
  geom_point() +
  geom_smooth(method='lm')+
  geom_hline(yintercept=0, linetype=2)+
  labs(x = "Linear predictor", y = "Deviance residual")+
  theme_bw() +
  theme(axis.text=element_text(size=14),
        axis.title=element_text(size=14),
        legend.text = element_text(size=13),
        legend.title= element_text(size=14),
        text = element_text(size = 13),
        strip.text.x = element_text(size = 20))

ggarrange(p1,p2, ncol=2) %>% 
  annotate_figure(top = text_grob(paste0(name), size = 16))
