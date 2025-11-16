## ---------------------------  ##
## 
## Name: 02_CS_MARINEFF_conformity.R
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
library(readxl)
library(writexl)
library(tidyr)
library(here)
library(irr)
library(MetBrewer)
library(stringr)
library(ggpubr)
library(vegan)
library(glue)
library(gmodels)
library(zoo)
library(lubridate)
library(DHARMa)
library(emmeans)
library(glmmTMB)
library(broom.mixed)
library(colorspace)
library(ggh4x)
library(irrCAC)
library(writexl)

# load data conformity ----------------------------------------------------
# load correspondance table for taxa labels consistency
corres_table=read_excel(here('data','raw','Morpho_corres_tot_2.xlsx'), sheet=2)

Taxa_name=read_excel(here('data','raw','Morpho_corres_english.xlsx'), sheet=2) %>% 
  mutate(sp_name_ok=ifelse(is.na(Latin), Sp_English, Latin),
         It=ifelse(is.na(Latin), F, T))

data_amb <- read_excel(here("BDD_2021_SP_VD.xlsx"), sheet=1) %>% filter(Quadrat=='RECIF') %>%
  gather(spe,val,-c(Date,Site,Diver,Quadrat,Code_replicat)) %>%
  filter(!spe %in% c('Eboulis_de_roche', 'Pierres_et_cailloutis','Sediment','Plateau','Tombant','Surplombs_Cavites')) %>% 
  mutate(val=as.numeric(replace(val,val %in% c('PONTE','ponte'),'1')))%>%
  left_join(corres_table, c('spe'='SP_Niv1')) %>% 
  mutate(val=as.numeric(val),
         Site=case_when(Site =='BI__'~ 'BI_',
                        Site =='VB__'~ 'VB_',
                        TRUE~ Site)) %>% 
  dplyr::select(Date,Site,Diver,Quadrat,Code_replicat,SP_commun,val)%>%
  group_by(Date,Site,Diver,Quadrat,Code_replicat,SP_commun) %>% 
  summarize(val=sum(val)) %>%  mutate(val=ifelse(val>0,1,0)) %>% 
  rename('spe'='SP_commun') %>% 
  # mutate(P=ifelse(is.na(P),0,1)) %>%
  spread(spe,val) %>% 
  dplyr::select(-'Hydraires_autres') %>% 
  filter(!Code_replicat =='CD_BI_16052021_RECIF') %>% 
  filter(!Site =='FE_') %>% 
  mutate(Site=str_sub(Site, end=-2))


# Other Prof scient. data
PS_divers=c('VD_','PC_','ZD_')


Data_PS_supp=data_amb %>% 
  filter(Diver %in% PS_divers)

# Removing PS quadrats
data_amb<- data_amb %>% 
  filter(!Diver %in% PS_divers)



data_val<- read_excel(here("BDD_2021_SP_VD.xlsx"), sheet=2) %>% filter(Quadrat=='RECIF') %>% na.omit() %>%
  pivot_longer(cols=-c(Date,Site,Diver,Quadrat,Code_replicat),
               names_to='spe',
               values_to='val') %>% 
  left_join(corres_table, c('spe'='SP_Niv1')) %>%
  mutate(val=as.numeric(val),
         Site=case_when(Site =='BI__'~ 'BI_',
                        Site =='VB__'~ 'VB_',
                        TRUE~ Site)) %>%   
  dplyr::select(Date,Site,Diver,Quadrat,Code_replicat,SP_commun,val)%>%
  group_by(Date,Site,Diver,Quadrat,Code_replicat,SP_commun) %>%
  mutate(val=as.numeric(val))%>%summarize(val=sum(val)) %>%  mutate(val=ifelse(val>0,1,0)) %>%
  rename('spe'='SP_commun') %>%
  pivot_wider(names_from=spe,
              values_from=val) %>% 
  filter(!Site =='FE_')%>% 
  mutate(Site=str_sub(Site, end=-2)) %>% 

  dplyr::select(-'Hydraires_autres') %>% 
  bind_rows(Data_PS_supp) 
  




# check for difference in column names --> ok
setdiff(colnames(data_amb),colnames(data_val))
setdiff(colnames(data_val),colnames(data_amb))

Common_col=intersect(colnames(data_amb),colnames(data_val))


# Merge 
Data_2021=rbind(data_amb,data_val)


data_amb2=data_amb[,Common_col]

# Table data available -----------------------------------------------------

data_amb_ok<-data_amb %>% ungroup() %>%   
  select(Date,Site,Diver) %>%  mutate(Comp='nPS')

data_val_ok<-data_val %>%ungroup() %>%   select(Date,Site,Diver)%>%  mutate(Comp='PS')

Table_divers=rbind(data_amb_ok,
                    data_val_ok) %>% 
  mutate(Date=as.Date(paste0(substr(Date, 5, 8),'-',substr(Date, 3, 4),'-', substr(Date, 1, 2)), format='%Y-%m-%d')) %>% 

  group_by(Comp) %>% 
  mutate(Diver_num = paste0('#', dense_rank(Diver))) %>%
    ungroup() %>% 
  mutate(Comp_Diver = paste(Comp, Diver_num, sep = "_")) %>%  
  mutate(value = 1) %>%                                   
  distinct(Date, Site, Comp_Diver, .keep_all = TRUE) %>%  # avoid duplicates
 select(-c(Diver,Diver_num,Comp)) %>% 
  pivot_wider(
    names_from = Comp_Diver,
    values_from = value,
    values_fill = 0)   %>%
  arrange(Date)



# write_xlsx(Table_divers, path=here('outputs','Table_divers.xlsx'))


#  comparison functions ---------------------------------------------------



## Percentage agreement computation
cross_table_PS=function(x) {
  if (sum(x$value_ref) == nrow(x) & sum(x$value_cand) == nrow(x)){
    GG= data.frame(x=c('0','1','0','1'), y=c('0','0','1','1'),
                   Freq=c(0,0,0,nrow(x)),
                   Comp=c('TA','FP','FA','TP'))  
    
  }else if (sum(x$value_ref) == 0 & sum(x$value_cand) == 0){
    GG= data.frame(x=c('0','1','0','1'), y=c('0','0','1','1'),
                   Freq=c(nrow(x),0,0,0),
                   Comp=c('TA','FP','FA','TP'))  
  }else if (sum(x$value_ref) == nrow(x) & sum(x$value_cand) == 0){
    GG= data.frame(x=c('0','1','0','1'), y=c('0','0','1','1'),
                   Freq=c(0,nrow(x),0,0),
                   Comp=c('TA','FP','FA','TP'))  
  }else if (sum(x$value_ref) == 0 & sum(x$value_cand) == nrow(x)){
    GG= data.frame(x=c('0','1','0','1'), y=c('0','0','1','1'),
                   Freq=c(0,0,nrow(x),0),
                   Comp=c('TA','FP','FA','TP'))  
  }else{
    g=CrossTable(x$value_ref  ,x$value_cand)
    GG=as.data.frame(g$t) %>% 
      mutate(Comp=case_when(x == 0 & y==0~'TA',
                            x==1 & y==0 ~'FP',
                            x==0 & y==1 ~'FA',
                            x==1 & y==1 ~'TP',
                            T ~'NA'))
  }
  # GG=GG %>% 
  #   mutate(diff_d)
  return(GG)
}


# Kappa computation

Kappa_agree_PS=function(x) {
  
  
  # Remove mirrored rows
  x <- x %>%
    group_by(Taxa) %>%
    slice(1) %>%
    ungroup()
  
  # 
  PS_PS=x %>%  dplyr::select(value_ref,value_cand)
  Agree=agree(PS_PS, tolerance=0)
  Kap=kappa2(PS_PS)
  
  ag_kap_df=data.frame(AG=Agree$value, Kappa=Kap$value, pval_kap=Kap$p.value)
  
  return(ag_kap_df)
}


# Gwet's AC1 computation
Gwet_AC1_agreement=function(x) {
  
  
  # Remove mirrored rows
  x <- x %>%
    group_by(Taxa) %>%
    slice(1) %>%
    ungroup()
  
  # 
  m=x %>%  dplyr::select(value_ref,value_cand)
  Agree=agree(m, tolerance=0)
  
  gwet=gwet.ac1.raw(data.frame(m$value_ref, m$value_cand))
  
  
  gwet_sel_df=gwet$est %>% dplyr::select(coeff.val, conf.int, p.value) %>% 
    rename('AC1'='coeff.val',
           'conf.int_gwet'='conf.int',
           'p.value_gwet'='p.value') %>% 
    bind_cols(Agree=Agree$value)
  
  return(gwet_sel_df)
}





fun_posthoc_comp<-function(x){
  
  
  em_values<-emmeans(x, pairwise ~ Comp | Site, cov.reduce = list(diff_d = mean,
                                                                  month=mean))
  em_contrast<-as.data.frame(em_values$contrasts) 
    
    return(em_contrast)
}


fun_posthoc_means_comp<-function(x){
  
  
  em_values<-emmeans(x, pairwise ~  Site, cov.reduce = list(diff_d = mean,
                                                                  month=mean))
  em_means<-as.data.frame(em_values$emmeans) 
  
  return(em_means)
}


# 1.1 Comparison CS & PS----------------------------------------------------------------
##  Comparison table -------------------------------------------------

data_amb_date=data_amb %>% ungroup %>% 
  mutate(Date=as.Date(paste0(substr(data_amb2$Date, 5, 8),'-',substr(data_amb2$Date, 3, 4),'-', substr(data_amb2$Date, 1, 2)), format='%Y-%m-%d'))

data_val_date=data_val %>% ungroup %>% 
  mutate(Date=as.Date(paste0(substr(data_val$Date, 5,8),'-',substr(data_val$Date, 3, 4),'-', substr(data_val$Date, 1, 2)), format='%Y-%m-%d'))




data_PS=data_val_date %>% 
  dplyr::select(Date,Site,Diver,Quadrat,Code_replicat)



data_CS=data_amb_date %>% 
  dplyr::select(Date,Site,Diver,Quadrat,Code_replicat)






#
DD_CS_PS_1 <- data_CS %>%
  # self-join to compare CS and PS
  inner_join(data_PS, by = character(), suffix = c("", "_2")) %>%
  filter(Site==Site_2) %>% 
  mutate(diff_d = as.numeric(abs(Date-Date_2))) %>%
  ungroup() %>% 
  # Remove comparison with dif days higher than 30
  filter(diff_d <=30) 

DD_CS_PS_2 <- data_PS %>%
  # self-join to compare CS and PS
  inner_join(data_CS, by = character(), suffix = c("", "_2")) %>%
  filter(Site==Site_2) %>% 
  mutate(diff_d = as.numeric(abs(Date-Date_2))) %>%
  ungroup() %>% 
  # Remove comparison with dif days higher than 30
  filter(diff_d <=30) 


DD_CS_PS=rbind(DD_CS_PS_1,DD_CS_PS_2)%>% 
  mutate(idx=seq(1,nrow(.),1)) %>% 
  # idx for Kappa
  mutate(
    fact_A = pmin(Code_replicat, Code_replicat_2),
    fact_B = pmax(Code_replicat, Code_replicat_2)) %>%
  group_by(fact_A, fact_B) %>% 
  mutate(idx2=cur_group_id())


levels(factor(DD_CS_PS$Code_replicat))




### PS-nPS formatting --------------------------------------------------------


data_both_date=rbind(data_amb_date,data_val_date)

data_amb_comb=data_both_date %>%
  pivot_longer(cols=-c(Date,Site,Diver,Quadrat,Code_replicat),
               names_to="Taxa",
               values_to="value_ref") %>% 
  left_join(DD_CS_PS, by=c('Date',
                           'Site',
                           'Diver',
                           'Quadrat',
                           'Code_replicat'))%>%
  distinct()


data_PS_comb=data_both_date %>% 
  rename('Date_2'='Date',
         'Site_2'='Site',
         'Quadrat_2'='Quadrat',
         'Code_replicat_2'='Code_replicat',
         'Diver_2'='Diver') %>%
  pivot_longer(cols=-c(Date_2,Site_2,Diver_2,Quadrat_2,Code_replicat_2),
               names_to="Taxa",
               values_to="value_cand") %>% 
  left_join(DD_CS_PS, by=c('Date_2',
                           'Site_2',
                           'Quadrat_2',
                           'Code_replicat_2',
                           'Diver_2'))%>%
  distinct() 


Data_merge=data_amb_comb %>% 
  left_join(data_PS_comb, by=c("Date" ,"Site",             
                               "Diver", "Quadrat" ,        
                               "Code_replicat","Taxa",
                               "Date_2","Site_2", 'Diver_2',
                               'Quadrat_2', 'Code_replicat_2', 'diff_d','idx','idx2',
                               'fact_A','fact_B')) %>% 
  drop_na()

#count number of raw without the Taxa and Taxa_2 column
test=Data_merge %>% 
  dplyr::select(-c(Taxa,Taxa, value_ref,value_cand)) %>% 
  distinct() 
  

DD2=as_tibble(DD_CS_PS) 

#check different rows between DD_CS_PS and test 
setdiff(DD2,test)
setdiff(test,DD2)




all=data_val %>%  ungroup() %>%  select(Code_replicat) %>%  distinct()
cons=DD2 %>% ungroup() %>% select(Code_replicat_2)%>%  distinct() %>%  
  rename('Code_replicat' = 'Code_replicat_2')

not_cons_inter=setdiff(all,cons)

### nPS-PS apply function ---------------------------------------------------


# # Split dataset by Diver
split_data_inter <- split(Data_merge, Data_merge$idx)
split_data_inter_idx2 <- split(Data_merge, Data_merge$idx2)




AC1_agree_table_inter <- lapply(split_data_inter_idx2, Gwet_AC1_agreement)
AC1_inter_df <- bind_rows(AC1_agree_table_inter, .id = "idx2") %>% 
  # mutate(idx=as.numeric(idx)) %>% 
  # left_join(DD_CS_PS, by='idx') %>% 
  mutate(idx2=as.numeric(idx2))






all_table_inter<-AC1_inter_df %>% 
mutate(idx2=as.numeric(idx2)) %>% 
  left_join(DD_CS_PS, by='idx2') %>% 
  mutate(Agree=Agree/100)



comparison_inter_df2=all_table_inter %>% 
  dplyr::select(Site,Agree, AC1) %>% 
  pivot_longer(cols=c(Agree,AC1),
               names_to='metric',
               values_to='values')


 
# plot test
ggplot(comparison_inter_df2, aes(x=Site, y=values, color=metric))+
  stat_summary(fun.data=mean_sdl, fun.args = list(mult=1), geom='pointrange', position=position_dodge(width=0.5), size=1)+
  stat_summary(fun='mean', geom='point', position=position_dodge(width=0.5), size=4)+
  theme_bw()+
  ylab('Metrics')+
  xlab('Site')+
  theme(axis.text=element_text(size=12),
        axis.title=element_text(size=15),
        legend.text = element_text(size=13),
        legend.title= element_text(size=14),
        text = element_text(size = 10),
        strip.text.x = element_text(size = 20))


# 1.2. Intra CS survey comparison ---------------------------------------------


data_nPS=data_amb_date %>% 
  dplyr::select(Date,Site,Diver,Quadrat,Code_replicat)




# Comparison table 
DD_intra_CS <- data_nPS %>%
  # self-join to compare each row with all others
  inner_join(data_nPS, by = character(), suffix = c("", "_other")) %>%
  filter(Code_replicat != Code_replicat_other) %>%
  filter(Diver != Diver_other) %>% 
  filter(Site==Site_other) %>% 
  mutate(diff_d = as.numeric(abs(Date-Date_other))) %>%
  group_by(Code_replicat, Diver_other) %>%
  ungroup() %>% 
  filter(diff_d <=30) %>% 
  # idx for Sens Spe
  mutate(idx=seq(1,nrow(.),1)) %>% 
  # idx for Kappa
  mutate(
    fact_A = pmin(Code_replicat, Code_replicat_other),
    fact_B = pmax(Code_replicat, Code_replicat_other)) %>%
  group_by(fact_A, fact_B) %>% 
  mutate(idx2=cur_group_id())





## nPS formatting -----------------------------------------------------


# Merge nPS1 & nPS2 for comparison
data_nPS_comb2_1=data_amb_date %>% 
  rename('Date'='Date',
         'Code_replicat'='Code_replicat',
         'Diver'='Diver') %>%
  pivot_longer(cols=-c(Date,Site,Diver,Quadrat,Code_replicat),
               names_to="Taxa",
               values_to="value_ref") %>% 
  left_join(DD_intra_CS, by=c('Date',
                              'Site',
                              'Diver',
                              'Quadrat',
                              'Code_replicat'))%>%
  distinct() %>% drop_na()


data_nPS_comb2_2=data_amb_date %>% 
  rename('Date_other'='Date',
         'Code_replicat_other'='Code_replicat',
         'PS_other'='Diver',
         'Quadrat_other'='Quadrat',
         'Site_other'='Site') %>%
  pivot_longer(cols=-c(Date_other,Site_other,PS_other,Quadrat_other,Code_replicat_other),
               names_to="Taxa",
               values_to="value_cand") %>% 

  distinct() 



Data_nPS_merge=data_nPS_comb2_1 %>% 
  left_join(data_nPS_comb2_2, by=c('Date_other',
                                  'Site_other',
                                  'Diver_other'='PS_other',
                                  'Quadrat_other', 
                                  'Code_replicat_other',
                                  'Taxa'))


## nPS  Apply function ---------------------------------------------------------



# # Split dataset by Diver (idx2 = removed mirrored rows for sens and spec)
split_data_nPS_idx <- split(Data_nPS_merge, Data_nPS_merge$idx)
split_data_nPS_idx2 <- split(Data_nPS_merge, Data_nPS_merge$idx2)


# Apply function with lapply
crosstables_nPS <- lapply(split_data_nPS_idx, cross_table_PS)

AC1_agree_table_nPS <- lapply(split_data_nPS_idx2, Gwet_AC1_agreement)





AC1_nPS_df <- bind_rows(AC1_agree_table_nPS, .id = "idx2") %>% 
  mutate(idx2=as.numeric(idx2)) %>% 
  left_join(DD_intra_CS, by='idx2')%>% 
  mutate(Agree=Agree/100)






all_table_nPS<-AC1_nPS_df %>% 
  mutate(idx2=as.numeric(idx2))



comparison_nPS_df2=all_table_nPS %>% 
  dplyr::select(Site,Agree, AC1) %>% 
  pivot_longer(cols=c(Agree, AC1),
               names_to='metric',
               values_to='values') 
  

# plot test
ggplot(comparison_nPS_df2, aes(x=Site, y=values, color=metric))+
  stat_summary(fun.data=mean_sdl, fun.args = list(mult=1), geom='pointrange', position=position_dodge(width=0.5), size=1)+
  stat_summary(fun='mean', geom='point', position=position_dodge(width=0.5), size=4)+
  theme_bw()+
  ylab('Metrics')+
  xlab('Site')+
  theme(axis.text=element_text(size=12),
        axis.title=element_text(size=15),
        legend.text = element_text(size=13),
        legend.title= element_text(size=14),
        text = element_text(size = 10),
        strip.text.x = element_text(size = 20))


# 1.3. Intra PS survey comparison ----------------------------------------------




data_PS=data_val_date %>% 
  dplyr::select(Date,Site,Diver,Quadrat,Code_replicat)




# Comparison table 
DD_intra_SP <- data_PS %>%
  # self-join to compare each row with all others
  inner_join(data_PS, by = character(), suffix = c("", "_other")) %>%
  filter(Code_replicat != Code_replicat_other) %>%
  filter(Diver != Diver_other) %>% 
  filter(Site==Site_other) %>% 
mutate(diff_d = as.numeric(abs(Date-Date_other))) %>%
  group_by(Code_replicat, Diver_other) %>%
  ungroup() %>% 
  filter(diff_d <=30) %>% 
  # idx for Sens Spe
  mutate(idx=seq(1,nrow(.),1)) %>% 
  # idx for Kappa
  mutate(
    fact_A = pmin(Code_replicat, Code_replicat_other),
    fact_B = pmax(Code_replicat, Code_replicat_other)) %>%
  group_by(fact_A, fact_B) %>% 
  mutate(idx2=cur_group_id())



all_PS=data_val %>%  ungroup() %>%  select(Code_replicat) %>%  distinct()
cons_PS=DD_intra_SP %>% ungroup() %>% select(Code_replicat)%>%  distinct() 


# Check if PS survey were considered in intra group comparison
not_cons_PS=setdiff(all_PS,cons_PS)

PS_not_cons=unique(rbind(not_cons_inter,not_cons_PS))

setwd(here('data'))
saveRDS(PS_not_cons, file='PS_not_considered.R')


## PS data formatting ---------------------
# Merge PS1 & PS2 for comparison
data_PS_comb2_1=data_val_date %>% 
  rename(
         'Code_replicat_PS'='Code_replicat') %>%
  pivot_longer(cols=-c(Date,Site,Diver,Quadrat,Code_replicat_PS),
               names_to="Taxa",
               values_to="value_ref") %>% 
  left_join(DD_intra_SP, by=c('Date',
                     'Site',
                     'Diver',
                     'Quadrat',
                     'Code_replicat_PS'='Code_replicat'))%>%
  distinct() %>% drop_na()


data_PS_comb2_2=data_val_date %>% 
  rename('Date_other'='Date',
         'Code_replicat_other'='Code_replicat',
         'PS_other'='Diver',
         'Quadrat_other'='Quadrat',
         'Site_other'='Site') %>%
  pivot_longer(cols=-c(Date_other,Site_other,PS_other,Quadrat_other,Code_replicat_other),
               names_to="Taxa",
               values_to="value_cand") %>% 
  # left_join(DD_intra_SP, by=c('Date_other',
  #                             'Site_other',
  #                             'PS_other'='Diver_other',
  #                             'Quadrat_other',
  #                             'Code_replicat_other'))%>%
  distinct() 



Data_PS_merge=data_PS_comb2_1 %>% 
  left_join(data_PS_comb2_2, by=c('Date_other',
                                  'Site_other',
                                  'Diver_other'='PS_other',
                                  'Quadrat_other', 
                                  'Code_replicat_other',
                                  'Taxa'))




## PS apply function ---------------------


# # Split dataset by Diver
split_data_PS_idx <- split(Data_PS_merge, Data_PS_merge$idx)
split_data_PS_idx2 <- split(Data_PS_merge, Data_PS_merge$idx2)


# Apply function with lapply
# crosstables_PS <- lapply(split_data_PS_idx, cross_table_PS)
# kappa_agree_table_PS <- lapply(split_data_PS_idx2, Kappa_agree_PS)
AC1_agree_table_PS <- lapply(split_data_PS_idx2, Gwet_AC1_agreement)



AC1_PS_df <- bind_rows(AC1_agree_table_PS, .id = "idx2") %>%
  mutate(idx2=as.numeric(idx2)) %>% 
  left_join(DD_intra_SP, by='idx2')%>% 
  mutate(Agree=Agree/100)



all_table_PS<-AC1_PS_df

comparison_PS_df2=all_table_PS %>% 
  dplyr::select(Site,Agree, AC1) %>% 
  pivot_longer(cols=c(Agree, AC1),
               names_to='metric',
               values_to='values')
 


Acc_kapp_PS_p=ggplot(comparison_PS_df2, aes(metric,values, color=Site))+
  stat_summary(fun.data= 'mean_sdl',fun.args = list(mult = 1), geom='pointrange', size=2,lwd=2, aes(group=metric))+
  facet_grid(.~Site)+
  geom_jitter(alpha=.3, size=3)+
 
  scale_color_manual(values=rev(met.brewer("Lakota", 3)))+
  theme_bw()+
  guides(color='none')+
  ylab('Accuracy/ Kappa values')+
  theme_bw()+
  theme(axis.text=element_text(size=14),
        axis.title=element_text(size=15),
        legend.text = element_text(size=13),
        legend.title= element_text(size=14),
        text = element_text(size = 13),
        strip.text.x = element_text(size = 20))




# 2.Grouping of CS_PS & PS_PS -----------------------------------

## 2.1. Models  -------------------------------------------------------------------------
# sel_var=c('Accuracy','Date', 'Diver','AC1', 'Sens', 'Spe', 'Site', 'diff_d',
# 'idx2')

sel_var=c('Agree','Date', 'Diver','AC1', 'Site', 'diff_d',
          'idx2')


all_comp=all_table_inter %>% 
  dplyr::select(sel_var) %>% 
  mutate(Comp= 'nPS-PS') %>% 
  
  bind_rows(all_table_nPS %>%
  dplyr::select(sel_var) %>%
  mutate(Comp= 'nPS-nPS')) %>%
  
  bind_rows(all_table_PS %>% 
              dplyr::select(sel_var)%>% 
              mutate(Comp= 'PS-PS')) %>% 
  mutate(month=lubridate::month(Date))




all_symm_comp=all_comp %>% 
  group_by(idx2, Comp) %>% 
  slice(1) 



## Compute models
responses_symm<-c("Agree", "AC1")


fits_symm <- lapply(responses_symm, function(responses_symm) {
  
  formula <- as.formula(
    paste(responses_symm, " ~ diff_d + month*Site  + Comp:Site + (1| Diver)")
  )
  glmmTMB(formula, family = gaussian(), data = all_symm_comp)
})



names(fits_symm)<- responses_symm



# Test
# data_test=all_symm_comp 
# 
# mod1=glmmTMB(AC1  ~ diff_d + month*Site + month:Comp + Comp:Site + (1| Diver), 
#         family = gaussian(), data = all_symm_comp)
# 
# mod2=glmmTMB(AC1 ~ diff_d + Site + Comp + Comp:Site + (1| Diver) + (1| month),
#              family = gaussian(), data = all_symm_comp)
# 
# 
# mod3=glmmTMB(AC1  ~ diff_d + month*Site  + Comp:Site + (1| Diver), 
#              family = gaussian(), data = all_symm_comp)
# 
# 
# mod4=glmmTMB(AC1  ~ diff_d + month*Site  + Comp*Site + (1| Diver), 
#              family = gaussian(), data = all_symm_comp)
# 
# 
# AIC(mod1,mod2,mod3,mod4)



emtrends_mod1 <- emtrends(mod1, ~ Comp, var = "month")

summary(emtrends_mod1, infer = c(TRUE, TRUE))
# pairs(emtrends_mod1)











fits<-c(fits_symm)



## Extract summary


summ_mods<-lapply(fits, function(x){
  
  broom.mixed::tidy(x)
  
})

summary_mod <- bind_rows(summ_mods, .id = "Model") %>% 
  mutate(p.value=as.numeric(p.value),
      p.value = case_when(
    p.value < 0.001 ~ "***",
    p.value < 0.01  ~ "**",
    p.value < 0.05  ~ "*",
    TRUE ~ as.character(p.value)))




## Extract paiwise comparison using emmeans
list_constrast<-lapply(fits,fun_posthoc_comp)
contrast_mod <- bind_rows(list_constrast, .id = "Model") %>% 
  mutate(p.value=as.numeric(p.value),
         p.value = case_when(
           p.value < 0.001 ~ "***",
           p.value < 0.01  ~ "**",
           p.value < 0.05  ~ "*",
           TRUE ~ as.character(p.value)))


contrast_mod



# extraction of estimated marginal mean
list_marg_means<-lapply(fits,fun_posthoc_means_comp)
contrast_marg_means <- bind_rows(list_marg_means, .id = "Model")



write_xlsx(contrast_mod, path=here('outputs','constrast_mod_summary.xlsx'))
write_xlsx(summary_mod, path=here('outputs','summary_mod_all_metrics.xlsx'))



## Check normality and 

Diag_plots <- mapply(function(mod, name) {  
  
  
  ## QQPlot
Resids <- residuals(mod, type = "pearson")

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
df_mod_age_size<- data.frame(log_fitted = fitted(mod),
                             residuals  = resid(mod, type = "deviance"))


p2 <- ggplot(df_mod_age_size, aes(x = log_fitted, y = residuals))+
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



}, fits, names(fits))

all_diagplots<-ggarrange(plotlist = Diag_plots, nrow = length(responses_symm))

ggsave(all_diagplots, file=here('outputs','Supp_info_diag_plots.png'),
       dpi=350,  width=10, height=16)


## 2.2. Plot Acc/ Kappa/ AC1 -------------------------------------------------

all_comp_long=all_comp %>% 
  pivot_longer(cols=-c(Site,Date,Diver,diff_d,Comp,month),
               names_to = 'metric',
               values_to='values') 
  
  
  
  
all_comp_acc_kappa=all_comp_long %>% 
  filter(metric %in% c('Agree','AC1')) %>%
  mutate(contrast_2=case_when(Comp=='nPS-PS' ~'Inter-groups',
                              Comp=='nPS-nPS' ~'Intra-nPS',
                              Comp=='PS-PS' ~'Intra-PS'))
  
  

all_comp_acc_kappa$metric <- factor(all_comp_acc_kappa$metric, levels = c("Agree", "AC1"))
all_comp_acc_kappa$contrast_2 <- factor(all_comp_acc_kappa$contrast_2, levels = c('Intra-nPS',
                                                                                  'Inter-groups',
                                                                                  'Intra-PS'))



# Extract significant differences for the plot
signif_table_acc_ac1=contrast_mod %>% 
  mutate(metric=case_when(Model =='Agree' ~'Agreement Percentage',
                   Model =='AC1' ~'AC1',
                   TRUE~ Model),
  
  xmin=case_when(metric== 'Agreement Percentage' & contrast== '(nPS-nPS) - (nPS-PS)' ~ 0.8,
                 metric== 'Agreement Percentage'& contrast== '(nPS-PS) - (PS-PS)' ~ 1,
                 metric== 'Agreement Percentage'& contrast== '(nPS-nPS) - (PS-PS)' ~ 0.8,
                 metric== 'AC1' & contrast== '(nPS-nPS) - (nPS-PS)' ~ 1.8,
                 metric== 'AC1'& contrast== '(nPS-PS) - (PS-PS)' ~ 2,
                 metric== 'AC1'& contrast== '(nPS-nPS) - (PS-PS)' ~ 1.8),
  
  xmax=case_when(metric== 'Agreement Percentage' & contrast== '(nPS-nPS) - (nPS-PS)' ~ 1,
                 metric== 'Agreement Percentage'& contrast== '(nPS-PS) - (PS-PS)' ~ 1.2,
                 metric== 'Agreement Percentage'& contrast== '(nPS-nPS) - (PS-PS)' ~ 1.2,
                 metric== 'AC1' & contrast== '(nPS-nPS) - (nPS-PS)' ~ 2,
                 metric== 'AC1'& contrast== '(nPS-PS) - (PS-PS)' ~ 2.2,
                 metric== 'AC1'& contrast== '(nPS-nPS) - (PS-PS)' ~ 2.2),
  ypos=0.99) %>% 
  filter(metric %in% c('Agreement Percentage','AC1')) %>% 
  filter(p.value<=0.05) 






Acc_agree_allcomp_p=ggplot(all_comp_acc_kappa, aes(metric,values, color=contrast_2))+
  stat_summary(fun.data= 'mean_sdl',fun.args = list(mult = 1), 
               geom='pointrange', size=1.5,lwd=1.5, aes(group=contrast_2),
               position=position_dodge(0.8))+
  facet_grid(.~Site)+
  geom_jitter(alpha=.1, size=2.5,  position=position_jitterdodge(jitter.width = 0.5,
                                                               dodge.width = 0.8))+
  scale_color_manual(values=rev(met.brewer("Johnson", 3)))+
  guides(color=guide_legend(title='Comparisons'))+
  
  
  geom_signif( data=signif_table_acc_ac1,
               aes(y_position = ypos, xmin =xmin, 
              xmax = xmax,annotations = p.value, group=Site),textsize=8,
               tip_length =0, manual=T,inherit.aes = FALSE, color='black')+
  theme_bw()+
  # guides(color='none')+
  ylab("Agreement percentage/ Gwet's AC1 values")+
  xlab('')+
  theme_bw()+
  theme(axis.text=element_text(size=14),
        axis.title=element_text(size=15),
        legend.text = element_text(size=13),
        legend.title= element_text(size=16),
        text = element_text(size = 13),
        strip.text.x = element_text(size = 20))+
  coord_cartesian(ylim=c(0,1.05))+
  scale_y_continuous(breaks=seq(0,1, by=0.2))


## 2.4.Save plots -------------------------------------------------------------



# 
ggsave(Acc_agree_allcomp_p, file=here('outputs','Fig_3_conformity_multiplot.png'),
       dpi=350,  width=12, height=8)




# 3. Permanova trend ---------------------------------------------------------

data_both=rbind(data_val_date,data_amb_date) %>% 
  filter(Site %in% c('BI','BU','VB')) %>% 
  mutate(Quat=as.yearqtr(Date),
         Month=month(Date)) %>% 
  mutate(Comp=case_when(Diver %in% c('VD_','PC_','ZD_') ~'PS',
                        TRUE ~'nPS'),
         Site2=paste0(Site,'_',Comp)) 




# Split by site
data_both_split=split(data_both, data_both$Site)


# apply function to perform PERMANOVA by site
Permanova_site=lapply(data_both_split, function(x){
  
  x_mat=x %>% 
    dplyr::select(-c(Date,Site,Diver,Quadrat,Code_replicat,
                     Quat, Comp,Month, Site2))
  
  x_mat_BC=vegdist(x_mat, method='bray')
  
  AD=adonis2(x_mat_BC ~ Comp, 
             strata=x$Month ,
             data=x)
  
  return(AD %>% broom::tidy())
  
})
  
  

Permanova_table<-bind_rows(Permanova_site, .id = "Site")
  


# export table
write_xlsx(Permanova_table, path=here('outputs',
                    'Permanova_table_comp.xlsx'))



## PCoA --------------------------------------------------------------------

data_both_mat=data_both %>% 
  dplyr::select(-c(Date,Site,Diver,Quadrat,Code_replicat,
                   Quat, Comp,Month, Site2))



SPtot_BC=vegdist(data_both_mat, method='bray')




pcoa_SP_tot<- cmdscale(d = as.matrix(SPtot_BC), eig = T) # vegan package
# Eigen Values

Eig=pcoa_SP_tot$eig[1:2]/sum(pcoa_SP_tot$eig)
Eig=Eig*100



# Data frame pour grp structure
PCOA_site_df=data.frame(x=pcoa_SP_tot$points[,1],y=pcoa_SP_tot$points[,2], 
                        site=data_both$Site,
                        month=data_both$Month,
                        quat=data_both$Quat,
                        comp=data_both$Comp,
                        site_comp=data_both$Site2)%>%
  group_by(comp,site, month) %>%
  mutate(Bar_x=mean(x), Bar_y=mean(y)) 



PCOA_site_mean_df=PCOA_site_df %>%  
  group_by(comp, site, month) %>% 
  
  # group_by(site_comp,month, site) %>% 
  summarize(Bar_x=mean(x), Bar_y=mean(y)) %>% 
  mutate(site_comp=paste0(site,'_',comp))





# palette color formatting
palette_light=lighten(rev(met.brewer("Lakota", 3)),amount=.4)
palette=darken(rev(met.brewer("Lakota", 3)),amount=.2)

Palette=data.frame(comp=c('nPS','nPS','nPS','PS', 'PS','PS'),
                   site=c('BI','BU','VB','BI','BU','VB'),
           palette=c(palette,palette_light)) %>% 
  arrange(factor(site, levels = c('BI','BU','VB')))





## plot
PCOA_p=ggplot(data=PCOA_site_df,aes(x,y))+
  geom_point(aes(col=site_comp), size=1.2, alpha=.2)+
  
  geom_point(data=PCOA_site_mean_df,
             aes(x=Bar_x, y=Bar_y,col=site_comp), size=4, alpha=.8)+
  geom_segment(aes(x=Bar_x, y=Bar_y, xend=x,yend=y,col=site_comp), alpha=.2)+
  scale_color_manual(values=Palette$palette,
                     " Site / group")+
  stat_ellipse(type = "norm",aes(color=site_comp), size=1.1, linetype=1)+

  theme_bw()+
  geom_hline(yintercept=0, linetype=2)+
  geom_vline(xintercept=0, linetype=2)+
  
  xlab(sprintf("PCo1 (%.02f %%)",Eig[1]))+
  ylab(sprintf("PCo2 (%.02f %%)",Eig[2]))+
  theme(axis.text=element_text(size=14),
        axis.title=element_text(size=15),
        legend.text = element_text(size=13),
        legend.title= element_text(size=17),
        text = element_text(size = 13),
        strip.text.y = element_text(size = 20))




ggsave(PCOA_p, file=here('outputs','Fig_4_PCOA_site_group.png'),
       dpi=350,  width=12, height=8)



# 4. Taxa level analysis ----------------------------------------------

Data_nPS_merge2=Data_nPS_merge %>% 
  filter(diff_d==0)



# Data split for nPS

split_data_nPS<-split(Data_nPS_merge2, list(Data_nPS_merge2$Taxa, Data_nPS_merge2$Date,
                                            Data_nPS_merge2$Site), drop=TRUE)


gwet_values_nPS=lapply(split_data_nPS, function(m){
  
  gwet=irrCAC::gwet.ac1.raw(data.frame(m$value_ref, m$value_cand))
  gwet_sel=gwet$est %>% dplyr::select(coeff.val,coeff.se, conf.int, p.value)
  
})
  

crosstables_nPS_taxa_df <- bind_rows(gwet_values_nPS, .id = "Taxa_date") %>% 
  mutate(Comp='nPS_nPS')




Data_merge2=Data_merge %>% 
  filter(diff_d==0) 


# Split of Data_merge by date and by taxa
split_data_taxa_date <- split(Data_merge2, list(Data_merge2$Taxa, Data_merge2$Date,
                                                Data_merge2$Site), drop=TRUE)

# # Apply function cross_table

gwet_values=lapply(split_data_taxa_date, function(m){
  
  gwet=irrCAC::gwet.ac1.raw(data.frame(m$value_ref, m$value_cand))
  gwet_sel=gwet$est %>% dplyr::select(coeff.val,coeff.se, conf.int, p.value)
  
})


crosstables_taxa_df <- bind_rows(gwet_values, .id = "Taxa_date") %>% 
   mutate(Comp='nPS_PS') %>% 
   bind_rows(crosstables_nPS_taxa_df)%>%
  mutate(Taxa=stringr::str_split_i(Taxa_date,'\\.', 1),
         Date=as.Date(stringr::str_split_i(Taxa_date,'\\.', 2)),
         Site=stringr::str_split_i(Taxa_date,'\\.', 3)) %>%
  left_join(Taxa_name, by=c('Taxa'='SP_Niv_good')) %>% 
  mutate(sp_name_ok2=ifelse(is.na(sp_name_ok), Taxa, sp_name_ok)) %>% 
  mutate(sp_name_ok3=if_else(It, glue("<i>{sp_name_ok2}</i>"),sp_name_ok2 ),
         sp_name_ok3=recode_factor(sp_name_ok3,
                                   '<i>Balanus spp.</i>' = '<i>Balanus</i> spp.', 
                                   '<i>Botryllus spp.</i>' = '<i>Botryllus</i> spp.',
                                   '<i>Ophiotrix spp.</i>'= '<i>Ophiotrix</i> spp.',
                                   '<i>Galathea spp.</i>' = '<i>Galathea</i> spp.')) %>% 
  mutate(Type_ok=
           case_when(Type=='Mobile'~ 'Mobile epibenthos',
                     Type == 'Fixed Epibenthos' ~'Fixed epibenthos',
                     Type == 'Fish' ~'Fishes',
                     TRUE ~Type))  %>% 
  mutate(Comparisons=case_when(Comp=='nPS_PS' ~'Inter-groups',
                              Comp=='nPS_nPS' ~'Intra-nPS'))






















### 4.2. Model taxa-dependant ------




mod_taxa_specific <- glmmTMB(
  coeff.val ~ Comp:sp_name_ok3 + (1 | Date) + (1 | Site),
  dispformula = ~ Comp  ,
  family = gaussian(),
  data = crosstables_taxa_df)



summary(mod_taxa_specific)
simulateResiduals(mod_taxa_specific) |> plot()




summary(mod_taxa_specific)

mod_taxa_specific_table=broom.mixed::tidy(mod_taxa_specific) %>% 
  mutate(p.value=as.numeric(p.value),
         p.value = case_when(
           p.value < 0.001 ~ "***",
           p.value < 0.01  ~ "**",
           p.value < 0.05  ~ "*",
           TRUE ~ as.character(p.value)))


# write_xlsx(mod_taxa_specific_table, path=here('outputs','summary_mod_taxa_specific.xlsx'))




## model check  
df_res_taxa_specific<- data.frame(log_fitted = fitted(mod_taxa_specific),
                             residuals  = resid(mod_taxa_specific, 'deviance'))


p1_taxamod<-ggplot(df_res_taxa_specific, aes(x = log_fitted, y = residuals))+
   geom_smooth(method='lm')+
   
  geom_point(alpha=.3) +
  geom_hline(yintercept=0, linetype=2)+
  labs(x = "Linear predictor", y = "Deviance residual")+
  theme_bw() +
  theme(axis.text=element_text(size=14),
        axis.title=element_text(size=14),
        legend.text = element_text(size=13),
        legend.title= element_text(size=14),
        text = element_text(size = 13),
        strip.text.x = element_text(size = 20))
 
 Resids <- residuals(mod_taxa_specific)
 
p2_taxamod<-ggplot(data.frame(resid = Resids), aes(sample = resid)) +
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


all_taxa_spec_diagplots<-ggarrange(p2_taxamod, p1_taxamod, nrow=1,ncol=2)
all_taxa_spec_diagplots

ggsave(all_taxa_spec_diagplots, file=here('outputs','all_taxa_spec_diagplots.png'),
       dpi=350,  width=8, height=6)



#



em_values <- emmeans(mod_taxa_specific, ~ Comp | sp_name_ok3)
contrast_values_bonferroni <- pairs(em_values, adjust = "bonferroni")
summary(contrast_values_bonferroni) %>% 
  filter(p.value<=0.05) %>% 
  arrange(sp_name_ok3)


contrast_taxon_agreement<-broom::tidy(contrast_values_bonferroni) %>% 
  mutate(p.value=as.numeric(p.value),
         p.value = case_when(
           p.value < 0.001 ~ "***",
           p.value < 0.01  ~ "**",
           p.value < 0.05  ~ "*",
           TRUE ~ as.character(p.value)))

write_xlsx(contrast_taxon_agreement, path=here('outputs','contrast_mod_taxa_specific.xlsx'))




## 4.3 Plot taxa specific agreement ----------




signif_table_taxa=contrast_values_bonferroni %>% 
  as.data.frame() %>% 
  filter(p.value<=0.05) %>% 
  mutate(p.value=as.numeric(p.value),
         p.value = case_when(
           p.value < 0.001 ~ "***",
           p.value < 0.01  ~ "**",
           p.value < 0.05  ~ "*",
           TRUE ~ as.character(p.value))) %>% 
  left_join( crosstables_taxa_df %>%  select(c(sp_name_ok3 ,Type_ok)),
             by=c('sp_name_ok3')) %>% 
  distinct()


# add bold name taxa if significant
crosstables_taxa_df=crosstables_taxa_df %>% 
  mutate(sp_name_ok3=case_when(sp_name_ok3 %in% signif_table_taxa$sp_name_ok3 ~  paste0("<b>", sp_name_ok3, "</b>"),
                            TRUE ~sp_name_ok3))

signif_table_taxa=signif_table_taxa %>% 
  mutate(sp_name_ok3=paste0("<b>", sp_name_ok3, "</b>"))




#  order the taxa by mean Gwet AC1 values
crosstables_taxa_df$sp_name_ok3 <- factor(crosstables_taxa_df$sp_name_ok3, 
                                          levels = crosstables_taxa_df %>%
                                            group_by(sp_name_ok3) %>%
                                            summarize(mean_gwet = mean(coeff.val, na.rm = TRUE)) %>%
                                            arrange(mean_gwet) %>%
                                            pull(sp_name_ok3))



crosstables_taxa_df$Comparisons <- factor(crosstables_taxa_df$Comparisons, 
                                          levels = c('Intra-nPS','Inter-groups'))

# Plot

Gwen_ac1_p=ggplot(crosstables_taxa_df)+
  scale_fill_identity() +
  theme_bw()+
  coord_flip(ylim=c(-0.5,1.2))+
  scale_y_continuous(breaks=seq(0,1, by=0.2))+
   facet_grid(Type_ok~., scale='free_y', space='free_y')+
   #facet_wrap(Type_ok~., scale='free_y',ncol=2)+
  
  
  stat_summary(data=crosstables_taxa_df, aes(sp_name_ok3,coeff.val, group=Comparisons,
                                             colour=Comparisons),
               
               fun.data= 'mean_se', geom='pointrange', size=0.5,lwd=0.5,
               position=position_dodge(0.5))+
  # stat_summary(data=crosstables_taxa_df, aes(sp_name_ok3,coeff.val, group=sp_name_ok3),
  #              
  #              fun.y= 'mean', geom='point', size=3, alpha=.7)+
  scale_color_manual(values=rev(met.brewer("Johnson", 2)))+
  geom_text( data=signif_table_taxa,
               aes(y=1.1, x=sp_name_ok3  ,label = p.value, group =Type_ok),size=8, color='black', vjust=0.5)+

  geom_hline(yintercept=c(0,0.2,0.4,0.6,0.8, 1), linetype=2, color='black', size=0.5)+
  ylab('Gwet AC1 values')+
  scale_y_continuous(expand = c(0, 0), breaks=c(-0.5,seq(0,1,by=0.2)))+
  xlab('')+
  theme(axis.text.y = ggtext::element_markdown(),
        axis.text=element_text(size=13),
        axis.title=element_text(size=16),
        legend.text = element_text(size=13),
        legend.title= element_text(size=15),
        text = element_text(size = 13),
        strip.text.x = element_text(size = 20),
        legend.position="top")


Gwen_ac1_p


# geom_density_ridge comparison between Type_ok 
x_marginalplot<-ggplot(crosstables_taxa_df, aes(coeff.val, '', color=Type_ok, fill=Type_ok))+
stat_summary(fun.data= 'mean_se', geom='pointrange', size=1,lwd=1,
             position=position_dodge(0.5), aes(group=Type_ok))+
  geom_jitter(aes(color=Type_ok), alpha=.3, size=2, position=position_jitterdodge(dodge.width = 0.5,
                                                                                  jitter.width=0.2,
                                                                                  jitter.height=0.2))+
  
  scale_color_manual(values=met.brewer("Hiroshige", 4))+
  theme_void()+
  guides(fill='none',color='none')+
  xlab('')+
  ylab('')+
  theme(axis.text=element_text(size=14),
        axis.title=element_text(size=15),
        legend.text = element_text(size=13),
        legend.title= element_text(size=14),
        text = element_text(size = 13),
        strip.text.x = element_text(size = 20))+
  coord_cartesian(xlim=c(-0.5,1.05))+
  scale_x_continuous(limits=c(-0.5,1.05), breaks=seq(0,1, by=0.2))+
  guides(color='none')





cowplot::plot_grid(x_marginalplot, Gwen_ac1_p, ncol = 1, align = "hv", 
          rel_widths = c(2, 1), rel_heights = c(1, 2))



## Reported values -----

taxa_gwet_table_phylum=crosstables_taxa_df %>% 
  group_by(Type_ok,sp_name_ok3, Comp) %>% 
  summarize(mean_gwet = mean(coeff.val, na.rm = TRUE),
            se_gwet=sd(coeff.val)/sqrt(n()))



taxa_gwet_table=crosstables_taxa_df %>% 
  group_by(sp_name_ok3,Type_ok, Comp) %>% 
  summarize(mean_gwet = mean(coeff.val, na.rm = TRUE),
            se_gwet=sd(coeff.val)/sqrt(n())) %>%
  arrange(mean_gwet) 





taxa_gwet_table_sel=taxa_gwet_table=crosstables_taxa_df %>% 
  group_by(sp_name_ok3,Type_ok, Comp) %>% 
  summarize(mean_gwet = mean(coeff.val, na.rm = TRUE),
            se_gwet=sd(coeff.val)/sqrt(n())) %>%
  arrange(mean_gwet)



taxa_gwet_table %>% 
  ungroup() %>% 
  group_by(Type_ok, Comp) %>% 
  summarize(mean=mean(mean_gwet),
            sd=sd(mean_gwet))


test=taxa_gwet_table_sel %>% 
  ungroup() %>% 
  select(-se_gwet) %>% 
  pivot_wider(names_from = Comp,
              values_from = c(mean_gwet)) %>% 
  filter(nPS_nPS > 0.4 & nPS_PS  > 0.4) 


nrow(test)/ 48



# 5. Plot description communities (Figure 1) ---------------------------------

Data_all_pivoted=Data_2021 %>% 
  mutate(Comp=case_when(Diver %in% c('VD_','PC_','ZD_') ~'PS',
                        TRUE ~'nPS')) %>% 
  pivot_longer(cols=-c(Date,Site,Diver,Quadrat,Code_replicat, Comp),
               names_to = 'Taxa',
               values_to='Abundance') %>% 
  # filter(! Comp == 'PS') %>%
  left_join(Taxa_name, by=c('Taxa'='SP_Niv_good')) %>% 
  mutate(sp_name_ok3=if_else(It, glue("<i>{sp_name_ok}</i>"),sp_name_ok),
         sp_name_ok3=recode_factor(sp_name_ok3,
                                   '<i>Balanus spp.</i>' = '<i>Balanus</i> spp.', 
                                   '<i>Botryllus spp.</i>' = '<i>Botryllus</i> spp.',
                                   '<i>Ophiotrix spp.</i>'= '<i>Ophiotrix</i> spp.',
                                   '<i>Galathea spp.</i>' = '<i>Galathea</i> spp.'))




  Data_all_pivoted$Taxa <- factor(Data_all_pivoted$Taxa, 
                                          levels = Data_all_pivoted %>%
                                     group_by(Taxa) %>% 
                                     summarize(Preval=sum(Abundance/n())) %>% 
                                     arrange(-Preval) %>% 
                                     pull(Taxa))


  
ggplot(Data_all_pivoted, aes(sp_name_ok3, Abundance, fill=Comp))+
  stat_summary(fun.data= 'mean_se',
               geom='pointrange', size=.3,lwd=.2, aes(color=Comp,group=Comp),
               position=position_dodge(0.8))+
  stat_summary(fun.y= 'mean', geom='point', size=2, alpha=.7,
               aes(group=1))+
  # facet_nested_wrap(vars(Type_2, Phylum), scales='free_x')
            # facet_grid(, scales='free', space='free')+
  facet_nested(Site~Type + Phylum, scale='free', space='free')+
 # scale angle 90° for labels
  xlab('')+
  ylab('Number of surveys')+
  theme_bw()+
  theme(axis.text.x = ggtext::element_markdown(angle = -45, vjust = 0.5, hjust=1),
    axis.text=element_text(size=14),
        axis.title=element_text(size=15),
        legend.text = element_text(size=13),
        legend.title= element_text(size=14),
        text = element_text(size = 13),
        strip.text.x = element_text(size = 11))



#### 5.2 Table prevalence (Appendix 2)------

Table_Preval<-Data_all_pivoted %>% 
  ungroup() %>% 
  distinct() %>% 
  group_by(sp_name_ok3,Site,Comp,Phylum, Type) %>% 
  mutate(tot=n()) %>% 
  ungroup() %>% 
  group_by(sp_name_ok,Site,Comp, tot,Phylum,Type) %>% 
  summarize(Preval=sum(Abundance)/tot) %>% 
  distinct() %>% 
  mutate(Comp_Site = paste(Comp, Site, sep = "_")) %>%  
  ungroup() %>% 
  select(-c(Comp,Site,tot)) %>% 
  pivot_wider(
    names_from = Comp_Site,
    values_from = Preval,
    values_fill = 0                                        # fill with 0 when missing
  ) %>% 
  arrange(Type)

write_xlsx(Table_Preval, path=here('outputs','Table_taxa_Preval.xlsx'))
