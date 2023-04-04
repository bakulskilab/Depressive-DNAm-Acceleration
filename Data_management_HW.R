library(tidyverse)
library(ggplot2)
library(dplyr)  
library(gtsummary)
library(rmarkdown)
library(knitr)
library(corrplot)
library(tidyselect)
library(labelled)
#install.packages("expss")
library(expss)

Sys.setenv(R_USER ="C:/Users/15340/Documents")
setwd("C:/Users/herongw/Desktop/Kelly_RA/Projects/depressive_DNAm_HRS/Brittant_original/code_dataset")
load("included_3915.RData") ## Brittany's analytic sample
load("HRS_data_VBSvars_2016 v3.RData")

hrs = HRS_data_VBSvars_2016
rm(HRS_data_VBSvars_2016)

load("randhrs1992_2018v1.rda")
rand = randhrs1992_2018v1
rm(randhrs1992_2018v1)


attr(rand$R13PROXY,"label") #Whether take the proxy interview in wave 13

##potential risk line 47-49 check the maximum of age variable
 ###checked
  ####summary(wave13_included$age)
  ####summary(wave13_included$BIRTHYR)
  ####table(wave13_included$RACE)
  ####table(wave13_included$HISPANIC)

##potential risk different code for race/Hispanic missing, maybe not just NA
 ###need double check the label of race/Hispanic

##check if any missingness for other three genetic clocks in included dataset
 ##checked #60-71

##check cross-table on R13SMOKV and R13SMOKN. each cell > 100
  ##checked #75-


#Bri's code 17-28
hrs_rand=rand %>% 
  #select(HHID,PN,R13PROXY) %>% 
  #right_join(hrs,by = c("HHID", "PN")) %>% 
  #filter (INW13 == 1 )


#df=rand %>% 
  #select(vars_select(names(.),starts_with("R13")))
#var_label(df)

#df1=rand %>% 
  #select(vars_select(names(.),starts_with("R13SMOK")))
#view(df1)
#table(df1$R13SMOKEV,df1$R13SMOKEN,useNA = "always")

## select subs with complete 5 genetic clocks and covariates
gene_clock = c("DNAMGRIMAGE", "LEVINE_DNAMAGE", "HORVATH_DNAMAGE", "MPOA", "HANNUM_DNAMAGE")
my_vars = c('HHID', 'PN','INW13', 'R13PROXY','BIRTHYR', 'GENDER', 'R13CESD', 'DNAMGRIMAGE', 'HORVATH_DNAMAGE',
            'RACE', 'HISPANIC', 'DEGREE', 'PMARST', 'R13CONDE',
            'R13SMOKEV', 'R13SMOKEN', 'R13DRINKD', 'PBASO', 'PEOS', 'PNEUT', 'PMONO', 'PLYMP')


#w13_clocks=hrs_rand %>% 
  #select(HHID, PN, gene_clock) %>% 
  #dplyr::filter (across(all_of(gene_clock), ~!is.na(.x)))

w13_clocks=hrs_rand %>% 
  select(HHID, PN, gene_clock, my_vars) %>% 
  dplyr::filter(if_all(gene_clock, ~!is.na(.x))) %>% 
  mutate(age = 2016 - BIRTHYR)

#w13_clocks %>% 
  #anti_join(wave13_2016_clocks, by = c("HHID", "PN")) ## complete Grim and Horvath also with complete other 3 clocks

#included_clocks <- wave13_included %>% 
  #inner_join(w13_clocks, by = c("HHID", "PN"))

#age acceleration based on GrimAge
w13_clocks <- w13_clocks %>% 
  mutate(age_diff_g = DNAMGRIMAGE - age)

w13_clocks$GRIM_residuals <- residuals(lm(DNAMGRIMAGE ~ age, data = w13_clocks, na.action = na.exclude))
summary(w13_clocks$GRIM_residuals)

w13_clocks <- w13_clocks %>% mutate(GRIMAccel = ifelse((GRIM_residuals > 0), 1, 0))
w13_clocks$GRIMAccel <- factor(w13_clocks$GRIMAccel, levels = c(0,1), labels = 
                                     c("Age decelerated", "Age accelerated"))
table(w13_clocks$GRIMAccel)

#age acceleration based on HORVATH
w13_clocks <- w13_clocks %>% 
  mutate(age_diff_h = HORVATH_DNAMAGE - age)

w13_clocks$HORVATH_residuals <- residuals(lm(HORVATH_DNAMAGE ~ age, data = w13_clocks, na.action = na.exclude))
summary(w13_clocks$HORVATH_residuals)

w13_clocks <- w13_clocks %>% mutate(HORVATHAccel = ifelse((HORVATH_residuals > 0), 1, 0))
w13_clocks$HORVATHAccel <- factor(w13_clocks$HORVATHAccel, levels = c(0,1), labels = 
                                 c("Age decelerated", "Age accelerated"))
table(w13_clocks$HORVATHAccel)

#age acceleration based on LEVINE
w13_clocks <- w13_clocks %>% 
  mutate(age_diff_l = LEVINE_DNAMAGE - age)

w13_clocks$LEVINE_residuals <- residuals(lm(LEVINE_DNAMAGE ~ age, data = w13_clocks, na.action = na.exclude))
summary(w13_clocks$LEVINE_residuals)

w13_clocks <- w13_clocks %>% mutate(LEVINEAccel = ifelse((LEVINE_residuals > 0), 1, 0))
w13_clocks$LEVINEAccel <- factor(w13_clocks$LEVINEAccel, levels = c(0,1), labels = 
                                    c("Age decelerated", "Age accelerated"))
table(w13_clocks$LEVINEAccel)

#age acceleration based on MPOA
w13_clocks <- w13_clocks %>% 
  mutate(age_diff_m = MPOA*age - age)

w13_clocks$MPOA_residuals <- residuals(lm(MPOA*age ~ age, data = w13_clocks, na.action = na.exclude))
summary(w13_clocks$MPOA_residuals)

w13_clocks <- w13_clocks %>% mutate(MPOAAccel = ifelse((MPOA_residuals > 0), 1, 0))
w13_clocks$MPOAAccel <- factor(w13_clocks$MPOAAccel, levels = c(0,1), labels = 
                                    c("Age decelerated", "Age accelerated"))
table(w13_clocks$MPOAAccel)

#full$age_diff_m = full$MPOA*full$age - full$age
#full$MPOA_residuals <- residuals(lm(MPOA*age ~ age, data = full, na.action = na.exclude))
#full <- full %>% mutate(MPOAAccel = ifelse((MPOA_residuals > 0), 1, 0))
#full$MPOAAccel <- factor(full$MPOAAccel, levels = c(0,1), labels = 
                                 #c("Age decelerated", "Age accelerated"))
#prop.table(table(full$MPOAAccel, useNA = "always"))
#mean(full$MPOA_residuals)
#sd(full$MPOA_residuals)

#complete = full %>% 
  #filter(complete.cases(.))
#table(complete$MPOAAccel, useNA = "always")
#prop.table(table(complete$MPOAAccel, useNA = "always"))
#mean(complete$MPOA_residuals)
#sd(complete$MPOA_residuals)

#exclude = full %>% 
  #anti_join(complete, by = c("HHID", "PN"))
#table(exclude$MPOAAccel, useNA = "always")
#prop.table(table(exclude$MPOAAccel, useNA = "always"))
#mean(exclude$MPOA_residuals)
#sd(exclude$MPOA_residuals)

#complete %>% 
  #group_by(GRIMAccel) %>% 
  #summarize(mean = mean(MPOA_residuals), sd = sd(MPOA_residuals))

#complete %>% 
  #group_by(depression) %>% 
  #summarize(mean = mean(MPOA_residuals), sd = sd(MPOA_residuals))

#age acceleration based on HANNUM
w13_clocks <- w13_clocks %>% 
  mutate(age_diff_ha = HANNUM_DNAMAGE - age)

w13_clocks$HANNUM_residuals <- residuals(lm(HANNUM_DNAMAGE ~ age, data = w13_clocks, na.action = na.exclude))
summary(w13_clocks$HANNUM_residuals)

w13_clocks <- w13_clocks %>% mutate(HANNUMAccel = ifelse((HANNUM_residuals > 0), 1, 0))
w13_clocks$HANNUMAccel <- factor(w13_clocks$HANNUMAccel, levels = c(0,1), labels = 
                                 c("Age decelerated", "Age accelerated"))
table(w13_clocks$HANNUMAccel)

# Create dichotomous depression variable cut off at 4

w13_clocks <- w13_clocks %>% mutate(depression = ifelse((R13CESD >= 4), 1, 0))

w13_clocks$depression <- factor(w13_clocks$depression, levels = c(0,1), labels=
                                   c("Low or no depressive symptoms", "Elevated depressive symptoms"))

table(w13_clocks$depression, useNA = "always")

# Create dichotomous depression variable cut off at 1- sensitivity analysis

w13_clocks <- w13_clocks %>% mutate(any_dep = ifelse((R13CESD > 0), 1, 0))

w13_clocks$any_dep <- factor(w13_clocks$any_dep, levels = c(0,1), labels = 
                                c("No depressive symptoms", "Has depressive symptoms"))
table(w13_clocks$any_dep, useNA = "always")

# Create granulocyte variable

w13_clocks <- w13_clocks %>% mutate(gran = PBASO + PEOS + PNEUT)
summary(w13_clocks$gran)
summary(w13_clocks$PBASO)
summary(w13_clocks$PEOS)
summary(w13_clocks$PNEUT)

# Create total cell type percentage variable

w13_clocks <- w13_clocks %>% mutate(cell_total = gran + PMONO + PLYMP)
summary(w13_clocks$cell_total)
summary(w13_clocks$PLYMP)
summary(w13_clocks$PMONO)

w13_clocks %>% 
  select (HHID, PN, PLYMP, PMONO, gran, cell_total) %>% 
  filter(PLYMP > 90)

# Make R13SHLT numeric

class(w13_clocks$R13SHLT)
w13_clocks$R13SHLT <- as.numeric(w13_clocks$R13SHLT)

####### Collapsing variable categories #######

# Collapse DEGREE to no degree, GED/HS diploma, degree unknown/some college
# Two or four year degree, master's/professional degree

# Categories in original variable are 0, 1, 2, 3, 4, 5, 6, 9- HRS file 
table(w13_clocks$DEGREE, useNA = "always")
w13_clocks <- w13_clocks %>% mutate(DEGREE_collapsed = case_when(
  DEGREE == 0 ~ 0,
  DEGREE == 1 ~ 1,
  DEGREE == 2 ~ 1,
  DEGREE == 3 ~ 3,
  DEGREE == 4 ~ 3, 
  DEGREE == 5 ~ 4, 
  DEGREE == 6 ~ 4, 
  DEGREE == 9 ~ 2))

table(w13_clocks$DEGREE_collapsed, useNA = "always")
w13_clocks$DEGREE_collapsed <- factor(w13_clocks$DEGREE_collapsed, 
                                       levels = c(0,1,2,3,4), 
                                       labels= c("No degree", "GED/HS diploma",
                                                 "Degree unknown/some college",
                                                 "Two or four year college degree",
                                                 "Master's or professional degree"))

table(w13_clocks$DEGREE_collapsed, w13_clocks$DEGREE) 

# Collapse marital status categories

table(w13_clocks$PMARST)
w13_clocks <- w13_clocks%>% mutate(PMARST_collapsed = case_when(
  PMARST == 1 ~ 1,
  PMARST == 2 ~ 2, 
  PMARST == 3 ~ 2, 
  PMARST == 4 ~ 2, 
  PMARST == 5 ~ NA_real_))
table(w13_clocks$PMARST_collapsed, useNA = "always")

w13_clocks$PMARST_collapsed <- factor(w13_clocks$PMARST_collapsed,
                                       levels = c(1,2),
                                       labels = c("Married", "Single"))

table(w13_clocks$PMARST_collapsed, w13_clocks$PMARST)

# Add labels to variables' values

w13_clocks$GENDER <- factor(w13_clocks$GENDER, levels = c(1, 2), labels=c("Male", "Female"))
w13_clocks$RACE <- factor(w13_clocks$RACE, levels = c(0, 1, 2, 7),
                           labels=c("Not obtained", "White", "Black", "Other"))

w13_clocks$HISPANIC <- factor(w13_clocks$HISPANIC, levels =c(0, 1, 2, 3, 5), 
                               labels=c("Not obtained", "Hispanic, Mexican", "Hispanic, other",
                                        "Hispanic, type unknown", "Non-Hispanic"))
w13_clocks$R13SMOKEN <- factor(w13_clocks$R13SMOKEN, levels = c(0,1), labels=c("Don't smoke now", "Smoke now"))
w13_clocks$R13SMOKEV <- factor(w13_clocks$R13SMOKEV, levels = c(0,1), labels=c("Never smoked", "Have smoked"))


# Setting reference levels for categorical variables 

w13_clocks$GENDER <- relevel(w13_clocks$GENDER, ref = "Male")
w13_clocks$RACE <- relevel(w13_clocks$RACE, ref = "White") 
w13_clocks$DEGREE_collapsed <- relevel(w13_clocks$DEGREE_collapsed, ref = "No degree") 
w13_clocks$PMARST_collapsed <- relevel(w13_clocks$PMARST_collapsed, ref = "Married") 
w13_clocks$depression <- relevel(w13_clocks$depression, ref = "Low or no depressive symptoms")
w13_clocks$any_dep <- relevel(w13_clocks$any_dep, ref = "No depressive symptoms")
w13_clocks$HISPANIC <- relevel(w13_clocks$HISPANIC, ref = "Non-Hispanic")

######## create complete analytic sample version one ######
w13_clocks %>% filter(rowSums(is.na(w13_clocks)) == 0)
w13_complete_1 <- w13_clocks %>% 
  filter(complete.cases(.)) %>% 
  filter(RACE != "Not obtained") %>%
  filter(PMARST_collapsed != "Marital status unknown") %>%
  filter(HISPANIC != "Not obtained")

######### Adding physical activity and financial situation before age 16

## cross-table on R13SMOKV and R13SMOKN. each cell > 100
table(wave13_included$R13SMOKEV,wave13_included$R13SMOKEN,useNA = "always")
cross_cases(wave13_included,R13SMOKEV,list(R13SMOKEN, total()))
cross_cases(wave13_included,R13SMOKEV,R13SMOKEN)

cross_cases(wave13_2016_clocks,R13SMOKEV,list(R13SMOKEN, total()))

## create three level smoking
table(w13_clocks$R13SMOKEN)
table(w13_clocks$R13SMOKEV)
w13_clocks_smok <- w13_clocks %>% 
  mutate(smoke = case_when(R13SMOKEV == "Have smoked" & R13SMOKEN == "Smoke now" ~ 2,
                           R13SMOKEV == "Have smoked" & R13SMOKEN == "Don't smoke now" ~ 1,
                           R13SMOKEV == "Never smoked" & R13SMOKEN == "Don't smoke now" ~ 0))
w13_clocks_smok$smoke <- factor(w13_clocks_smok$smoke, levels = c(0,1,2), labels = c("Never", "Former", "Current"))

table(w13_clocks_smok$smoke, useNA = "always")
w13_clocks_smok$smoke <- relevel(w13_clocks_smok$smoke, ref = "Never")

## add covariate physical activities
physical_activity = c("R13VGACTX", "R13MDACTX", "R13LTACTX")
w13_clocks_phy=rand %>% 
  select(HHID,PN, physical_activity) %>% 
  right_join(w13_clocks_smok,by = c("HHID", "PN"))


table(w13_clocks_phy$R13LTACTX, useNA = "always")
table(w13_clocks_phy$R13MDACTX, useNA = "always")  
table(w13_clocks_phy$R13VGACTX, useNA = "always") 


###reverse the coding of physical activity
class(w13_clocks_phy$R13LTACTX)
class(w13_clocks_phy$R13MDACTX)
class(w13_clocks_phy$R13VGACTX)

w13_clocks_phy_re = w13_clocks_phy%>% 
  mutate (R13LTACTR = case_when(R13LTACTX == 1 ~ 4,
                               R13LTACTX == 2 ~ 3,
                               R13LTACTX == 3 ~ 2,
                               R13LTACTX == 4 ~ 1,
                               R13LTACTX == 5 ~ 0,
                               is.na(R13LTACTX) ~ NA_real_))

table(w13_clocks_phy_re$R13LTACTR, w13_clocks_phy_re$R13LTACTX, useNA = "always")

w13_clocks_phy_re = w13_clocks_phy_re %>% 
  mutate (R13MDACTR = case_when(R13MDACTX == 1 ~ 4,
                               R13MDACTX == 2 ~ 3,
                               R13MDACTX == 3 ~ 2,
                               R13MDACTX == 4 ~ 1,
                               R13MDACTX == 5 ~ 0,
                               is.na(R13MDACTX) ~ NA_real_))

table(w13_clocks_phy_re$R13MDACTR, w13_clocks_phy_re$R13MDACTX, useNA = "always")

w13_clocks_phy_re = w13_clocks_phy_re %>% 
  mutate (R13VGACTR = case_when(R13VGACTX == 1 ~ 4,
                               R13VGACTX == 2 ~ 3,
                               R13VGACTX == 3 ~ 2,
                               R13VGACTX == 4 ~ 1,
                               R13VGACTX == 5 ~ 0,
                               is.na(R13VGACTX) ~ NA_real_))

table(w13_clocks_phy_re$R13VGACTR, w13_clocks_phy_re$R13VGACTX, useNA = "always")

###created the weighted variable
w13_clocks_phy_we = w13_clocks_phy_re %>% 
  mutate (R13LTACTW = R13LTACTR * 2.5,
          R13MDACTW = R13MDACTR * 4,
          R13VGACTW = R13VGACTR * 8,
          R13ACT = R13LTACTW + R13MDACTW + R13VGACTW)

table(w13_clocks_phy_we$R13LTACTR, w13_clocks_phy_we$R13LTACTW, useNA = "always")
table(w13_clocks_phy_we$R13MDACTR, w13_clocks_phy_we$R13MDACTW, useNA = "always")
table(w13_clocks_phy_we$R13VGACTR, w13_clocks_phy_we$R13VGACTW, useNA = "always")

w13_smk_phy <- w13_clocks_phy_we
  

### QC of the new created physical activity
summary(w13_clocks_phy_we$R13ACT)
#hist(w13_clocks_phy_act_we$R13ACT, xlab = "new weighted physical activity", 
    #ylab = "frequency", main = "Histogram of the Weighted Physical Activity")

w13_clocks_phy_we %>% 
  ggplot(aes(x=R13ACT)) + 
  geom_histogram(aes(y=..density..), colour="black", fill="grey", binwidth = 2)+
  geom_density(alpha=.1, fill="#FF6666")+
  labs(title = "Density of the Weighted Physical Activity",
       x = "weighted physical activity" )

table(w13_clocks_phy_we$R13ACT, useNA = "always")
w13_clocks_phy_we %>% 
  select(HHID, PN, R13LTACTR, R13LTACTW, R13MDACTR, R13MDACTW, R13VGACTR, R13VGACTW, R13ACT) %>% 
  filter (R13ACT == 58)


plot(w13_clocks_phy_we$R13LTACTR, w13_clocks_phy_we$R13ACT)

w13_clocks_phy_we = w13_clocks_phy_we %>% 
  group_by(R13LTACTR) %>% 
  mutate(R13LTACTM = mean(R13ACT, na.rm = T)) %>% 
  ungroup %>% 
  group_by(R13MDACTR) %>% 
  mutate(R13MDACTM = mean(R13ACT,na.rm= T)) %>% 
  ungroup %>% 
  group_by(R13VGACTR) %>% 
  mutate(R13VGACTM = mean(R13ACT,na.rm= T)) %>% 
  ungroup


w13_clocks_phy_we %>% 
  ggplot() +
  geom_jitter(aes(y = R13ACT, x = R13LTACTR))+
  geom_line(aes(y = R13LTACTM, x = R13LTACTR), size = 1, color = 'red')+
  labs(title = "Scatter plot of original light vs. weighted physical activity",
       x = "original light physical activity",
       y = "weight physical activity")

w13_clocks_phy_we %>% 
  ggplot() +
  geom_violin(aes(x = factor(R13LTACTR), y = R13ACT))+
  labs(title = "violin plot of original light vs. weighted physical activity",
       x = "original light physical activity",
       y = "weight physical activity")

w13_clocks_phy_we %>% 
  ggplot() +
  geom_jitter(aes(y = R13ACT, x = R13MDACTR))+
  geom_line(aes(y = R13MDACTM, x = R13MDACTR), size = 1, color = 'red')+
  labs(title = "Scatter plot of original moderate vs. weighted physical activity",
       x = "original moderate physical activity",
       y = "weight physical activity")
w13_clocks_phy_we %>% 
  ggplot() +
  geom_violin(aes(x = factor(R13MDACTR), y = R13ACT))+
  labs(title = "violin plot of original moderate vs. weighted physical activity",
       x = "original moderate physical activity",
       y = "weight physical activity")

w13_clocks_phy_we %>% 
  ggplot() +
  geom_jitter(aes(y = R13ACT, x = R13VGACTR))+
  geom_line(aes(y = R13VGACTM, x = R13VGACTR), size = 1, color = 'red')+
  labs(title = "Scatter plot of original vigorous vs. weighted physical activity",
       x = "original vigorous physical activity",
       y = "weight physical activity")
w13_clocks_phy_we %>% 
  ggplot() +
  geom_violin(aes(x = factor(R13VGACTR), y = R13ACT))+
  labs(title = "violin plot of original vigorous vs. weighted physical activity",
       x = "original vigorous physical activity",
       y = "weight physical activity")

plot(w13_clocks_phy_we$R13LTACTR,w13_clocks_phy_we$R13ACT)


##adding covariate financial situation before 16

library(haven)
library(rlang)
w2016 <- read_dta ("C:/Users/15340/Desktop/Kelly_RA/database/HRS/Public_Core/w2016/H16B_R.dta")
w2014 <- read_dta ("C:/Users/15340/Desktop/Kelly_RA/database/HRS/Public_Core/w2016/H14B_R.dta")
w2012 <- read_dta ("C:/Users/15340/Desktop/Kelly_RA/database/HRS/Public_Core/w2016/H12B_R.dta")
w2010 <- read_dta ("C:/Users/15340/Desktop/Kelly_RA/database/HRS/Public_Core/w2016/H10B_R.dta")
w2008 <- read_dta ("C:/Users/15340/Desktop/Kelly_RA/database/HRS/Public_Core/w2016/H08B_R.dta")
w2006 <- read_dta ("C:/Users/15340/Desktop/Kelly_RA/database/HRS/Public_Core/w2016/H06B_R.dta")
w2004 <- read_dta ("C:/Users/15340/Desktop/Kelly_RA/database/HRS/Public_Core/w2016/H04B_R.dta")
w2002 <- read_dta ("C:/Users/15340/Desktop/Kelly_RA/database/HRS/Public_Core/w2016/H02B_R.dta")
w2000 <- read_dta ("C:/Users/15340/Desktop/Kelly_RA/database/HRS/Public_Core/w2016/H00A_R.dta")
w1998 <- read_dta ("C:/Users/15340/Desktop/Kelly_RA/database/HRS/Public_Core/w2016/H98A_R.dta")

w16 = w2016 %>% 
  select(HHID, PN, PB020)
w14 = w2014 %>% 
  select(HHID, PN,OB020)
w12 = w2012 %>% 
  select(HHID, PN,NB020)
w10 = w2010 %>% 
  select(HHID, PN,MB020)
w08 = w2008 %>% 
  select(HHID, PN,LB020)
w06 = w2006 %>% 
  select(HHID, PN,KB020)
w04 = w2004 %>% 
  select(HHID, PN,JB020)
w02 = w2002 %>% 
  select(HHID, PN,HB020)
w00 = w2000 %>% 
  select(HHID, PN,G1080)
w98 = w1998 %>% 
  select(HHID, PN,F993)

rm(w2016, w2014, w2012, w2010, w2008, w2006, w2004, w2002, w2000, w1998)
df_list = list(w16, w14, w12, w10, w08, w06, w04, w02, w00, w98)
var_list = c("PB020", "OB020", 'NB020', 'MB020', 'LB020', 'KB020', 'JB020', 'HB020', 'G1080','F993')

w_merge = df_list %>% 
  reduce(full_join, by = c('HHID', "PN"))

head(w_merge)
rm(w16, w14, w12, w10, w08, w06, w04, w02, w00, w98)


##calculate the number of missing per row
w_merge$count_na <- rowSums(is.na(w_merge))
table(w_merge$count_na, useNA = "always")

w_pov = w_merge %>% 
  filter(count_na == 9) %>% 
  mutate(child_poverty = rowMeans(.[,c(3:12)], na.rm = T))

table(w_pov$child_poverty, useNA = "always")


w_pov = w_pov%>% 
  mutate(poverty = if_else(child_poverty == 6 | child_poverty == 8 | child_poverty == 9, NA_real_, child_poverty))
table(w_pov$poverty, useNA = "always")

w13_smk_phy_pov = w_pov %>% 
  select(HHID, PN, poverty) %>% 
  right_join(w13_smk_phy, by = c('HHID', "PN")) 
table(w13_smk_phy_pov$poverty, useNA = "always")

w13_smk_phy_pov$poverty <- factor(w13_smk_phy_pov$poverty, levels = c(1, 3, 5), labels = c("Well off", "Average", "Poor"))
w13_smk_phy_pov$poverty <- relevel(w13_smk_phy_pov$poverty, ref = "Well off")


######## create complete analytic sample #########

complete <- w13_smk_phy_pov %>% 
  filter(complete.cases(.)) %>% 
  filter(RACE != "Not obtained") %>%
  filter(PMARST_collapsed != "Marital status unknown") %>%
  filter(HISPANIC != "Not obtained") %>% 
  mutate(included = 1)


full <- complete %>% 
  select(HHID, PN, included) %>% 
  full_join(w13_smk_phy_pov, by = c("HHID", "PN")) %>% 
  mutate(included = case_when(included == 1~ "Included",
                             TRUE ~ "Excluded"))
table(full$included, useNA = "always")

save(full, file = "full_4018.RData")
save(complete, file = "analytic_sample_3790.RData")


### change not obtained in the race; hispanic; marital status into missing ###

load("full_4018.RData")
load("analytic_sample_3790.RData")

full <- full %>% 
  mutate(RACE = if_else(RACE == "Not obtained", NA_real_, as.numeric(RACE)))
table(full$RACE)
full$RACE <- factor(full$RACE, levels = c(1, 3, 4),
                          labels=c( "White", "Black", "Other"))
table(full$RACE, useNA = "always")
full$RACE <- relevel(full$RACE, ref = "White") 

full <- full %>% 
  mutate(HISPANIC = if_else(HISPANIC == "Not obtained", NA_real_, as.numeric(HISPANIC)))
table(full$HISPANIC)
full$HISPANIC <- factor(full$HISPANIC, levels =c( 1, 3, 4, 5), 
                              labels=c("Non-Hispanic", "Hispanic, Mexican",
                                       "Hispanic, other", "Hispanic, type unknown"))
table(full$HISPANIC, useNA = "always")

full <- full %>% 
  mutate(PMARST_collapsed = if_else(PMARST_collapsed == "Marital status unknown",
         NA_real_, as.numeric(PMARST_collapsed)))
table(full$PMARST_collapsed)
full$PMARST_collapsed <- factor(full$PMARST_collapsed, levels = c(1,2),
                                labels = c("Married", "Single"))
table(full$PMARST_collapsed)


### collapse race and hispanic into race/ethnicity ###
table(full$RACE, full$HISPANIC, useNA = "always")

full2 <- full %>% 
  mutate(race_ethnicity = case_when(RACE == "White" & HISPANIC == "Non-Hispanic" ~ "Non-Hispanic White",
                                    RACE == "Black" & HISPANIC == "Non-Hispanic" ~ "Non-Hispanic Black",
                                    RACE == "Other" & HISPANIC == "Non-Hispanic" ~ "Non-Hispanic Other",
                                    str_detect(HISPANIC, "Hispanic") & !is.na(RACE) ~ "Hispanic"))
table(full2$race_ethnicity, useNA = "always")

full <- full2

### collapse education ###

table(full$DEGREE_collapsed, useNA = "always")
full3 <- full %>% 
  mutate(DEGREE_collapsed = case_when(DEGREE_collapsed == "No degree" ~ "No degree",
                                      DEGREE_collapsed == "GED/HS diploma" ~ "GED/HS diploma",
                                      TRUE ~ "College and more than college"))
table(full3$DEGREE_collapsed)
full <- full3
  
# Setting reference levels for categorical variables 

full$race_ethnicity <- relevel(factor(full$race_ethnicity), ref = "Non-Hispanic White") 
full$DEGREE_collapsed <- relevel(factor(full$DEGREE_collapsed), ref = "No degree") 
full$PMARST_collapsed <- relevel(full$PMARST_collapsed, ref = "Married") 

full$race_ethnicity <- factor(full$race_ethnicity, levels = c("Non-Hispanic White", "Non-Hispanic Black", 
                                                              "Non-Hispanic Other", "Hispanic"))

full$DEGREE_collapsed <- factor(full$DEGREE_collapsed, levels = c("No degree", "GED/HS diploma", "College and more than college"))


### filter the missing in covariate set one by one ###
nomiss_1 = full %>% 
  filter(!is.na(R13CESD))
loss_depression = nrow(full) - nrow(nomiss_1)

nomiss_2 = nomiss_1 %>% 
  filter(!is.na(GENDER))
loss_gender = nrow(nomiss_1) - nrow(nomiss_2)

nomiss_3 = nomiss_2 %>% 
  filter(!is.na(race_ethnicity))
loss_RaceEthnicity = nrow(nomiss_2) - nrow(nomiss_3)

nomiss_4 = nomiss_3 %>% 
  filter(!is.na(PMARST_collapsed))
loss_martial = nrow(nomiss_3) - nrow(nomiss_4)

nomiss_5 = nomiss_4 %>% 
  filter(!is.na(DEGREE_collapsed))
loss_education = nrow(nomiss_4) - nrow(nomiss_5)

nomiss_6 = nomiss_5 %>% 
  filter(!is.na(R13CONDE))
loss_chronic = nrow(nomiss_5) - nrow(nomiss_6)

nomiss_7 = nomiss_6 %>% 
  filter(!is.na(smoke))
loss_smoke = nrow(nomiss_6) - nrow(nomiss_7)

nomiss_8 = nomiss_7 %>% 
  filter(!is.na(R13DRINKD))
loss_drink = nrow(nomiss_7) - nrow(nomiss_8)

nomiss_9 = nomiss_8 %>% 
  filter(!is.na(R13ACT))
loss_act = nrow(nomiss_8) - nrow(nomiss_9)

nomiss_10 = nomiss_9 %>% 
  filter(!is.na(poverty))
loss_poverty = nrow(nomiss_9) - nrow(nomiss_10)

nomiss_11 = nomiss_10 %>% 
  filter(!is.na(gran))
loss_gran = nrow(nomiss_10) - nrow(nomiss_11)

nomiss_12 = nomiss_11 %>% 
  filter(!is.na(PMONO))
loss_mono = nrow(nomiss_11) - nrow(nomiss_12)

nomiss_13 = nomiss_12 %>% 
  filter(!is.na(PLYMP))
loss_LYMP = nrow(nomiss_12) - nrow(nomiss_13)


#### SAVE NEW COMPLETE AND FULL DATASET ####
full <- full %>% 
  select(-included)

complete <- full %>% 
  filter(complete.cases(.)) %>% 
  mutate(included = 1)

full <- complete %>% 
  select(HHID, PN, included) %>% 
  full_join(full, by = c("HHID", "PN")) %>% 
  mutate(included = case_when(included == 1 ~ "Included",
                              TRUE ~ "Excluded"))

table(full$included, useNA = "always")

save(full, file = "full_4018.RData")
save(complete, file = "analytic_sample_3793.RData")

