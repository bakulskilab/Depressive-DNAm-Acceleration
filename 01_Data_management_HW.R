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


setwd("C:/Users/herongw/Desktop/Kelly_RA/Projects/HW_depressive_DNAm_HRS/Submission/Psychology_aging/Revision_20240806_final/code")

load("HRS_data_VBSvars_2016 v3.RData")

hrs = HRS_data_VBSvars_2016
rm(HRS_data_VBSvars_2016)

load("randhrs1992_2018v1.rda")
rand = randhrs1992_2018v1
rm(randhrs1992_2018v1)


## select subs with complete 5 genetic clocks and covariates
gene_clock = c("DNAMGRIMAGE", "LEVINE_DNAMAGE", "HORVATH_DNAMAGE", "MPOA", "HANNUM_DNAMAGE")
my_vars = c('HHID', 'PN', 'BIRTHYR', 'GENDER', 'R13CESD', 'DNAMGRIMAGE', 'HORVATH_DNAMAGE',
            'RACE', 'HISPANIC', 'DEGREE', 'PMARST', 'R13CONDE',
            'R13SMOKEV', 'R13SMOKEN', 'R13DRINKD', 'PBASO', 'PEOS', 'PNEUT', 'PMONO', 'PLYMP')



w13_clocks=hrs %>% 
  select(HHID, PN, all_of(gene_clock), all_of(my_vars)) %>% 
  dplyr::filter(if_all(gene_clock, ~!is.na(.x))) %>% 
  mutate(age = 2016 - BIRTHYR)


#age acceleration based on GrimAge
w13_clocks$GRIM_residuals <- residuals(lm(DNAMGRIMAGE ~ age, data = w13_clocks, na.action = na.exclude))
summary(w13_clocks$GRIM_residuals)

w13_clocks <- w13_clocks %>% mutate(GRIMAccel = ifelse((GRIM_residuals > 0), 1, 0))
w13_clocks$GRIMAccel <- factor(w13_clocks$GRIMAccel, levels = c(0,1), labels = 
                                     c("Age decelerated", "Age accelerated"))
table(w13_clocks$GRIMAccel)

#age acceleration based on HORVATH
w13_clocks$HORVATH_residuals <- residuals(lm(HORVATH_DNAMAGE ~ age, data = w13_clocks, na.action = na.exclude))
summary(w13_clocks$HORVATH_residuals)

w13_clocks <- w13_clocks %>% mutate(HORVATHAccel = ifelse((HORVATH_residuals > 0), 1, 0))
w13_clocks$HORVATHAccel <- factor(w13_clocks$HORVATHAccel, levels = c(0,1), labels = 
                                 c("Age decelerated", "Age accelerated"))
table(w13_clocks$HORVATHAccel)

#age acceleration based on LEVINE
w13_clocks$LEVINE_residuals <- residuals(lm(LEVINE_DNAMAGE ~ age, data = w13_clocks, na.action = na.exclude))
summary(w13_clocks$LEVINE_residuals)

w13_clocks <- w13_clocks %>% mutate(LEVINEAccel = ifelse((LEVINE_residuals > 0), 1, 0))
w13_clocks$LEVINEAccel <- factor(w13_clocks$LEVINEAccel, levels = c(0,1), labels = 
                                    c("Age decelerated", "Age accelerated"))
table(w13_clocks$LEVINEAccel)

#age acceleration based on MPOA
w13_clocks$MPOA_residuals <- residuals(lm(MPOA*age ~ age, data = w13_clocks, na.action = na.exclude))
summary(w13_clocks$MPOA_residuals)

w13_clocks <- w13_clocks %>% mutate(MPOAAccel = ifelse((MPOA_residuals > 0), 1, 0))
w13_clocks$MPOAAccel <- factor(w13_clocks$MPOAAccel, levels = c(0,1), labels = 
                                    c("Age decelerated", "Age accelerated"))
table(w13_clocks$MPOAAccel)


#age acceleration based on HANNUM

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

w13_complete_1 <- w13_clocks %>% 
  filter(complete.cases(.)) %>% 
  filter(RACE != "Not obtained") %>%
  filter(PMARST_collapsed != "Marital status unknown") %>%
  filter(HISPANIC != "Not obtained")


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

### change not obtained in the race; hispanic; marital status into missing ###

table(w13_clocks_smok$RACE, useNA = "always")

full <- w13_clocks_smok %>% 
  mutate(RACE = if_else(RACE == "Not obtained", NA_real_, as.numeric(RACE)))
table(full$RACE, useNA = "always")

full$RACE <- factor(full$RACE, levels = c(1, 3, 4),
                          labels=c( "White", "Black", "Other"))
table(full$RACE, useNA = "always")

full$RACE <- relevel(full$RACE, ref = "White") 

table(full$HISPANIC, useNA = "always")

full <- full %>% 
  mutate(HISPANIC = if_else(HISPANIC == "Not obtained", NA_real_, as.numeric(HISPANIC)))

table(full$HISPANIC, useNA = "always")
full$HISPANIC <- factor(full$HISPANIC, levels =c(1, 3, 4, 5), 
                              labels=c("Non-Hispanic", "Hispanic, Mexican",
                                       "Hispanic, other", "Hispanic, type unknown"))
table(full$HISPANIC, useNA = "always")

table(full$PMARST_collapsed, useNA = "always")
full <- full %>% 
  mutate(PMARST_collapsed = if_else(PMARST_collapsed == "Marital status unknown",
         NA_real_, as.numeric(PMARST_collapsed)))
table(full$PMARST_collapsed, useNA = "always")
full$PMARST_collapsed <- factor(full$PMARST_collapsed, levels = c(1,2),
                                labels = c("Married", "Single"))
table(full$PMARST_collapsed, useNA = "always")


### collapse race and hispanic into race/ethnicity ###
table(full$RACE, full$HISPANIC, useNA = "always")

full2 <- full %>% 
  mutate(race_ethnicity = case_when(RACE == "White" & HISPANIC == "Non-Hispanic" ~ "Non-Hispanic White",
                                    RACE == "Black" & HISPANIC == "Non-Hispanic" ~ "Non-Hispanic Black",
                                    RACE == "Other" & HISPANIC == "Non-Hispanic" ~ "Non-Hispanic Other",
                                    str_detect(HISPANIC, "Hispanic") & !is.na(RACE) ~ "Hispanic"))
table(full2$race_ethnicity, full2$RACE, useNA = "always")
table(full2$race_ethnicity, full2$HISPANIC, useNA = "always")


### collapse education ###

table(full2$DEGREE_collapsed, useNA = "always")
full3 <- full2 %>% 
  mutate(DEGREE_collapsed = case_when(DEGREE_collapsed == "No degree" ~ "No degree",
                                      DEGREE_collapsed == "GED/HS diploma" ~ "GED/HS diploma",
                                      TRUE ~ "College and more than college"))
table(full3$DEGREE_collapsed)
  
# Setting reference levels for categorical variables 

full3$race_ethnicity <- relevel(factor(full3$race_ethnicity), ref = "Non-Hispanic White") 
full3$DEGREE_collapsed <- relevel(factor(full3$DEGREE_collapsed), ref = "No degree") 
full3$PMARST_collapsed <- relevel(full3$PMARST_collapsed, ref = "Married") 

full3$race_ethnicity <- factor(full3$race_ethnicity, levels = c("Non-Hispanic White", "Non-Hispanic Black", 
                                                              "Non-Hispanic Other", "Hispanic"))

full3$DEGREE_collapsed <- factor(full3$DEGREE_collapsed, levels = c("No degree", "GED/HS diploma", "College and more than college"))


######### Adding physical activity and financial situation before age 16

## add covariate physical activities
physical_activity = c("R13VGACTX", "R13MDACTX", "R13LTACTX")

w13_clocks_phy=rand %>% 
  select(HHID,PN, all_of(physical_activity)) %>% 
  right_join(full3,by = c("HHID", "PN"))

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

table(w13_clocks_phy_we$R13ACT, useNA = "always")

full4 <- w13_clocks_phy_we

### filter the missing in covariate set one by one ###
nomiss_1 = full4 %>% 
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
  filter(!is.na(gran))
loss_gran = nrow(nomiss_8) - nrow(nomiss_9)

nomiss_10 = nomiss_9 %>% 
  filter(!is.na(PMONO))
loss_mono = nrow(nomiss_9) - nrow(nomiss_10)

nomiss_11 = nomiss_10 %>% 
  filter(!is.na(PLYMP))
loss_LYMP = nrow(nomiss_10) - nrow(nomiss_11)

nomiss_12 = nomiss_11 %>% 
  filter(!is.na(R13ACT))
loss_physical = nrow(nomiss_11) - nrow(nomiss_12)


######## create complete analytic sample #########

full_bo = full4 %>% 
  mutate(race_ethnicity = if_else(race_ethnicity == "Non-Hispanic Black"|race_ethnicity == "Non-Hispanic Other", "NH-Black or Other", as.character(race_ethnicity)))
table(full_bo$race_ethnicity)
full_bo$race_ethnicity = as.factor(full_bo$race_ethnicity)
full_bo$race_ethnicity = relevel(full_bo$race_ethnicity, ref = "Non-Hispanic White")


complete <- full_bo %>% 
  filter(complete.cases(.)) %>% 
  filter(RACE != "Not obtained") %>%
  filter(PMARST_collapsed != "Marital status unknown") %>%
  filter(HISPANIC != "Not obtained") %>% 
  mutate(included = 1)


full <- complete %>% 
  select(HHID, PN, included) %>% 
  full_join(full4, by = c("HHID", "PN")) %>% 
  mutate(included = case_when(included == 1~ "Included",
                              TRUE ~ "Excluded"))
table(full$included, useNA = "always")

#### SAVE NEW COMPLETE AND FULL DATASET ####

save(full, file = "full_4018.RData")
save(complete, file = "analytic_3882.RData")


