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

df_list = list(w16, w14, w12, w10, w08, w06, w04, w02, w00, w98)
var_list = c("PB020", "OB020", 'NB020', 'MB020', 'LB020', 'KB020', 'JB020', 'HB020', 'G1080','F993')

w_merge = df_list %>% 
  reduce(full_join, by = c('HHID', "PN"))

head(w_merge)

w13_merge2 = w13_merge1 %>% 
  left_join(w_merge, by = c('HHID', "PN")) %>% 
  select(HHID,PN, .dots = var_list)

head(w13_merge2)

w13_merge3 = w13_merge2 %>% 
  mutate(child_poverty = rowMeans(.[,c(3:12)], na.rm = T))

table(w13_merge3$child_poverty, useNA = "always")

class(w13_merge3$child_poverty)

w13_merge3 %>% 
  filter(child_poverty == 2 | child_poverty == 4 |child_poverty == 4.5)
