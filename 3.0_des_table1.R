# Project PFAS-ome pilot, metabolomics
# Sister Study variables Table 1

# Table 1 ----
library(table1)
library(tidyverse)
library(Hmisc)

SSdat2 <- read_csv("clean_data/PFAS_ome_upload/SS.csv")

tableone <- SSdat2 %>%
  select(SS_age_T0, SS_race,  SS_edu, SS_income_T0, SS_menostatus_T0, SS_BMI_T0) %>%
  mutate_at(vars(SS_race,  SS_edu, SS_income_T0, SS_menostatus_T0), as.factor) 


tableone2 <- tableone %>%
  mutate(SS_race = factor(SS_race, levels = c(0, 1, 2), 
                          labels = c("Non-Hisp White", "Non-Hisp Black", "Hispanic")),
         SS_edu = factor(SS_edu, levels = c(0, 1, 2), 
                         labels = c("High school or less", "Some college", "College and above")),
         SS_income_T0 = factor(SS_income_T0, levels = c(1, 2, 3, 4, 5), 
                               labels = c("<20,000", "20,000-49,999",
                                          "50,000-99,999", "100,000-200,000",">200,000")), 
         SS_menostatus_T0 = factor(SS_menostatus_T0, levels = c("0", "1"), 
                                      labels = c("Premenopause", "Postmenopause")))


label(tableone2$SS_age_T0) <- "Age (years)"
label(tableone2$SS_income_T0) <- "Annual income at baseline"
label(tableone2$SS_edu) <- "Educational attainmentat baseline"
label(tableone2$SS_BMI_T0) <- "BMI at baseline (kg/m2)"
label(tableone2$SS_race) <- "Race and ethnicity"
label(tableone2$SS_menostatus_T0) <- "Menopausal status at baseline"


(tab1 = table1(~., data = tableone2,
               render.continuous=c(.="Median [IQR]")))




