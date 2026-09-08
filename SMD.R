setwd("S:/HI-SPEED projects/HS26_02 GLP1 pancreatitis/")
library(tableone)
library(haven)
library(dplyr)
library(survey)
library(stringr)
library(ggplot2)
library(tibble)

# p_load(Hmisc, tidyverse, tableone, haven, kableExtra, ggplot2, 
#                survivalAnalysis, skimr, survminer, broom, haven, sjPlot,
#                dplyr, epitools, survival, lubridate, janitor, writexl, table1,
#                magrittr, survey, twang)

data.w=read_dta("data/ana_data_iptw.dta")


ps.var <- c("age_index", "age_g", "gender", "indexyr", "bc", "region_g", "education3", "fam_panc", "t2donsetyr_g",
	"bmi_6y_g", "bp_sys_6y_g", "bp_dias_6y_g", "gfr_6y_g", "hba1c_6y_g", "triglyceride_6y_g", "ldl_6y_g",
	"smokinghabit_6y_g", "pre_a10a", "pre_a10ba02", "pre_a10bb", "pre_a10bg", "pre_c02", "pre_c03", "pre_c07",
  "pre_c08", "pre_c09", "pre_c10", "pre_alco_dis", "pre_cereb_dis", "pre_cholelithiasis", "pre_diab_nephro",
  "pre_diab_retino", "pre_heart_fail", "pre_isc_heart_dis", "pre_jxx00")
ps.var.cat  <- c("age_g", "gender", "indexyr", "bc", "region_g", "education3", "fam_panc", "t2donsetyr_g",
                 "bmi_6y_g", "bp_sys_6y_g", "bp_dias_6y_g", "gfr_6y_g", "hba1c_6y_g", "triglyceride_6y_g", "ldl_6y_g",
                 "smokinghabit_6y_g", "pre_a10a", "pre_a10ba02", "pre_a10bb", "pre_a10bg", "pre_c02", "pre_c03", "pre_c07",
                 "pre_c08", "pre_c09", "pre_c10", "pre_alco_dis", "pre_cereb_dis", "pre_cholelithiasis", "pre_diab_nephro",
                 "pre_diab_retino", "pre_heart_fail", "pre_isc_heart_dis", "pre_jxx00")

################ GVD part ##############
data.gvd <- data.w %>% filter(expo != "sglt2")
tab1.gvd <- as.data.frame(print(CreateTableOne(vars = ps.var, data = data.gvd, factorVars = ps.var.cat, strata="exposure"), 
                            smd=TRUE, printToggle = FALSE)) %>% 
  rownames_to_column("row_name") %>% dplyr::select(-c(p,test))
colnames(tab1.gvd) <- c("var", "dpp4_before", "glp1_before", "SMD_before")

design.gvd.ATT   <- svydesign(ids=~1, weights=data.gvd$iptw_smr, data=data.gvd)
tab1.w.gvd <- as.data.frame(print(svyCreateTableOne(vars = ps.var, data = design.gvd.ATT, factorVars = ps.var.cat, 
                                                strata="exposure"), smd=TRUE, printToggle = FALSE)) %>% 
  rownames_to_column("row_name") %>% dplyr::select(-c(p,test))
colnames(tab1.w.gvd) <- c("var", "dpp4_after", "glp1_after", "SMD_after")

gvdtab <- tab1.gvd %>% full_join(tab1.w.gvd, by="var")


################ GVS part ##############
data.gvs <- data.w %>% filter(expo != "dpp4")
tab1.gvs <- as.data.frame(print(CreateTableOne(vars = ps.var, data = data.gvs, factorVars = ps.var.cat, strata="exposure"), 
                            smd=TRUE, printToggle = FALSE)) %>% 
  rownames_to_column("row_name") %>% dplyr::select(-c(p,test))
colnames(tab1.gvs) <- c("var","glp1_before", "sglt2_before", "SMD_before")

design.gvs.ATT   <- svydesign(ids=~1, weights=data.gvs$iptw_smr, data=data.gvs)
tab1.w.gvs <- as.data.frame(print(svyCreateTableOne(vars = ps.var, data = design.gvs.ATT, factorVars = ps.var.cat, 
                                                strata="exposure"), smd=TRUE, printToggle = FALSE)) %>% 
  rownames_to_column("row_name") %>% dplyr::select(-c(p,test))
colnames(tab1.w.gvs) <- c("var", "glp1_after", "sglt2_after", "SMD_after")

gvstab <- tab1.gvs %>% full_join(tab1.w.gvs, by="var")


######  output results
writexl::write_xlsx(
  list("GLP1_DPP4" = gvdtab,
       "GLP1_SGLT2" = gvstab), 
  path = "output/20260825/SMD calcu.xlsx")
