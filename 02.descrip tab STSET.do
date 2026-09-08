*** 02 Descriptive tables

cd "S:\HI-SPEED projects\HS26_02 GLP1 pancreatitis\data\"
global resultpath "S:\HI-SPEED projects\HS26_02 GLP1 pancreatitis\output\20260825"

use study_data, clear
keep if age_index>=18 & t2d_indication_ndr==1 & pre3y_3drug==0 & pre1y_panc_wk852==0 & preindex_biliarycancer==0 & preindex_pancreacancer==0
encode expo, gen(exposure)
label define exposure 1 "DPP-4" 2 "GLP-1" 3 "SGLT-2", modify
gen bc=2
replace bc=1 if fodelselandgrupp=="Sverige"
label define bc 1 "Sverige" 2 "Other"
label values bc bc
encode healthcare_region, gen(region_g)
replace region_g=99 if region_g==.
label define region_g 99 "missing", modify

destring r_yearofonset, gen(yearofonset)
gen t2dyr=year(t2ddate)
egen t2donsetyr=rowmin(yearofonset t2dyr)
egen t2donsetyr_g = cut(t2donsetyr), at(0, 2015, 2017, 2019, 2020, 2021, 2022, 2023, 2024, 2025, 2026) label
egen age_g=cut(age_index), at(0,18, 40,60,70,80,190) label
gen indexyr=year(indexdate)

label variable cout1_all "Appendicitis"
label variable cout2_all "Osteoarthritis"
label variable cout3_all "Herpes zoster"

foreach i in 1 2 3 6{
	destring bmi_`i'y, gen(bmi_`i'y_num)
	egen bmi_`i'y_g=cut(bmi_`i'y_num), at(0,18.5,25,30,999) label
	replace bmi_`i'y_g=. if bmi_`i'y_g==0
}

foreach i in 1 2 3 6{
	destring bp_sys_`i'y, gen(bp_sys_`i'y_num)
	egen bp_sys_`i'y_g=cut(bp_sys_`i'y_num), at(0,120,130,140,999) label
	replace bp_sys_`i'y_g=. if bp_sys_`i'y_g==0
}

foreach i in 1 2 3 6{
	destring bp_dias_`i'y, gen(bp_dias_`i'y_num)
	egen bp_dias_`i'y_g=cut(bp_dias_`i'y_num), at(0,60,80,90,100,110,999) label
	replace bp_dias_`i'y_g=. if bp_dias_`i'y_g==0
}

foreach i in 1 2 3 6{
	destring gfr_`i'y, gen(gfr_`i'y_num)
	egen gfr_`i'y_g=cut(gfr_`i'y_num), at(0,15,30,60,90,999) label
	replace gfr_`i'y_g=. if gfr_`i'y_g==0
}

foreach i in 1 2 3 6{
	destring hba1c_`i'y, gen(hba1c_`i'y_num)
	egen hba1c_`i'y_g=cut(hba1c_`i'y_num), at(0,42,47,999) label
	replace hba1c_`i'y_g=. if hba1c_`i'y_g==0
}

foreach i in 1 2 3 6{
	destring triglyceride_`i'y, gen(triglyceride_`i'y_num)
	egen triglyceride_`i'y_g=cut(triglyceride_`i'y_num), at(0,1.7,2.2,5.6,999) label
	replace triglyceride_`i'y_g=. if triglyceride_`i'y_g==0
}

foreach i in 1 2 3 6{
	destring ldl_`i'y, gen(ldl_`i'y_num)
	egen ldl_`i'y_g=cut(ldl_`i'y_num), at(0,2.6, 3.3, 4.1, 4.9,999) label
	replace ldl_`i'y_g=. if ldl_`i'y_g==0
}

destring smokinghabit_6y, gen(smokinghabit_6y_g)

foreach i in 1 2 3 6{
	egen covarmiss_`i'y=rowmiss(age_index gender indexyr bc region_g education3 fam_panc t2donsetyr_g bmi_`i'y_num bp_sys_`i'y_num bp_dias_`i'y_num gfr_`i'y_num hba1c_`i'y_num triglyceride_`i'y_num ldl_`i'y_num smokinghabit_6y_g pre_a10a pre_a10ba02 pre_a10bb pre_a10bg pre_c02 pre_c03 pre_c07 pre_c08 pre_c09 pre_c10 pre_alco_dis pre_cereb_dis pre_cholelithiasis pre_diab_nephro pre_diab_retino pre_heart_fail pre_isc_heart_dis pre_jxx00)
	gen complete_`i'y=(covarmiss_`i'y==0)
}
label define complete 0 "No" 1 "Yes"
label values complete_* complete


recode smokinghabit_6y_g (3=2) (4=3) (.=99)
label define smokinghabit_6y_g 1 "Never smoker" 2 "Current smoker" 3 "Former smoker" 99 "missing", modify
label values smokinghabit_6y_g smokinghabit_6y_g

replace education3 = 99 if education3==.
label define education3 1 "primary" 2 "secondary" 3 "tertiary" 99 "missing", modify
label values education3 education3

foreach var in bmi bp_sys bp_dias gfr hba1c triglyceride ldl{
	foreach i in 1 2 3 6{
		replace `var'_`i'y_g=99 if `var'_`i'y_g==.
		label define `var'_`i'y_g 99 "missing", modify
	}
}

save ana_data, replace


****** Descriptive Table 1.
use ana_data, clear
dtable age_index i.age_g i.gender i.indexyr i.bc i.region_g i.education3 i.fam_panc i.t2donsetyr_g ///
	bmi_1y_num i.bmi_1y_g bmi_2y_num i.bmi_2y_g bmi_3y_num i.bmi_3y_g bmi_6y_num i.bmi_6y_g ///
	bp_sys_1y_num i.bp_sys_1y_g bp_sys_2y_num i.bp_sys_2y_g bp_sys_3y_num i.bp_sys_3y_g bp_sys_6y_num i.bp_sys_6y_g ///
	bp_dias_1y_num i.bp_dias_1y_g bp_dias_2y_num i.bp_dias_2y_g bp_dias_3y_num i.bp_dias_3y_g bp_dias_6y_num i.bp_dias_6y_g ///
	gfr_1y_num i.gfr_1y_g gfr_2y_num i.gfr_2y_g gfr_3y_num i.gfr_3y_g gfr_6y_num i.gfr_6y_g ///
	hba1c_1y_num i.hba1c_1y_g hba1c_2y_num i.hba1c_2y_g hba1c_3y_num i.hba1c_3y_g hba1c_6y_num i.hba1c_6y_g ///
	triglyceride_1y_num i.triglyceride_1y_g triglyceride_2y_num i.triglyceride_2y_g triglyceride_3y_num i.triglyceride_3y_g triglyceride_6y_num i.triglyceride_6y_g ///
	ldl_1y_num i.ldl_1y_g ldl_2y_num i.ldl_2y_g ldl_3y_num i.ldl_3y_g ldl_6y_num i.ldl_6y_g ///
	i.complete_1y i.complete_2y i.complete_3y i.complete_6y ///
	i.smokinghabit_6y_g i.pre_a10a i.pre_a10ba02 i.pre_a10bb i.pre_a10bg i.pre_c02 i.pre_c03 i.pre_c07 i.pre_c08 i.pre_c09 i.pre_c10 i.pre_alco_dis i.pre_cereb_dis i.pre_cholelithiasis i.pre_diab_nephro i.pre_diab_retino i.pre_heart_fail i.pre_isc_heart_dis i.pre_jxx00 i.panc_all_wk852 i.panc_hosp_wk852 i.panc_cod_wk852 i.panc_all_wok852 i.panc_hosp_wok852 i.panc_cod_wok852, by(exposure) export("${resultpath}\GLP1 pancreatitis T1.xlsx", modify)

collect clear
table (glp1_atc) (exposure) if exposure==2 & t2d_indication_ndr==1
table (dpp4_atc) (exposure) if exposure==1 & t2d_indication_ndr==1, append
table (sglt2_atc) (exposure) if exposure==3 & t2d_indication_ndr==1, append
collect layout (glp1_atc dpp4_atc sglt2_atc) (exposure)
collect export "${resultpath}\GLP1 pancreatitis T1.xlsx", sheet(drug_class, replace) modify