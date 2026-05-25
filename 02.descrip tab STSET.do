*** 02 Descriptive tables

cd "S:\HI-SPEED projects\HS26_02 GLP1 pancreatitis\data\"

use study_data, clear
keep if t2d_indication_ndr==1 & pre3y_3drug==0 & pre1y_pancreatitis==0 & preindex_biliarycancer==0 & preindex_pancreacancer==0 & healthcare_region!=""
encode expo, gen(exposure)
gen bc=2
replace bc=1 if fodelselandgrupp=="Sverige"
label define bc 1 "Sverige" 2 "Other"
label values bc bc
encode healthcare_region, gen(region_g)
destring r_yearofonset, gen(yearofonset)
gen t2dyr=year(t2ddate)
egen t2donsetyr=rowmin(yearofonset t2dyr)
egen t2donsetyr_g = cut(t2donsetyr), at(0, 2015, 2019, 2026) label
egen age_g=cut(age), at(0,40, 60,70,80,190) label
gen indexyr=year(indexdate)

destring bmi, gen(bmi_num)
egen bmi_g=cut(bmi_num), at(0,18.5,25,30,999) label
replace bmi_g=. if bmi_g==0

destring smokinghabit, gen(smokinghabit_g)
recode smokinghabit_g (3=2) (4=3)
label define smokinghabit_g 1 "Never smoker" 2 "Current smoker" 3 "Former smoker", modify
label values smokinghabit_g smokinghabit_g

destring hba1c, gen(hba1c_num)
egen hba1c_g=cut(hba1c_num), group(4) label

destring triglyceride, gen(triglyceride_num)
egen triglyceride_g=cut(triglyceride_num), group(4) label

egen covarmiss=rowmiss(age gender indexyr bc region_g education3 fam_panc t2donsetyr_g bmi_g smokinghabit_g hba1c_num preindex_lipidlowerdrug preindex_alcoholdisorder preindex_gallstones preindex_abdomsurgery)
gen complete=(covarmiss==0)

replace education3 = 99 if education3==.
label define education3 1 "primary" 2 "secondary" 3 "tertiary" 99 "missing", modify
label values education3 education3
foreach v in bmi smokinghabit hba1c triglyceride{
	replace `v'_g=99 if `v'_g==.
	label define `v'_g 99 "missing", modify
}



save ana_data, replace

*dtable age i.age_g i.gender i.indexyr i.bc i.region_g i.education3 i.fam_panc i.t2donsetyr_g bmi_num i.bmi_g i.smokinghabit_g hba1c_num i.hba1c_g /*triglyceride_num i.triglyceride_g*/ i.preindex_lipidlowerdrug i.preindex_alcoholdisorder i.preindex_gallstones i.preindex_abdomsurgery i.complete i.pancreatitis_all i.pancreatitis_hosp i.pancreatitis_cod, by(exposure) export("S:\HI-SPEED projects\HS26_02 GLP1 pancreatitis\output\GLP1 pancreatitis T1.xlsx", replace)


dtable age i.age_g i.gender i.indexyr i.bc i.region_g i.education3 i.fam_panc i.t2donsetyr_g bmi_num i.bmi_g i.smokinghabit_g hba1c_num i.hba1c_g /*triglyceride_num i.triglyceride_g*/ i.preindex_lipidlowerdrug i.preindex_alcoholdisorder i.preindex_gallstones i.preindex_biliarycancer i.preindex_pancreacancer i.preindex_abdomsurgery i.complete i.pancreatitis_all i.pancreatitis_hosp i.pancreatitis_cod if t2d_indication_ndr==1, by(exposure) export("S:\HI-SPEED projects\HS26_02 GLP1 pancreatitis\output\GLP1 pancreatitis T1_NDRcaseonly.xlsx", replace)

collect clear
table (glp1_atc) (exposure) if exposure==2 & t2d_indication_ndr==1
table (dpp4_atc) (exposure) if exposure==1 & t2d_indication_ndr==1, append
table (sglt2_atc) (exposure) if exposure==3 & t2d_indication_ndr==1, append
collect layout (glp1_atc dpp4_atc sglt2_atc) (exposure)
collect export "S:\HI-SPEED projects\HS26_02 GLP1 pancreatitis\output\GLP1 pancreatitis T1_NDRcaseonly.xlsx", sheet(drug_class, replace) modify