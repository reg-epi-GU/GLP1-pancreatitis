options nofmterr validvarname=any compress=yes spool;

libname im "E:\scifi\data\2.0.intermediate\2025-12";
libname sos "E:\scifi\data\1.0.original\2025-12\sos";
libname scb "E:\scifi\data\1.0.original\2025-12\scb";
libname proj "S:\HI-SPEED projects\HS26_02 GLP1 pancreatitis\data";

proc delete lib=work data=_all_; run;

/* Run this part only once;
data proj.a10 proj.c10 proj.c02 proj.c03 proj.c07 proj.c08 proj.c09;
set im.im_010_clean_lmr(keep=lopnr edatum atc);
if atc in: ("A10") then output proj.a10; 
else if atc in: ("C02") then output proj.c02;
else if atc in: ("C03") then output proj.c03;
else if atc in: ("C07") then output proj.c07;
else if atc in: ("C08") then output proj.c08;
else if atc in: ("C09") then output proj.c09;
else if atc in: ("C10") then output proj.c10;
run;

data proj.pancreatitis_hosp;
set im.im_020_clean_par(keep=lopnr indatum sv dia1-dia30);
where sv=1;
array adia dia1-dia30;
do over adia;
if adia in: ("K850" "K851" "K852" "K853" "K858" "K859") then k=1;
end;
if k;
drop k;
rename indatum=pancreatitis_hosp_date;
run;

data proj.pancreatitis_cod;
set im.Im_060_deathinfo;
array cod ulorsak morsak1-morsak48;
do over cod;
if cod in: ("K850" "K851" "K852" "K853" "K858" "K859") then k=1;
end;
if k;
rename deathdate=pancreatitis_cod_date;
drop k;
run;

data proj.r1_cod(keep=lopnr r1_cod_date);
set im.Im_060_deathinfo;
array cod ulorsak morsak1-morsak48;
do over cod;
if cod in: ("R1") then k=1;
end;
if k;
rename deathdate=r1_cod_date;
run;

data jxx00;
set im.im_020_clean_par(keep=lopnr indatum op1-op30 sv);
array aop op1-op30;
do over aop;
if prxmatch('/^J[A-Za-z]{2}/', aop) then jxx00=1;
end;
if jxx00;
drop op1-op30;
run;

data diab_nephro diab_retino isc_heart_dis cereb_dis heart_fail alco_dis cholelithiasis;
merge study_pop(keep=lopnr indexdate in=a)
	im.im_020_clean_par(keep=lopnr indatum dia1-dia30 sv);
by lopnr; if a;
if indexdate-365.25*6<=indatum<=indexdate;
array adia dia1-dia30;
do over adia;
if adia in: ("E112") then diab_nephro=1;
if adia in: ("E113") then diab_retino=1;
if adia in: ("I20" "I21" "I22" "I23" "I24" "I25") then isc_heart_dis=1;
if adia in: ("I60" "I61" "I62" "I63" "I64" "I65" "I66" "I67" "I68" "I69") then cereb_dis=1;
if adia in: ("I110" "I130" "I132" "I50") then heart_fail=1;
if adia in: ("F10") then alco_dis=1;
if adia in: ("K80") then cholelithiasis=1;
end;
if diab_nephro then output diab_nephro;
if diab_retino then output diab_retino;
if isc_heart_dis then output isc_heart_dis;
if cereb_dis then output cereb_dis;
if heart_fail then output heart_fail;
if alco_dis then output alco_dis;
if cholelithiasis then output cholelithiasis;
drop dia1-dia30;
run;

%macro v(v);
data proj.&v;
set &v;
keep lopnr indatum &v sv;
run;
%mend;

%v(jxx00);
%v(diab_nephro);
%v(diab_retino);
%v(isc_heart_dis);
%v(cereb_dis);
%v(heart_fail);
%v(alco_dis);
%v(cholelithiasis);

data proj.ndr;
set sos.Sos_ndr_2020_04085_2025_12_09;
proc sort; by lopnr; run;

data proj.c2223 proj.c25;
set sos.Ut_r_can_26256_2025(keep=lopnr diadat icdo3);
if substr(diadat,1,4)>2014;
dia_dat=input(diadat, yymmdd8.); format dia_dat yymmdd10.;
if icdo3 in: ("C22" "C23") then output proj.c2223;
if icdo3 in: ("C25") then output proj.c25;
run;
*/

proc copy in=proj out=work;
select a10 c10 c02 c03 c07 c08 c09 pancreatitis_hosp pancreatitis_cod c2223 c25 
	jxx00 ndr diab_nephro diab_retino isc_heart_dis cereb_dis heart_fail alco_dis cholelithiasis;
run;

data glp1 dpp4 sglt2 other_a10;
set a10;
if atc in: ("A10BJ" "A10BX16") then output glp1; else
if atc in: ("A10BH") then output dpp4; else
if atc in: ("A10BK") then output sglt2; else
output other_a10;
run;

data A10A A10BA02 A10BB A10BG;
set other_a10;
if atc in: ("A10A") then output A10A; else
if atc in: ("A10BA02") then output A10BA02; else
if atc in: ("A10BB") then output A10BB; else
if atc in: ("A10BG") then output A10BG; 
run;


%macro init(d);
data &d._1;
set &d;
where edatum>="1jan2021"d;
by lopnr; 
retain fd;
if first.lopnr then fd=edatum; else drop=(fd<edatum); 
if drop=1 then delete;
rename edatum=&d._initdate
	atc = &d._initatc;
drop fd drop;
proc sort nodup; by lopnr &d._initatc; run;

proc transpose data=&d._1 out=&d._2(drop=_name_ _label_) prefix=&d._initatc;
var &d._initatc;
by lopnr &d._initdate;
run;

data &d._3;
set &d._2;
&d._atc=catx("_", of &d._initatc:);
drop &d._initatc:;
run;
%mend;

%init(glp1);
%init(dpp4);
%init(sglt2);

data expodrug;
merge glp1_3 dpp4_3 sglt2_3;
by lopnr; 
if glp1_initdate ne . and 
(glp1_initdate<dpp4_initdate or dpp4_initdate=.) and 
(glp1_initdate<sglt2_initdate or sglt2_initdate=.) then do;
	indexdate=glp1_initdate; expo="glp1 "; end; 
else if dpp4_initdate ne . and 
(dpp4_initdate<glp1_initdate or glp1_initdate=.) and 
(dpp4_initdate<sglt2_initdate or sglt2_initdate=.) then do;
	indexdate=dpp4_initdate; expo="dpp4 "; end;
else if sglt2_initdate ne . and 
(sglt2_initdate<glp1_initdate or glp1_initdate=.) and 
(sglt2_initdate<dpp4_initdate or dpp4_initdate=.) then do;
	indexdate=sglt2_initdate; expo="sglt2"; end;
if indexdate ne .; format indexdate yymmdd10.;
run;

data pre3y_drug;
set glp1 dpp4 sglt2;
proc sort; by lopnr edatum; run;

data pre3y_drug;
merge pre3y_drug expodrug(keep=lopnr indexdate in=a);
by lopnr; if a;
if indexdate-3*365<= edatum < indexdate;
pre3y_3drug=1;
keep lopnr pre3y_3drug;
proc sort nodup; by lopnr; run;

*** drug discontinuation;
%macro discont(drug);
data &drug._1;
merge pop_&drug.(in=a) &drug.;
by lopnr; if a;
if edatum>=indexdate;
run; 

data &drug._2;
set &drug._1(keep=lopnr edatum atc);
retain epi pre_edatum;
by lopnr;
if first.lopnr then do;
	epi=1; pre_edatum=edatum; end;
else do;
	if edatum<=pre_edatum+180 then do; 
		pre_edatum=edatum; end;
	else do; epi=epi+1; pre_edatum=edatum; end;
end; 
format pre_edatum yymmdd10.;
run;

proc sql;
create table discon_&drug. as select
lopnr, max(edatum)+180 as discon_date
from &drug._2 where epi=1
group by lopnr;
quit;
%mend;

%discont(glp1);
%discont(dpp4);
%discont(sglt2);

data discont;
set discon_:;
proc sort; by lopnr; run;
proc sort nodupkey; by lopnr; run;

*** DM diagnosis from NDR;
data ndr_preindex;
merge ndr(keep=lopnr r_contactdate r_diabetestype) expodrug(in=a keep=lopnr indexdate);
by lopnr; if a ;
ndr_date=input(r_contactdate, yymmdd10.);
if indexdate-365.25*6<=ndr_date<=indexdate;
format ndr_date yymmdd10.;
*drop r_contactdate;
run;

**** Keep only those with consistent/uniform T2D records;
proc sql;
create table ndr_indi as
select 
    lopnr,
    sum(r_diabetestype = "1") as N_t1d,
    sum(r_diabetestype = "2") as N_t2d
from ndr_preindex
group by lopnr;
quit;

data ndr_date;
merge ndr_indi(in=a where=(N_t1d=0)) ndr_preindex(keep=lopnr ndr_date);
by lopnr; if a;
keep lopnr ndr_date;
proc sort; by lopnr ndr_date; run;

data t2d_indication_ndr;
set ndr_date( rename=(ndr_date=t2ddate_ndr));
proc sort nodupkey; by lopnr; run;

data study_pop;
merge expodrug(in=a) pre3y_drug discont t2d_indication_ndr;
by lopnr; if a;
t2d_indication_ndr=(t2ddate_ndr ne . and t2ddate_ndr<=indexdate);
if pre3y_3drug=. then pre3y_3drug=0;
run;

*** outcome - acute pancreatitis AP;
*** pancreatitis overall, even K852;
data pancreatitis_all_wK852;
set pancreatitis_hosp pancreatitis_cod;
panc_date=coalesce(pancreatitis_hosp_date, pancreatitis_cod_date); format panc_date yymmdd10.;
keep lopnr panc_date sv;
proc sort; by lopnr panc_date; run;

data panc_data_wK852;
merge study_pop(in=a keep=lopnr indexdate) pancreatitis_all_wK852;
by lopnr; if a;
if indexdate-365<=panc_date<indexdate then pre1y_panc=1;
run;

proc sql;
create table pre1y_panc_wK852 as select
lopnr, max(pre1y_panc) as pre1y_panc_wK852
from panc_data_wK852 group by lopnr;

create table out_panc_all_wK852 as select
lopnr, min(panc_date) as panc_all_wK852_date format yymmdd10., 1 as panc_all_wK852
from panc_data_wK852 where panc_date>=indexdate group by lopnr;

create table out_panc_hosp_wK852 as select
lopnr, min(panc_date) as panc_hosp_wK852_date format yymmdd10., 1 as panc_hosp_wK852
from panc_data_wK852 where panc_date>=indexdate and sv=1 group by lopnr;

create table out_panc_cod_wK852 as select
lopnr, min(panc_date) as panc_cod_wK852_date format yymmdd10., 1 as panc_cod_wK852
from panc_data_wK852 where panc_date>=indexdate and sv=. group by lopnr;
quit;

*** pancreatitis , without K852 (named as pancreatitis_all to align with the first plan;
data panc_hosp_woK852;
set pancreatitis_hosp;
array adia dia1-dia30;
do over adia;
if adia in: ("K850" "K851" "K853" "K858" "K859") then k=1;
end;
if k;
drop k;
run;

data panc_cod_woK852;
set pancreatitis_cod;
array cod ulorsak morsak1-morsak48;
do over cod;
if cod in: ("K850" "K851" "K853" "K858" "K859") then k=1;
end;
if k;
drop k;
run;

data pancreatitis_all_woK852;
set panc_hosp_woK852 panc_cod_woK852;
panc_date=coalesce(panc_hosp_date, panc_cod_date); format panc_date yymmdd10.;
keep lopnr panc_date sv;
proc sort; by lopnr panc_date; run;

data panc_data_woK852;
merge study_pop(in=a keep=lopnr indexdate) pancreatitis_all_woK852;
by lopnr; if a;
if indexdate-365<=panc_date<indexdate then pre1y_panc=1;
run;

proc sql;
create table pre1y_panc_woK852 as select
lopnr, max(pre1y_panc) as pre1y_panc_woK852
from panc_data_woK852 group by lopnr;

create table out_panc_all_woK852 as select
lopnr, min(panc_date) as panc_all_woK852_date format yymmdd10., 1 as panc_all_woK852
from panc_data_woK852 where panc_date>=indexdate group by lopnr;

create table out_panc_hosp_woK852 as select
lopnr, min(panc_date) as panc_hosp_woK852_date format yymmdd10., 1 as panc_hosp_woK852
from panc_data_woK852 where panc_date>=indexdate and sv=1 group by lopnr;

create table out_panc_cod_woK852 as select
lopnr, min(panc_date) as panc_cod_woK852_date format yymmdd10., 1 as panc_cod_woK852
from panc_data_woK852 where panc_date>=indexdate and sv=. group by lopnr;
quit;



*******************   BASIC study data is ready till here. Below are for covariates;


/*For the list and timewindow for covariates, 
please refer to Supplementary table S2. Definitions of exposures, outcomes, and covariates.*/

data edu;
merge study_pop(keep=lopnr indexdate in=a)
im.im_141_ses_2015(keep=lopnr education3_: healthcare_region_:)
im.im_141_ses_2016(keep=lopnr education3_: healthcare_region_:)
im.im_141_ses_2017(keep=lopnr education3_: healthcare_region_:)
im.im_141_ses_2018(keep=lopnr education3_: healthcare_region_:)
im.im_141_ses_2019(keep=lopnr education3_: healthcare_region_:)
im.im_141_ses_2020(keep=lopnr education3_: healthcare_region_:)
im.im_141_ses_2021(keep=lopnr education3_: healthcare_region_:)
im.im_141_ses_2022(keep=lopnr education3_: healthcare_region_:)
im.im_141_ses_2023(keep=lopnr education3_: healthcare_region_:)
im.im_141_ses_2024(keep=lopnr education3_: healthcare_region_:);
by lopnr; if a;
education3_2024=education3_2023;
indexyr=year(indexdate);
    array aedu education3_: ;   /* include all edu year variables */
	array hcr healthcare_region_: ;
    education3 = aedu[indexyr - 2015];      /* adjust base year */
	healthcare_region = hcr[indexyr - 2015];
if education3=. then education3=max(of education3_2015-education3_2024);
keep lopnr education3 healthcare_region;
run;

*** familiy history of pancreatitis;
	data foraldrar;
	set scb.Fn19_lev_foraldrar(keep=lopnr lopnr_mor lopnr_far);
	proc sort nodupkey; by lopnr; run;

	data fam;
	merge study_pop(keep=lopnr in=a)
		foraldrar;
	by lopnr; if a;
	run;

	proc sql;
	create table far_panc as select
	fam.lopnr, fam.lopnr_far, min(p.pancreatitis_date) as far_panc format yymmdd10.
	from fam left join pancreatitis_all_wK852 p on fam.lopnr_far=p.lopnr
	group by fam.lopnr;

	create table mor_panc as select
	fam.lopnr, fam.lopnr_mor, min(p.pancreatitis_date) as mor_panc format yymmdd10.
	from fam left join pancreatitis_all_wK852 p on fam.lopnr_mor=p.lopnr
	group by fam.lopnr;
	quit;

	proc sort data=far_panc nodupkey; by lopnr; run;
	proc sort data=mor_panc nodupkey; by lopnr; run;

data fam_panc;
merge study_pop(keep=lopnr indexdate in=a)
far_panc mor_panc;
by lopnr; if a;
fam_panc=((far_panc ne . and far_panc<indexdate) or (mor_panc ne . and mor_panc<indexdate));
run;

/*
Diabetic nephropathy	National Patient Register, 2015 to 2025	Any records of hospitalisation or specialised outpatient care visit in the preceding 6 years before index date	ICD-10 code E11.2
Diabetic retinopathy	National Patient Register, 2015 to 2025	Any records of hospitalisation or specialised outpatient care visit in the preceding 6 years before index date	ICD-10 code E11.3
Ischemic heart disease	National Patient Register, 2015 to 2025	Any records of hospitalisation or specialised outpatient care visit in the preceding 6 years before index date	ICD-10 codes I20, I21, I22, I23, I24, I25
Cerebrovascular disease	National Patient Register, 2015 to 2025	Any records of hospitalisation or specialised outpatient care visit in the preceding 6 years before index date	ICD-10 codes I60, I61, I62, I63, I64, I65, I66, I67, I68, I69
Heart failure	National Patient Register, 2015 to 2025	Any records of hospitalisation or specialised outpatient care visit in the preceding 6 years before index date	ICD-10 codes I11.0, I13.0, I13.2, I50
Alcohol use disorder	National Patient Register, 2015 to 2025	Any records of hospitalisation or specialised outpatient care visit in the preceding 6 years before index date	ICD-10 code F10
Cholelithiasis	National Patient Register, 2015 to 2025	Any records of hospitalisation or specialised outpatient care visit in the preceding 6 years before index date	ICD-10 code K80
Abdominal surgery	National Patient Register, 2015 to 2025	Any records of abdominal surgery in the preceding 6 years before index date	KKÅ code J
Biliary cancer	Cancer Register, 1958 to 2025	Any records of diagnosis in the preceding data before index date	ICD-10 codes C22, C23
Pancreatic cancer	Cancer Register, 1958 to 2025	Any records of diagnosis in the preceding data before index date	ICD-10 code C25
Lipid-lowering medications	Prescribed Drug Register, 2018 to 2025	Any dispensation in the preceding 1 year before index date	ATC code C10
Insulin	Prescribed Drug Register, 2018 to 2025	Any dispensation in the preceding 1 year before index date	ATC code A10A
Metformin	Prescribed Drug Register, 2018 to 2025	Any dispensation in the preceding 1 year before index date	ATC code A10BA02
Sulfonylureas	Prescribed Drug Register, 2018 to 2025	Any dispensation in the preceding 1 year before index date	ATC code A10BB
Thiazolidinediones	Prescribed Drug Register, 2018 to 2025	Any dispensation in the preceding 1 year before index date	ATC code A10BG
Antihypertensives	Prescribed Drug Register, 2018 to 2025	Any dispensation in the preceding 1 year before index date	ATC code C02
Diuretics	Prescribed Drug Register, 2018 to 2025	Any dispensation in the preceding 1 year before index date	ATC code C03
Beta-blocking agents	Prescribed Drug Register, 2018 to 2025	Any dispensation in the preceding 1 year before index date	ATC code C07
Calcium-channel blockers	Prescribed Drug Register, 2018 to 2025	Any dispensation in the preceding 1 year before index date	ATC code C08
Agents acting on the renin-angiotensin system	Prescribed Drug Register, 2018 to 2025	Any dispensation in the preceding 1 year before index date	ATC code C09
Systolic blood pressure	National Diabetes Register, 2015 to 2025	Most recent value in the preceding 1 year before index date	Measured, mmHg, quartiles
Diastolic blood pressure	National Diabetes Register, 2015 to 2025	Most recent value in the preceding 1 year before index date	Measured, mmHg, quartiles
Body mass index	National Diabetes Register, 2015 to 2025	Most recent value in the preceding 1 year before index date	Measured, kg/m2, WHO categories
HbA1c 	National Diabetes Register, 2015 to 2025	Most recent value in the preceding 1 year before index date	Measured, mmol/mol, quartiles
Estimated glomerular filtration rate	National Diabetes Register, 2015 to 2025	Most recent value in the preceding 1 year before index date	Estimated, mL/min/1.73m2, quartiles
Triglycerides	National Diabetes Register, 2015 to 2025	Most recent value in the preceding 1 year before index date	Measured, mmol/l, quartiles
Low-density lipoprotein cholesterol	National Diabetes Register, 2015 to 2025	Most recent value in the preceding 1 year before index date	Measured, mmol/l, quartiles
Smoking status	National Diabetes Register, 2015 to 2025	Most recent value in the preceding 6 years before index date	Self-reported, never/current/former
*/


*** covars from NDR
Systolic blood pressure
Diastolic blood pressure
Body mass index
HbA1c 
Estimated glomerular filtration rate
Triglycerides
Low-density lipoprotein cholesterol
Smoking status
;
proc sort data=ndr; by lopnr; run;
%macro infondr;
%do i=1 %to 3;
data ndr_&i.y;
merge study_pop(keep=lopnr indexdate in=a) 
	ndr(keep=lopnr r_contactdate r_bpsystolic r_bpdiastolic r_bmi 
	r_gfr r_hba1c r_triglyceride r_ldl);
by lopnr; if a;
date=input(r_contactdate, yymmdd10.); format date yymmdd10.;
if indexdate-365*&i.<=date<=indexdate and date ne .;
proc sort; by lopnr date; run;

data infondr_&i.y;
set ndr_&i.y;
by lopnr;
retain bmi_&i.y bp_sys_&i.y bp_dias_&i.y gfr_&i.y hba1c_&i.y triglyceride_&i.y ldl_&i.y;
if first.lopnr then do; bmi_&i.y=r_bmi;  bp_sys_&i.y=r_bpsystolic; bp_dias_&i.y=r_bpdiastolic;
	gfr_&i.y=r_gfr; hba1c_&i.y=r_hba1c; triglyceride_&i.y=r_triglyceride; ldl_&i.y=r_ldl; 
	end;
else do;
	if r_bmi="" then r_bmi=bmi_&i.y; else bmi_&i.y=r_bmi;
	if r_bpsystolic="" then r_bpsystolic=bp_sys_&i.y; else bp_sys_&i.y=r_bpsystolic;
	if r_bpdiastolic="" then r_bpdiastolic=bp_dias_&i.y; else bp_dias_&i.y=r_bpdiastolic;
	if r_gfr="" then r_gfr=gfr_&i.y; else gfr_&i.y=r_gfr;
	if r_hba1c="" then r_hba1c=hba1c_&i.y; else hba1c_&i.y=r_hba1c;
	if r_triglyceride="" then r_triglyceride=triglyceride_&i.y; else triglyceride_&i.y=r_triglyceride;
	if r_ldl="" then r_ldl=ldl_&i.y; else ldl_&i.y=r_ldl;
end;
if last.lopnr then output;
keep lopnr bmi_&i.y bp_sys_&i.y bp_dias_&i.y gfr_&i.y hba1c_&i.y triglyceride_&i.y ldl_&i.y;
run;
%end;
%mend;

%infondr;

data ndr_6y;
merge study_pop(keep=lopnr indexdate in=a) 
	ndr(keep=lopnr r_contactdate r_bpsystolic r_bpdiastolic r_bmi r_smokinghabit
	r_gfr r_hba1c r_triglyceride r_ldl);
by lopnr; if a;
date=input(r_contactdate, yymmdd10.); format date yymmdd10.;
if indexdate-365*6<=date<=indexdate and date ne .;
proc sort; by lopnr date; run;

data infondr_6y;
set ndr_6y;
by lopnr;
retain bmi_6y bp_sys_6y bp_dias_6y gfr_6y hba1c_6y triglyceride_6y ldl_6y smokinghabit_6y;
if first.lopnr then do; bmi_6y=r_bmi;  bp_sys_6y=r_bpsystolic; bp_dias_6y=r_bpdiastolic;
	gfr_6y=r_gfr; hba1c_6y=r_hba1c; triglyceride_6y=r_triglyceride; ldl_6y=r_ldl; smokinghabit_6y=r_smokinghabit;
	end;
else do;
	if r_bmi="" then r_bmi=bmi_6y; else bmi_6y=r_bmi;
	if r_bpsystolic="" then r_bpsystolic=bp_sys_6y; else bp_sys_6y=r_bpsystolic;
	if r_bpdiastolic="" then r_bpdiastolic=bp_dias_6y; else bp_dias_6y=r_bpdiastolic;
	if r_gfr="" then r_gfr=gfr_6y; else gfr_6y=r_gfr;
	if r_hba1c="" then r_hba1c=hba1c_6y; else hba1c_6y=r_hba1c;
	if r_triglyceride="" then r_triglyceride=triglyceride_6y; else triglyceride_6y=r_triglyceride;
	if r_ldl="" then r_ldl=ldl_6y; else ldl_6y=r_ldl;
	if r_smokinghabit="" then r_smokinghabit=smokinghabit_6y; else smokinghabit_6y=r_smokinghabit;
end;
if last.lopnr then output;
keep lopnr bmi_6y bp_sys_6y bp_dias_6y gfr_6y hba1c_6y triglyceride_6y ldl_6y smokinghabit_6y;
run;

data yr_onset;
merge study_pop(keep=lopnr in=a) 
ndr(keep=lopnr r_yearofonset);
by lopnr; if a and not missing(r_yearofonset);
r_yearofonset=r_yearofonset*1;
proc sort; by lopnr r_yearofonset; run;
proc sort nodupkey; by lopnr; run;
 

data infondr;
merge infondr_1y infondr_2y infondr_3y infondr_6y yr_onset;
by lopnr; 
run;


*** comorbidity;
%macro pre(code, date, washyear);
proc sort data=&code.; by lopnr &date.; run;
data pre_&code.;
merge study_pop(keep=lopnr indexdate in=a) &code.(in=b);
by lopnr; if a and b;
if indexdate - &washyear.*365 <= &date.<=indexdate;
if last.lopnr;
pre_&code.=1;
keep lopnr pre_&code.;
run;
%mend;

%pre(diab_nephro, indatum, 6);
%pre(diab_retino, indatum, 6);
%pre(isc_heart_dis, indatum, 6);
%pre(cereb_dis, indatum, 6);
%pre(heart_fail, indatum, 6);
%pre(alco_dis, indatum, 6);
%pre(cholelithiasis, indatum, 6);
%pre(jxx00, indatum, 1);
%pre(c2223,dia_dat,99);
%pre(c25,dia_dat,99);
%pre(c10,edatum,1);
%pre(a10a,edatum,1);
%pre(a10ba02,edatum,1);
%pre(a10bb,edatum,1);
%pre(a10bg,edatum,1);
%pre(c02,edatum,1);
%pre(c03,edatum,1);
%pre(c07,edatum,1);
%pre(c08,edatum,1);
%pre(c09,edatum,1);


** Merge to get all data;

data proj.study_data;
merge study_pop(in=a) out_panc_all_wk852 out_panc_all_wok852 
					  out_panc_hosp_wk852 out_panc_hosp_wok852
					  out_panc_cod_wk852 out_panc_cod_wok852
					  pre1y_panc_wk852 pre1y_panc_wok852 
edu fam_panc infondr pre_:
im.im_140_demographics(keep=lopnr fodelselandgrupp birthdate deathdate emig_date immig_date gender);
by lopnr; if a;
array change pre1y_pancreatitis fam_panc pre_:;
do over change; if change=. then change=0; end;
if deathdate ne . and deathdate<indexdate then delete;
age_index=(indexdate-birthdate)/365.25;
run;





***************  For Positive/Negative control outcomes;
/*********************************
K35, K36, K37, K38 for Appendicitis
M15, M16, M17, M18, M19 for Osteoarthritis
B02 for Herpes zoster
*********************************/

data pop;
set proj.study_data(keep=lopnr indexdate);
run;


data cout1 cout2 cout3;
set im.im_020_clean_par(keep=lopnr indatum sv dia1-dia30);
where sv=1;
array adia dia1-dia30;
do over adia;
if adia in: ("K35" "K36" "K37" "K38") then cout1=1;
if adia in: ("M15" "M16" "M17" "M18" "M19") then cout2=1;
if adia in: ("B02") then cout3=1;
end;
if cout1 then output cout1;
if cout2 then output cout2;
if cout3 then output cout3;
drop dia1-dia30;
run;

data cout1_cod cout2_cod cout3_cod;
merge pop(in=a) im.Im_060_deathinfo;
by lopnr; if a; 
array cod ulorsak morsak1-morsak48;
do over cod;
if cod in: ("K35" "K36" "K37" "K38") then cout1=1;
if cod in: ("M15" "M16" "M17" "M18" "M19") then cout2=1;
if cod in: ("B02") then cout3=1;
end;
if cout1 then output cout1_cod;
if cout2 then output cout2_cod;
if cout3 then output cout3_cod;
keep lopnr deathdate;
run;

%macro contrlout;
%do i=1 %to 3;
data cout&i._hosp;
merge pop(in=a) cout&i.(in=b);
by lopnr; if a and b;
keep lopnr indatum sv;
rename indatum=cout&i._hosp_date;
run;

data cout&i._all;
set cout&i._hosp cout&i._cod;
cout&i._date=coalesce(cout&i._hosp_date, deathdate); format cout&i._date yymmdd10.;
keep lopnr cout&i._date sv;
proc sort; by lopnr cout&i._date; run;

data cout&i._data;
merge pop(in=a keep=lopnr indexdate) cout&i._all;
by lopnr; if a;
if indexdate-365<=cout&i._date<indexdate then pre1y_cout&i.=1;
run;

proc sql;
create table out_pre1y_cout&i. as select
lopnr, max(pre1y_cout&i.) as pre1y_cout&i.
from cout&i._data group by lopnr;

create table out_cout&i._all as select
lopnr, min(cout&i._date) as cout&i._all_date format yymmdd10., 1 as cout&i._all
from cout&i._data where cout&i._date>=indexdate group by lopnr;

create table out_cout&i._hosp as select
lopnr, min(cout&i._date) as cout&i._hosp_date format yymmdd10., 1 as cout&i._hosp
from cout&i._data where cout&i._date>=indexdate and sv=1 group by lopnr;

create table out_cout&i._cod as select
lopnr, min(cout&i._date) as cout&i._cod_date format yymmdd10., 1 as cout&i._cod
from cout&i._data where cout&i._date>=indexdate and sv=. group by lopnr;
quit;
%end;
%mend;

%contrlout;

data proj.control_outcomes;
merge out_cout: out_pre1y_cout:;
by lopnr;
proc sort; by lopnr; run;



**** To get AP subtypes;
data pop;
set proj.study_data;
if pancreatitis_all_wk852 =1 ;
if age_index>=18 & t2d_indication_ndr=1 & pre3y_3drug=0 & pre1y_pancreatitis_wk852=. & 
pre_c2223 =0 & pre_c25=0;
keep lopnr pancreatitis: ;
run;

data pop_panhosp;
merge pop(in=a keep=lopnr) proj.pancreatitis_hosp(in=b);
by lopnr; if a and b;
if first.lopnr;
array adia dia:;
do over adia;
if adia not in: ("K850" "K851" "K852" "K853" "K858" "K859") then adia="";
end;
panhosp_dia=catx("_", of dia1-dia30);
if panhosp_dia="K858_K858" then panhosp_dia="K858";
if panhosp_dia="K859_K859" then panhosp_dia="K859";
run;

data pop_pancod;
merge pop(in=a keep=lopnr) proj.pancreatitis_cod(in=b);
by lopnr; if a and b;
if first.lopnr;
array adia ulorsak morsak1-morsak48;
do over adia;
if adia not in: ("K850" "K851" "K852" "K853" "K858" "K859") then adia="";
end;
pancod_dia=catx("_", of ulorsak morsak1-morsak48);
if pancod_dia="K858_K858" then pancod_dia="K858";
if pancod_dia="K859_K859" then pancod_dia="K859";
run;

data pop_pan;
merge pop_panhosp(keep=lopnr pan:)  pop_pancod(keep=lopnr pan:);
by lopnr;
pantype=coalescec(panhosp_dia, pancod_dia);
keep lopnr pantype;
run;

data proj.study_data;
merge proj.study_data pop_pan;
by lopnr;
run;

proc freq data=pop_pan;
table pantype;
run;
