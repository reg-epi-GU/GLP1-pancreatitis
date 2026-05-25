options nofmterr validvarname=any compress=yes spool;

libname im "E:\scifi\data\2.0.intermediate\2025-12";
libname sos "E:\scifi\data\1.0.original\2025-12\sos";
libname scb "E:\scifi\data\1.0.original\2025-12\scb";
libname proj "S:\HI-SPEED projects\HS26_02 GLP1 pancreatitis\data";

proc delete lib=work data=_all_; run;

/* Run this part only once;
data proj.a10 proj.c10;
set im.im_010_clean_lmr(keep=lopnr edatum atc);
if atc in: ("A10") then output proj.a10; 
else if atc in: ("C10") then output proj.c10;
run;

data proj.pancreatitis_hosp;
set im.im_020_clean_par(keep=lopnr indatum sv dia1-dia30);
where sv=1;
array adia dia1-dia30;
do over adia;
if adia in: ("K850" "K853" "K858" "K859") then k=1;
end;
if k;
drop k;
rename indatum=pancreatitis_hosp_date;
run;

data proj.pancreatitis_cod;
set im.Im_060_deathinfo;
array cod ulorsak morsak1-morsak48;
do over cod;
if cod in: ("K850" "K853" "K858" "K859") then k=1;
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
set im.im_020_clean_par(keep=lopnr indatum op1-op30);
array aop op1-op30;
do over aop;
if prxmatch('/^J[A-Za-z]{2}/', aop) then jxx00=1;
end;
if jxx00;
drop op1-op30;
run;

data f10 k80 r1 e11;
set im.im_020_clean_par(keep=lopnr indatum sv dia1-dia30);
array adia dia1-dia30;
do over adia;
if adia in: ("F10") then f10=1;
if adia in: ("K80") then k80=1;
if adia in: ("R1") then r1=1;
if adia in: ("E11") then e11=1;
end;
if f10 then output f10;
if k80 then output k80;
if r1 then output r1;
if e11 then output e11;
drop dia1-dia30;
run;

%macro v(v);
data proj.&v;
set &v;
keep lopnr indatum &v sv;
run;
%mend;

%v(f10);
%v(k80);
%v(r1);
%v(e11);
%v(jxx00);

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
select a10 c10 pancreatitis_hosp pancreatitis_cod f10 k80 r1 r1_cod c2223 c25 e11 jxx00 ndr;
run;

data glp1 dpp4 sglt2;
set a10;
if atc in: ("A10BJ" "A10BX16") then output glp1; else
if atc in: ("A10BH") then output dpp4; else
if atc in: ("A10BK") then output sglt2;
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

**** Keep only those with consistent/uniform T2D records;
proc sql;
create table ndr_indi as
select 
    lopnr,
    sum(r_diabetestype = "1") as N_t1d,
    sum(r_diabetestype = "2") as N_t2d
from ndr
group by lopnr;
quit;

data ndr_date;
merge ndr_indi(in=a where=(N_t1d=0)) ndr(keep=lopnr r_contactdate);
by lopnr; if a;
ndr_date=input(r_contactdate, yymmdd10.);
format ndr_date yymmdd10.;
keep lopnr ndr_date;
proc sort; by lopnr ndr_date; run;

data other_a10;
set a10;
where atc not in: ("A10BJ" "A10BX16" "A10BH" "A10BK");
run;

data t2d_indication_ndr;
set ndr_date( rename=(ndr_date=t2ddate_ndr));
proc sort nodupkey; by lopnr; run;

data study_pop;
merge expodrug(in=a) pre3y_drug t2d_indication_ndr;
by lopnr; if a;
t2d_indication_ndr=(t2ddate_ndr ne . and t2ddate_ndr<=indexdate);
if pre3y_3drug=. then pre3y_3drug=0;
run;

*** outcome - acute pancreatitis AP;
data pancreatitis_all;
set pancreatitis_hosp pancreatitis_cod;
pancreatitis_date=coalesce(pancreatitis_hosp_date, pancreatitis_cod_date); format pancreatitis_date yymmdd10.;
keep lopnr pancreatitis_date sv;
proc sort; by lopnr pancreatitis_date; run;

data panc_data;
merge study_pop(in=a keep=lopnr indexdate) pancreatitis_all;
by lopnr; if a;
if indexdate-365<=pancreatitis_date<indexdate then pre1y_pancreatitis=1;
run;

proc sql;
create table pre1y_panc as select
lopnr, max(pre1y_pancreatitis) as pre1y_pancreatitis
from panc_data group by lopnr;

create table out_panc_all as select
lopnr, min(pancreatitis_date) as pancreatitis_all_date format yymmdd10., 1 as pancreatitis_all
from panc_data where pancreatitis_date>=indexdate group by lopnr;

create table out_panc_hosp as select
lopnr, min(pancreatitis_date) as pancreatitis_hosp_date format yymmdd10., 1 as pancreatitis_hosp
from panc_data where pancreatitis_date>=indexdate and sv=1 group by lopnr;

create table out_panc_cod as select
lopnr, min(pancreatitis_date) as pancreatitis_cod_date format yymmdd10., 1 as pancreatitis_cod
from panc_data where pancreatitis_date>=indexdate and sv=. group by lopnr;
quit;


*******************   BASIC study data is ready till here. Below are for covariates;


/*We adjusted our analysis for measured covariates which we assumed to 
be causally related to dispensation of the study drugs and incident 
pancreatitis (ie, confounders). These included age, sex, country of birth,
geographic region, education level, family history of pancreatitis, 
year of T2D onset, body mass index, smoking status, HbA1c levels,
triglycerides levels/dispensation of lipid-lowering medications,
alcohol use disorder, gallstones, biliary cancer, 
and pancreatic cancer (all at index date or most recent availability)
, as well as abdominal surgery in the 30 days prior to index date*/

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
	from fam left join pancreatitis_all p on fam.lopnr_far=p.lopnr
	group by fam.lopnr;

	create table mor_panc as select
	fam.lopnr, fam.lopnr_mor, min(p.pancreatitis_date) as mor_panc format yymmdd10.
	from fam left join pancreatitis_all p on fam.lopnr_mor=p.lopnr
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

*** covars from NDR;
	proc sort data=ndr; by lopnr; run;
	data infondr;
	merge study_pop(keep=lopnr indexdate in=a) 
	ndr(keep=lopnr r_contactdate r_yearofonset r_bmi r_smoker 
		r_smokinghabit r_smokingendyear r_hba1c r_triglyceride);
	by lopnr; if a;
	date=input(r_contactdate, yymmdd10.); format date yymmdd10.;
	if date<=indexdate and date ne .;
	proc sort; by lopnr date; run;

data infondr1;
set infondr;
by lopnr;
retain bmi smoker smokinghabit smokingendyear hba1c triglyceride;
if first.lopnr then do; bmi=r_bmi; 
	smoker=r_smoker; smokinghabit=r_smokinghabit; smokingendyear=r_smokingendyear; 
	hba1c=r_hba1c; triglyceride=r_triglyceride; end;
else do;
	if r_bmi="" then r_bmi=bmi; else bmi=r_bmi;
	if r_smoker="" then r_smoker=smoker; else smoker=r_smoker;
	if r_smokinghabit="" then r_smokinghabit=smokinghabit; else smokinghabit=r_smokinghabit;
	if r_smokingendyear="" then r_smokingendyear=smokingendyear; else smokingendyear=r_smokingendyear;
	if r_hba1c="" then r_hba1c=hba1c; else hba1c=r_hba1c;
	if r_triglyceride="" then r_triglyceride=triglyceride; else triglyceride=r_triglyceride;
end;
if last.lopnr then output;
keep lopnr r_yearofonset bmi smoker smokinghabit smokingendyear hba1c triglyceride;
run;

*** comorbidity;
%macro pre(code, date, washyear);
proc sort data=&code.; by lopnr &date.; run;
data pre_&code.;
merge study_pop(keep=lopnr indexdate in=a) &code.(in=b);
by lopnr; if a and b;
if indexdate - &washyear.*365 < &date.<=indexdate;
if last.lopnr;
keep lopnr &date.;
rename &date.=last_&code._preindex;
run;
%mend;

%pre(c10,edatum,3);
%pre(f10,indatum,6);
%pre(k80,indatum,6);
%pre(c2223,dia_dat,6);
%pre(c25,dia_dat,6);

proc sort data=jxx00; by lopnr indatum; run;
data pre_jxx00;
merge study_pop(keep=lopnr indexdate in=a) jxx00(in=b);
by lopnr; if a and b;
if indexdate-30<=indatum<=indexdate;
if last.lopnr;
keep lopnr indatum;
rename indatum=last_jxx00_preindex;
run;


** Merge to get all data;

data proj.study_data;
merge study_pop(in=a) out_panc_all out_panc_hosp out_panc_cod pre1y_panc 
out_gisymp_all out_gisymp_hosp out_gisymp_cod pre1y_gisymp
out_c19_all out_c19_hosp out_c19_cod pre1y_c19
edu fam_panc infondr1 pre_c10 pre_f10 pre_k80 pre_c2223 pre_c25 pre_jxx00
im.im_140_demographics(keep=lopnr fodelselandgrupp birthdate deathdate emig_date immig_date gender);
by lopnr; if a;
preindex_lipidlowerdrug=(last_c10_preindex ne . and last_c10_preindex<=indexdate);
preindex_alcoholdisorder=(last_f10_preindex ne . and last_f10_preindex<=indexdate);
preindex_gallstones=(last_k80_preindex ne . and last_k80_preindex<=indexdate);
preindex_biliarycancer=(last_c2223_preindex ne . and last_c2223_preindex<=indexdate);
preindex_pancreacancer=(last_c25_preindex ne . and last_c25_preindex<=indexdate);
preindex_abdomsurgery=(last_jxx00_preindex ne . and last_jxx00_preindex<=indexdate);
array change pre1y_pancreatitis gisymp_all gisymp_hosp gisymp_cod pre1y_gisymp c19_all c19_hosp c19_cod pre1y_c19 far_panc mor_panc fam_panc;
do over change; if change=. then change=0; end;
if deathdate ne . and deathdate<indexdate then delete;
run;





***************  For Positive/Negative control outcomes;
/*********************************
M10 for GOUT
B02 for herpes zoster
*********************************/

data pop;
set proj.study_data(keep=lopnr indexdate);
run;


data cout1 cout2;
set im.im_020_clean_par(keep=lopnr indatum sv dia1-dia30);
where sv=1;
array adia dia1-dia30;
do over adia;
if adia in: ("M10") then gout=1;
if adia in: ("B02") then herpeszoster=1;
end;
if gout then output cout1;
if herpeszoster then output cout2;
drop dia1-dia30;
run;

data cout1_cod cout2_cod;
merge pop(in=a) im.Im_060_deathinfo;
by lopnr; if a; 
array cod ulorsak morsak1-morsak48;
do over cod;
if cod in: ("M10") then gout=1;
if cod in: ("B02") then herpeszoster=1;
end;
if gout then output cout1_cod;
if herpeszoster then output cout2_cod;
keep lopnr deathdate;
run;

%macro contrlout;
%do i=1 %to 2;
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

data control_outcomes;
merge out_:;
by lopnr;
proc sort; by lopnr; run;

proc export data=control_outcomes
outfile="S:\HI-SPEED projects\HS26_02 GLP1 pancreatitis\data\controloutcomes.dta"
dbms=dta replace;
run;

**** To get AP subtypes;
data pop;
set proj.study_data;
if pancreatitis_all =1 ;
if t2d_indication_ndr=1 & pre3y_3drug=0 & pre1y_pancreatitis=0 & preindex_biliarycancer=0 & preindex_pancreacancer=0 & healthcare_region ne "";
keep lopnr pancreatitis: ;
run;

data pop_panhosp;
merge pop(in=a keep=lopnr) proj.pancreatitis_hosp(in=b);
by lopnr; if a and b;
if first.lopnr;
array adia dia:;
do over adia;
if adia not in: ("K850" "K853" "K858" "K859") then adia="";
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
if adia not in: ("K850" "K853" "K858" "K859") then adia="";
end;
pancod_dia=catx("_", of ulorsak morsak1-morsak48);
if pancod_dia="K858_K858" then pancod_dia="K858";
if pancod_dia="K859_K859" then pancod_dia="K859";
run;

data pop_pan;
merge pop_panhosp(keep=lopnr pan:)  pop_pancod(keep=lopnr pan:);
by lopnr;
pantype=coalescec(panhosp_dia, pancod_dia);
run;

proc freq data=pop_pan;
table pantype;
run;
