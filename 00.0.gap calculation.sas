options nofmterr validvarname=any compress=yes spool;

libname im "E:\scifi\data\2.0.intermediate\2025-12";
libname imnew "E:\scifi\data\2.0.intermediate\2026-03";
libname sos "E:\scifi\data\1.0.original\2025-12\sos";
libname scb "E:\scifi\data\1.0.original\2025-12\scb";
libname proj "S:\HI-SPEED projects\HS26_02 GLP1 pancreatitis\data";

data pop_glp1 pop_dpp4 pop_sglt2;
set proj.study_data;
if age_index>=18 and t2d_indication_ndr=1 and pre3y_3drug=0 and pre1y_pancreatitis=0 and pre_c2223=0 and pre_c25=0;
if expo="glp1" then output pop_glp1;
if expo="dpp4" then output pop_dpp4;
if expo="sglt2" then output pop_sglt2;
keep lopnr expo indexdate;
run;

%macro gap(d);
data &d._gap0;
merge pop_&d.(in=a) &d.;
by lopnr; if a and edatum>= indexdate;
run;

data &d._gap;
set &d._gap0;
by lopnr;
prev_edatum=lag(edatum);
format prev_edatum yymmdd10.;
gap=edatum-prev_edatum;
if first.lopnr then do; prev_edatum=.; gap=.; end;
keep lopnr edatum prev_edatum gap;
run;

data &d._gap;
set &d._gap;
gap_m=floor(gap/30);
run;

proc freq data=&d._gap;
table gap_m;
run;

proc univariate data=&d._gap noprint;
histogram gap / barwidth=15;
run;
%mend;

%gap(glp1);
%gap(dpp4);
%gap(sglt2);


***** Define discontinuation;

/*
data glp1 dpp4 sglt2;
set imnew.im_010_clean_lmr_enddate;
if atc in: ("A10BJ" "A10BX16") then output glp1; else
if atc in: ("A10BH") then output dpp4; else
if atc in: ("A10BK") then output sglt2; 
run;
*/

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
