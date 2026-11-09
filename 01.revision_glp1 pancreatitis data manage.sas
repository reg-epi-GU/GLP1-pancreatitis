options nofmterr validvarname=any compress=yes spool;

libname im "E:\scifi\data\2.0.intermediate\2026-03";
libname sos "E:\scifi\data\1.0.original\2026-03\sos";
libname scb "E:\scifi\data\1.0.original\2026-03\scb";
libname proj "S:\HI-SPEED projects\HS26_02 GLP1 pancreatitis\data";

proc delete lib=work data=_all_; run;

data pancreatitis_all_wK852;
set proj.pancreatitis_hosp proj.pancreatitis_cod;
panc_date=coalesce(pancreatitis_hosp_date, pancreatitis_cod_date); format panc_date yymmdd10.;
keep lopnr panc_date sv;
proc sort; by lopnr panc_date; run;

data panc_data_wK852;
merge proj.study_data(in=a keep=lopnr indexdate) pancreatitis_all_wK852;
by lopnr; if a;
if indexdate-6*365<=panc_date<indexdate then pre6y_panc=1;
run;

proc sql;
create table pre6y_panc_wK852 as select
lopnr, max(pre6y_panc) as pre6y_panc_wK852
from panc_data_wK852 group by lopnr;
quit;

data proj.pre6y_panc_wk852;
set pre6y_panc_wK852;
run;

*** search for ERCP and cholecystectomy;
data proj.biliary;
set im.im_020_clean_par(keep=lopnr indatum op1-op30);
array aop op1-op30;
do over aop;
if aop in: ("UJK02" "UJK05" "JKA20" "JKA21") then otherbiliary=1;
end;
if otherbiliary;
run;

