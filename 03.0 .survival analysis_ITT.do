*** 03 Survival analysis

cd "S:\HI-SPEED projects\HS26_02 GLP1 pancreatitis\data\"
global resultpath "S:\HI-SPEED projects\HS26_02 GLP1 pancreatitis\output\20260825"

***************************** We decided to use t2d_indication_ndr==1 as the inclusion...
/* RUN once
use ana_data, clear

replace emig_date=. if emig_date<=indexdate
gen fuend=mdy(12,31,2025)
gen datein=indexdate-1

recode exposure (1=0) (2=1) (3=.), gen(expo_gvd)
label variable expo_gvd "GLP1 (vs DPP-4)"
recode exposure (3=0) (2=1) (1=.), gen(expo_gvs)
label variable expo_gvs "GLP1 (vs SGLT-2)" 

cap rename (panc_all_date_wk852 panc_all_date_wok852 panc_cod_date_wk852 panc_cod_date_wok852 panc_hosp_date_wk852 panc_hosp_date_wok852) (panc_all_wk852_date panc_all_wok852_date panc_cod_wk852_date panc_cod_wok852_date panc_hosp_wk852_date panc_hosp_wok852_date)

mkspline agesp = age_index, cubic
*mkspline bmisp = bmi_num, cubic //not possible to do due to missings. Use categories instead
*mkspline hba1csp = hba1c_num, cubic //not possible to do due to missings. Use categories instead

save, replace
*/


/* RUN ONCE
use ana_data, clear
***** IPTW, SMR weighting. 
gen iptw_smr=1 if exposure==2
forvalues i=1/2{		
	local v=cond(`i'==1, "gvd", "gvs")
	local c=cond(`i'==1, "DPP-4", "SGLT-2")
	quietly logit expo_`v' agesp* i.gender i.indexyr i.bc i.region_g i.education3 i.fam_panc i.t2donsetyr_g ///
	i.bmi_6y_g i.bp_sys_6y_g i.bp_dias_6y_g i.gfr_6y_g i.hba1c_6y_g i.triglyceride_6y_g i.ldl_6y_g ///
	i.smokinghabit_6y_g i.pre_a10a i.pre_a10ba02 i.pre_a10bb i.pre_a10bg i.pre_c02 i.pre_c03 i.pre_c07 i.pre_c08 i.pre_c09 i.pre_c10 i.pre_alco_dis i.pre_cereb_dis i.pre_cholelithiasis i.pre_diab_nephro i.pre_diab_retino i.pre_heart_fail i.pre_isc_heart_dis i.pre_jxx00 if expo_`v'!=.
	predict p_expo_`v'
	replace iptw_smr=p_expo_`v'/(1-p_expo_`v') if expo_`v'==0
	twoway (kdensity p_expo_`v' if expo_`v' == 1, color(navy%40)) ///
		   (kdensity p_expo_`v' if expo_`v' == 0, color(maroon%40)), ///
           legend(order(1 "GLP1" 2 "`c'")) xtitle("") ytitle("Density")
	graph export "${resultpath}\probability of treatment GLP1 vs `c'.svg", replace
}
replace iptw_smr=10 if iptw_smr>10
save ana_data_iptw, replace
*/

**** check balance afterwards.
use ana_data_iptw, clear
dtable age_index i.age_g i.gender i.indexyr i.bc i.region_g i.education3 i.fam_panc i.t2donsetyr_g ///
	bmi_1y_num i.bmi_1y_g bmi_2y_num i.bmi_2y_g bmi_3y_num i.bmi_3y_g bmi_6y_num i.bmi_6y_g ///
	bp_sys_1y_num i.bp_sys_1y_g bp_sys_2y_num i.bp_sys_2y_g bp_sys_3y_num i.bp_sys_3y_g bp_sys_6y_num i.bp_sys_6y_g ///
	bp_dias_1y_num i.bp_dias_1y_g bp_dias_2y_num i.bp_dias_2y_g bp_dias_3y_num i.bp_dias_3y_g bp_dias_6y_num i.bp_dias_6y_g ///
	gfr_1y_num i.gfr_1y_g gfr_2y_num i.gfr_2y_g gfr_3y_num i.gfr_3y_g gfr_6y_num i.gfr_6y_g ///
	hba1c_1y_num i.hba1c_1y_g hba1c_2y_num i.hba1c_2y_g hba1c_3y_num i.hba1c_3y_g hba1c_6y_num i.hba1c_6y_g ///
	triglyceride_1y_num i.triglyceride_1y_g triglyceride_2y_num i.triglyceride_2y_g triglyceride_3y_num i.triglyceride_3y_g triglyceride_6y_num i.triglyceride_6y_g ///
	ldl_1y_num i.ldl_1y_g ldl_2y_num i.ldl_2y_g ldl_3y_num i.ldl_3y_g ldl_6y_num i.ldl_6y_g ///
	i.complete_1y i.complete_2y i.complete_3y i.complete_6y ///
	i.smokinghabit_6y_g i.pre_a10a i.pre_a10ba02 i.pre_a10bb i.pre_a10bg i.pre_c02 i.pre_c03 i.pre_c07 i.pre_c08 i.pre_c09 i.pre_c10 i.pre_alco_dis i.pre_cereb_dis i.pre_cholelithiasis i.pre_diab_nephro i.pre_diab_retino i.pre_heart_fail i.pre_isc_heart_dis i.pre_jxx00 i.panc_all_wk852 i.panc_hosp_wk852 i.panc_cod_wk852 i.panc_all_wok852 i.panc_hosp_wok852 i.panc_cod_wok852 [pweight=iptw_smr], by(exposure) export("${resultpath}\GLP1 pancreatitis T1_iptw.xlsx", modify)


foreach o in panc_all_wk852 panc_hosp_wk852 panc_all_wok852 panc_hosp_wok852 cout1_all cout2_all cout3_all{
	egen dateout_`o'=rowmin(`o'_date emig_date deathdate fuend)
	format datein dateout_`o' fuend %d
	gen `o'_out1=(dateout_`o'==`o'_date)

	stset dateout_`o' [pweight=iptw_smr], fail(`o'_out1) enter(datein) origin(datein) scale(365.25)
	gen outyear_`o'=year(dateout_`o') 
	gen out_satt_`o'=_d
	replace out_satt_`o'=2 if out_satt_`o'==0 & dateout_`o'==deathdate
	replace out_satt_`o'=3 if out_satt_`o'==0 & dateout_`o'==emig_date
	replace out_satt_`o'=4 if out_satt_`o'==0 & dateout_`o'==fuend
	label define out_satt_`o' 1 "Event" 2 "Death" 3 "Emigration" 4 "EndFollowUp"
	label values out_satt_`o' out_satt_`o'
	table (exposure out_satt_`o') (outyear_`o') [pw=iptw_smr]
	collect export "${resultpath}\GLP1 outcome freq_iptw.xlsx", sheet(`o') modify
	
	if "`o'"=="panc_all_wk852" {
		table pantype [pw=iptw_smr]
		collect export "${resultpath}\GLP1 outcome freq_iptw.xlsx", sheet(`o') cell(J5) modify
	}
	
	putexcel set "${resultpath}\GLP1 outcome freq_iptw.xlsx", sheet(`o') modify
		quietly su _t [weight=iptw_smr], d 
		local fu_50=round(r(p50),.1)
		local fu_25=round(r(p25),.1)
		local fu_75=round(r(p75),.1)
		local fu_100=round(r(max),.1)

		quietly total iptw_smr if _d==1
		matrix t=e(b)
		local event=el(t,1,1)
		quietly total iptw_smr
		matrix t=e(b)
		local totalpop=el(t,1,1)
		local event_prop=round(`=100*`event'/`totalpop'', .01)

		local event_name=cond("`o'"=="panc_all_wk852", "an acute pancreatitis event (including K85.2) hospitalisation or death", ///
				cond("`o'"=="panc_hosp_wk852", "a pancreatitis (including K85.2) hospitalisation", ///
				cond("`o'"=="panc_all_wok852", "a pancreatitis (excluding K85.2) hospitalisation or death", ///
				cond("`o'"=="panc_hosp_wok852", "a pancreatitis (excluding K85.2) hospitalisation", ///
				cond("`o'"=="cout1_all", "Appendicitis hospitalisation or death", ///
				cond("`o'"=="cout2_all", "Osteoarthritis hospitalisation or death", ///				
				cond("`o'"=="cout3_all", "Herpes zoster hospitalisation or death", "")))))))

		quietly total iptw_smr if _d==1 & exposure==2
		matrix t=e(b)
		local event_glp1=el(t,1,1)
		quietly total iptw_smr if _d==1 & exposure==1
		matrix t=e(b)
		local event_dpp4=el(t,1,1)
		quietly total iptw_smr if _d==1 & exposure==3
		matrix t=e(b)
		local event_sglt2=el(t,1,1)

		cap gen age_event=round((dateout_`o'-birthdate)/365.25, .1) if _d==1
		quietly su age_event [weight=iptw_smr] if _d==1 & exposure==2, d
		local age_glp1="`=round(r(p50),.1)' (`=round(r(p25),.1)' to `=round(r(p75),.1)')"
		quietly su age_event [weight=iptw_smr] if _d==1 & exposure==1, d
		local age_dpp4="`=round(r(p50),.1)' (`=round(r(p25),.1)' to `=round(r(p75),.1)')"
		quietly su age_event [weight=iptw_smr] if _d==1 & exposure==3, d
		local age_sglt2="`=substr("`=round(r(p50),.1)'", 1,4)' (`=round(r(p25),.1)' to `=round(r(p75),.1)')"
		cap drop age_event

		putexcel I2= "During a median (interquartile range) follow-up of `fu_50' (`fu_25' to `fu_75') years, (maximum `fu_100' years) a total of `event' participants (`event_prop'%) experienced `event_name'. The events were distributed as: `event_glp1' among GLP-1 initiators, `event_dpp4' among DPP-4 initiators, and `event_sglt2' among SGLT-2 initiators. The median (interquartile range) age of pancreatitis event among GLP-1, DPP-4, and SGLT-2 initiators was `age_glp1', `age_dpp4', and `age_sglt2' years, respectively."
		putexcel close

		sts graph, by(exposure) yscale(r(0.99 1)) ylabel(#6, format(%4.3f)) ymtick(##5) xsize(8) xtitle("     ") legend(position(7) ring(0))
		graph export "${resultpath}\km_`o'_iptw.png", replace
}


	* Final analysis
		************* TVC
use ana_data_iptw, clear
range time 0 5 101 
foreach o in panc_all_wk852 panc_hosp_wk852 panc_all_wok852 panc_hosp_wok852 cout1_all cout2_all cout3_all{
		egen dateout_`o'=rowmin(`o'_date emig_date deathdate fuend)
		format datein dateout_`o' fuend %d
		gen `o'_out1=(dateout_`o'==`o'_date)

		stset dateout_`o' [pweight=iptw_smr], fail(`o'_out1) enter(datein) origin(datein) scale(365.25)
		
	forvalues i=1/2{		
		local v=cond(`i'==1, "gvd", "gvs")
		local c=cond(`i'==1, "DPP-4", "SGLT-2")
		
		collect clear
		cap collect:stpm3 i.expo_`v' if expo_`v'!=.,  knots(5 35 65 95, percentile) tvc(expo_`v') knotstvc(10 50 90, percentile) scale(lncumhazard) eform 
		if _rc != 0 {
				di as error "Model failed. Skipping the analysis for outcome `o' and comparison between GLP-1 vs `c'"
				putexcel set "${resultpath}\TableS3_ITT_tvc.xlsx", sheet(`o'_`v') modify
				putexcel A1="Model failed. Skipping the analysis for outcome `o' and comparison between GLP-1 vs `c'"
				putexcel close
				continue
		} 
		collect layout (colname) (result[_r_b _r_p _r_lb _r_ub])
		collect export "${resultpath}\TableS3_ITT_tvc.xlsx", sheet(`o'_`v', replace) cell(A15) modify
		
		***Standardised incidence
		standsurv if expo_`v'!=., failure at1(expo_`v' 0) at2(expo_`v' 1) timevar(time) ///
		atvar(`v'_std glp_`v'_std) contrast(difference) contrastvars(riskdiff_glp_`v') ci frame(msurv, replace)  //indweights(iptw_smr)    weights should not be used again in prediction.

		***Graph across time 
		quietly: mylabels 0(0.1)1, local(labels) myscale(@/100) suffix("%") format(%9.1f)
		frame msurv: tw (line glp_`v'_std time, lc(blue%80)) ///
		(line `v'_std time, lc(red%80)) ///
		(rarea glp_`v'_std_lci glp_`v'_std_uci time, color(blue%30)) ///
		(rarea `v'_std_lci `v'_std_uci time, color(red%30)), ///
		scheme(tab2) ///
		ylabel(`labels', labsize(small)) ///
		xtitle("{bf:Years since treatment initiation}", size(*1) margin(0 0 0 3)) ///
		xlabel(0(0.5)3, nogrid format(%9.1f) labsize(small)) ///
		xscale(range(0 3)) ///
		ytitle("{bf:Standardised cumulative incidence (95% CI)}", size(*1) margin(0 0 0 3)) ///
		name("`v'", replace) ///
		legend(order(3 "GLP-1 initiation" 4 "`c' initiation") col(1) size(*1.5) title("`o'") ring(0) pos(11) yoffset(-11) xoffset(2) region(fcolor(none)))

		graph save "${resultpath}\Fig `o'_`v'_ITT_tvc.gph", replace
		graph export "${resultpath}\Fig `o'_`v'_ITT_tvc.png", replace
				
		***Derive incidences and differences Table
		foreach var in `v'_std `v'_std_lci `v'_std_uci glp_`v'_std glp_`v'_std_lci glp_`v'_std_uci riskdiff_glp_`v' riskdiff_glp_`v'_lci riskdiff_glp_`v'_uci {
			frame msurv{ 
				cap drop `var'_100
				gen `var'_100 = `var'*100
			}
		}

		frame msurv: export excel time *_100 using "${resultpath}\TableS3_ITT_tvc.xlsx" if inlist(time, 0.5, 1, 1.5, 2, 2.5, 3, 3.5, 4, 4.5, 5), sheet(`o'_`v', modify)  firstrow(variables)
		
		***Standardised hazard
		standsurv if expo_`v'!=., hazard at1(expo_`v' 0) at2(expo_`v' 1) timevar(time) ///
		atvar(`v'_std glp_`v'_std) contrast(ratio) contrastvars(hr_glp_`v') ci frame(mhaza, replace)  //indweights(iptw_smr)    weights should not be used again in prediction.

		frame mhaza: export excel time * using "${resultpath}\TableS3_ITT_tvc_hazaHR.xlsx" if inlist(time, 0.5, 1, 1.5, 2, 2.5, 3, 3.5, 4, 4.5, 5), sheet(`o'_`v', modify)  firstrow(variables)
		
		***Standardised incidence and RR
		standsurv if expo_`v'!=., failure at1(expo_`v' 0) at2(expo_`v' 1) timevar(time) ///
		atvar(`v'_std glp_`v'_std) contrast(ratio) contrastvars(riskratio_glp_`v') ci frame(msurvr, replace)  //indweights(iptw_smr)    weights should not be used again in prediction.

		foreach var in `v'_std `v'_std_lci `v'_std_uci glp_`v'_std glp_`v'_std_lci glp_`v'_std_uci {
			frame msurvr{ 
				cap drop `var'_100
				gen `var'_100 = `var'*100
			}
		}
		
		frame msurvr: export excel using "${resultpath}\TableS3_ITT_tvc_failratio.xlsx" if inlist(time, 0.5, 1, 1.5, 2, 2.5, 3, 3.5, 4, 4.5, 5), sheet(`o'_`v', modify)  firstrow(variables)
	}
}	
