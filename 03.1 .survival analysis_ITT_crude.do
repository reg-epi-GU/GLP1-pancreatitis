*** 03 Survival analysis

cd "S:\HI-SPEED projects\HS26_02 GLP1 pancreatitis\data\"
global resultpath "S:\HI-SPEED projects\HS26_02 GLP1 pancreatitis\output\20260825"

use ana_data, clear
foreach o in panc_all_wk852 panc_hosp_wk852 panc_all_wok852 panc_hosp_wok852 cout1_all cout2_all cout3_all{
	egen dateout_`o'=rowmin(`o'_date emig_date deathdate fuend)
	format datein dateout_`o' fuend %d
	gen `o'_out1=(dateout_`o'==`o'_date)

	stset dateout_`o', fail(`o'_out1) enter(datein) origin(datein) scale(365.25)
	gen outyear_`o'=year(dateout_`o') 
	gen out_satt_`o'=_d
	replace out_satt_`o'=2 if out_satt_`o'==0 & dateout_`o'==deathdate
	replace out_satt_`o'=3 if out_satt_`o'==0 & dateout_`o'==emig_date
	replace out_satt_`o'=4 if out_satt_`o'==0 & dateout_`o'==fuend
	label define out_satt_`o' 1 "Event" 2 "Death" 3 "Emigration" 4 "EndFollowUp"
	label values out_satt_`o' out_satt_`o'
	table (exposure out_satt_`o') (outyear_`o')
	collect export "${resultpath}\GLP1 outcome freq_crude.xlsx", sheet(`o') modify
	
	if "`o'"=="panc_all_wk852" {
		table pantype
		collect export "${resultpath}\GLP1 outcome freq_crude.xlsx", sheet(`o') cell(J5) modify
	}
	
	putexcel set "${resultpath}\GLP1 outcome freq_crude.xlsx", sheet(`o') modify
		quietly su _t, d 
		local fu_50=round(r(p50),.1)
		local fu_25=round(r(p25),.1)
		local fu_75=round(r(p75),.1)
		local fu_100=round(r(max),.1)

		quietly count if _d==1
		local event=r(N)
		quietly count 
		local totalpop=r(N)
		local event_prop=round(`=100*`event'/`totalpop'', .01)

		local event_name=cond("`o'"=="panc_all_wk852", "an acute pancreatitis event (including K85.2) hospitalisation or death", ///
				cond("`o'"=="panc_hosp_wk852", "a pancreatitis (including K85.2) hospitalisation", ///
				cond("`o'"=="panc_all_wok852", "a pancreatitis (excluding K85.2) hospitalisation or death", ///
				cond("`o'"=="panc_hosp_wok852", "a pancreatitis (excluding K85.2) hospitalisation", ///
				cond("`o'"=="cout1_all", "Appendicitis hospitalisation or death", ///
				cond("`o'"=="cout2_all", "Osteoarthritis hospitalisation or death", ///				
				cond("`o'"=="cout3_all", "Herpes zoster hospitalisation or death", "")))))))

		quietly count if _d==1 & exposure==2
		local event_glp1=r(N)
		quietly count if _d==1 & exposure==1
		local event_dpp4=r(N)
		quietly count if _d==1 & exposure==3
		local event_sglt2=r(N)

		cap gen age_event=round((dateout_`o'-birthdate)/365.25, .1) if _d==1
		quietly su age_event if _d==1 & exposure==2, d
		local age_glp1="`=round(r(p50),.1)' (`=round(r(p25),.1)' to `=round(r(p75),.1)')"
		quietly su age_event if _d==1 & exposure==1, d
		local age_dpp4="`=round(r(p50),.1)' (`=round(r(p25),.1)' to `=round(r(p75),.1)')"
		quietly su age_event if _d==1 & exposure==3, d
		local age_sglt2="`=substr("`=round(r(p50),.1)'", 1,4)' (`=round(r(p25),.1)' to `=round(r(p75),.1)')"
		cap drop age_event

		putexcel I2= "During a median (interquartile range) follow-up of `fu_50' (`fu_25' to `fu_75') years, (maximum `fu_100' years) a total of `event' participants (`event_prop'%) experienced `event_name'. The events were distributed as: `event_glp1' among GLP-1 initiators, `event_dpp4' among DPP-4 initiators, and `event_sglt2' among SGLT-2 initiators. The median (interquartile range) age of pancreatitis event among GLP-1, DPP-4, and SGLT-2 initiators was `age_glp1', `age_dpp4', and `age_sglt2' years, respectively."
		putexcel close
}


	* Final analysis
		************* TVC
use ana_data, clear
range time 0 5 101 
foreach o in panc_all_wk852 panc_hosp_wk852 panc_all_wok852 panc_hosp_wok852 cout1_all cout2_all cout3_all{
		egen dateout_`o'=rowmin(`o'_date emig_date deathdate fuend)
		format datein dateout_`o' fuend %d
		gen `o'_out1=(dateout_`o'==`o'_date)

		stset dateout_`o', fail(`o'_out1) enter(datein) origin(datein) scale(365.25)
		
	forvalues i=1/2{		
		local v=cond(`i'==1, "gvd", "gvs")
		local c=cond(`i'==1, "DPP-4", "SGLT-2")
		
		collect clear
		cap collect:stpm3 i.expo_`v' if expo_`v'!=.,  knots(5 35 65 95, percentile) tvc(expo_`v') knotstvc(10 50 90, percentile) scale(lncumhazard) eform 
		if _rc != 0 {
				di as error "Model failed. Skipping the analysis for outcome `o' and comparison between GLP-1 vs `c'"
				putexcel set "${resultpath}\TableS3_ITT_tvc_crude.xlsx", sheet(`o'_`v') modify
				putexcel A1="Model failed. Skipping the analysis for outcome `o' and comparison between GLP-1 vs `c'"
				putexcel close
				continue
		} 
		collect layout (colname) (result[_r_b _r_p _r_lb _r_ub])
		collect export "${resultpath}\TableS3_ITT_tvc_crude.xlsx", sheet(`o'_`v', replace) cell(A15) modify
		
		***Standardised incidence
		standsurv if expo_`v'!=., failure at1(expo_`v' 0) at2(expo_`v' 1) timevar(time) ///
		atvar(`v'_std glp_`v'_std) contrast(difference) contrastvars(riskdiff_glp_`v') ci frame(msurv, replace)  //indweights(iptw_smr)    weights should not be used again in prediction.

// 		***Graph across time 
// 		quietly: mylabels 0(0.1)1, local(labels) myscale(@/100) suffix("%") format(%9.1f)
// 		frame msurv: tw (line glp_`v'_std time, lc(blue%80)) ///
// 		(line `v'_std time, lc(red%80)) ///
// 		(rarea glp_`v'_std_lci glp_`v'_std_uci time, color(blue%30)) ///
// 		(rarea `v'_std_lci `v'_std_uci time, color(red%30)), ///
// 		scheme(tab2) ///
// 		ylabel(`labels', labsize(small)) ///
// 		xtitle("{bf:Years since treatment initiation}", size(*1) margin(0 0 0 3)) ///
// 		xlabel(0(0.5)3, nogrid format(%9.1f) labsize(small)) ///
// 		xscale(range(0 3)) ///
// 		ytitle("{bf:Standardised cumulative incidence (95% CI)}", size(*1) margin(0 0 0 3)) ///
// 		name("`v'", replace) ///
// 		legend(order(3 "GLP-1 initiation" 4 "`c' initiation") col(1) size(*1.5) title("`o'") ring(0) pos(11) yoffset(-11) xoffset(2) region(fcolor(none)))
//
// 		graph save "${resultpath}\Fig `o'_`v'_ITT_tvc.gph", replace
// 		graph export "${resultpath}\Fig `o'_`v'_ITT_tvc.png", replace
				
		***Derive incidences and differences Table
		foreach var in `v'_std `v'_std_lci `v'_std_uci glp_`v'_std glp_`v'_std_lci glp_`v'_std_uci riskdiff_glp_`v' riskdiff_glp_`v'_lci riskdiff_glp_`v'_uci {
			frame msurv{ 
				cap drop `var'_100
				gen `var'_100 = `var'*100
			}
		}

		frame msurv: export excel time *_100 using "${resultpath}\TableS3_ITT_tvc_crude.xlsx" if inlist(time, 0.5, 1, 1.5, 2, 2.5, 3, 3.5, 4, 4.5, 5), sheet(`o'_`v', modify)  firstrow(variables)
		
		***Standardised hazard
		standsurv if expo_`v'!=., hazard at1(expo_`v' 0) at2(expo_`v' 1) timevar(time) ///
		atvar(`v'_std glp_`v'_std) contrast(ratio) contrastvars(hr_glp_`v') ci frame(mhaza, replace)  //indweights(iptw_smr)    weights should not be used again in prediction.

		frame mhaza: export excel time * using "${resultpath}\TableS3_ITT_tvc_crude_hazaHR.xlsx" if inlist(time, 0.5, 1, 1.5, 2, 2.5, 3, 3.5, 4, 4.5, 5), sheet(`o'_`v', modify)  firstrow(variables)
		
		***Standardised incidence and RR
		standsurv if expo_`v'!=., failure at1(expo_`v' 0) at2(expo_`v' 1) timevar(time) ///
		atvar(`v'_std glp_`v'_std) contrast(ratio) contrastvars(riskratio_glp_`v') ci frame(msurvr, replace)  //indweights(iptw_smr)    weights should not be used again in prediction.

		foreach var in `v'_std `v'_std_lci `v'_std_uci glp_`v'_std glp_`v'_std_lci glp_`v'_std_uci {
			frame msurvr{ 
				cap drop `var'_100
				gen `var'_100 = `var'*100
			}
		}
		
		frame msurvr: export excel using "${resultpath}\TableS3_ITT_tvc_crude_failratio.xlsx" if inlist(time, 0.5, 1, 1.5, 2, 2.5, 3, 3.5, 4, 4.5, 5), sheet(`o'_`v', modify)  firstrow(variables)
	}
}	
