*** 06 Survival analysis Per protocol -- Naive analysis.

cd "S:\HI-SPEED projects\HS26_02 GLP1 pancreatitis\data\"
global resultpath "S:\HI-SPEED projects\HS26_02 GLP1 pancreatitis\output\20260825\PPnaive"
cap mkdir "${resultpath}"

/* RUN once
use ana_data_iptw, clear
keep lopnr exposure birthdate panc_all_wk852_date panc_hosp_wk852_date panc_all_wok852_date panc_hosp_wok852_date cout1_all_date cout2_all_date cout3_all_date *_initdate expo emig_date deathdate fuend discon_date datein expo_gvd expo_gvs agesp* gender indexyr bc region_g education3 fam_panc t2donsetyr_g bmi_6y_g bp_sys_6y_g bp_dias_6y_g gfr_6y_g hba1c_6y_g triglyceride_6y_g ldl_6y_g smokinghabit_6y_g pre_a10a pre_a10ba02 pre_a10bb pre_a10bg pre_c02 pre_c03 pre_c07 pre_c08 pre_c09 pre_c10 pre_alco_dis pre_cereb_dis pre_cholelithiasis pre_diab_nephro pre_diab_retino pre_heart_fail pre_isc_heart_dis pre_jxx00 iptw_smr
egen drug_switch_date_g=rowmin(dpp4_initdate sglt2_initdate) if expo=="glp1"
egen drug_switch_date_d=rowmin(glp1_initdate sglt2_initdate) if expo=="dpp4"
egen drug_switch_date_s=rowmin(dpp4_initdate glp1_initdate) if expo=="sglt2"
egen drug_switch_date=rowmin(drug_switch_date_g drug_switch_date_d drug_switch_date_s)
format drug_switch_date %d
drop drug_switch_date_g drug_switch_date_d drug_switch_date_s
save ana_data_iptw_pp, replace
*/

***** Descriptives
foreach o in panc_all_wk852 panc_all_wok852{
	use ana_data_iptw_pp, clear
	
	egen info_censor_date=rowmin(discon_date drug_switch_date)
	replace info_censor_date = info_censor_date-1 if info_censor_date==`o'_date
	gen info_censor_date_psuedo = datein + 30*(ceil((info_censor_date - datein)/30))
	format info_censor_date_psuedo %d
	
	egen dateout_`o'=rowmin(`o'_date emig_date deathdate fuend info_censor_date_psuedo)
	format datein dateout_`o' fuend %d
	gen `o'_out1=(dateout_`o'==`o'_date)

	stset dateout_`o' [pw=iptw_smr], fail(`o'_out1) enter(datein) origin(datein) scale(365.25)
	gen outyear_`o'=year(dateout_`o') 
	gen out_satt_`o'=_d
	replace out_satt_`o'=2 if out_satt_`o'==0 & dateout_`o'==deathdate
	replace out_satt_`o'=3 if out_satt_`o'==0 & dateout_`o'==emig_date
	replace out_satt_`o'=4 if out_satt_`o'==0 & dateout_`o'==fuend
	replace out_satt_`o'=5 if out_satt_`o'==0 & dateout_`o'==info_censor_date_psuedo
	label define out_satt_`o' 1 "Event" 2 "Death" 3 "Emigration" 4 "EndFollowUp" 5 "Discontinue/Switch"
	label values out_satt_`o' out_satt_`o'
	table (exposure out_satt_`o') (outyear_`o') [pw=iptw_smr]
	collect export "${resultpath}\GLP1 outcome freq IPTW_PPanalys.xlsx", sheet(`o') modify
	
	putexcel set "${resultpath}\GLP1 outcome freq IPTW_PPanalys.xlsx", sheet(`o') modify
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
				cond("`o'"=="panc_all_wok852", "a pancreatitis (excluding K85.2) hospitalisation or death", ""))

		quietly total iptw_smr if _d==1 & exposure==2
		matrix t=e(b)
		local event_glp1=el(t,1,1)
		quietly total iptw_smr if _d==1 & exposure==1
		matrix t=e(b)
		local event_dpp4=el(t,1,1)
		quietly total iptw_smr if _d==1 & exposure==3
		matrix t=e(b)
		local event_sglt2=el(t,1,1)

		gen age_event=round((dateout_`o'-birthdate)/365.25, .1) if _d==1
		quietly su age_event [weight=iptw_smr] if _d==1 & exposure==2, d
		local age_glp1="`=round(r(p50),.1)' (`=round(r(p25),.1)' to `=round(r(p75),.1)')"
		quietly su age_event [weight=iptw_smr] if _d==1 & exposure==1, d
		local age_dpp4="`=round(r(p50),.1)' (`=round(r(p25),.1)' to `=round(r(p75),.1)')"
		quietly su age_event [weight=iptw_smr] if _d==1 & exposure==3, d
		local age_sglt2="`=substr("`=round(r(p50),.1)'", 1,4)' (`=round(r(p25),.1)' to `=round(r(p75),.1)')"

		putexcel I2= "During a median (interquartile range) follow-up of `fu_50' (`fu_25' to `fu_75') years, (maximum `fu_100' years) a total of `event' participants (`event_prop'%) experienced `event_name'. The events were distributed as: `event_glp1' among GLP-1 initiators, `event_dpp4' among DPP-4 initiators, and `event_sglt2' among SGLT-2 initiators. The median (interquartile range) age of pancreatitis event among GLP-1, DPP-4, and SGLT-2 initiators was `age_glp1', `age_dpp4', and `age_sglt2' years, respectively."
	putexcel close

	sts graph, by(exposure) yscale(r(0.99 1)) ylabel(#6, format(%4.3f)) ymtick(##5) xsize(8) xtitle("     ") legend(position(7) ring(0))
	graph export "${resultpath}\km_crude_IPTW_PPnaive`o'.png", replace

}


	* Final analysis
		************* Non-proportional (tvc)
use ana_data_iptw_pp, clear
range time 0 5 101 
foreach o in panc_all_wk852 panc_all_wok852{	
		egen dateout_`o'=rowmin(`o'_date emig_date deathdate fuend drug_switch_date)
		format datein dateout_`o' fuend %d
		gen `o'_out1=(dateout_`o'==`o'_date)

		stset dateout_`o' [pweight=iptw_smr], fail(`o'_out1) enter(datein) origin(datein) scale(365.25)
		
	forvalues i=1/2{		
		local v=cond(`i'==1, "gvd", "gvs")
		local c=cond(`i'==1, "DPP4", "SGLT2")
		
		collect clear
		cap collect:stpm3 i.expo_`v' if expo_`v'!=., knots(5 35 65 95, percentile) tvc(expo_`v') knotstvc(10 50 90, percentile) scale(lncumhazard) eform 	
		if _rc!=0 {
			di as error "Model failed. Skipping the analysis for outcome `o' and comparison between GLP-1 vs `c'"
			putexcel set "${resultpath}\TableS3_PPnaive_tvc.xlsx", sheet(`o'_`v') modify
			putexcel A1="Model failed. Skipping the analysis for outcome `o' and comparison between GLP-1 vs `c'"
			putexcel close
			continue
		}
		collect layout (colname) (result[_r_b _r_p _r_lb _r_ub])
		collect export "${resultpath}\TableS3_PPnaive_tvc.xlsx", sheet(`o'_`v', replace) cell(A15) modify
		
		***Standardised incidence
		standsurv if expo_`v'!=., failure at1(expo_`v' 0) at2(expo_`v' 1) timevar(time) ///
		atvar(`c'_std glp_`c'_std) contrast(difference) contrastvars(riskdiff_glp_`c') ci frame(msurv, replace)

		***Graph across time 
		quietly: mylabels 0(0.1)1, local(labels) myscale(@/100) suffix("%") format(%9.1f)
		frame msurv: tw (line glp_`c'_std time, lc(blue%80)) ///
		(line `c'_std time, lc(red%80)) ///
		(rarea glp_`c'_std_lci glp_`c'_std_uci time, color(blue%30)) ///
		(rarea `c'_std_lci `c'_std_uci time, color(red%30)), ///
		scheme(tab2) ///
		ylabel(`labels', labsize(small)) ///
		xtitle("{bf:Years since treatment initiation}", size(*1) margin(0 0 0 3)) ///
		xlabel(0(0.5)3, nogrid format(%9.1f) labsize(small)) ///
		xscale(range(0 3)) ///
		ytitle("{bf:Standardised cumulative incidence (95% CI)}", size(*1) margin(0 0 0 3)) ///
		name("`v'", replace) ///
		legend(order(3 "GLP-1 initiation" 4 "`c' initiation") col(1) size(*1.5) title("") ring(0) pos(11) yoffset(-11) xoffset(2) region(fcolor(none)))

		graph save "${resultpath}\Fig `o'_`v'_PPnaive_tvc.gph", replace
		graph export "${resultpath}\Fig `o'_`v'_PPnaive_tvc.png", replace
				
		***Derive incidences and differences Table
		foreach var in `c'_std `c'_std_lci `c'_std_uci glp_`c'_std glp_`c'_std_lci glp_`c'_std_uci riskdiff_glp_`c' riskdiff_glp_`c'_lci riskdiff_glp_`c'_uci {
			frame msurv{ 
				cap drop `var'_100
				gen `var'_100 = `var'*100
			}
		}

		frame msurv: export excel time *_100 using "${resultpath}\TableS3_PPnaive_tvc.xlsx" if inlist(time, 0.5, 1, 1.5, 2, 2.5, 3, 3.5, 4, 4.5, 5), sheet(`o'_`v', modify)  firstrow(variables)
		
		***Standardised hazard
		standsurv if expo_`v'!=., hazard at1(expo_`v' 0) at2(expo_`v' 1) timevar(time) ///
		atvar(`v'_std glp_`v'_std) contrast(ratio) contrastvars(hr_glp_`v') ci frame(mhaza, replace)  

		frame mhaza: export excel time * using "${resultpath}\TableS3_PPnaive_tvc_hazaHR.xlsx" if inlist(time, 0.5, 1, 1.5, 2, 2.5, 3, 3.5, 4, 4.5, 5), sheet(`o'_`v', modify)  firstrow(variables)
		
		***Standardised incidence and RR
		standsurv if expo_`v'!=., failure at1(expo_`v' 0) at2(expo_`v' 1) timevar(time) ///
		atvar(`v'_std glp_`v'_std) contrast(ratio) contrastvars(riskratio_glp_`v') ci frame(msurvr, replace) 

		foreach var in `v'_std `v'_std_lci `v'_std_uci glp_`v'_std glp_`v'_std_lci glp_`v'_std_uci {
			frame msurvr{ 
				cap drop `var'_100
				gen `var'_100 = `var'*100
			}
		}
		
		frame msurvr: export excel using "${resultpath}\TableS3_PPnaive_tvc_failratio.xlsx" if inlist(time, 0.5, 1, 1.5, 2, 2.5, 3, 3.5, 4, 4.5, 5), sheet(`o'_`v', modify)  firstrow(variables)
	}
}	
