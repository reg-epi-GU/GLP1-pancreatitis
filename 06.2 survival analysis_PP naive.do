*** 06 Survival analysis Per protocol -- Naive analysis.

cd "S:\HI-SPEED projects\HS26_02 GLP1 pancreatitis\data\"
global resultpath "S:\HI-SPEED projects\HS26_02 GLP1 pancreatitis\output\20260422\PP_naive"
cap mkdir "${resultpath}"


***** Descriptives
use ana_data_perprotocol, clear
foreach o in pancreatitis_all pancreatitis_hosp cout1_all cout2_all{
	egen dateout_`o'=rowmin(`o'_date emig_date deathdate fuend drug_switch_date)
	format datein dateout_`o' fuend %d
	gen `o'_out1=(dateout_`o'==`o'_date)

	stset dateout_`o', fail(`o'_out1) enter(datein) origin(datein) scale(365.25)
	gen outyear_`o'=year(dateout_`o') 
	gen out_satt_`o'=_d
	replace out_satt_`o'=2 if out_satt_`o'==0 & dateout_`o'==deathdate
	replace out_satt_`o'=3 if out_satt_`o'==0 & dateout_`o'==emig_date
	replace out_satt_`o'=4 if out_satt_`o'==0 & dateout_`o'==fuend
	replace out_satt_`o'=5 if out_satt_`o'==0 & dateout_`o'==drug_switch_date
	label define out_satt_`o' 1 "Event" 2 "Death" 3 "Emigration" 4 "EndFollowUp" 5 "DrugSwitch"
	label values out_satt_`o' out_satt_`o'
	table (exposure out_satt_`o') (outyear_`o') 
	collect export "${resultpath}\GLP1 outcome freq PPanalys.xlsx", sheet(`o') modify
	
	putexcel set "${resultpath}\GLP1 outcome freq PPanalys.xlsx", sheet(`o') modify
		quietly su _t, d 
		local fu_50=round(r(p50),.1)
		local fu_25=round(r(p25),.1)
		local fu_75=round(r(p75),.1)
		local fu_100=round(r(max),.1)

		quietly count if _d==1
		local event=r(N)
		quietly count
		local totalpop=r(N)
		local event_prop=round(100*`event'/`totalpop', .01)

		local event_name=cond("`o'"=="pancreatitis_all", "an acute pancreatitis event, of which 570 were hospitalisations and 10 deaths", ///
				cond("`o'"=="pancreatitis_hosp", "a pancreatitis hospitalisation", ///
				cond("`o'"=="gisymp_hosp", "a GI symptom hospitalisation", "a severe COVID-19")))

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

	sts graph, by(exposure) yscale(r(0.99 1)) ylabel(#6, format(%4.3f)) ymtick(##5) xsize(8) xtitle("     ") legend(position(7) ring(0))
	graph export "${resultpath}\km_crude_PPnaive`o'.png", replace

}


	* Final analysis
		************* Proportional Assumption (no tvc)
use ana_data_perprotocol, clear
range time 0 3 49 
foreach o in pancreatitis_all pancreatitis_hosp cout1_all cout2_all{	
		egen dateout_`o'=rowmin(`o'_date emig_date deathdate fuend drug_switch_date)
		format datein dateout_`o' fuend %d
		gen `o'_out1=(dateout_`o'==`o'_date)

		stset dateout_`o', fail(`o'_out1) enter(datein) origin(datein) scale(365.25)
		
	forvalues i=1/2{		
		local v=cond(`i'==1, "gvd", "gvs")
		local c=cond(`i'==1, "DPP4", "SGLT2")
		
		collect clear
		cap collect:stpm3 i.expo_`v' agesp* i.gender i.indexyr ib2.region_g i.bc i.education3 i.fam_panc i.t2donsetyr_g ib2.bmi_g i.smokinghabit_g i.hba1c_g i.preindex_lipidlowerdrug i.preindex_alcoholdisorder i.preindex_gallstones i.preindex_abdomsurgery if expo_`v'!=., knots(5 35 65 95, percentile) scale(lncumhazard) eform 	
		if _rc==0 | ("`o'"=="cout2_all")  {
			di as error "Full PH model for `o' `v' failed, running model without abdomsurgery"
			collect clear
			collect:stpm3 i.expo_`v' agesp* i.gender i.indexyr ib2.region_g i.bc i.education3 i.fam_panc i.t2donsetyr_g ib2.bmi_g i.smokinghabit_g i.hba1c_g i.preindex_lipidlowerdrug i.preindex_alcoholdisorder i.preindex_gallstones if expo_`v'!=., knots(5 35 65 95, percentile) scale(lncumhazard) eform 	
		}
		collect layout (colname) (result[_r_b _r_p _r_lb _r_ub])
		collect export "${resultpath}\TableS3_PP_naive_ph.xlsx", sheet(`o'_`v', replace) cell(A9) modify
		
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

		graph save "${resultpath}\Fig `o'_`v'_PP_naive_ph.gph", replace
				
		***Derive incidences and differences Table
		foreach var in `c'_std `c'_std_lci `c'_std_uci glp_`c'_std glp_`c'_std_lci glp_`c'_std_uci riskdiff_glp_`c' riskdiff_glp_`c'_lci riskdiff_glp_`c'_uci {
			frame msurv{ 
				cap drop `var'_100
				gen `var'_100 = `var'*100
			}
		}

		frame msurv: export excel time *_100 using "${resultpath}\TableS3_PP_naive_ph.xlsx" if inlist(time, 0.5, 1, 1.5, 2, 2.5, 3), sheet(`o'_`v', modify)  firstrow(variables)		
	}		
}	

	* Final analysis
		************* Non-proportional (tvc)
use ana_data_perprotocol, clear
range time 0 3 49 
foreach o in pancreatitis_all pancreatitis_hosp cout1_all /*cout2_all*/{	
		egen dateout_`o'=rowmin(`o'_date emig_date deathdate fuend drug_switch_date)
		format datein dateout_`o' fuend %d
		gen `o'_out1=(dateout_`o'==`o'_date)

		stset dateout_`o', fail(`o'_out1) enter(datein) origin(datein) scale(365.25)
		
	forvalues i=1/2{		
		local v=cond(`i'==1, "gvd", "gvs")
		local c=cond(`i'==1, "DPP4", "SGLT2")
		
		collect clear
		cap collect:stpm3 i.expo_`v' agesp* i.gender i.indexyr ib2.region_g i.bc i.education3 i.fam_panc i.t2donsetyr_g ib2.bmi_g i.smokinghabit_g i.hba1c_g i.preindex_lipidlowerdrug i.preindex_alcoholdisorder i.preindex_gallstones i.preindex_abdomsurgery if expo_`v'!=., knots(5 35 65 95, percentile) tvc(expo_`v') knotstvc(10 50 90, percentile) scale(lncumhazard) eform 	
		if _rc==0 | ("`o'"=="cout2_all")  {
			di as error "Full PH model for `o' `v' failed, running model without abdomsurgery"
			collect clear
			collect:stpm3 i.expo_`v' agesp* i.gender i.indexyr ib2.region_g i.bc i.education3 i.fam_panc i.t2donsetyr_g ib2.bmi_g i.smokinghabit_g i.hba1c_g i.preindex_lipidlowerdrug i.preindex_alcoholdisorder i.preindex_gallstones if expo_`v'!=., knots(5 35 65 95, percentile) scale(lncumhazard) eform 	
		}
		collect layout (colname) (result[_r_b _r_p _r_lb _r_ub])
		collect export "${resultpath}\TableS3_PP_naive_tvc.xlsx", sheet(`o'_`v', replace) cell(A9) modify
		
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

		graph save "${resultpath}\Fig `o'_`v'_PP_naive_tvc.gph", replace
				
		***Derive incidences and differences Table
		foreach var in `c'_std `c'_std_lci `c'_std_uci glp_`c'_std glp_`c'_std_lci glp_`c'_std_uci riskdiff_glp_`c' riskdiff_glp_`c'_lci riskdiff_glp_`c'_uci {
			frame msurv{ 
				cap drop `var'_100
				gen `var'_100 = `var'*100
			}
		}

		frame msurv: export excel time *_100 using "${resultpath}\TableS3_PP_naive_tvc.xlsx" if inlist(time, 0.5, 1, 1.5, 2, 2.5, 3), sheet(`o'_`v', modify)  firstrow(variables)
		
	}
}	
