cd "S:\HI-SPEED projects\HS26_02 GLP1 pancreatitis\data\"
global resultpath "S:\HI-SPEED projects\HS26_02 GLP1 pancreatitis\output\20261015\"
cap mkdir "${resultpath}"

******** weighted non-parametric curve (KM), overall, and 2021~2022, 2023~2025
** If a better graph is needed, we must use STATA19 with new option "altrisktable" for sts graph.
use ana_data_iptw, clear
gen yr_in=year(datein)
egen yr_in_g=cut(yr_in), at(2021,2023,2026) label
foreach o in panc_all_wk852 {
	egen dateout_`o'=rowmin(`o'_date emig_date deathdate fuend)
	format datein dateout_`o' fuend %d
	gen `o'_out1=(dateout_`o'==`o'_date)

	stset dateout_`o' [pweight=iptw_smr], fail(`o'_out1) enter(datein) origin(datein) scale(365.25)
	
**** Overall figure.
		sts graph, by(exposure) yscale(r(0.99 1)) ylabel(#6, format(%4.3f)) ymtick(##5) xscale(r(0 5)) xsize(12) xtitle("     ") xlabel(0(1)5) legend(position(7) ring(0)) title("KM curve for those who initiated in year 2021-2025")
		graph export "${resultpath}\km_`o'_iptw_allyears.png", replace	
		
		**** calculate risktable
		capture frame drop risktable
		frame create risktable

		tempfile riskdata
		tempname handle

		postfile `handle' double exposure time nrisk nfail using `riskdata'

		forvalues g = 1/3 {
			forvalues t = 0/5 {
				quietly summarize iptw_smr if ///
					_st == 1 & exposure == `g' & ///
					_t0 <= `t' & _t >= `t', meanonly
				local nr = r(sum)

				quietly summarize iptw_smr if ///
					_st == 1 & exposure == `g' & _d == 1 & ///
					_t >= `t' & _t < `t' + 1, meanonly
				local nf = r(sum)

				post `handle' (`g') (`t') (round(`nr')) (round(`nf'))
			}
		}

		postclose `handle'
		frame risktable: use `riskdata', clear
		frame risktable: gen str40 risklabel = strtrim(string(nrisk, "%12.0fc")) + " (" + strtrim(string(nfail, "%12.0fc")) + ")"
		frame risktable: label define exposure 1 "DPP-4" 2 "GLP-1 RA" 3 "SGLT-2"
		frame risktable: label values exposure exposure
		frame risktable: drop nrisk nfail
		frame risktable: reshape wide risklabel,i(exposure) j(time)
		frame risktable: export excel using "${resultpath}\risktable.xlsx", sheet("overall") replace

	

	forvalues yg=0/1{
		preserve 
		keep if yr_in_g==`yg'
		
		local lab=cond(`yg'==0, "2021-2022", "2023-2025")
		sts graph, by(exposure) yscale(r(0.99 1)) ylabel(#6, format(%4.3f)) ymtick(##5) xscale(r(0 5)) xsize(12) xtitle("     ") xlabel(0(1)`=5-2*`yg'') legend(position(7) ring(0)) title("KM curve for those who initiated in year `lab'")
		graph export "${resultpath}\km_`o'_iptw_`lab'.png", replace		
		
		**** calculate risktable
		capture frame drop risktable
		frame create risktable

		tempfile riskdata
		tempname handle

		postfile `handle' double exposure time nrisk nfail using `riskdata'
		
		forvalues g = 1/3 {
			forvalues t = 0/5 {
				quietly summarize iptw_smr if ///
					_st == 1 & exposure == `g' & ///
					_t0 <= `t' & _t >= `t', meanonly
				local nr = r(sum)

				quietly summarize iptw_smr if ///
					_st == 1 & exposure == `g' & _d == 1 & ///
					_t >= `t' & _t < `t' + 1, meanonly
				local nf = r(sum)

				post `handle' (`g') (`t') (round(`nr')) (round(`nf'))
			}
		}

		postclose `handle'
		frame risktable: use `riskdata', clear
		frame risktable: gen str40 risklabel = strtrim(string(nrisk, "%12.0fc")) + " (" + strtrim(string(nfail, "%12.0fc")) + ")"
		frame risktable: label define exposure 1 "DPP-4" 2 "GLP-1 RA" 3 "SGLT-2"
		frame risktable: label values exposure exposure
		frame risktable: drop nrisk nfail
		frame risktable: reshape wide risklabel,i(exposure) j(time)
		frame risktable: export excel using "${resultpath}\risktable.xlsx", sheet("`lab'")
		
		restore
	}
}



*********** Calculate effective sample size ESS for the comparison group
use ana_data, clear

putexcel set "${resultpath}\revision_results.xlsx", sheet(effectivesamplesize, replace) modify
putexcel B1 = "Actual N"
putexcel C1 = "ESS (capped at 10)"
putexcel D1 = "ESS (uncapped)"


***** IPTW, SMR weighting. 
gen iptw_smr_uncap=1 if exposure==2
gen iptw_smr_cap=1 if exposure==2
gen wuncap2=1 if exposure==2
gen wcap2=1 if exposure==2
forvalues i=1/2{		
	local v=cond(`i'==1, "gvd", "gvs")
	local c=cond(`i'==1, "DPP-4", "SGLT-2")
	quietly logit expo_`v' agesp* i.gender i.indexyr i.bc i.region_g i.education3 i.fam_panc i.t2donsetyr_g ///
	i.bmi_6y_g i.bp_sys_6y_g i.bp_dias_6y_g i.gfr_6y_g i.hba1c_6y_g i.triglyceride_6y_g i.ldl_6y_g ///
	i.smokinghabit_6y_g i.pre_a10a i.pre_a10ba02 i.pre_a10bb i.pre_a10bg i.pre_c02 i.pre_c03 i.pre_c07 i.pre_c08 i.pre_c09 i.pre_c10 i.pre_alco_dis i.pre_cereb_dis i.pre_cholelithiasis i.pre_diab_nephro i.pre_diab_retino i.pre_heart_fail i.pre_isc_heart_dis i.pre_jxx00 if expo_`v'!=.
	predict p_expo_`v'
	replace iptw_smr_uncap=p_expo_`v'/(1-p_expo_`v') if expo_`v'==0
	replace wuncap2=iptw_smr_uncap^2
	
	replace iptw_smr_cap=cond((iptw_smr_uncap>10 & iptw_smr_uncap!=.), 10, iptw_smr_uncap) if expo_`v'==0
	replace wcap2=iptw_smr_cap^2
	
	quietly count if expo_`v'==0
	local cnt = r(N)
	
	quietly summarize iptw_smr_cap if expo_`v' == 0
	local sum_wcap = r(sum)
	quietly summarize wcap2 if expo_`v' == 0
	local sum_wcap2 = r(sum)
	
	quietly summarize iptw_smr_uncap if expo_`v' == 0
	local sum_wuncap = r(sum)
	quietly summarize wuncap2 if expo_`v' == 0
	local sum_wuncap2 = r(sum)
	
	local k=`i'+1
	putexcel A`k' = "Comparison group `c'"
	putexcel B`k' = `cnt'
	putexcel C`k' = (`sum_wcap'^2 / `sum_wcap2')
	putexcel D`k' = (`sum_wuncap'^2 / `sum_wuncap2')
}
putexcel close
table () (exposure) [pw=iptw_smr_cap]
collect export "${resultpath}\revision_results.xlsx", sheet(effectivesamplesize) cell(A6) modify
table () (exposure) [pw=iptw_smr_uncap]
collect export "${resultpath}\revision_results.xlsx", sheet(effectivesamplesize) cell(A16) modify


twoway (kdensity iptw_smr_cap if expo_gvd == 0, color(navy%40)) ///
	   (kdensity iptw_smr_cap if expo_gvs == 0, color(maroon%40)), ///
	   legend(order(1 "DPP-4" 2 "SGLT-2")) xtitle("") ytitle("Density")
graph export "${resultpath}\IPTW SMR weighting capped distribution.png", replace

twoway (kdensity iptw_smr_uncap if expo_gvd == 0, color(navy%40)) ///
	   (kdensity iptw_smr_uncap if expo_gvs == 0, color(maroon%40)), ///
	   legend(order(1 "DPP-4" 2 "SGLT-2")) xtitle("") ytitle("Density")
graph export "${resultpath}\IPTW SMR weighting uncapped distribution.png", replace


range time 0 5 101 
foreach o in panc_all_wk852{
	egen dateout_`o'=rowmin(`o'_date emig_date deathdate fuend)
	format datein dateout_`o' fuend %d
	gen `o'_out1=(dateout_`o'==`o'_date)

	stset dateout_`o' [pweight=iptw_smr_uncap], fail(`o'_out1) enter(datein) origin(datein) scale(365.25)
		
	forvalues i=1/2{		
		local v=cond(`i'==1, "gvd", "gvs")
		local c=cond(`i'==1, "DPP-4", "SGLT-2")
		
		collect clear
		cap collect:stpm3 i.expo_`v' if expo_`v'!=.,  knots(5 35 65 95, percentile) tvc(expo_`v') knotstvc(10 50 90, percentile) scale(lncumhazard) eform 
		if _rc != 0 {
				di as error "Model failed. Skipping the analysis for outcome `o' and comparison between GLP-1 vs `c'"
				putexcel set "${resultpath}\revision_results.xlsx", sheet(uncapped_`o'_`v') modify
				putexcel A1="Model failed. Skipping the analysis for outcome `o' and comparison between GLP-1 vs `c'"
				putexcel close
				continue
		} 
		collect layout (colname) (result[_r_b _r_p _r_lb _r_ub])
		collect export "${resultpath}\revision_results.xlsx", sheet(uncapped_`o'_`v', replace) cell(A45) modify
		
		***Standardised incidence
		standsurv if expo_`v'!=., failure at1(expo_`v' 0) at2(expo_`v' 1) timevar(time) ///
		atvar(`v'_std glp_`v'_std) contrast(difference) contrastvars(riskdiff_glp_`v') ci frame(msurv, replace)  //indweights(iptw_smr)    weights should not be used again in prediction.
				
		***Derive incidences and differences Table
		foreach var in `v'_std `v'_std_lci `v'_std_uci glp_`v'_std glp_`v'_std_lci glp_`v'_std_uci riskdiff_glp_`v' riskdiff_glp_`v'_lci riskdiff_glp_`v'_uci {
			frame msurv{ 
				cap drop `var'_100
				gen `var'_100 = `var'*100
			}
		}

		frame msurv: export excel time *_100 using "${resultpath}\revision_results.xlsx" if inlist(time, 0.5, 1, 1.5, 2, 2.5, 3, 3.5, 4, 4.5, 5), sheet(uncapped_`o'_`v', modify)  firstrow(variables)
		
		***Standardised hazard
		standsurv if expo_`v'!=., hazard at1(expo_`v' 0) at2(expo_`v' 1) timevar(time) ///
		atvar(`v'_std glp_`v'_std) contrast(ratio) contrastvars(hr_glp_`v') ci frame(mhaza, replace)  //indweights(iptw_smr)    weights should not be used again in prediction.

		frame mhaza: export excel time * using "${resultpath}\revision_results.xlsx" if inlist(time, 0.5, 1, 1.5, 2, 2.5, 3, 3.5, 4, 4.5, 5), sheet(uncapped_`o'_`v', modify) cell(A15) firstrow(variables)
		
		***Standardised incidence and RR
		standsurv if expo_`v'!=., failure at1(expo_`v' 0) at2(expo_`v' 1) timevar(time) ///
		atvar(`v'_std glp_`v'_std) contrast(ratio) contrastvars(riskratio_glp_`v') ci frame(msurvr, replace)  //indweights(iptw_smr)    weights should not be used again in prediction.

		foreach var in `v'_std `v'_std_lci `v'_std_uci glp_`v'_std glp_`v'_std_lci glp_`v'_std_uci {
			frame msurvr{ 
				cap drop `var'_100
				gen `var'_100 = `var'*100
			}
		}
		
		frame msurvr: export excel using "${resultpath}\revision_results.xlsx" if inlist(time, 0.5, 1, 1.5, 2, 2.5, 3, 3.5, 4, 4.5, 5), sheet(uncapped_`o'_`v', modify) cell(A30) firstrow(variables)
	}
}	




************* more groups in BMI grouping
use ana_data, clear
egen bmi_6y_gnew=cut(bmi_6y_num), at(0,18.5,25,30,35,40,999) label
replace bmi_6y_gnew=. if bmi_6y_gnew==0
replace bmi_6y_gnew=99 if bmi_6y_gnew==.

***** IPTW, SMR weighting. 
gen iptw_smr1=1 if exposure==2
forvalues i=1/2{		
	local v=cond(`i'==1, "gvd", "gvs")
	local c=cond(`i'==1, "DPP-4", "SGLT-2")
	quietly logit expo_`v' agesp* i.gender i.indexyr i.bc i.region_g i.education3 i.fam_panc i.t2donsetyr_g ///
	i.bmi_6y_gnew i.bp_sys_6y_g i.bp_dias_6y_g i.gfr_6y_g i.hba1c_6y_g i.triglyceride_6y_g i.ldl_6y_g ///
	i.smokinghabit_6y_g i.pre_a10a i.pre_a10ba02 i.pre_a10bb i.pre_a10bg i.pre_c02 i.pre_c03 i.pre_c07 i.pre_c08 i.pre_c09 i.pre_c10 i.pre_alco_dis i.pre_cereb_dis i.pre_cholelithiasis i.pre_diab_nephro i.pre_diab_retino i.pre_heart_fail i.pre_isc_heart_dis i.pre_jxx00 if expo_`v'!=.
	predict p_expo_`v'
	replace iptw_smr1=p_expo_`v'/(1-p_expo_`v') if expo_`v'==0
	replace iptw_smr1=cond((iptw_smr1>10 & iptw_smr1!=.), 10, iptw_smr1) if expo_`v'==0
}

range time 0 5 101 
foreach o in panc_all_wk852{
	egen dateout_`o'=rowmin(`o'_date emig_date deathdate fuend)
	format datein dateout_`o' fuend %d
	gen `o'_out1=(dateout_`o'==`o'_date)

	stset dateout_`o' [pweight=iptw_smr1], fail(`o'_out1) enter(datein) origin(datein) scale(365.25)
		
	forvalues i=1/2{		
		local v=cond(`i'==1, "gvd", "gvs")
		local c=cond(`i'==1, "DPP-4", "SGLT-2")
		
		collect clear
		cap collect:stpm3 i.expo_`v' if expo_`v'!=.,  knots(5 35 65 95, percentile) tvc(expo_`v') knotstvc(10 50 90, percentile) scale(lncumhazard) eform 
		if _rc != 0 {
				di as error "Model failed. Skipping the analysis for outcome `o' and comparison between GLP-1 vs `c'"
				putexcel set "${resultpath}\revision_results.xlsx", sheet(bmi_gnew_`o'_`v') modify
				putexcel A1="Model failed. Skipping the analysis for outcome `o' and comparison between GLP-1 vs `c'"
				putexcel close
				continue
		} 
		collect layout (colname) (result[_r_b _r_p _r_lb _r_ub])
		collect export "${resultpath}\revision_results.xlsx", sheet(bmi_gnew_`o'_`v', replace) cell(A45) modify
		
		***Standardised incidence
		standsurv if expo_`v'!=., failure at1(expo_`v' 0) at2(expo_`v' 1) timevar(time) ///
		atvar(`v'_std glp_`v'_std) contrast(difference) contrastvars(riskdiff_glp_`v') ci frame(msurv, replace)  //indweights(iptw_smr)    weights should not be used again in prediction.
				
		***Derive incidences and differences Table
		foreach var in `v'_std `v'_std_lci `v'_std_uci glp_`v'_std glp_`v'_std_lci glp_`v'_std_uci riskdiff_glp_`v' riskdiff_glp_`v'_lci riskdiff_glp_`v'_uci {
			frame msurv{ 
				cap drop `var'_100
				gen `var'_100 = `var'*100
			}
		}

		frame msurv: export excel time *_100 using "${resultpath}\revision_results.xlsx" if inlist(time, 0.5, 1, 1.5, 2, 2.5, 3, 3.5, 4, 4.5, 5), sheet(bmi_gnew_`o'_`v', modify)  firstrow(variables)
		
		***Standardised hazard
		standsurv if expo_`v'!=., hazard at1(expo_`v' 0) at2(expo_`v' 1) timevar(time) ///
		atvar(`v'_std glp_`v'_std) contrast(ratio) contrastvars(hr_glp_`v') ci frame(mhaza, replace)  //indweights(iptw_smr)    weights should not be used again in prediction.

		frame mhaza: export excel using "${resultpath}\revision_results.xlsx" if inlist(time, 0.5, 1, 1.5, 2, 2.5, 3, 3.5, 4, 4.5, 5), sheet(bmi_gnew_`o'_`v', modify) cell(A15) firstrow(variables)
		
		***Standardised incidence and RR
		standsurv if expo_`v'!=., failure at1(expo_`v' 0) at2(expo_`v' 1) timevar(time) ///
		atvar(`v'_std glp_`v'_std) contrast(ratio) contrastvars(riskratio_glp_`v') ci frame(msurvr, replace)  //indweights(iptw_smr)    weights should not be used again in prediction.

		foreach var in `v'_std `v'_std_lci `v'_std_uci glp_`v'_std glp_`v'_std_lci glp_`v'_std_uci {
			frame msurvr{ 
				cap drop `var'_100
				gen `var'_100 = `var'*100
			}
		}
		
		frame msurvr: export excel using "${resultpath}\revision_results.xlsx" if inlist(time, 0.5, 1, 1.5, 2, 2.5, 3, 3.5, 4, 4.5, 5), sheet(bmi_gnew_`o'_`v', modify) cell(A30) firstrow(variables)
	}
}	




******  what if wash-out for pancreatitis is 6 years.
import sas "pre6y_panc_wK852.sas7bdat", clear
rename *, lower
replace pre6y_panc_wk852=0 if pre6y_panc_wk852==.
merge 1:1 lopnr using ana_data, keep(2 3) nogen
dtable i.pre6y_panc_wk852, export("${resultpath}\revision_results.xlsx", sheet(panc_6y_washout_panc_all_wk852) modify)

keep if pre6y_panc_wk852==0

***** IPTW, SMR weighting. 
gen iptw_smr1=1 if exposure==2
forvalues i=1/2{		
	local v=cond(`i'==1, "gvd", "gvs")
	local c=cond(`i'==1, "DPP-4", "SGLT-2")
	quietly logit expo_`v' agesp* i.gender i.indexyr i.bc i.region_g i.education3 i.fam_panc i.t2donsetyr_g ///
	i.bmi_6y_g i.bp_sys_6y_g i.bp_dias_6y_g i.gfr_6y_g i.hba1c_6y_g i.triglyceride_6y_g i.ldl_6y_g ///
	i.smokinghabit_6y_g i.pre_a10a i.pre_a10ba02 i.pre_a10bb i.pre_a10bg i.pre_c02 i.pre_c03 i.pre_c07 i.pre_c08 i.pre_c09 i.pre_c10 i.pre_alco_dis i.pre_cereb_dis i.pre_cholelithiasis i.pre_diab_nephro i.pre_diab_retino i.pre_heart_fail i.pre_isc_heart_dis i.pre_jxx00 if expo_`v'!=.
	predict p_expo_`v'
	replace iptw_smr1=p_expo_`v'/(1-p_expo_`v') if expo_`v'==0
	replace iptw_smr1=cond((iptw_smr1>10 & iptw_smr1!=.), 10, iptw_smr1) if expo_`v'==0
}

table () (exposure) () [pw=iptw_smr1]

range time 0 5 101 
foreach o in panc_all_wk852{
	egen dateout_`o'=rowmin(`o'_date emig_date deathdate fuend)
	format datein dateout_`o' fuend %d
	gen `o'_out1=(dateout_`o'==`o'_date)

	stset dateout_`o' [pweight=iptw_smr1], fail(`o'_out1) enter(datein) origin(datein) scale(365.25)
		
	forvalues i=1/2{		
		local v=cond(`i'==1, "gvd", "gvs")
		local c=cond(`i'==1, "DPP-4", "SGLT-2")
		
		collect clear
		cap collect:stpm3 i.expo_`v' if expo_`v'!=.,  knots(5 35 65 95, percentile) tvc(expo_`v') knotstvc(10 50 90, percentile) scale(lncumhazard) eform 
		if _rc != 0 {
				di as error "Model failed. Skipping the analysis for outcome `o' and comparison between GLP-1 vs `c'"
				putexcel set "${resultpath}\revision_results.xlsx", sheet(p6ywash_`o'_`v') modify
				putexcel A1="Model failed. Skipping the analysis for outcome `o' and comparison between GLP-1 vs `c'"
				putexcel close
				continue
		} 
		collect layout (colname) (result[_r_b _r_p _r_lb _r_ub])
		collect export "${resultpath}\revision_results.xlsx", sheet(p6ywash_`o'_`v', replace) cell(A45) modify
		
		***Standardised incidence
		standsurv if expo_`v'!=., failure at1(expo_`v' 0) at2(expo_`v' 1) timevar(time) ///
		atvar(`v'_std glp_`v'_std) contrast(difference) contrastvars(riskdiff_glp_`v') ci frame(msurv, replace)  //indweights(iptw_smr)    weights should not be used again in prediction.
				
		***Derive incidences and differences Table
		foreach var in `v'_std `v'_std_lci `v'_std_uci glp_`v'_std glp_`v'_std_lci glp_`v'_std_uci riskdiff_glp_`v' riskdiff_glp_`v'_lci riskdiff_glp_`v'_uci {
			frame msurv{ 
				cap drop `var'_100
				gen `var'_100 = `var'*100
			}
		}

		frame msurv: export excel time *_100 using "${resultpath}\revision_results.xlsx" if inlist(time, 0.5, 1, 1.5, 2, 2.5, 3, 3.5, 4, 4.5, 5), sheet(p6ywash_`o'_`v', modify)  firstrow(variables)
		
		***Standardised hazard
		standsurv if expo_`v'!=., hazard at1(expo_`v' 0) at2(expo_`v' 1) timevar(time) ///
		atvar(`v'_std glp_`v'_std) contrast(ratio) contrastvars(hr_glp_`v') ci frame(mhaza, replace)  //indweights(iptw_smr)    weights should not be used again in prediction.

		frame mhaza: export excel using "${resultpath}\revision_results.xlsx" if inlist(time, 0.5, 1, 1.5, 2, 2.5, 3, 3.5, 4, 4.5, 5), sheet(p6ywash_`o'_`v', modify) cell(A15) firstrow(variables)
		
		***Standardised incidence and RR
		standsurv if expo_`v'!=., failure at1(expo_`v' 0) at2(expo_`v' 1) timevar(time) ///
		atvar(`v'_std glp_`v'_std) contrast(ratio) contrastvars(riskratio_glp_`v') ci frame(msurvr, replace)  //indweights(iptw_smr)    weights should not be used again in prediction.

		foreach var in `v'_std `v'_std_lci `v'_std_uci glp_`v'_std glp_`v'_std_lci glp_`v'_std_uci {
			frame msurvr{ 
				cap drop `var'_100
				gen `var'_100 = `var'*100
			}
		}
		
		frame msurvr: export excel using "${resultpath}\revision_results.xlsx" if inlist(time, 0.5, 1, 1.5, 2, 2.5, 3, 3.5, 4, 4.5, 5), sheet(p6ywash_`o'_`v', modify) cell(A30) firstrow(variables)
	}
}	


*********** subtype of AP by treatment group
use ana_data_iptw, clear
		table (pantype) (exposure) [pw=iptw_smr]
		collect export "${resultpath}\revision_results.xlsx", sheet(AP_types, replace) modify


**** Combining graphs now ITT into one and PP into another
*** separate gvd and gvs, do one above wide panel for curves, and below two plots for rd and rr respectively
graph use "${resultpath}\figure\fig2A_curv"
*** do edits. 
*** remove title
*** edit y-axis title "Weighted risk (95% CI)"
*** change size, xsize(10)
*** change GLP-1 to "GLP-1 receptor agonist" and SGLT-2 to "SGLT-2 inhibitor" and DPP-4 to "DPP-4 inhibitor"
graph use "${resultpath}\figure\fig3A_curv"
graph use "${resultpath}\figure\fig2B_curv"
graph use "${resultpath}\figure\fig3B_curv"

graph use "${resultpath}\figure\fig2A_rd"  
*** change color to black%30 and black 
*** remove the percentage point from the "Risk difference, percentage point (95% CI)" so it just says "Risk difference (95% CI). 
*** Add % to the Y-axis
graph use "${resultpath}\figure\fig3A_rd" 
graph use "${resultpath}\figure\fig2B_rd" 
graph use "${resultpath}\figure\fig3B_rd" 



*** Now draw the RR fig. and combine with RD and curve fig.
forvalues i=1/2{
// 	local v=cond(`i'==2, "gvd", "gvs")
// 	local c=cond(`i'==2, "DPP-4", "SGLT-2")
//	
// 	import excel using "S:\HI-SPEED projects\HS26_02 GLP1 pancreatitis\output\20260825\TableS3_ITT_tvc_failratio.xlsx", sheet("panc_all_wk852_`v'") firstrow clear
//	
// 	tw (rcap riskratio_glp_`v'_lci riskratio_glp_`v'_uci time if inlist(time, 0, 0.5, 1, 1.5, 2, 2.5, 3, 3.5, 4, 4.5, 5), color(black%30)) ///
// 		(scatter riskratio_glp_`v' time if inlist(time, 0, 0.5, 1, 1.5, 2, 2.5, 3, 3.5, 4, 4.5, 5), color(black)), scheme(tab2) ///
// 		ylabel(0.5(0.5)2, nogrid format(%4.1f) labsize(small)) ///
// 		xtitle("{bf:Years since treatment initiation}", size(*1) margin(0 0 0 0)) ///
// 		xlabel(0(1)5, nogrid format(%9.1f) labsize(small)) ///
// 		xscale(range(0 5)) ///
// 		ytitle("{bf:Risk ratio} {bf:(95% CI)}", size(*1) margin(0 0 0 0)) ///
// 		yscale(range(0.5 2)  log) ///
// 		yline(1, lp(-)) ///
// 		legend(off) ///
// 		name("fig`=`i'+1'A_rr", replace)
//		
// 	graph save "${resultpath}\figure\fig`=`i'+1'A_rr", replace
		
	graph combine "${resultpath}\figure\fig`=`i'+1'A_rd" "${resultpath}\figure\fig`=`i'+1'A_rr", cols(2) name(bottom, replace)
	graph combine "${resultpath}\figure\fig`=`i'+1'A_curv" bottom, cols(1) ysize(8)
	
	graph save "${resultpath}\figure\fig`=`i'+1'A_final", replace
	graph export "${resultpath}\figure\fig`=`i'+1'A_final.svg", replace
}

forvalues i=1/2{
	graph combine "${resultpath}\figure\fig`=`i'+1'B_rd" "${resultpath}\figure\fig`=`i'+1'B_rr", cols(2) name(bottom, replace)
	graph combine "${resultpath}\figure\fig`=`i'+1'B_curv" bottom, cols(1) ysize(8)
	
	graph save "${resultpath}\figure\fig`=`i'+1'B_final", replace
	graph export "${resultpath}\figure\fig`=`i'+1'B_final.svg", replace
}


************* other biliary outcomes
import sas using biliary.sas7bdat, clear
rename *, lower
drop op1-op30
format %d indatum 
keep if indatum<=mdy(11,30,2025)
gsort lopnr indatum
save otherbiliary, replace

use ana_data_iptw, clear
egen dateout_panc_all_wk852=rowmin(panc_all_wk852_date emig_date deathdate fuend)
keep lopnr exposure indexdate dateout_panc_all_wk852 iptw_smr 
merge 1:m lopnr using otherbiliary, nogen keep(1 3)
replace indatum=. if indatum<=indexdate | indatum>dateout_panc_all_wk852
replace otherbiliary=0 if otherbiliary==. | indatum==.
bys lopnr: egen max_otherbiliary=max(otherbiliary)
egen x=tag(lopnr)

		table (max_otherbiliary) (exposure) [pw=iptw_smr] if x==1
		collect export "${resultpath}\revision_results.xlsx", sheet(otherbiliary, replace) modify


************ competeing risk analysis (CIF, crude risk)
use ana_data_iptw, clear
range time  0 5 101 
foreach o in panc_all_wk852{
	egen dateout_`o'=rowmin(`o'_date emig_date deathdate fuend)
	format datein dateout_`o' fuend %d
	gen `o'_out1=(dateout_`o'==`o'_date)
	
	stset dateout_`o' [pweight=iptw_smr], fail(`o'_out1) enter(datein) origin(datein) scale(365.25)
	
	gen out_satt_`o'=_d
	replace out_satt_`o'=2 if out_satt_`o'==0 & dateout_`o'==deathdate
	replace out_satt_`o'=3 if out_satt_`o'==0 & dateout_`o'==emig_date
	replace out_satt_`o'=4 if out_satt_`o'==0 & dateout_`o'==fuend
	label define out_satt_`o' 1 "Event" 2 "Death" 3 "Emigration" 4 "EndFollowUp"
	label values out_satt_`o' out_satt_`o'
	
	forvalues i=1/2{		
		local v=cond(`i'==1, "gvd", "gvs")
		local c=cond(`i'==1, "DPP-4", "SGLT-2")
		
		stset dateout_`o' [pweight=iptw_smr], fail(out_satt_`o'==1) enter(datein) origin(datein) scale(365.25)
		collect clear
		quietly stpm3 i.expo_`v' if expo_`v'!=.,  knots(5 35 65 95, percentile) tvc(expo_`v') knotstvc(10 50 90, percentile) scale(lncumhazard) eform 
		est store panc_`v'
		
		stset dateout_`o' [pweight=iptw_smr], fail(out_satt_`o'==2) enter(datein) origin(datein) scale(365.25)
		collect clear
		quietly stpm3 i.expo_`v' if expo_`v'!=.,  knots(5 35 65 95, percentile) tvc(expo_`v') knotstvc(10 50 90, percentile) scale(lncumhazard) eform 
		est store death_`v'
		
		
		***CIF using standsurv 
		standsurv if expo_`v'!=., crmodels(panc_`v' death_`v') cif timevar(time) verbose at1(expo_`v' 0) at2(expo_`v' 1) atvar(`v'_std glp_`v'_std) contrast(difference) contrastvars(cif_rd_`v') ci frame(cifrd, replace)  
		frame cifrd: export excel using "${resultpath}\revision_results.xlsx" if inlist(time, 0.5, 1, 1.5, 2, 2.5, 3, 3.5, 4, 4.5, 5), sheet(CIF_`o'_`v', modify) firstrow(variables)
		
		standsurv if expo_`v'!=., crmodels(panc_`v' death_`v') cif timevar(time) verbose at1(expo_`v' 0) at2(expo_`v' 1) atvar(`v'_std glp_`v'_std) contrast(ratio) contrastvars(cif_rr_`v') ci frame(cifrr, replace) 
		frame cifrr: export excel using "${resultpath}\revision_results.xlsx" if inlist(time, 0.5, 1, 1.5, 2, 2.5, 3, 3.5, 4, 4.5, 5), sheet(CIF_`o'_`v', modify) cell(A15) firstrow(variables)
		
	}
}	


************ diabetes duration, instead of onsetyear, and use splines.
use ana_data, clear
gen dm_dur=indexyr-t2donsetyr
mkspline dmdur_sp = dm_dur, cubic

gen iptw_smr1=1 if exposure==2
forvalues i=1/2{		
	local v=cond(`i'==1, "gvd", "gvs")
	local c=cond(`i'==1, "DPP-4", "SGLT-2")
	quietly logit expo_`v' agesp* i.gender i.indexyr i.bc i.region_g i.education3 i.fam_panc dmdur_sp* ///
	i.bmi_6y_g i.bp_sys_6y_g i.bp_dias_6y_g i.gfr_6y_g i.hba1c_6y_g i.triglyceride_6y_g i.ldl_6y_g ///
	i.smokinghabit_6y_g i.pre_a10a i.pre_a10ba02 i.pre_a10bb i.pre_a10bg i.pre_c02 i.pre_c03 i.pre_c07 i.pre_c08 i.pre_c09 i.pre_c10 i.pre_alco_dis i.pre_cereb_dis i.pre_cholelithiasis i.pre_diab_nephro i.pre_diab_retino i.pre_heart_fail i.pre_isc_heart_dis i.pre_jxx00 if expo_`v'!=.
	predict p_expo_`v'
	replace iptw_smr1=p_expo_`v'/(1-p_expo_`v') if expo_`v'==0
	replace iptw_smr1=cond((iptw_smr1>10 & iptw_smr1!=.), 10, iptw_smr1) if expo_`v'==0
}

range time 0 5 101 
foreach o in panc_all_wk852{
	egen dateout_`o'=rowmin(`o'_date emig_date deathdate fuend)
	format datein dateout_`o' fuend %d
	gen `o'_out1=(dateout_`o'==`o'_date)

	stset dateout_`o' [pweight=iptw_smr1], fail(`o'_out1) enter(datein) origin(datein) scale(365.25)
		
	forvalues i=1/2{		
		local v=cond(`i'==1, "gvd", "gvs")
		local c=cond(`i'==1, "DPP-4", "SGLT-2")
		
		collect clear
		cap collect:stpm3 i.expo_`v' if expo_`v'!=.,  knots(5 35 65 95, percentile) tvc(expo_`v') knotstvc(10 50 90, percentile) scale(lncumhazard) eform 
		if _rc != 0 {
				di as error "Model failed. Skipping the analysis for outcome `o' and comparison between GLP-1 vs `c'"
				putexcel set "${resultpath}\revision_results.xlsx", sheet(dmduration_`o'_`v') modify
				putexcel A1="Model failed. Skipping the analysis for outcome `o' and comparison between GLP-1 vs `c'"
				putexcel close
				continue
		} 
		collect layout (colname) (result[_r_b _r_p _r_lb _r_ub])
		collect export "${resultpath}\revision_results.xlsx", sheet(dmduration_`o'_`v', replace) cell(A45) modify
		
		***Standardised incidence
		standsurv if expo_`v'!=., failure at1(expo_`v' 0) at2(expo_`v' 1) timevar(time) ///
		atvar(`v'_std glp_`v'_std) contrast(difference) contrastvars(riskdiff_glp_`v') ci frame(msurv, replace)  //indweights(iptw_smr)    weights should not be used again in prediction.
				
		***Derive incidences and differences Table
		foreach var in `v'_std `v'_std_lci `v'_std_uci glp_`v'_std glp_`v'_std_lci glp_`v'_std_uci riskdiff_glp_`v' riskdiff_glp_`v'_lci riskdiff_glp_`v'_uci {
			frame msurv{ 
				cap drop `var'_100
				gen `var'_100 = `var'*100
			}
		}

		frame msurv: export excel time *_100 using "${resultpath}\revision_results.xlsx" if inlist(time, 0.5, 1, 1.5, 2, 2.5, 3, 3.5, 4, 4.5, 5), sheet(dmduration_`o'_`v', modify)  firstrow(variables)
		
		***Standardised hazard
		standsurv if expo_`v'!=., hazard at1(expo_`v' 0) at2(expo_`v' 1) timevar(time) ///
		atvar(`v'_std glp_`v'_std) contrast(ratio) contrastvars(hr_glp_`v') ci frame(mhaza, replace)  //indweights(iptw_smr)    weights should not be used again in prediction.

		frame mhaza: export excel using "${resultpath}\revision_results.xlsx" if inlist(time, 0.5, 1, 1.5, 2, 2.5, 3, 3.5, 4, 4.5, 5), sheet(dmduration_`o'_`v', modify) cell(A15) firstrow(variables)
		
		***Standardised incidence and RR
		standsurv if expo_`v'!=., failure at1(expo_`v' 0) at2(expo_`v' 1) timevar(time) ///
		atvar(`v'_std glp_`v'_std) contrast(ratio) contrastvars(riskratio_glp_`v') ci frame(msurvr, replace)  //indweights(iptw_smr)    weights should not be used again in prediction.

		foreach var in `v'_std `v'_std_lci `v'_std_uci glp_`v'_std glp_`v'_std_lci glp_`v'_std_uci {
			frame msurvr{ 
				cap drop `var'_100
				gen `var'_100 = `var'*100
			}
		}
		
		frame msurvr: export excel using "${resultpath}\revision_results.xlsx" if inlist(time, 0.5, 1, 1.5, 2, 2.5, 3, 3.5, 4, 4.5, 5), sheet(dmduration_`o'_`v', modify) cell(A30) firstrow(variables)
	}
}	
