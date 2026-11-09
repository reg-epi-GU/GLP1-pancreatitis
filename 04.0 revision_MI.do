cd "S:\HI-SPEED projects\HS26_02 GLP1 pancreatitis\data\"
global resultpath "S:\HI-SPEED projects\HS26_02 GLP1 pancreatitis\output\20261015"

cap mkdir "$resultpath"
cap mkdir "$resultpath\MI_results\"

use ana_data, clear
drop if region_g==99

vl create outcome = (panc_all_wk852)

* Create survival information for imputation
foreach o of global outcome {
    egen dateout_`o' = rowmin(`o'_date emig_date deathdate fuend)
    gen `o'_out1 = !missing(`o'_date) & dateout_`o' == `o'_date
    gen double fu_`o' = (dateout_`o' - datein)/365.25
    stset fu_`o', failure(`o'_out1)

    * Nelson-Aalen cumulative hazard
    sts gen H_`o' = na
}

global mi_survival

foreach o of global outcome {
    global mi_survival ///
        $mi_survival `o'_out1 H_`o'
}

vl create covar =  (gender indexyr bc region_g education3 fam_panc t2donsetyr_g bmi_6y_g bp_sys_6y_g bp_dias_6y_g gfr_6y_g hba1c_6y_g triglyceride_6y_g ldl_6y_g smokinghabit_6y_g pre_a10a pre_a10ba02 pre_a10bb pre_a10bg pre_c02 pre_c03 pre_c07 pre_c08 pre_c09 pre_c10 pre_alco_dis pre_cereb_dis pre_cholelithiasis pre_diab_nephro pre_diab_retino pre_heart_fail pre_isc_heart_dis pre_jxx00)
vl substitute covar_i = i.covar

foreach var of varlist $covar {
	di "`var'"
	replace `var'=. if `var'==99
}


** Multiple imputation
misstable summarize $covar
misstable patterns $covar

vl create missing=(smokinghabit_6y_g education3 bmi_6y_g bp_sys_6y_g bp_dias_6y_g gfr_6y_g triglyceride_6y_g ldl_6y_g hba1c_6y_g)
vl create nonmiss= covar - missing
vl substitute nonmiss_i = i.nonmiss

mi set mlong
mi register imputed $missing

set maxiter 200

mi impute chained ///
    (mlogit, augment) smokinghabit_6y_g ///
    (ologit, augment) education3 bmi_6y_g bp_sys_6y_g bp_dias_6y_g gfr_6y_g triglyceride_6y_g ldl_6y_g hba1c_6y_g ///
    = i.exposure $mi_survival agesp* $nonmiss, ///
    add(50) burnin(20) rseed(12345)

compress
save mi_master, replace
	
	
	
* =====================================================
* In each MI: IPTW (SMR) + stpm3 + standsurv
* Outcome: panc_all_wk852
* =====================================================	
forvalues m = 1/50 {
	di as text "=================================="
    di as result "Processing imputation `m' of 50"
    di as text "=================================="
	
	use mi_master, clear
	mi extract `m', clear
	
	***** IPTW, SMR weighting. 
	gen iptw_smr=1 if exposure==2
	forvalues i=1/2{		
		local v=cond(`i'==1, "gvd", "gvs")
		local c=cond(`i'==1, "DPP-4", "SGLT-2")
		quietly logit expo_`v' agesp* $covar_i if expo_`v'!=.
		predict p_expo_`v'
		replace iptw_smr=p_expo_`v'/(1-p_expo_`v') if expo_`v'==0
		replace iptw_smr=10 if iptw_smr>10 & iptw_smr<.
	}

	cap drop time
	range time 0 5 101 
	foreach o in panc_all_wk852 {
		forvalues i=1/2{		
			local v=cond(`i'==1, "gvd", "gvs")
			local c=cond(`i'==1, "DPP-4", "SGLT-2")
			
			cap drop dateout_`o' `o'_out1
			egen dateout_`o'=rowmin(`o'_date emig_date deathdate fuend)
			format datein dateout_`o' fuend %d
			gen `o'_out1=(dateout_`o'==`o'_date)

			stset dateout_`o' [pweight=iptw_smr], fail(`o'_out1) enter(datein) origin(datein) scale(365.25)
			
			collect clear
			cap collect:stpm3 i.expo_`v' if expo_`v'!=.,  knots(5 35 65 95, percentile) tvc(expo_`v') knotstvc(10 50 90, percentile) scale(lncumhazard) eform 
			if _rc != 0 {
				di as error "Model failed in Imputation `m': comparison `v'. Skipping the analysis for outcome `o' and comparison between GLP-1 vs `c'"
				putexcel set "${resultpath}\TableS3_ITT_tvc.xlsx", sheet(`o'_`v') modify
				putexcel A1="Model failed in Imputation `m': comparison `v'. Skipping the analysis for outcome `o' and comparison between GLP-1 vs `c'"
				putexcel close
				continue
			} 
			
				
			***Standardised incidence
			standsurv if expo_`v'!=., failure at1(expo_`v' 0) at2(expo_`v' 1) timevar(time) ///
			atvar(`v'_std glp_`v'_std) contrast(difference) contrastvars(riskdiff_glp_`v') ci frame(msurv, replace)  //indweights(iptw_smr)    weights should not be used again in prediction.
			frame msurv: save "${resultpath}\MI_results\inc_`v'_m`m'.dta", replace

			***Standardised hazard
			standsurv if expo_`v'!=., hazard at1(expo_`v' 0) at2(expo_`v' 1) timevar(time) ///
			atvar(`v'_std glp_`v'_std) contrast(ratio) contrastvars(hr_glp_`v') ci frame(mhaza, replace)  //indweights(iptw_smr)    weights should not be used again in prediction.
			frame mhaza: save "${resultpath}\MI_results\hr_`v'_m`m'.dta", replace				
			
			***Standardised incidence and RR
			standsurv if expo_`v'!=., failure at1(expo_`v' 0) at2(expo_`v' 1) timevar(time) ///
			atvar(`v'_std glp_`v'_std) contrast(ratio) contrastvars(riskratio_glp_`v') ci frame(msurvr, replace)  //indweights(iptw_smr)    weights should not be used again in prediction.
			frame msurvr: save "${resultpath}\MI_results\rr_`v'_m`m'.dta", replace
		}
	}
}


* Pool risk differences across 50 imputations
forvalues i=1/2{
	local v=cond(`i'==1, "gvd", "gvs")
	
	use "${resultpath}\MI_results\rr_`v'_m1.dta", clear
	gen imputation = 1
	
	foreach n in glp_`v' glp_`v'_lci glp_`v'_uci{
		gen lnrr_`n' = ln(riskratio_`n')
	}
	keep imputation time lnrr_*	
	
	merge 1:1 time using "${resultpath}\MI_results\inc_`v'_m1.dta", nogen
	save combined_`v', replace
	
	forvalues m = 2/50 {
		use "${resultpath}\MI_results\rr_`v'_m`m'.dta", clear
		gen imputation = `m'
		
		foreach n in glp_`v' glp_`v'_lci glp_`v'_uci{
			gen lnrr_`n' = ln(riskratio_`n')
		}
		keep imputation time lnrr_*	
		
		merge 1:1 time using "${resultpath}\MI_results\inc_`v'_m`m'.dta", nogen
		append using combined_`v'
		save combined_`v', replace
	}
}


merge 1:1 imputation time using combined_gvd, nogen
save combined, replace



*** Final pooling using Rubin's rules
use combined, clear
****** do this for 7 estimates (each with 95CI): risk estimate (3 groups) risk difference (2 comparison) risk ratio(2 comparison)
foreach var in glp_gvd_std gvd_std riskdiff_glp_gvd lnrr_glp_gvd glp_gvs_std gvs_std riskdiff_glp_gvs lnrr_glp_gvs{

    * Within-imputation variance
    gen se_`var' = (`var'_uci - `var'_lci) / (2*1.96)
	gen within_`var' = se_`var'^2

	* Rubin's pooling components
	bysort time: egen pooled_`var' = mean(`var')
	bysort time: egen W_`var' = mean(within_`var')
	bysort time: egen Bsd_`var' = sd(`var')
	bysort time: egen M_`var' = count(`var')

	gen B_`var' = Bsd_`var'^2

	gen totalvar_`var' = W_`var' + (1 + 1/M_`var')*B_`var'
	gen pooledse_`var' = sqrt(totalvar_`var')

	* Degrees of freedom
	gen df_`var' = (M_`var'-1) * (1 + W_`var'/((1+1/M_`var')*B_`var'))^2 if B_`var' > 0

	gen crit_`var' = invttail(df_`var', .025) if B_`var' > 0
	replace crit_`var' = invnormal(.975) if B_`var' == 0

	* Confidence intervals
	gen lci_`var' = pooled_`var' - crit_`var'*pooledse_`var'
	gen uci_`var' = pooled_`var' + crit_`var'*pooledse_`var'
		
    if strpos("`var'", "lnrr_") == 1 {
        gen str50 result_`var' = strtrim(string(exp(pooled_`var'), "%9.3f")) + " (" + strtrim(string(exp(lci_`var'), "%9.3f")) + " - " + strtrim(string(exp(uci_`var'), "%9.3f")) + ")"
    }
    else {
		gen str50 result_`var' = strtrim(string(pooled_`var'*100, "%9.3f")) + " (" + strtrim(string(lci_`var'*100, "%9.3f")) + " - " + strtrim(string(uci_`var'*100, "%9.3f")) + ")"
    }
}

rename (result_lnrr_glp_gvd result_lnrr_glp_gvs result_gvd_std result_gvs_std) (result_riskratio_glp_gvd result_riskratio_glp_gvs result_dpp4_std result_sglt2_std)
* Keep one observation per time point
egen tag = tag(time)
keep if tag

save "${resultpath}\MI_pooled_all.dta", replace


*** result cleanning
keep time result_*
keep if inlist(time,  1, 2, 3, 4, 5)
export excel using "${resultpath}\MI_pooled_all.xlsx", replace firstrow(variables)
		