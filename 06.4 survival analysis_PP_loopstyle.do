cd "S:\HI-SPEED projects\HS26_02 GLP1 pancreatitis\data\"
global resultpath "S:\HI-SPEED projects\HS26_02 GLP1 pancreatitis\output\20260825\PP_IPCW"
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

foreach o in panc_all_wk852 panc_all_wok852{   //UPDATE NEED
	forvalues i =1/2{		//UPDATE NEED
		local v=cond(`i'==1, "gvd", "gvs")
		local c=cond(`i'==1, "DPP-4", "SGLT-2")
		
		foreach fil in maxweight10 maxweight10_hazaHR maxweight10_failratio{
			putexcel set "${resultpath}\PP_IPCW_tvc_`fil'_`o'_`v'.xlsx", sheet(`o'_`v', replace) modify
			putexcel B1 = "HR"
			putexcel C1 = "Time"
			putexcel D1 = "`v'_std_100"
			putexcel E1 = "glp_`v'_std_100" 
			putexcel F1 = "riskdiff_glp_`v'_100"
			putexcel close 
		}
	
		forvalues k=1/300{   //UPDATE NEED
			di "Outcome `o', bootstrapping iteration `k'/300"
			use ana_data_iptw_pp, clear 
			keep if expo_`v' != .		

			if `k' > 1 {
				set seed `=`k'*10'
				bsample
			}
			
			egen info_censor_date=rowmin(discon_date drug_switch_date)
			replace info_censor_date = info_censor_date-1 if info_censor_date==`o'_date
			gen info_censor_date_psuedo = datein + 30*(ceil((info_censor_date - datein)/30))
			format info_censor_date_psuedo %d
			
			egen dateout=rowmin(`o'_date emig_date deathdate fuend info_censor_date_psuedo)
			format datein dateout fuend %d
			gen out1=(dateout==`o'_date)
			
			gen treat_vio = (info_censor_date_psuedo==dateout)
			
			gen fup = dateout - datein 
			gen time_event = `o'_date - datein 
			replace fup = 1 if fup==0
			replace time_event = 1 if time_event ==0
			
			gen fakeid=_n
			
			stset fup, fail(treat_vio) id(fakeid)
			quietly stpm3 i.expo_`v' agesp* i.gender i.indexyr i.bc i.region_g i.education3 i.fam_panc i.t2donsetyr_g ///
	i.bmi_6y_g i.bp_sys_6y_g i.bp_dias_6y_g i.gfr_6y_g i.hba1c_6y_g i.triglyceride_6y_g i.ldl_6y_g ///
	i.smokinghabit_6y_g i.pre_a10a i.pre_a10ba02 i.pre_a10bb i.pre_a10bg i.pre_c02 i.pre_c03 i.pre_c07 i.pre_c08 i.pre_c09 i.pre_c10 i.pre_alco_dis i.pre_cereb_dis i.pre_cholelithiasis i.pre_diab_nephro i.pre_diab_retino i.pre_heart_fail i.pre_isc_heart_dis i.pre_jxx00, scale(lncumhazard) df(5) eform //nolog		

			stsplit, at(failures)
			gen t_enter=_t0
			gen t_out=_t
			replace out1 = (t_out == time_event)			
							
			predict p, survival
			gen ipcw = 1 / p	
			
			gen total_weight=iptw_smr*ipcw
			replace total_weight=10 if total_weight>10
			
			stset t_out [pweight=total_weight], fail(out1) enter(t_enter) origin(time 0) scale(365.25) 
			
			***** Final Analysis
			***** TVC option
			cap drop time
			cap stpm3 i.expo_`v' if expo_`v'!=., knots(5 35 65 95, percentile) tvc(expo_`v') knotstvc(10 50 90, percentile) scale(lncumhazard) eform 
			
			if _rc!=0 {
				di as error "Model failed. Skipping the analysis for outcome `o' and comparison between GLP-1 vs `c'"
				foreach fil in maxweight10 maxweight10_hazaHR maxweight10_failratio{
					putexcel set "${resultpath}\PP_IPCW_tvc_`fil'_`o'_`v'.xlsx", sheet(`o'_`v') modify
					putexcel A`=`k'*11-9'="Run `k': Model failed. Skipping the analysis for outcome `o' and comparison between GLP-1 vs `c'"
					putexcel close
				}
				continue 	
			} 
			else {
				foreach fil in maxweight10 maxweight10_hazaHR maxweight10_failratio{
					putexcel set "${resultpath}\PP_IPCW_tvc_`fil'_`o'_`v'.xlsx", sheet(`o'_`v') modify
					putexcel A`=`k'*11-9':A`=`k'*11' = "Run `k'"
					putexcel B`=`k'*11-9' = _b[1.expo_`v']
					putexcel close
				}
			}

					
			use ana_data_iptw_pp, clear 
			keep if expo_`v' != .	
			range time 0 5 101 	

			***Standardised incidence
			standsurv, failure at1(expo_`v' 0) at2(expo_`v' 1) timevar(time) ///
			atvar(`v'_std glp_`v'_std) contrast(difference) contrastvars(riskdiff_glp_`v') frame(msurv, replace)
			
			***Derive incidences and differences Table
			foreach var in `v'_std glp_`v'_std riskdiff_glp_`v' {
				frame msurv{ 
					cap drop `var'_100
					gen `var'_100 = `var'*100
				}
			}

			frame msurv: export excel time *_100 using "${resultpath}\PP_IPCW_tvc_maxweight10_`o'_`v'.xlsx" if inlist(time, 0.5, 1, 1.5, 2, 2.5, 3, 3.5, 4, 4.5, 5), sheet(`o'_`v', modify) cell(C`=`k'*11-9') 
		
			***Standardised hazard
			standsurv, hazard at1(expo_`v' 0) at2(expo_`v' 1) timevar(time) ///
			atvar(`v'_std glp_`v'_std) contrast(ratio) contrastvars(hr_glp_`v') ci frame(mhaza, replace)  

			frame mhaza: export excel time * using "${resultpath}\PP_IPCW_tvc_maxweight10_hazaHR_`o'_`v'.xlsx" if inlist(time, 0.5, 1, 1.5, 2, 2.5, 3, 3.5, 4, 4.5, 5), sheet(`o'_`v', modify)  cell(C`=`k'*11-9') firstrow(variables)
			
			***Standardised incidence and RR
			standsurv, failure at1(expo_`v' 0) at2(expo_`v' 1) timevar(time) ///
			atvar(`v'_std glp_`v'_std) contrast(ratio) contrastvars(riskratio_glp_`v') ci frame(msurvr, replace) 

			foreach var in `v'_std `v'_std_lci `v'_std_uci glp_`v'_std glp_`v'_std_lci glp_`v'_std_uci {
				frame msurvr{ 
					cap drop `var'_100
					gen `var'_100 = `var'*100
				}
			}
			
			frame msurvr: export excel using "${resultpath}\PP_IPCW_tvc_maxweight10_failratio_`o'_`v'.xlsx" if inlist(time, 0.5, 1, 1.5, 2, 2.5, 3, 3.5, 4, 4.5, 5), sheet(`o'_`v', modify)  cell(C`=`k'*11-9') firstrow(variables)
		}
	}
}
