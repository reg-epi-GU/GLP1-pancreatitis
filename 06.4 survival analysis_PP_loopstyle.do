cd "S:\HI-SPEED projects\HS26_02 GLP1 pancreatitis\data\"
global resultpath "S:\HI-SPEED projects\HS26_02 GLP1 pancreatitis\output\20260422\PP_IPCW"
cap mkdir "${resultpath}"


foreach o in pancreatitis_all /*pancreatitis_hosp cout1_all cout2_all*/ {    //UPDATE NEED
	forvalues i =1/2{		//UPDATE NEED
		local v=cond(`i'==1, "gvd", "gvs")
		local c=cond(`i'==1, "DPP4", "SGLT2")
		
		putexcel set "${resultpath}\TableS3_PP_IPCW_ph_new.xlsx", sheet(`o'_`v', replace) modify
		putexcel B1 = "HR"
		putexcel C1 = "Time"
		putexcel D1 = "`c'_std_100"
		putexcel E1 = "glp_`c'_std_100" 
		putexcel F1 = "riskdiff_glp_`c'_100"
		putexcel close 
	
		forvalues k=1/301{   //UPDATE NEED
			di "Outcome `o', bootstrapping iteration `k'/301"
			use ana_data_perprotocol, clear 
			keep `o'_date expo emig_date deathdate fuend drug_switch_date datein expo_`v' agesp* gender indexyr bc region_g education3 fam_panc t2donsetyr_g bmi_g smokinghabit_g hba1c_g preindex_lipidlowerdrug preindex_alcoholdisorder preindex_gallstones preindex_abdomsurgery
			rename preindex_* p_*
			keep if expo=="glp1" | expo==lower("`c'")		

			if `k' > 1 {
				set seed `=`k'*10'
				bsample
			}
			
			replace drug_switch_date = drug_switch_date-1 if drug_switch_date==`o'_date
			gen drug_switch_date_psuedo = datein + 30*(ceil((drug_switch_date - datein)/30))
			format drug_switch_date_psuedo %d
			
			egen dateout=rowmin(`o'_date emig_date deathdate fuend drug_switch_date_psuedo)
			format datein dateout fuend %d
			gen out1=(dateout==`o'_date)
			
			gen treat_vio = (drug_switch_date_psuedo==dateout)
			
			gen fup = dateout - datein 
			gen time_event = `o'_date - datein 
			replace fup = 1 if fup==0
			replace time_event = 1 if time_event ==0
			
			gen fakeid=_n
			
			stset fup, fail(treat_vio) id(fakeid)
			quietly stpm3 i.expo_`v' agesp* i.gender i.indexyr i.bc ib2.region_g i.education3 i.fam_panc i.t2donsetyr_g ib2.bmi_g i.smokinghabit_g i.hba1c_g i.p_lipidlowerdrug i.p_alcoholdisorder i.p_gallstones i.p_abdomsurgery, scale(lncumhazard) df(5) eform //nolog		
						
			stsplit, at(failures)
			gen t_enter=_t0
			gen t_out=_t
			replace out1 = (t_out == time_event)			

			predict p, survival
			gen ipcw = 1 / p
			
			stset t_out [pweight=ipcw], fail(out1) enter(t_enter) origin(time 0) scale(365.25) 
			
			***** Final Analysis
			***** PH assumption
			cap drop time
			cap stpm3 i.expo_`v' agesp* i.gender i.indexyr ib2.region_g i.bc i.education3 i.fam_panc i.t2donsetyr_g ib2.bmi_g i.smokinghabit_g i.hba1c_g i.p_lipidlowerdrug i.p_alcoholdisorder i.p_gallstones i.p_abdomsurgery if expo_`v'!=., knots(5 35 65 95, percentile) scale(lncumhazard) eform 	
			if _rc!=0 {
				di as error "Full PH model for `o' `c' failed, running model without abdomsurgery"
				cap stpm3 i.expo_`v' agesp* i.gender i.indexyr ib2.region_g i.bc i.education3 i.fam_panc i.t2donsetyr_g ib2.bmi_g i.smokinghabit_g i.hba1c_g i.p_lipidlowerdrug i.p_alcoholdisorder i.p_gallstones if expo_`v'!=., knots(5 35 65 95, percentile) scale(lncumhazard) eform 	
			}
			if _rc!=0 {
				di as error "Model without abdomsurgery failed. Skipping this iteration"
				continue
			}
			putexcel set "${resultpath}\TableS3_PP_IPCW_ph_new.xlsx", sheet(`o'_`v') modify
			putexcel A`=`k'*7-5':A`=`k'*7' = "Run `k'"
			putexcel B`=`k'*7-5' = _b[1.expo_`v']
			putexcel close
			
			
			use ana_data_perprotocol, clear 
			keep `o'_date expo emig_date deathdate fuend drug_switch_date datein expo_`v' agesp* gender indexyr bc region_g education3 fam_panc t2donsetyr_g bmi_g smokinghabit_g hba1c_g preindex_lipidlowerdrug preindex_alcoholdisorder preindex_gallstones preindex_abdomsurgery
			rename preindex_* p_*
			keep if expo=="glp1" | expo==lower("`c'")
			range time 0 3 49 	

			***Standardised incidence
			standsurv, failure at1(expo_`v' 0) at2(expo_`v' 1) timevar(time) ///
			atvar(`c'_std glp_`c'_std) contrast(difference) contrastvars(riskdiff_glp_`c') frame(msurv, replace)
			
			***Derive incidences and differences Table
			foreach var in `c'_std glp_`c'_std riskdiff_glp_`c' {
				frame msurv{ 
					cap drop `var'_100
					gen `var'_100 = `var'*100
				}
			}

			frame msurv: export excel time *_100 using "${resultpath}\TableS3_PP_IPCW_ph_new.xlsx" if inlist(time, 0.5, 1, 1.5, 2, 2.5, 3), sheet(`o'_`v', modify) cell(C`=`k'*7-5') 
		}
	}
}
	
	

foreach o in pancreatitis_all pancreatitis_hosp cout1_all cout2_all {    //UPDATE NEED
	forvalues i =1/2{		//UPDATE NEED
		local v=cond(`i'==1, "gvd", "gvs")
		local c=cond(`i'==1, "DPP4", "SGLT2")
		
		putexcel set "${resultpath}\TableS3_PP_IPCW_tvc.xlsx", sheet(`o'_`v', replace) modify
		putexcel B1 = "HR"
		putexcel C1 = "Time"
		putexcel D1 = "`c'_std_100"
		putexcel E1 = "glp_`c'_std_100" 
		putexcel F1 = "riskdiff_glp_`c'_100"
		putexcel close 
		
		forvalues k=1/201{   //UPDATE NEED
			use ana_data_perprotocol, clear 
			keep `o'_date expo emig_date deathdate fuend drug_switch_date datein expo_`v' agesp* gender indexyr bc region_g education3 fam_panc t2donsetyr_g bmi_g smokinghabit_g hba1c_g preindex_lipidlowerdrug preindex_alcoholdisorder preindex_gallstones preindex_abdomsurgery
			rename preindex_* p_*
			keep if expo=="glp1" | expo==lower("`c'")		

			if `k' > 1 {
				set seed `=`k'*10'
				bsample
			}
			
			replace drug_switch_date = drug_switch_date-1 if drug_switch_date==`o'_date
			gen drug_switch_date_psuedo = datein + 30*(ceil((drug_switch_date - datein)/30))
			format drug_switch_date_psuedo %d
			
			egen dateout=rowmin(`o'_date emig_date deathdate fuend drug_switch_date_psuedo)
			format datein dateout fuend %d
			gen out1=(dateout==`o'_date)
			
			gen treat_vio = (drug_switch_date_psuedo==dateout)
			
			gen fup = dateout - datein 
			gen time_event = `o'_date - datein 
			replace fup = 1 if fup==0
			replace time_event = 1 if time_event ==0
			
			gen fakeid=_n
			
			stset fup, fail(treat_vio) id(fakeid)
			stpm3 i.expo_`v' agesp* i.gender i.indexyr i.bc ib2.region_g i.education3 i.fam_panc i.t2donsetyr_g ib2.bmi_g i.smokinghabit_g i.hba1c_g i.p_lipidlowerdrug i.p_alcoholdisorder i.p_gallstones i.p_abdomsurgery, scale(lncumhazard) df(5) eform //nolog		
						
			stsplit, at(failures)
			gen t_enter=_t0
			gen t_out=_t
			replace out1 = (t_out == time_event)			

			predict p, survival
			gen ipcw = 1 / p
			
			stset t_out [pweight=ipcw], fail(out1) enter(t_enter) origin(time 0) scale(365.25) 
			
			***** Final Analysis
			***** TVC option
			cap drop time
			cap stpm3 i.expo_`v' agesp* i.gender i.indexyr ib2.region_g i.bc i.education3 i.fam_panc i.t2donsetyr_g ib2.bmi_g i.smokinghabit_g i.hba1c_g i.p_lipidlowerdrug i.p_alcoholdisorder i.p_gallstones i.p_abdomsurgery if expo_`v'!=., knots(5 35 65 95, percentile) tvc(expo_`v') knotstvc(10 50 90, percentile) scale(lncumhazard) eform 	
			if _rc!=0 {
				di as error "Full PH model for `o' `c' failed, running model without abdomsurgery"
				stpm3 i.expo_`v' agesp* i.gender i.indexyr ib2.region_g i.bc i.education3 i.fam_panc i.t2donsetyr_g ib2.bmi_g i.smokinghabit_g i.hba1c_g i.p_lipidlowerdrug i.p_alcoholdisorder i.p_gallstones i.p_abdomsurgery if expo_`v'!=., knots(5 35 65 95, percentile) tvc(expo_`v') knotstvc(10 50 90, percentile) scale(lncumhazard) eform 	
			}
			putexcel set "${resultpath}\TableS3_PP_IPCW_tvc.xlsx", sheet(`o'_`v') modify
			putexcel A`=`k'*7-5':A`=`k'*7' = "Run `k'"
			putexcel B`=`k'*7-5' = _b[1.expo_`v']
			putexcel close
			
			
			use ana_data_perprotocol, clear 
			keep `o'_date expo emig_date deathdate fuend drug_switch_date datein expo_`v' agesp* gender indexyr bc region_g education3 fam_panc t2donsetyr_g bmi_g smokinghabit_g hba1c_g preindex_lipidlowerdrug preindex_alcoholdisorder preindex_gallstones preindex_abdomsurgery
			rename preindex_* p_*
			keep if expo=="glp1" | expo==lower("`c'")
			range time 0 3 49 	

			***Standardised incidence
			standsurv, failure at1(expo_`v' 0) at2(expo_`v' 1) timevar(time) ///
			atvar(`c'_std glp_`c'_std) contrast(difference) contrastvars(riskdiff_glp_`c') frame(msurv, replace)
			
			***Derive incidences and differences Table
			foreach var in `c'_std glp_`c'_std riskdiff_glp_`c' {
				frame msurv{ 
					cap drop `var'_100
					gen `var'_100 = `var'*100
				}
			}

			frame msurv: export excel time *_100 using "${resultpath}\TableS3_PP_IPCW_tvc.xlsx" if inlist(time, 0.5, 1, 1.5, 2, 2.5, 3), sheet(`o'_`v', modify) cell(C`=`k'*7-5') 
		}
	}
}
