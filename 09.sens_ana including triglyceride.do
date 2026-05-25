*** 09 sensitivity analysis of including triglyceride in the models

cd "S:\HI-SPEED projects\HS26_02 GLP1 pancreatitis\data\"
global resultpath "S:\HI-SPEED projects\HS26_02 GLP1 pancreatitis\output\20260422"
cap mkdir "${resultpath}"


**********************************************************************ITT
	* Final analysis
		************* Proportional Assumption (no tvc)
use ana_data, clear
range time 0 3 49 
foreach o in pancreatitis_all{	
		egen dateout_`o'=rowmin(`o'_date emig_date deathdate fuend)
		format datein dateout_`o' fuend %d
		gen `o'_out1=(dateout_`o'==`o'_date)

		stset dateout_`o', fail(`o'_out1) enter(datein) origin(datein) scale(365.25)
		
	forvalues i=1/2{		
		local v=cond(`i'==1, "gvd", "gvs")
		local c=cond(`i'==1, "DPP4", "SGLT2")
		
		collect clear
		collect:stpm3 i.expo_`v' agesp* i.gender i.indexyr ib2.region_g i.bc i.education3 i.fam_panc i.t2donsetyr_g ib2.bmi_g i.smokinghabit_g i.hba1c_g i.triglyceride_g i.preindex_lipidlowerdrug i.preindex_alcoholdisorder i.preindex_gallstones i.preindex_abdomsurgery if expo_`v'!=., knots(5 35 65 95, percentile) scale(lncumhazard) eform 	

		collect layout (colname) (result[_r_b _r_p _r_lb _r_ub])
		collect export "${resultpath}\TableS4_sens_triglyceride.xlsx", sheet(ITTph_`v', replace) cell(A9) modify
		
		***Standardised incidence
		standsurv if expo_`v'!=., failure at1(expo_`v' 0) at2(expo_`v' 1) timevar(time) ///
		atvar(`c'_std glp_`c'_std) contrast(difference) contrastvars(riskdiff_glp_`c') ci frame(msurv, replace)

		***Derive incidences and differences Table
		foreach var in `c'_std `c'_std_lci `c'_std_uci glp_`c'_std glp_`c'_std_lci glp_`c'_std_uci riskdiff_glp_`c' riskdiff_glp_`c'_lci riskdiff_glp_`c'_uci {
			frame msurv{ 
				cap drop `var'_100
				gen `var'_100 = `var'*100
			}
		}

		frame msurv: export excel time *_100 using "${resultpath}\TableS4_sens_triglyceride.xlsx" if inlist(time, 0.5, 1, 1.5, 2, 2.5, 3), sheet(ITTph_`v', modify)  firstrow(variables)
		
	}		
}	

use ana_data_perprotocol, clear
range time 0 3 49 
foreach o in pancreatitis_all{	
		egen dateout_`o'=rowmin(`o'_date emig_date deathdate fuend drug_switch_date)
		format datein dateout_`o' fuend %d
		gen `o'_out1=(dateout_`o'==`o'_date)

		stset dateout_`o', fail(`o'_out1) enter(datein) origin(datein) scale(365.25)
		
	forvalues i=1/2{		
		local v=cond(`i'==1, "gvd", "gvs")
		local c=cond(`i'==1, "DPP4", "SGLT2")
		
		collect clear
		collect:stpm3 i.expo_`v' agesp* i.gender i.indexyr ib2.region_g i.bc i.education3 i.fam_panc i.t2donsetyr_g ib2.bmi_g i.smokinghabit_g i.hba1c_g i.triglyceride_g i.preindex_lipidlowerdrug i.preindex_alcoholdisorder i.preindex_gallstones i.preindex_abdomsurgery if expo_`v'!=., knots(5 35 65 95, percentile) scale(lncumhazard) eform 	

		collect layout (colname) (result[_r_b _r_p _r_lb _r_ub])
		collect export "${resultpath}\TableS4_sens_triglyceride.xlsx", sheet(PPphnaive_`v', replace) cell(A9) modify
		
		***Standardised incidence
		standsurv if expo_`v'!=., failure at1(expo_`v' 0) at2(expo_`v' 1) timevar(time) ///
		atvar(`c'_std glp_`c'_std) contrast(difference) contrastvars(riskdiff_glp_`c') ci frame(msurv, replace)

		***Derive incidences and differences Table
		foreach var in `c'_std `c'_std_lci `c'_std_uci glp_`c'_std glp_`c'_std_lci glp_`c'_std_uci riskdiff_glp_`c' riskdiff_glp_`c'_lci riskdiff_glp_`c'_uci {
			frame msurv{ 
				cap drop `var'_100
				gen `var'_100 = `var'*100
			}
		}

		frame msurv: export excel time *_100 using "${resultpath}\TableS4_sens_triglyceride.xlsx" if inlist(time, 0.5, 1, 1.5, 2, 2.5, 3), sheet(PPphnaive_`v', modify)  firstrow(variables)
		
	}		
}	