*** 03 Survival analysis

cd "S:\HI-SPEED projects\HS26_02 GLP1 pancreatitis\data\"
global resultpath "S:\HI-SPEED projects\HS26_02 GLP1 pancreatitis\output\20260422"
cap mkdir "${resultpath}"

* Final analysis with IPMW
		************* Proportional Assumption (no tvc)

foreach o in pancreatitis_all pancreatitis_hosp cout1_all cout2_all{			
	forvalues i=1/2{		
		local v=cond(`i'==1, "gvd", "gvs")
		local c=cond(`i'==1, "DPP4", "SGLT2")
		
		use ana_data, clear
		keep `o'_date expo emig_date deathdate fuend datein expo_`v' agesp* gender indexyr bc region_g education3 fam_panc t2donsetyr_g bmi_g smokinghabit_g hba1c_g preindex_lipidlowerdrug preindex_alcoholdisorder preindex_gallstones preindex_abdomsurgery
		keep if expo=="glp1" | expo==lower("`c'")		
		
		foreach var in gender indexyr region_g bc education3 fam_panc t2donsetyr_g bmi_g smokinghabit_g hba1c_g preindex_lipidlowerdrug preindex_alcoholdisorder preindex_gallstones{
			di "`var'"
			replace `var'=. if `var'==99
		}
		egen complete = rowmiss(education3 bmi_g smokinghabit_g hba1c_g)
		recode complete (0=1) (else=0)
		
		logit complete i.expo_`v' agesp* i.gender i.indexyr ib2.region_g i.bc i.fam_panc i.t2donsetyr_g i.preindex_lipidlowerdrug i.preindex_alcoholdisorder i.preindex_gallstones
		predict pr_nomiss
		gen ipmw=1/pr_nomiss
		
		keep if complete==1
		
		egen dateout_`o'=rowmin(`o'_date emig_date deathdate fuend)
		format datein dateout_`o' fuend %d
		gen `o'_out1=(dateout_`o'==`o'_date)

		stset dateout_`o' [pweight=ipmw], fail(`o'_out1) enter(datein) origin(datein) scale(365.25)
		range time 0 3 49 

		collect clear
		cap collect: stpm3 i.expo_`v' agesp* i.gender i.indexyr ib2.region_g i.bc i.education3 i.fam_panc i.t2donsetyr_g ib2.bmi_g i.smokinghabit_g i.hba1c_g i.preindex_lipidlowerdrug i.preindex_alcoholdisorder i.preindex_gallstones i.preindex_abdomsurgery, knots(5 35 65 95, percentile) scale(lncumhazard) eform 
		
		if _rc != 0 {
			di as error "Full PH model for `o' `c' failed, running model without abdomsurgery"
			collect clear
			collect:stpm3 i.expo_`v' agesp* i.gender i.indexyr ib2.region_g i.bc i.education3 i.fam_panc i.t2donsetyr_g ib2.bmi_g i.smokinghabit_g i.hba1c_g i.preindex_lipidlowerdrug i.preindex_alcoholdisorder i.preindex_gallstones, knots(5 35 65 95, percentile) scale(lncumhazard) eform 	
		} 
		
		collect layout (colname) (result[_r_b _r_p _r_lb _r_ub])
		collect export "${resultpath}\TableS3_ITT_ph_ipmw.xlsx", sheet(`o'_`v', replace) cell(A9) modify
		
		***Standardised incidence
		standsurv, failure at1(expo_`v' 0) at2(expo_`v' 1) timevar(time) ///
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

		graph save "${resultpath}\Fig `o'_`v'_ITT_ph_ipmw.gph", replace
				
		***Derive incidences and differences Table
		foreach var in `c'_std `c'_std_lci `c'_std_uci glp_`c'_std glp_`c'_std_lci glp_`c'_std_uci riskdiff_glp_`c' riskdiff_glp_`c'_lci riskdiff_glp_`c'_uci {
			frame msurv{ 
				cap drop `var'_100
				gen `var'_100 = `var'*100
			}
		}

		frame msurv: export excel time *_100 using "${resultpath}\TableS3_ITT_ph_ipmw.xlsx" if inlist(time, 0.5, 1, 1.5, 2, 2.5, 3), sheet(`o'_`v', modify)  firstrow(variables)

		
		************* Non-proportional (tvc)		
		collect clear
		cap collect: stpm3 i.expo_`v' agesp* i.gender i.indexyr ib2.region_g i.bc i.education3 i.fam_panc i.t2donsetyr_g ib2.bmi_g i.smokinghabit_g i.hba1c_g i.preindex_lipidlowerdrug i.preindex_alcoholdisorder i.preindex_gallstones i.preindex_abdomsurgery, knots(5 35 65 95, percentile) tvc(expo_`v') knotstvc(10 50 90, percentile) scale(lncumhazard) eform 
		
		if _rc != 0 {
			di as error "Full TVC model for `o' `c' failed, running model without abdomsurgery"
			collect clear
			collect:stpm3 i.expo_`v' agesp* i.gender i.indexyr ib2.region_g i.bc i.education3 i.fam_panc i.t2donsetyr_g ib2.bmi_g i.smokinghabit_g i.hba1c_g i.preindex_lipidlowerdrug i.preindex_alcoholdisorder i.preindex_gallstones, knots(5 35 65 95, percentile) tvc(expo_`v') knotstvc(10 50 90, percentile) scale(lncumhazard) eform 	
		} 
		
		collect layout (colname) (result[_r_b _r_p _r_lb _r_ub])
		collect export "${resultpath}\TableS3_ITT_tvc_ipmw.xlsx", sheet(`o'_`v', replace) cell(A9) modify
		
		***Standardised incidence
		standsurv, failure at1(expo_`v' 0) at2(expo_`v' 1) timevar(time) ///
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

		graph save "${resultpath}\Fig `o'_`v'_ITT_tvc_ipmw.gph", replace
				
		***Derive incidences and differences Table
		foreach var in `c'_std `c'_std_lci `c'_std_uci glp_`c'_std glp_`c'_std_lci glp_`c'_std_uci riskdiff_glp_`c' riskdiff_glp_`c'_lci riskdiff_glp_`c'_uci {
			frame msurv{ 
				cap drop `var'_100
				gen `var'_100 = `var'*100
			}
		}

		frame msurv: export excel time *_100 using "${resultpath}\TableS3_ITT_tvc_ipmw.xlsx" if inlist(time, 0.5, 1, 1.5, 2, 2.5, 3), sheet(`o'_`v', modify)  firstrow(variables)
	}		
}	
