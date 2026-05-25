cd "S:\HI-SPEED projects\HS26_02 GLP1 pancreatitis\output\fig for pub"


use "S:\HI-SPEED projects\HS26_02 GLP1 pancreatitis\data\ana_data_perprotocol.dta", clear
range time 0 3 49
quietly: mylabels 0(0.1)1, local(labels) myscale(@/100) suffix("%") format(%4.1f)
***** Descriptives
foreach o in pancreatitis_all{
	
	egen dateout_`o'=rowmin(`o'_date emig_date deathdate fuend)
	format datein dateout_`o' fuend %d
	gen `o'_out1=(dateout_`o'==`o'_date)

	egen dateout_`o'_pp=rowmin(`o'_date emig_date deathdate fuend drug_switch_date)
	gen `o'_out1_pp=(dateout_`o'_pp==`o'_date)
	
	forvalues i=1/2{	
		
	******* ITT graph
		stset dateout_`o', fail(`o'_out1) enter(datein) origin(datein) scale(365.25)
	
		local v=cond(`i'==1, "gvs", "gvd")
		local c=cond(`i'==1, "SGLT2", "DPP4")
		local nam=cond(`i'==1, "SGLT-2", "DPP-4")
		

		cap stpm3 i.expo_`v' agesp* i.gender i.indexyr ib2.region_g i.bc i.education3 i.fam_panc i.t2donsetyr_g ib2.bmi_g i.smokinghabit_g i.hba1c_g i.preindex_lipidlowerdrug i.preindex_alcoholdisorder i.preindex_gallstones i.preindex_abdomsurgery if expo_`v'!=., knots(5 35 65 95, percentile) scale(lncumhazard) eform 	
		if _rc != 0 {
			stpm3 i.expo_`v' agesp* i.gender i.indexyr ib2.region_g i.bc i.education3 i.fam_panc i.t2donsetyr_g ib2.bmi_g i.smokinghabit_g i.hba1c_g i.preindex_lipidlowerdrug i.preindex_alcoholdisorder i.preindex_gallstones if expo_`v'!=., knots(5 35 65 95, percentile) scale(lncumhazard) eform 	
		}
		if _rc != 0{
			di as error "Both Model failed, check how to continue"
			break
		}
		cap frame drop msurv
		standsurv if expo_`v'!=., failure at1(expo_`v' 0) at2(expo_`v' 1) timevar(time) ///
		atvar(`c'_std glp_`c'_std) contrast(difference) contrastvars(riskdiff_glp_`c') ci frame(msurv, replace)
		
		foreach var in `c'_std `c'_std_lci `c'_std_uci glp_`c'_std glp_`c'_std_lci glp_`c'_std_uci riskdiff_glp_`c' riskdiff_glp_`c'_lci riskdiff_glp_`c'_uci {
			frame msurv{ 
				cap drop `var'_100
				gen `var'_100 = `var'*100
			}
		}

		***Graph across time 
		frame msurv: tw (line glp_`c'_std time, lc(blue%80)) ///
		(line `c'_std time, lc(red%80)) ///
		(rline glp_`c'_std_lci glp_`c'_std_uci time, color(blue%30) lpattern(dash)) ///
		(rline `c'_std_lci `c'_std_uci time, color(red%30) lpattern(dash)), ///
		scheme(tab2) ///
		title("{bf:Intention-to-treat}") ///
		ylabel(`labels', labsize(small)) ///
		xtitle("{bf:Years since treatment initiation}", size(*1) margin(0 0 0 0)) ///
		xlabel(0(0.5)3, nogrid format(%9.1f) labsize(small)) ///
		xscale(range(0 3)) ///
		ytitle("{bf:Standardised cumulative incidence}" "{bf:% (95% CI)}", size(*1) margin(0 0 0 0)) ///
		name("fig`=`i'+1'A_curv", replace) ///
		legend(order(1 "GLP-1" 2 "`nam'") col(1) size(*1.5) title("") ring(0) pos(11) yoffset(-11) xoffset(2) region(fcolor(none)))

		graph save "fig`=`i'+1'A_curv", replace
		
		frame msurv: tw (rcap riskdiff_glp_`c'_lci_100 riskdiff_glp_`c'_uci_100 time if inlist(time, 0, 0.5, 1, 1.5, 2, 2.5, 3), color(blue%30)) ///
		(scatter riskdiff_glp_`c'_100 time if inlist(time, 0, 0.5, 1, 1.5, 2, 2.5, 3)), scheme(tab2) ///
		ylabel(-0.2(0.05)0.2, nogrid format(%04.2f) labsize(small)) ///
		xtitle("{bf:Years since treatment initiation}", size(*1) margin(0 0 0 0)) ///
		xlabel(0(0.5)3, nogrid format(%9.1f) labsize(small)) ///
		xscale(range(0 3)) ///
		ytitle("{bf:Risk difference}" "{bf:percentage point (95% CI)}", size(*1) margin(0 0 0 0)) ///
		yscale(range(-0.2 0.2)) ///
		yline(0, lp(-)) ///
		legend(off) ///
		name("fig`=`i'+1'A_rd", replace)
		
		graph save "fig`=`i'+1'A_rd", replace
		
		
	******* PP graph
		stset dateout_`o'_pp, fail(`o'_out1_pp) enter(datein) origin(datein) scale(365.25)
	
		cap stpm3 i.expo_`v' agesp* i.gender i.indexyr ib2.region_g i.bc i.education3 i.fam_panc i.t2donsetyr_g ib2.bmi_g i.smokinghabit_g i.hba1c_g i.preindex_lipidlowerdrug i.preindex_alcoholdisorder i.preindex_gallstones i.preindex_abdomsurgery if expo_`v'!=., knots(5 35 65 95, percentile) scale(lncumhazard) eform 	
		if _rc != 0 {
			stpm3 i.expo_`v' agesp* i.gender i.indexyr ib2.region_g i.bc i.education3 i.fam_panc i.t2donsetyr_g ib2.bmi_g i.smokinghabit_g i.hba1c_g i.preindex_lipidlowerdrug i.preindex_alcoholdisorder i.preindex_gallstones if expo_`v'!=., knots(5 35 65 95, percentile) scale(lncumhazard) eform 	
		}
		if _rc != 0{
			di as error "Both Model failed, check how to continue"
			break
		}
		cap frame drop msurv
		standsurv if expo_`v'!=., failure at1(expo_`v' 0) at2(expo_`v' 1) timevar(time) ///
		atvar(`c'_std glp_`c'_std) contrast(difference) contrastvars(riskdiff_glp_`c') ci frame(msurv, replace)
		
		foreach var in `c'_std `c'_std_lci `c'_std_uci glp_`c'_std glp_`c'_std_lci glp_`c'_std_uci riskdiff_glp_`c' riskdiff_glp_`c'_lci riskdiff_glp_`c'_uci {
			frame msurv{ 
				cap drop `var'_100
				gen `var'_100 = `var'*100
			}
		}
		
		***Graph across time 
		frame msurv: tw (line glp_`c'_std time, lc(blue%80)) ///
		(line `c'_std time, lc(red%80)) ///
		(rline glp_`c'_std_lci glp_`c'_std_uci time, color(blue%30) lpattern(dash)) ///
		(rline `c'_std_lci `c'_std_uci time, color(red%30) lpattern(dash)), ///
		scheme(tab2) ///
		title("{bf:Per-protocol}") ///
		ylabel(`labels', labsize(small)) ///
		xtitle("{bf:Years since treatment initiation}", size(*1) margin(0 0 0 0)) ///
		xlabel(0(0.5)3, nogrid format(%9.1f) labsize(small)) ///
		xscale(range(0 3)) ///
		name("fig`=`i'+1'B_curv", replace) ///
		legend(order(1 "GLP-1" 2 "`nam'") col(1) size(*1.5) title("") ring(0) pos(11) yoffset(-11) xoffset(2) region(fcolor(none)))

		graph save "fig`=`i'+1'B_curv", replace	
		
		frame msurv: tw (rcap riskdiff_glp_`c'_lci_100 riskdiff_glp_`c'_uci_100 time if inlist(time, 0, 0.5, 1, 1.5, 2, 2.5, 3), color(blue%30)) ///
		(scatter riskdiff_glp_`c'_100 time if inlist(time, 0, 0.5, 1, 1.5, 2, 2.5, 3)), scheme(tab2) ///
		ylabel(-0.2(0.05)0.2, nogrid format(%04.2f) labsize(small)) ///
		xtitle("{bf:Years since treatment initiation}", size(*1) margin(0 0 0 0)) ///
		xlabel(0(0.5)3, nogrid format(%9.1f) labsize(small)) ///
		xscale(range(0 3)) ///
		yscale(range(-0.2 0.2)) ///
		yline(0, lp(-)) ///
		legend(off) ///
		name("fig`=`i'+1'B_rd", replace)
		
		graph save "fig`=`i'+1'B_rd", replace
	}	
}


forvalues i=1/2{	
		**** Combining graphs
		graph combine "fig`=`i'+1'A_curv" "fig`=`i'+1'B_curv" "fig`=`i'+1'A_rd" "fig`=`i'+1'B_rd", col(2) imargin(2 2 10 2) xsize(18) ysize(16) xcommon
		graph save "fig`=`i'+1'", replace
		graph export "fig`=`i'+1'.svg", replace
}

graph use "fig2A_rd"

graph combine "fig2A_curv" "fig2A_rd", col(1) xcommon

graph combine "fig2A_curv" "fig2A_rd", col(1)
