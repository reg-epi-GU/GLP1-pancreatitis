*** PPgraph

cd "S:\HI-SPEED projects\HS26_02 GLP1 pancreatitis\output\20261015\figure"
global path "S:\HI-SPEED projects\HS26_02 GLP1 pancreatitis\output\20261015\PP_IPCW\with excel cal"


forvalues i=1/2{
	quietly: mylabels 0(0.2)1.6, local(labels) myscale(@) suffix("%") format(%4.1f)

	local v=cond(`i'==1, "gvs", "gvd")
	local c=cond(`i'==1, "sglt2", "dpp4")
	local nam=cond(`i'==1, "SGLT-2", "DPP-4")
	
	import excel using "${path}\PP_IPCW_tvc_maxweight10_panc_all_wk852_`v'", cellrange(M1:V11) firstrow clear
	rename *, lower	

	tw (line glp_`c'_std_100 time, lc(blue%80)) ///
		(line `c'_std_100 time, lc(red%80)) ///
		(rline glp_`c'_std_100_ci_low glp_`c'_std_100_ci_high time, color(blue%30) lpattern(dash)) ///
		(rline `c'_std_100_ci_low `c'_std_100_ci_high time, color(red%30) lpattern(dash)), ///
		scheme(tab2) ///
		title("{bf:Per protocol}") ///
		ylabel(`labels', labsize(small)) ///
		xtitle("{bf:Years since treatment initiation}", size(*1) margin(0 0 0 0)) ///
		xlabel(0(1)5, nogrid format(%9.1f) labsize(small)) ///
		xscale(range(0 5)) ///
		ytitle("{bf:Weighted risk}" "{bf:% (95% CI)}", size(*1) margin(0 0 0 0)) ///
		legend(order(1 "GLP-1" 2 "`nam'") col(1) size(*1.5) title("") ring(0) pos(11) yoffset(-5) xoffset(2) region(fcolor(none))) 
	graph save "fig`=`i'+1'B_curv", replace
	
	tw (rcap riskdiff_glp_`c'_100_ci_low riskdiff_glp_`c'_100_ci_high time if inlist(time, 0, 0.5, 1, 1.5, 2, 2.5, 3, 3.5, 4, 4.5, 5), color(black%30)) ///
		(scatter riskdiff_glp_`c'_100 time if inlist(time, 0, 0.5, 1, 1.5, 2, 2.5, 3, 3.5, 4, 4.5, 5), color(black)), scheme(tab2) ///
		ylabel(-0.4(0.1)0.4, nogrid format(%04.2f) labsize(small)) ///
		xtitle("{bf:Years since treatment initiation}", size(*1) margin(0 0 0 0)) ///
		xlabel(0(1)5, nogrid format(%9.1f) labsize(small)) ///
		xscale(range(0 5)) ///
		ytitle("{bf:Risk difference}" "{bf:percentage point (95% CI)}", size(*1) margin(0 0 0 0)) ///
		yscale(range(-0.2 0.2)) ///
		yline(0, lp(-)) ///
		legend(off)
	graph save "fig`=`i'+1'B_rd", replace
		

	import excel using "${path}\PP_IPCW_tvc_maxweight10_failratio_panc_all_wk852_`v'", cellrange(Z1:AI11) firstrow clear
	rename *, lower

	tw (rcap failratio_glp_`c'_ci_low failratio_glp_`c'_ci_high time if inlist(time, 0, 0.5, 1, 1.5, 2, 2.5, 3, 3.5, 4, 4.5, 5), color(black%30)) ///
		(scatter failratio_glp_`c' time if inlist(time, 0, 0.5, 1, 1.5, 2, 2.5, 3, 3.5, 4, 4.5, 5), color(black)), scheme(tab2) ///
		ylabel(0.5(0.5)2, nogrid format(%4.1f) labsize(small)) ///
		xtitle("{bf:Years since treatment initiation}", size(*1) margin(0 0 0 0)) ///
		xlabel(0(1)5, nogrid format(%9.1f) labsize(small)) ///
		xscale(range(0 5)) ///
		ytitle("{bf:Risk ratio} {bf:(95% CI)}", size(*1) margin(0 0 0 0)) ///
		yscale(range(0.5 2)  log) ///
		yline(1, lp(-)) ///
		legend(off)
	graph save "fig`=`i'+1'B_rr", replace
}