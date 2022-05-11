********************************************************************************
* File: Do_master_file
* Author: Kecy Wu
* Date: 20211215
* Purpose: WMS and VC: management quality (long term impact)
********************************************************************************

cap log close
clear all
set maxvar 5000
set matsize 10000
set more off, permanently
set type double, permanently

**Global path
global rawpath "C:\Users\kecyw\Dropbox\Graduate\Research\Venture_Capital_Paper\Empirical\rawdata"
global intpath "C:\Users\kecyw\Dropbox\Graduate\Research\Venture_Capital_Paper\Empirical\intdata"
global finalpath "C:\Users\kecyw\Dropbox\Graduate\Research\Venture_Capital_Paper\Empirical\finaldata"
global output "C:\Users\kecyw\Dropbox\Graduate\Research\Venture_Capital_Paper\Empirical\output"

log using "$output\wms_vc_121521", replace
use "$finalpath/merged_vc_wms_cleaned.dta", clear
tab vcfunded
summarize management operations monitor target people firmage org_a

estpost summarize management operations monitor target people firmage org_a,listwise
esttab using "$output\wmstats.tex", cells("count mean sd min max") nomtitle nonumber

foreach var of varlist management operations monitor target people firmage org_a post_g i_comptenure i_posttenure {
	tab vcfunded, summarize(`var')
}

*operations: lean production process, automation, flexible manpower
*monitor: track and review KPI, enact and follow up changes 
*target: set short and long term goals, follow up, debrief
*people: attract and keep talent, set bonus and promotion

*regression - management and people

foreach var of varlist management operations monitor target people {
	reg `var' vcfunded, robust
}

foreach var of varlist management operations monitor target people {
	reg `var' vcfunded firmage org_a, robust
}

* don't do the following as no. of observations is few
*foreach var of varlist management operations monitor target people {
*	reghdfe `var' vcfunded firmage org_a, absorb(sic_wms) vce(robust)
*}

*foreach var of varlist management operations monitor target people {
*	reg `var' vcfunded firmage org_a post_a reliability, robust
*}

*compare matched vc firms to rest of vc sample
use "$finalpath/SDC_Vxpert_ALL_vcid_summary.dta",clear
summarize totalamt totalround totalinvestor highqualityinvestor earlystage
merge 1:m vcid using "$finalpath/merged_vc_wms_cleaned.dta", gen(merge)
drop if merge == 2
duplicates drop vcid, force
gen wms = 1 if merge == 3
replace wms = 0 if wms == . /* 80 unique vc firms */
foreach var of varlist totalamt totalround totalinvestor highqualityinvestor earlystage {
	tab wms, summarize(`var')
}

*regression - high and low quality VC firms
use "$finalpath/SDC_Vxpert_ALL_vcid_summary.dta",clear
merge 1:m vcid using "$finalpath/merged_vc_wms_cleaned.dta", gen(merge)
keep if merge == 3

foreach var of varlist management operations monitor target people {
	tab highqualityinvestor, summarize(`var')
}

foreach var of varlist management operations monitor target people {
	reg `var' highqualityinvestor, robust
}

foreach var of varlist management operations monitor target people {
	reg `var' highqualityinvestor firmage org_a, robust
}

log close