********************************************************************************
* File: Do_master_file
* Author: Kecy Wu
* Date: 20200921
* Purpose: Creat a matched group of firms
*		   Regression
********************************************************************************

cap log close
clear all
set maxvar 5000
set more off, permanently
set type double, permanently

**Global path
global rawpath "H:\Dissertation_data\rawdata\"
global intpath "H:\Dissertation_data\intdata"
global finalpath "H:\Dissertation_data\finaldata"
global output "H:\Dissertation_data\output"

*get first year of VC funding 
use "$finalpath/merged_uniquelinks_id.dta", clear
merge 1:m vcid using "$finalpath/vxpert_expanded_us_forfailed_edit.dta"
keep if _merge == 3
keep vcid roundnum roundyear dunsno
bysort vcid: gen firstfundyear = roundyear[1]
duplicates drop vcid firstfundyear, force
keep dunsno firstfundyear
save "$finalpath\merged_links_fundyear.dta", replace

*reduce size of dataset, ignore sales as in Crane and Decker (2020)
use "H:\Dissertation_data\rawdata\nets_all_expanded.dta", clear
drop yearstarted state zipcode age1 sales
merge m:1 dunsno using "$finalpath\merged_links_fundyear.dta" /* all matched */
drop _merge
save "$rawpath\nets_all_expanded_addyear.dta", replace

*psmatch regression
gen lnemp = log(emp)
teffects psmatch (lnemp) (vcfunded firstyear age statenum parentstatus sic1)

*regression: control year 0, test year 3, year 5, year 10 (LT)
*cluster standard errors by industry
use "$rawpath\nets_all_expanded_addyear.dta", clear
sort dunsno year
drop sales_gr firstfundyear
bysort dunsno: gen emp3 = emp if age == 3
bysort dunsno: gen emp5 = emp if age == 5
bysort dunsno: gen emp10 = emp if age == 10
bysort dunsno: gen emp0 = emp[1]
drop emp firstyear age emp_gr
save "$rawpath\nets_empreg.dta", replace

keep if !missing(emp3)
drop emp5 emp10
gen lnemp0 = log(emp0)
gen lnemp3 = log(emp3)

eststo model1: reghdfe lnemp3 vcfunded lnemp0, absorb(sic1 parentstatus statenum year) vce(cluster sic1)
*esttab model1 using "$output\temp.tex", ///
*		      r2 se label nonumbers mtitle("model 1")

use "$rawpath\nets_empreg.dta", clear
keep if !missing(emp5)
drop emp3 emp10
gen lnemp0 = log(emp0)
gen lnemp5 = log(emp5)

eststo model2: reghdfe lnemp5 vcfunded lnemp0, absorb(sic1 parentstatus statenum year) vce(cluster sic1)

use "$rawpath\nets_empreg.dta", clear
keep if !missing(emp10)
drop emp3 emp5
gen lnemp0 = log(emp0)
gen lnemp10 = log(emp10)

eststo model3: reghdfe lnemp10 vcfunded lnemp0, absorb(sic1 parentstatus statenum year) vce(cluster sic1)

esttab model* using "$output\empreg.tex", r2 se label nonumbers
