********************************************************************************
* File: Do_master_file
* Author: Kecy Wu
* Date: 20200919
* Purpose: Plot graph of vc vs non-vc firms, vc in overall number and employment
*	   unconditional exit rates comparison
*	   exit times for vc-funded firms
********************************************************************************
cap log close
clear all
set maxvar 5000
set more off, permanently
set type double, permanently

**Global path
global rawpath "/home/kecywu/Dropbox/Graduate/Research/Venture_Capital_Paper/Empirical/rawdata"
global intpath "/home/kecywu/Dropbox/Graduate/Research/Venture_Capital_Paper/Empirical/intdata"
global finalpath "/home/kecywu/Dropbox/Graduate/Research/Venture_Capital_Paper/Empirical/finaldata"
global output "/home/kecywu/Dropbox/Graduate/Research/Venture_Capital_Paper/Empirical/output"

use "$finalpath/nets_sumstats.dta", clear
sort vcfunded age
keep if age <= 10

gen empu = emp + empse*1.96
gen empl = emp - empse*1.96
gen empgru = emp_gr + empgrse*1.96
gen empgrl = emp_gr - empgrse*1.96

*employment and growth
*check year 5 and 6, large std for vc funded firms
*check employment growth for non-vc fudned firms
set scheme s2color
graph twoway (line emp age if vcfunded == 1) (rcap empu empl age if vcfunded == 1) ///
	(line emp age if vcfunded == 0) (rcap empu empl age if vcfunded == 0), ///
	graphregion(fcolor(white)) legend(label(1 VC) label(2 VC C.I.) label(3 non-VC) label(4 non-VC C.I.)) ///
	ylabel(0(20)100) xlabel(0(1)10) ///
	xline(5, lp(dash) lc(black)) ///
	xtitle("age") ytitle("employment") title("Employment by Age")

*xaxis(1 2) xla(1940 "start", axis(2) grid glcolor(red)) xtitle("", axis(2))
	
graph twoway (line emp_gr age if vcfunded == 1) (rcap empgru empgrl age if vcfunded == 1) ///
	(line emp_gr age if vcfunded == 0) (rcap empgru empgrl age if vcfunded == 0), ///
	graphregion(fcolor(white)) legend(label(1 VC) label(2 VC C.I.) label(3 non-VC) label(4 non-VC C.I.)) ///
	ylabel(0(0.1)0.4) xlabel(0(1)10) ///
	xline(5, lp(dash) lc(black)) ///
	xtitle("age") ytitle("employment growth") title("Employment Growth by Age")

graph twoway (line emp age if vcfunded == 1) (line emp age if vcfunded == 0), ///
	graphregion(fcolor(white)) legend(label(1 VC) label(2 non-VC)) ///
	ylabel(0(20)100) xlabel(0(1)10) ///
	xtitle("age") ytitle("employment") title("Employment by Age")	

graph twoway (line emp_gr age if vcfunded == 1) (line emp_gr age if vcfunded == 0), ///
	graphregion(fcolor(white)) legend(label(1 VC) label(2 non-VC)) ///
	ylabel(0(0.1)0.4) xlabel(0(1)10) ///
	xtitle("age") ytitle("employment") title("Employment Growth by Age")	

*presentation uses this
use "$finalpath/nets_sumstats_2.dta", clear
sort vcfunded age
keep if age <= 10 & age >= 0

set scheme s2color
graph twoway (line emp age if vcfunded == 1) (line emp age if vcfunded == 0), ///
	graphregion(fcolor(white)) legend(label(1 VC) label(2 non-VC)) ///
	ylabel(0(20)100) xlabel(0(1)10) ///
	xtitle("age") ytitle("employment") title("Employment by Age")
	
graph twoway (line emp_gr age if vcfunded == 1 & age > 0) (line emp_gr age if vcfunded == 0 & age > 0), ///
	graphregion(fcolor(white)) legend(label(1 VC) label(2 non-VC)) ///
	ylabel(0(0.1)0.5) xlabel(0(1)10) ///
	xtitle("age") ytitle("employment growth") title("Employment Growth by Age")

	
*sales and growth
*plot later, use thousands, and at 2000$
graph twoway (line sales age if vcfunded == 1) (line sales age if vcfunded == 0), ///
	graphregion(fcolor(white)) legend(label(1 VC) label(2 non-VC)) ///
	ylabel(0(20)100) xlabel(0(1)10) ///
	xtitle("age") ytitle("sales") title("Average Sales")
	
graph twoway (line sales_gr age if vcfunded == 1) (line sales_gr age if vcfunded == 0), ///
	graphregion(fcolor(white)) legend(label(1 VC) label(2 non-VC)) ///
	ylabel(0(0.1)0.4) xlabel(0(1)10) ///
	xtitle("age") ytitle("sales growth") title("Average Sales Growth")

*exit rates, unconditional
gen N1 = 68481845
gen N2 = 33991
gen failedpct = failed_age/N1*100 if vcfunded == 0
replace failedpct = failed_age/N2*100 if vcfunded == 1
keep age vcfunded failedpct
save "$finalpath/exitbyage_nonvc.dta", replace /* this dataset only non-vc firms */

append using "$finalpath/exitbyage_vc.dta"

*didn't use age_2 to calculate failed pct
graph twoway (line failedpct age if vcfunded == 1) (line failedpct age if vcfunded == 0), ///
	graphregion(fcolor(white)) legend(label(1 VC) label(2 non-VC)) ///
	ylabel(0(5)15) xlabel(0(1)10) ///
	xtitle("age") ytitle("percentage") title("Exit Rates by Age")
		*xline(5, lp(dash) lc(black)) ///

use "$finalpath/vxpert_expanded_us_forfailed_edit2.dta", clear
gen age = roundyear - foundyear
gen failed_age = failed*lastround
gen success_age = success*lastround
collapse (sum) failed_age success_age, by(age)
keep if age>=0 & age <=10
save "$finalpath/exitbyage.dta", replace

use "$finalpath/exitbyage.dta", clear
gen N = 33991
gen failedpct = failed_age/N*100
gen successpct = success_age/N*100

set scheme s2color
graph twoway (line failedpct age) (line successpct age), ///
	graphregion(fcolor(white)) legend(label(1 failed) label(2 success)) ///
	ylabel(0(1)5) xlabel(0(1)10) ///
	xtitle("age") ytitle("percentage") title("Percentage of Failures/Successes by Age")

gen vcfunded = 1
keep vcfunded age failedpct
save "$finalpath/exitbyage_vc.dta", replace

*overall pct of vc-funded firms
use "$finalpath/nets_vcpct.dta", clear
twoway line vc_num_pct vc_emp_pct year, graphregion(fcolor(white)) ///
	ylabel(0(0.05)0.3) xlabel(1985(5)2015) legend(label(1 number) label(2 employment)) ///
	xtitle("Year") ytitle("Percentage") title("VC Funded Firms in the Economy") 


*average exit year
use "$finalpath/vxpert_expanded_us_forfailed_edit2.dta", clear
sort vcid roundyear
keep vcid roundyear rounddate
bysort vcid: gen length = (rounddate[_N] - rounddate[1])/30
keep if length > 0
sum length, detail
*10% 13.2, 50% 51.7, mean 61.8, 90% 124


*Length of relationship
use "$finalpath/vxpert_expanded_us_investor_stage.dta", clear
sort vcid rounddate
keep vcid rounddate lastround
duplicates drop vcid rounddate, force
bysort vcid: gen nextrounddate = rounddate[_n+1]
format nextrounddate %td
replace nextrounddate = rounddate if lastround == 1
merge 1:m vcid rounddate lastround using "$finalpath/vxpert_expanded_us_investor_stage.dta"
keep if _merge == 3
drop _merge

sort investorid vcid rounddate
order investorid vcid rounddate
gen roundlength = (nextrounddate - rounddate)/30
*26% of the data have 0 roundlength, don't know lastround to exit times
drop if roundlength == 0

bysort investorid vcid: egen totallength = total(roundlength)
bysort investorid: egen avelength = mean(totallength)
duplicates drop investorid, force
keep investorid avelength
sum avelength, detail
hist avelength, graphregion(fcolor(white)) ///
	xtitle("Months") ytitle("Density") title("Distribution of Average Length of Relationship") 
*mean 28 months, medium 24 months
