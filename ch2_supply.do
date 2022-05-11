********************************************************************************
* File: Do_master_file
* Author: Kecy Wu
* Date: 20220225
* Purpose: VC investor description
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

*number of rounds, number of investors 
use "$intpath/SDC_Vxpert_ALL_investor_vcid.dta", clear
gen year = year(date)
save "$intpath/SDC_Vxpert_ALL_investor_vcid.dta", replace

keep if year >= 1980 & year <= 2015
bysort year: egen numrounds = count(vcid)
drop if investor == "Undisclosed Fund"
duplicates drop investor year, force
bysort year: egen numinvestor = count(investor)
duplicates drop year, force

set scheme s2color
graph twoway line numrounds year, graphregion(fcolor(white)) ///
	ylabel(0(10000)40000) xlabel(1980(5)2015) ///
	xtitle("year") ytitle("number of rounds") 

graph twoway line numinvestor year, graphregion(fcolor(white)) ///
	ylabel(0(1000)7000) xlabel(1980(5)2015) ///
	xtitle("year") ytitle("number of active investors") 

use "$intpath/SDC_Vxpert_ALL_investor_vcid.dta", clear
gen year = year(date)
keep if year >= 1980 & year <= 2015
drop if investor == "Undisclosed Fund"
bysort investor year: egen indrounds = count(vcid)
bysort year: egen averounds = mean(indrounds)
* supply tracks demand
bysort investor: egen rounds = mean(indrounds)
duplicates drop investor, force
sum rounds, detail

*sum stats

-------------------------------------------------------------
      Percentiles      Smallest
 1%            1              1
 5%            1              1
10%            1              1       Obs              28,793
25%            1              1       Sum of Wgt.      28,793

50%     1.666667                      Mean            2.96837
                        Largest       Std. Dev.      4.191002
75%            3       84.82016
90%     6.272727       86.26536       Variance       17.56449
95%       9.6875       114.7868       Skewness       7.408905
99%     19.61314       150.0009       Kurtosis       124.0809

*investment staff from VentureSource
*VS has less coverage of investors not accurate
use "C:\Users\kecyw\Dropbox\Graduate\Research\Venture_Capital_Paper\Empirical\rawdata\vs_roundinvestor_aum_emp.dta", clear
* staff number not changing, 47% non-missing firm staff no.
duplicates drop entityid year, force
gen tag = 0
replace tag = 1 if missing(invfirmstaff)
tab tag
bysort year: egen totalstaff = total(invfirmstaff)
keep if year >= 1980 & year <= 2015
duplicates drop year, force

set scheme s2color
graph twoway line totalstaff year, graphregion(fcolor(white)) ///
	ylabel(0(10000)70000) xlabel(1980(5)2015) ///
	xtitle("year") ytitle("number of investment staff") 

*plot early stage and late stage firms trend
use "$intpath/SDC_Vxpert_ALL_roundinfo.dta", clear
gen year = year(date)
keep if year >= 1980 & year <= 2015
duplicates drop vcid year, force
merge 1:m vcid year using "$intpath/SDC_Vxpert_ALL_investor_vcid.dta"
keep if _merge == 3

gen earlystage = 1 if stage == "Early Stage" | stage == "Startup/Seed" | stage == "Expansion"
gen earlystage = 1 if stage == "Early Stage" | stage == "Startup/Seed" /* different results */
replace earlystage = 0 if earlystage == .
gen latestage = 1-earlystage
bysort year: egen early = total(earlystage)
bysort year: egen late = total(latestage)
duplicates drop year, force

set scheme s2color
graph twoway line early late year

*entry and exit
duplicates drop investor year, force
bysort investor: gen new = 1 if year == year[1]
replace new = 0 if new == .
*exit means inactive 5 years
bysort investor: gen out = 1 if year[_N] + 5 <= 2015 & year == year[_N]
replace out = 0 if out == .
bysort year: egen numinvestor = count(investor)
bysort year: egen entry = total(new)
bysort year: egen exit = total(out)
duplicates drop year, force
keep if year >= 1985 & year <= 2010
gen entry_rate = entry/numinvestor
gen exit_rate = exit/numinvestor

set scheme s2color
graph twoway line entry_rate exit_rate year
graph twoway line entry exit year