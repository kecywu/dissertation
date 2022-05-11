********************************************************************************
* File: Do_master_file
* Author: Kecy Wu
* Date: 20211215
* Purpose: VC dataset with stage, total rounds, total funding amount, investor experience
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

**create names_only datasets for reclink
use "$rawpath/SDC_Vxpert_ALL.dta",clear
drop CompanyCurrentSituationDate_3 Company6digitCUSIP_3 CompanyCompetitors_3
renvarlab, postdrop(2)
*take 0604 to get vcid
drop if CompanyName == "Undisclosed Company"
drop if missing(CompanyName)
duplicates drop CompanyName, force
replace CompanyName = "ALPHABET INC" if CompanyName == "Google, Inc."
replace CompanyName = "STARBUCKS" if CompanyName == "Starbucks Coffee & Company"
merge m:1 CompanyName using "$intpath/VXpert_names_all_standardized.dta"
keep CompanyName CompanyFoundingDate RoundDates RoundAmtDisclosed000 RoundNumbers NumberofInvestorseaRound RoundInfoDateDisclosedAmt CompanyStageLevel1ateachRo vcid
*keep if CompanyNationCode == "US", in line with merged compustat 
save "$intpath/SDC_Vxpert_ALL_cleaned.dta", replace

*parse variables 
use "$intpath/SDC_Vxpert_ALL_cleaned.dta", clear
local round_vars = "RoundDates RoundAmtDisclosed000 RoundNumbers NumberofInvestorseaRound CompanyStageLevel1ateachRo"
foreach var of local round_vars {
	split `var', parse("*")
	drop `var'
}
reshape long `round_vars', i(vcid) j(j) string
drop j
drop if missing(RoundDates)

*correct dates
gen date = date(CompanyFoundingDate,"DMY")
format date %td
gen foundyear = year(date)
drop date
gen date = date(RoundDates,"MDY")
format date %td

*destring variables
local destring_vars = "RoundAmtDisclosed000 RoundNumbers NumberofInvestorseaRound"
foreach var of local destring_vars{
	replace `var' = subinstr(`var',",","",.)
}
destring `destring_vars', replace
ren RoundAmtDisclosed000 roundamt
ren RoundNumbers roundnum
ren NumberofInvestorseaRound numinvestor
ren CompanyStageLevel1ateachRo stage
keep vcid CompanyName RoundInfoDateDisclosedAmt roundamt roundnum numinvestor stage date foundyear
sort vcid roundnum
save "$intpath/SDC_Vxpert_ALL_roundinfo.dta", replace
keep vcid roundnum numinvestor RoundInfoDateDisclosedAmt
save "$intpath/SDC_Vxpert_ALL_investorinfo.dta", replace

*parse investor variable
use "$intpath/SDC_Vxpert_ALL_cleaned.dta", clear
keep vcid RoundInfoDateDisclosedAmt
split RoundInfoDateDisclosedAmt, parse("*")
drop RoundInfoDateDisclosedAmt
reshape long RoundInfoDateDisclosedAmt, i(vcid) j(j)
keep if !missing(RoundInfoDateDisclosedAmt)
replace RoundInfoDateDisclosedAmt = subinstr(RoundInfoDateDisclosedAmt, "  ", "*", 2)
split RoundInfoDateDisclosedAmt, parse("*") limit(3)
keep vcid RoundInfoDateDisclosedAmt1 RoundInfoDateDisclosedAmt3
ren RoundInfoDateDisclosedAmt1 RoundDates
ren RoundInfoDateDisclosedAmt3 investor
gen date = date(RoundDates,"MDY")
format date %td
drop RoundDates
save "$intpath/SDC_Vxpert_ALL_investor_vcid.dta", replace

**define high quality investor - number of rounds participated is 90th percentile
use "$intpath/SDC_Vxpert_ALL_investor_vcid.dta", clear
drop if investor == "Undisclosed Fund"
drop if missing(investor)
gen round = _n
bysort investor: egen numdeals = count(round)
duplicates drop investor, force 
sort numdeals
xtile pct = numdeals, nq(100)
keep investor numdeals pct 
gen high_quality = 0
replace high_quality = 1 if pct >= 90 /* 30 rounds of investing */
save "$finalpath/SDC_Vxpert_ALL_investor_pct.dta", replace

use "$intpath/SDC_Vxpert_ALL_investor_vcid.dta", clear
merge m:1 investor using "$finalpath/SDC_Vxpert_ALL_investor_pct.dta"
drop if missing(investor)
sort vcid
bysort vcid: egen numhighquality = total(high_quality)
duplicates drop vcid, force
keep vcid numhighquality
save "$finalpath/SDC_Vxpert_ALL_vc_quality.dta", replace

**calculate total number of rounds, funding, investors, stage, whether high quality for each vcid
use "$intpath/SDC_Vxpert_ALL_roundinfo.dta",clear
keep vcid CompanyName roundamt roundnum numinvestor stage foundyear
*all started in first round
bysort vcid: gen firststage = stage[1]
drop stage
bysort vcid: egen totalamt = total(roundamt)
bysort vcid: egen totalround = count(roundnum)
bysort vcid: egen totalinvestor = total(numinvestor)
drop roundamt roundnum numinvestor
duplicates drop vcid, force
merge 1:1 vcid using "$finalpath/SDC_Vxpert_ALL_vc_quality.dta"
keep if _merge == 3
drop _merge
gen highqualityinvstor = 1 if numhighquality/totalinvestor >= 0.5
replace highqualityinvstor = 0 if highqualityinvstor == .
drop if totalamt < 0
tab firststage
gen earlystage = 1 if firststage == "Early Stage" | firststage == "Startup/Seed" | firststage == "Expansion"
replace earlystage = 0 if earlystage == .
save "$finalpath/SDC_Vxpert_ALL_vcid_summary.dta", replace


estpost summarize totalamt totalround totalinvestor earlystage, listwise
esttab using "$output\vcstats.tex", cells("count mean sd min max") nomtitle nonumber