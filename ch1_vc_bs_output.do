********************************************************************************
* File: Do_master_file
* Author: Kecy Wu
* Date: 20211221
* Purpose: Compustat and VC: organizational capital (quantity)
* Robust: 10%-50% of depreciation rate, drop first 5 years of data
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

*clean data
*sample from 1970 - 2014, 16913 firms, 4176 vc-funded (25%)
*log using "$output\bs_vc_121521", replace
use "$finalpath\merged_bs.dta", clear
*keep fiscal year ending in December
gen month = month(datadate)
keep if month == 12
*drop duplicates
gen temp = -xsga 
sort compustatid year temp
duplicates drop compustatid year, force
drop temp
*get year since IPO
sort compustatid year
egen ipoyear = tag(compustatid)
by compustatid: gen yeartoipo = _n-1
keep year CompanyName at cogs emp revt xrd xsga mkvalt sic state vcid compustatid vcfunded yeartoipo
destring sic, replace
*treat missing SGA as zero
replace xrd = 0 if xrd == .
replace xsga = 0 if xsga == .
*xsga - xrd
gen temp = xrd < cogs & xrd > xsga
gen xsga_net = xsga - xrd if temp == 0
replace xsga_net = xsga if temp == 1
*drop financials and regulated utilities
drop if inrange(sic, 6000, 6799)
drop if inrange(sic, 4900, 4999)
*merge to get cpi index
merge m:1 year using "$rawpath\cpi_1970_2014.dta"
drop _merge

*calculate real growth rate of SG&A in sample - 0.11
replace xsga = xsga / (cpi/100)
replace xsga_net = xsga_net / (cpi/100)
sort compustatid year
bysort compustatid: gen g = xsga / xsga[_n-1] - 1 if xsga[_n-1] != 0
bysort compustatid: gen g_net = xsga_net / xsga_net[_n-1] - 1 if xsga_net[_n-1] != 0
summarize g g_net /* too large, 0.85*/
winsor g, gen(Wg) p(0.05)
winsor g_net, gen(Wgnet) p(0.05)
summarize Wg Wgnet 

*assign industry by Fama French classification
gen ind=.
label var ind "FF 17 Ind. Class."

replace ind=1 if inrange(sic,0100,0199) | inrange(sic,0200,0299) | inrange(sic,0700,0799) | inrange(sic,0900,0999) | ///
inrange(sic,2000,2087) | inrange(sic,2090,2092) |inrange(sic,2095,2099) | inrange(sic,5140,5159) | inrange(sic,5180,5191)

replace ind=2 if inrange(sic,1000,1049) | inrange(sic,1060,1069) | inrange(sic,1080,1099) | inrange(sic,1200,1299) | ///
inrange(sic,1400,1499) | inrange(sic,5050,5052)

replace ind=3 if inrange(sic,1300,1329) | inrange(sic,1380,1382) | inrange(sic,1389,1389) | inrange(sic,2900,2912) | ///
inrange(sic,5170,5172)

replace ind=4 if inrange(sic,2200,2284) | inrange(sic,2290,2399) | inrange(sic,3020,3021) | inrange(sic,3100, 3111) | ///
inrange(sic,3130,3131) | inrange(sic,3140,3151) | inrange(sic,3963,3965) | inrange(sic,5130,5139)

replace ind=5 if inrange(sic,2510,2519) | inrange(sic,2590,2599) | inrange(sic,3060,3099) | inrange(sic,3630,3639) | ///
inrange(sic,3650,3652) | inrange(sic,3860,3861) | inrange(sic,3870,3873) | inrange(sic,3910,3911) | inrange(sic,3914,3915)| ///
inrange(sic,3930,3931) | inrange(sic,3940,3949) | inrange(sic,3960,3962) | inrange(sic,5020,5023) | inrange(sic,5064,5064)| ///
inrange(sic,5094,5094) | inrange(sic,5099,5099)

replace ind=6 if inrange(sic,2800,2829) | inrange(sic,2860,2879) | inrange(sic,2890,2899) | inrange(sic,5160,5169)

replace ind=7 if inrange(sic,2100,2199) | inrange(sic,2830,2834) | inrange(sic,2840,2844) | inrange(sic,5120,5122) | inrange(sic,5194,5194)

replace ind=8 if inrange(sic,0800,0899) | inrange(sic,1500,1511) | inrange(sic,1520,1549) | inrange(sic,1600,1699) | ///
inrange(sic,1700,1799) | inrange(sic,2400,2459) | inrange(sic,2490,2499) | inrange(sic,2850,2859) | inrange(sic,2950,2952) | ///
inrange(sic,3200,3200) | inrange(sic,3210,3211) | inrange(sic,3240,3241) | inrange(sic,3250,3259) | inrange(sic,3261,3261) | ///
inrange(sic,3270,3275) | inrange(sic,3280,3281) | inrange(sic,3290,3293) | inrange(sic,3420,3433) | inrange(sic,3440,3442) | ///
inrange(sic,3446,3446) | inrange(sic,3448,3452) | inrange(sic,5030,5039) | inrange(sic,5070,5078) | inrange(sic,5198,5198) | ///
inrange(sic,5210,5211) | inrange(sic,5230,5231) | inrange(sic,5250,5251)

replace ind=9 if inrange(sic,3300,3300) | inrange(sic,3310,3317) | inrange(sic,3320,3325) | inrange(sic,3330,3341) | ///
inrange(sic,3350,3357) | inrange(sic,3360,3369) | inrange(sic,3390,3399)

replace ind=10 if inrange(sic,3410,3412) | inrange(sic,3443,3444) | inrange(sic,3460,3499)

replace ind=11 if inrange(sic,3510,3536) | inrange(sic,3540,3582) | inrange(sic,3585,3586) | inrange(sic,3589,3600) | ///
inrange(sic,3610,3613) | inrange(sic,3620,3629) | inrange(sic,3670,3695) | inrange(sic,3699,3699) | inrange(sic,3810,3812) | ///
inrange(sic,3820,3839) | inrange(sic,3950,3955) | inrange(sic,5060,5060) | inrange(sic,5063,5063) | inrange(sic,5065,5065) | ///
inrange(sic,5080,5081)

replace ind=12 if inrange(sic,3710,3711) | inrange(sic,3714,3714) | inrange(sic,3716,3716) | inrange(sic,3750,3751) | ///
inrange(sic,3792,3792) | inrange(sic,5010,5015) | inrange(sic,5510,5521) | inrange(sic,5530,5531) | inrange(sic,5560,5561) | ///
inrange(sic,5570,5571) | inrange(sic,5590,5599)

replace ind=13 if inrange(sic,3713,3713) | inrange(sic,3715,3715) | inrange(sic,3720,3721) | inrange(sic,3724,3725) | ///
inrange(sic,3728,3728) | inrange(sic,3730,3732) | inrange(sic,3740,3743) | inrange(sic,3760,3769) | inrange(sic,3790,3790) | ///
inrange(sic,3795,3795) | inrange(sic,3799,3799) | inrange(sic,4000,4013) | inrange(sic,4100,4100) | inrange(sic,4110,4121) | ///
inrange(sic,4130,4131) | inrange(sic,4140,4142) | inrange(sic,4150,4151) | inrange(sic,4170,4173) | inrange(sic,4190,4200) | ///
inrange(sic,4210,4231) | inrange(sic,4400,4499) | inrange(sic,4500,4599) | inrange(sic,4600,4699) | inrange(sic,4700,4700) | ///
inrange(sic,4710,4712) | inrange(sic,4720,4742) | inrange(sic,4780,4780) | inrange(sic,4783,4783) | inrange(sic,4785,4785) | ///
inrange(sic,4789,4789)

replace ind=15 if inrange(sic,5260,5261) | inrange(sic,5270,5271) | inrange(sic,5300,5300) | inrange(sic,5310,5311) | ///
inrange(sic,5320,5320) | inrange(sic,5330,5331) | inrange(sic,5334,5334) | inrange(sic,5390,5400) | inrange(sic,5410,5412) | ///
inrange(sic,5420,5421) | inrange(sic,5430,5431) | inrange(sic,5440,5441) | inrange(sic,5450,5451) | inrange(sic,5460,5461) | ///
inrange(sic,5490,5499) | inrange(sic,5540,5541) | inrange(sic,5550,5551) | inrange(sic,5600,5700) | inrange(sic,5710,5722) | ///
inrange(sic,5730,5736) | inrange(sic,5750,5750) | inrange(sic,5800,5813) | inrange(sic,5890,5890) | inrange(sic,5900,5900) | ///
inrange(sic,5910,5912) | inrange(sic,5920,5921) | inrange(sic,5930,5932) | inrange(sic,5940,4963) | inrange(sic,5980,5995) | ///
inrange(sic,5999,5999)

replace ind=17 if ind == .

label define industries 1 "Food" 2 "Mines" ///
3 "Oil" 4 "Clothes" 5 "Durables" 6 "Chems" 7 "Consumption" ///
8 "Construction" 9 "Steel" 10 "Fabrics" 11 "Machine" 12 "Cars" ///
13 "Transportation" 15 "Retail" 17 "Others"

label value ind industries 

*calculate organizational capital - use delta = 0.15, g = 0.11
sort compustatid year
gen OC = 0
gen OCnet = 0
*initial stock
bysort compustatid: gen firstyear = year[1]
bysort compustatid: replace OC = xsga / (0.15 + 0.11) if year == firstyear
bysort compustatid: replace OCnet = xsga_net / (0.15 + 0.11) if year == firstyear
*cumulative
bysort compustatid: replace OC = (1-0.15)*OC[_n-1] + xsga if year != firstyear
bysort compustatid: replace OCnet = (1-0.15)*OCnet[_n-1] + xsga_net if year != firstyear
*O/K ratio
gen OKratio = OC/at 
gen OKratio2 = OCnet/at
drop if missing(OKratio) | missing(OKratio2)

order compustatid year yeartoipo OKratio OKratio2 OC OCnet at ind emp vcid vcfunded 

save "$finalpath\oc_compustat_vc.dta",replace

use "$finalpath\oc_compustat_vc.dta", clear 
*plot graphs of vc vs non-vc funded firms' OC and O/K
collapse (mean) OC OCnet OKratio OKratio2, by(year vcfunded)

graph twoway (line OC year if vcfunded==0) (line OC year if vcfunded==1), ///
	graphregion(fcolor(white)) legend(label(1 non-VC) label(2 VC)) ///
	ylabel(0(500)3000) xlabel(1970(5)2015) ///
	xtitle("year") ytitle("organizational capital") 

graph twoway (line OCnet year if vcfunded==0) (line OCnet year if vcfunded==1), ///
	graphregion(fcolor(white)) legend(label(1 non-VC) label(2 VC)) ///
	ylabel(0(500)1500) xlabel(1970(5)2015) ///
	xtitle("year") ytitle("organizational capital") 

graph twoway (line OKratio year if vcfunded==0) (line OKratio year if vcfunded==1), ///
	graphregion(fcolor(white)) legend(label(1 non-VC) label(2 VC)) ///
	ylabel(0(20)80) xlabel(1970(5)2015) ///
	xtitle("year") ytitle("O/K ratio") 

graph twoway (line OKratio2 year if vcfunded==0) (line OKratio2 year if vcfunded==1), ///
	graphregion(fcolor(white)) legend(label(1 non-VC) label(2 VC)) ///
	ylabel(0(10)40) xlabel(1970(5)2015) ///
	xtitle("year") ytitle("O/K ratio") 

log using "$output\bs_vc_121521", replace
use "$finalpath\oc_compustat_vc.dta", clear 
gen emplog = ln(emp)

reg OKratio2 vcfunded emplog
reg OKratio2 vcfunded##c.yeartoipo emplog
reghdfe OKratio2 vcfunded##c.yeartoipo emplog, absorb(year ind) vce(cluster ind)
*reghdfe OKratio2 vcfunded emplog, absorb(year sic) vce(cluster sic)

xtset compustatid year
gen rev_gr = (revt - l1.revt)/revt
gen emp_gr = (emp - l1.emp)/emp
xtset, clear

reghdfe rev_gr vcfunded##c.yeartoipo emplog, absorb(year ind) vce(cluster ind)
reghdfe emp_gr vcfunded##c.yeartoipo emplog, absorb(year ind) vce(cluster ind)

keep if vcfunded == 1
merge m:1 vcid using "$finalpath\SDC_Vxpert_ALL_vcid_summary.dta"
keep if _merge == 3
gen lnfunding = ln(totalamt)

reghdfe OKratio2 highqualityinvestor lnfunding emplog, absorb(year ind) vce(cluster ind)

reg OKratio2 highqualityinvestor yeartoipo lnfunding emplog, robust
reghdfe OKratio2 highqualityinvestor yeartoipo lnfunding emplog, absorb(year ind) vce(cluster ind)

log close