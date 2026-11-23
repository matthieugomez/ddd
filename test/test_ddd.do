/***************************************************************************************************
test_ddd.do

Purpose:
	Tests for ddd: the accounting identities, the closed forms of Table 1 for the
	built-in statistics, the agreement of user-supplied ado programs (numerical drift) with
	the built-ins (analytic drift), weights, by(), horizons, lagged weights, the CES aggregates,
	rank bands, and the mass parts of demography.

	Run from the package folder: do test/test_ddd.do
***************************************************************************************************/

adopath ++ "`c(pwd)'"
discard

capture program drop _sim
program define _sim
	* Unbalanced panel of positive positions with entry and exit
	args n T seed
	clear
	set seed `seed'
	set obs `n'
	gen id = _n
	gen double x0 = exp(rnormal())
	expand `T'
	bys id: gen year = 2000 + _n - 1
	bys id (year): gen double x = x0 * exp(sum(0.02 + 0.3 * rnormal()))
	drop if runiform() < 0.1
	drop x0
	tsset id year
end

capture program drop _check
program define _check
	* Accounting identities on the output of ddd
	assert abs(total - (drift + displacement + demography)) < 1e-10 if !missing(drift)
	assert abs(demography - (entry + exit)) < 1e-10 if !missing(drift)
	assert abs(entry - (entry_mass + entry_selection)) < 1e-10 if !missing(drift)
	assert abs(exit - (exit_mass + exit_selection)) < 1e-10 if !missing(drift)
end

/***************************************************************************************************
Test 1: Variance. Drift = 2Cov(x, dx) and displacement = Var(dx) among continuing units
***************************************************************************************************/
di _newline "=== Test 1: variance ==="
_sim 300 3 1
tempfile panel
save `panel'
gen double dx = F.x - x
qui sum dx if year == 2000
local var = r(Var) * (r(N) - 1) / r(N)
qui corr x dx if year == 2000 & !missing(dx), cov
local cov = r(cov_12) * (r(N) - 1) / r(N)
ddd x, statistic(var) clear
_check
assert abs(drift[1] - 2 * `cov') < 1e-10
assert abs(displacement[1] - `var') < 1e-10
di "Test 1 passed"

/***************************************************************************************************
Test 2: User-supplied ado programs give the built-in decomposition
***************************************************************************************************/
di _newline "=== Test 2: ado programs ==="
capture program drop myvar
program define myvar, rclass
	syntax varname [iw]
	tempvar d
	qui sum `varlist' [`weight'`exp'], meanonly
	qui gen double `d' = (`varlist' - r(mean))^2
	qui sum `d' [`weight'`exp'], meanonly
	return scalar nu = r(mean)
end
capture program drop mygini
program define mygini, rclass
	* Gini as the mean absolute difference over twice the mean (unweighted)
	syntax varname [iw]
	tempvar r
	sort `varlist'
	qui gen double `r' = (2 * _n - 1 - _N) * `varlist'
	qui sum `r', meanonly
	local num = r(sum)
	qui sum `varlist', meanonly
	return scalar nu = `num' / (_N * r(sum))
end
capture program drop mytheil
program define mytheil, rclass
	syntax varname [iw]
	tempvar k
	qui sum `varlist' [`weight'`exp'], meanonly
	qui gen double `k' = (`varlist' / r(mean)) * ln(`varlist' / r(mean))
	qui sum `k' [`weight'`exp'], meanonly
	return scalar nu = r(mean)
end

foreach spec in "var|myvar|" "ge 1|mytheil|normalize" "gini|mygini|normalize"{
	tokenize "`spec'", parse("|")
	local a `1'
	local b `3'
	local o `5'
	use `panel', clear
	ddd x, statistic(`a') saving(`"`c(tmpdir)'/dd_a.dta"', replace) `o'
	use `panel', clear
	ddd x, statistic(`b') clear `o'
	_check
	rename (drift displacement demography entry_mass exit_mass) b_=
	merge 1:1 year0 year1 using `"`c(tmpdir)'/dd_a.dta"', assert(match) nogen
	foreach v in drift displacement demography entry_mass exit_mass{
		assert reldif(`v', b_`v') < 1e-6
	}
	di "  `a' = `b': ok"
}
* Masses reach the program: a weighted ado statistic equals the weighted built-in
use `panel', clear
gen w = 1 + mod(id, 3)
ddd x [fw = w], statistic(var) saving(`"`c(tmpdir)'/dd_a.dta"', replace)
use `panel', clear
gen w = 1 + mod(id, 3)
ddd x [fw = w], statistic(myvar) clear
_check
rename (drift displacement demography entry_mass exit_mass) b_=
merge 1:1 year0 year1 using `"`c(tmpdir)'/dd_a.dta"', assert(match) nogen
foreach v in drift displacement demography entry_mass exit_mass{
	assert reldif(`v', b_`v') < 1e-6
}
di "Test 2 passed"

/***************************************************************************************************
Test 3: Closed forms of Table 1 (normalized positions)
***************************************************************************************************/
di _newline "=== Test 3: closed forms ==="
* Continuing units between 2000 and 2001, positions normalized by their mean at each date
use `panel', clear
gen double x1 = F.x
keep if year == 2000 & !missing(x1)
rename x x0
foreach t in 0 1{
	qui sum x`t', meanonly
	gen double n`t' = x`t' / r(mean)
}
* Theil: displacement = E[n1 ln(n1/n0)]
gen double k = n1 * ln(n1 / n0)
qui sum k, meanonly
local theil = r(mean)
* Exponential mean of log positions, gamma = 2: displacement = E[s ln(s/s')] / gamma
gen double e0 = exp(2 * ln(x0))
gen double e1 = exp(2 * ln(x1))
qui sum e0, meanonly
gen double s0 = e0 / r(mean)
qui sum e1, meanonly
gen double c1 = e1 / r(mean)
gen double kl = s0 * ln(s0 / c1)
qui sum kl, meanonly
local expm = r(mean) / 2
* Gini: displacement = 2 E[n1 (R1 - R0)], midrank percentile ranks among continuing units
foreach t in 0 1{
	sort x`t'
	gen double R`t' = (_n - 0.5) / _N
}
gen double gk = 2 * n1 * (R1 - R0)
qui sum gk, meanonly
local gini = r(mean)

use `panel', clear
ddd x, statistic(theil) clear normalize
assert abs(displacement[1] - `theil') < 1e-10
use `panel', clear
gen double lx = ln(x)
ddd lx, statistic(expmean 2) clear
assert abs(displacement[1] - `expm') < 1e-10
use `panel', clear
ddd x, statistic(gini) clear normalize
assert abs(displacement[1] - `gini') < 1e-10
di "Test 3 passed"

/***************************************************************************************************
Test 4: Every built-in satisfies the identities; analytic drift = numerical drift of the same statistic
***************************************************************************************************/
di _newline "=== Test 4: built-ins ==="
* Subprograms of an ado-file are private to it; loading the file with run makes _ddd_builtin visible to
* the program bi below
foreach p in ddd _ddd_saving _ddd_nu _ddd_drift _ddd_correction _ddd_builtin _ddd_rankweight _ddd_ado{
	capture program drop `p'
}
run "`c(pwd)'/ddd.ado"
capture program drop bi
program define bi, rclass
	* A built-in statistic called through the ado interface, so that drift is computed numerically
	syntax varname [iw/], name(string) [par(string) par2(string)]
	if "`par'" == "" local par .
	if "`par2'" == "" local par2 .
	tempvar cell nu
	gen byte `cell' = 1
	_ddd_builtin `varlist', mass(`exp') cell(`cell') nu(`nu') name(`name') par(`par') par2(`par2')
	return scalar nu = `nu'[1]
end
foreach spec in "mean|" "var|" "powermean|0.5" "powermean|0" "powersum|2" "expmean|0.5" ///
	"expsum|0.5" "ge|2" "ge|0" "ge|1" "atkinson|0.5" "atkinson|1" "hhi|" "gini|" "share|0.9 1" ///
	"share|0.5 0.9" "average|0.5 0.9" "average|0.9 1"{
	gettoken a p : spec, parse("|")
	local p = substr("`p'", 2, .)
	local q
	if `: word count `p'' == 1 local q par(`p')
	if `: word count `p'' == 2 local q par(`: word 1 of `p'') par2(`: word 2 of `p'')
	use `panel', clear
	ddd x [fw = 1 + mod(id, 3)], statistic(`a' `p') saving(`"`c(tmpdir)'/dd_a.dta"', replace)
	use `panel', clear
	gen w = 1 + mod(id, 3)
	ddd x [fw = w], statistic(bi name(`a') `q') clear
	_check
	rename (drift displacement) b_=
	merge 1:1 year0 year1 using `"`c(tmpdir)'/dd_a.dta"', assert(match) nogen
	_check
	* Rank statistics are piecewise linear along the path, so numerical drift is exact unless units
	* cross within the step
	assert reldif(drift, b_drift) < 1e-8
	di "  `a' `p': ok"
}
di "Test 4 passed"

/***************************************************************************************************
Test 5: Frequency weights are equivalent to duplicated units
***************************************************************************************************/
di _newline "=== Test 5: weights ==="
foreach a in hhi gini "share 0.9 1" "average 0.9 1" "powersum 2" var{
	use `panel', clear
	gen w = 1 + mod(id, 3)
	ddd x [fw = w], statistic(`a') saving(`"`c(tmpdir)'/dd_a.dta"', replace)
	use `panel', clear
	gen w = 1 + mod(id, 3)
	expand w
	bys id year: replace id = id * 10 + _n
	tsset id year
	ddd x, statistic(`a') clear
	rename (drift displacement entry exit entry_mass exit_mass) b_=
	merge 1:1 year0 year1 using `"`c(tmpdir)'/dd_a.dta"', assert(match) nogen
	foreach v in drift displacement entry exit entry_mass exit_mass{
		assert abs(`v' - b_`v') < 1e-10
	}
	di "  `a': ok"
}
di "Test 5 passed"

/***************************************************************************************************
Test 6: by(), horizons, and lagged weights
***************************************************************************************************/
di _newline "=== Test 6: by(), horizon(), lag() ==="
* by(): same as running each group separately
use `panel', clear
gen str1 grp = cond(mod(id, 2), "a", "b")
ddd x, statistic(gini) by(grp) horizon(2) saving(`"`c(tmpdir)'/dd_a.dta"', replace)
use `panel', clear
gen str1 grp = cond(mod(id, 2), "a", "b")
keep if grp == "b"
ddd x, statistic(gini) horizon(2) clear
gen grp = "b"
merge 1:1 grp year0 year1 using `"`c(tmpdir)'/dd_a.dta"', keep(match) nogen
assert _N == 1
_check
* Lagged weights on a balanced panel: drift = 2Cov(x at t-1, dx). The first interval has no lagged date and
* is dropped. chi_lag recovers drift and displacement with weights at the start.
_sim 300 3 2
bys id: keep if _N == 3
tsset id year
gen double dx = F.x - x
qui corr L.x dx if year == 2001, cov
local cov = r(cov_12) * (r(N) - 1) / r(N)
ddd x, statistic(var) saving(`"`c(tmpdir)'/dd_a.dta"', replace)
ddd x, statistic(var) lag(1) clear
assert _N == 1 & year0[1] == 2001
_check
assert abs(drift[1] - 2 * `cov') < 1e-10
rename (drift displacement) l_=
merge 1:1 year0 year1 using `"`c(tmpdir)'/dd_a.dta"', keep(match) nogen
assert _N == 1
assert abs(drift + displacement - l_drift - l_displacement) < 1e-10
assert reldif(chi_lag, (l_drift - drift) / displacement) < 1e-8
assert reldif(displacement, l_displacement / (1 - chi_lag)) < 1e-8
assert reldif(drift, l_drift - chi_lag * l_displacement / (1 - chi_lag)) < 1e-8
di "Test 6 passed"

/***************************************************************************************************
Test 7: Two variables per unit (sales-weighted average productivity)
***************************************************************************************************/
di _newline "=== Test 7: several variables ==="
capture program drop wavg
program define wavg, rclass
	syntax varlist(min=2 max=2) [iw]
	gettoken x s : varlist
	tempvar sx
	qui gen double `sx' = `s' * `x'
	qui sum `sx' [`weight'`exp'], meanonly
	local num = r(mean)
	qui sum `s' [`weight'`exp'], meanonly
	return scalar nu = `num' / r(mean)
end
_sim 300 3 3
gen double s = exp(rnormal())
ddd x s, statistic(wavg) clear
_check
* Drift of a ratio of sums: E[s dx]/E[s] + E[(x - nu) ds]/E[s]
_sim 300 3 3
gen double s = exp(rnormal())
gen double dx = F.x - x
gen double ds = F.s - s
keep if year == 2000 & !missing(dx)
qui sum s, meanonly
local Es = r(mean)
gen double sx = s * x
qui sum sx, meanonly
local nu = r(mean) / `Es'
gen double k = (s * dx + (x - `nu') * ds) / `Es'
qui sum k, meanonly
local d = r(mean)
_sim 300 3 3
gen double s = exp(rnormal())
ddd x s, statistic(wavg) clear
assert reldif(drift[1], `d') < 1e-8
di "Test 7 passed"

/***************************************************************************************************
Test 8: CES aggregate. Expsum is expmean plus ln(N)/g: same drift and displacement, the mass parts are the
love of variety. Powermean 1 is the mean.
***************************************************************************************************/
di _newline "=== Test 8: expsum, powermean ==="
use `panel', clear
gen double lx = ln(x)
ddd lx, statistic(expmean 3) clear
rename (drift displacement entry_mass entry_selection exit_mass exit_selection) e_=
tempfile em
save `em'
use `panel', clear
gen double lx = ln(x)
ddd lx, statistic(expsum 3) clear
_check
merge 1:1 year0 year1 using `em', assert(match) nogen
assert abs(drift - e_drift) < 1e-10
assert abs(displacement - e_displacement) < 1e-10
assert abs(entry_mass - ln(N1 / NC) / 3) < 1e-10
assert abs(exit_mass - ln(NC / N0) / 3) < 1e-10
assert abs(e_entry_mass) < 1e-10 & abs(e_exit_mass) < 1e-10
assert abs(entry_selection - e_entry_selection) < 1e-10
assert abs(exit_selection - e_exit_selection) < 1e-10
* Demography is ln of the continuing units' share of exp(3 lx) at t over that at t + 1, divided by 3
use `panel', clear
gen double e = x^3
gen byte s0 = !missing(F.x)
gen byte c1 = !missing(L.x)
qui sum e if year == 2000, meanonly
local E0 = r(sum)
qui sum e if year == 2000 & s0, meanonly
local share0 = r(sum) / `E0'
qui sum e if year == 2001, meanonly
local E1 = r(sum)
qui sum e if year == 2001 & c1, meanonly
local share1 = r(sum) / `E1'
use `panel', clear
gen double lx = ln(x)
ddd lx, statistic(expsum 3) clear
assert abs(demography[1] - ln(`share0' / `share1') / 3) < 1e-10
use `panel', clear
ddd x, statistic(mean) saving(`"`c(tmpdir)'/dd_a.dta"', replace)
use `panel', clear
ddd x, statistic(powermean 1) clear
rename (nu0 nu1 drift displacement entry exit) b_=
merge 1:1 year0 year1 using `"`c(tmpdir)'/dd_a.dta"', assert(match) nogen
foreach v in nu0 nu1 drift displacement entry exit{
	assert reldif(`v', b_`v') < 1e-12
}
di "Test 8 passed"

/***************************************************************************************************
Test 9: Rank bands. Each component of the decile shares sums to zero, since the shares sum to one.
average 0 1 is the mean, and each component of the decile averages, weighted by 1/10, sums to that of
the mean, with and without lagged weights.
***************************************************************************************************/
di _newline "=== Test 9: rank bands ==="
forvalues d = 1/10{
	use `panel', clear
	ddd x [fw = 1 + mod(id, 3)], statistic(share `=(`d' - 1) / 10' `=`d' / 10') clear
	_check
	gen byte decile = `d'
	tempfile dec`d'
	save `dec`d''
}
clear
forvalues d = 1/10{
	append using `dec`d''
}
collapse (sum) nu0 nu1 total drift displacement demography, by(year0 year1)
assert abs(nu0 - 1) < 1e-10 & abs(nu1 - 1) < 1e-10
foreach v in total drift displacement demography{
	assert abs(`v') < 1e-10
}
foreach l in 0 1{
	use `panel', clear
	ddd x [fw = 1 + mod(id, 3)], statistic(mean) lag(`l') saving(`"`c(tmpdir)'/dd_a.dta"', replace)
	use `panel', clear
	ddd x [fw = 1 + mod(id, 3)], statistic(average 0 1) lag(`l') clear
	rename (nu0 nu1 drift displacement entry exit) b_=
	merge 1:1 year0 year1 using `"`c(tmpdir)'/dd_a.dta"', assert(match) nogen
	foreach v in nu0 nu1 drift displacement entry exit{
		assert abs(`v' - b_`v') < 1e-8 * abs(nu0)
	}
	forvalues d = 1/10{
		use `panel', clear
		ddd x [fw = 1 + mod(id, 3)], statistic(average `=(`d' - 1) / 10' `=`d' / 10') lag(`l') clear
		_check
		gen byte decile = `d'
		tempfile dec`d'
		save `dec`d''
	}
	clear
	forvalues d = 1/10{
		append using `dec`d''
	}
	foreach v in nu0 nu1 drift displacement demography{
		replace `v' = `v' / 10
	}
	collapse (sum) nu0 nu1 drift displacement demography, by(year0 year1)
	rename (nu0 nu1 drift displacement demography) s_=
	merge 1:1 year0 year1 using `"`c(tmpdir)'/dd_a.dta"', assert(match) nogen
	foreach v in nu0 nu1 drift displacement demography{
		assert abs(s_`v' - `v') < 1e-8 * abs(nu0)
	}
}
di "Test 9 passed"

/***************************************************************************************************
Test 10: Lagged weights for a rank band. Its correction is that of the upper bands [lo, 1] and [hi, 1]:
for share, their difference; for average, their difference weighted by the bands' widths.
***************************************************************************************************/
di _newline "=== Test 10: lagged weights for rank bands ==="
foreach s in share average{
	local j 0
	foreach b in "0.5 1" "0.9 1" "0.5 0.9"{
		local ++j
		use `panel', clear
		ddd x [fw = 1 + mod(id, 3)], statistic(`s' `b') saving(`"`c(tmpdir)'/dd_a.dta"', replace)
		ddd x [fw = 1 + mod(id, 3)], statistic(`s' `b') lag(1) clear
		rename drift l
		merge 1:1 year0 year1 using `"`c(tmpdir)'/dd_a.dta"', keep(match) keepusing(drift) nogen
		gen double c`j' = l - drift
		rename drift d`j'
		keep year0 c`j' d`j'
		tempfile band`j'
		save `band`j''
	}
	use `band1', clear
	merge 1:1 year0 using `band2', assert(match) nogen
	merge 1:1 year0 using `band3', assert(match) nogen
	assert _N == 1
	if "`s'" == "share"{
		assert abs(c3 - (c1 - c2)) < 1e-12
		assert abs(d3 - (d1 - d2)) < 1e-12
	}
	else{
		assert reldif(c3, (0.5 * c1 - 0.1 * c2) / 0.4) < 1e-10
		assert reldif(d3, (0.5 * d1 - 0.1 * d2) / 0.4) < 1e-10
	}
	di "  `s': ok"
}
di "Test 10 passed"

/***************************************************************************************************
Test 11: Mass parts of demography: adding entrants, or removing exiters, distributed like the continuing
units. Replicating the population k times divides the HHI by k and multiplies powersum a by k^(1/a),
and leaves the average of a rank band unchanged. normalize requires a scale-invariant statistic.
***************************************************************************************************/
di _newline "=== Test 11: mass ==="
use `panel', clear
ddd x, statistic(hhi) clear
assert reldif(exit_mass, (nu0 + exit) * (1 - NC / N0)) < 1e-10
assert reldif(entry_mass, (nu1 - entry) * (NC / N1 - 1)) < 1e-10
use `panel', clear
ddd x, statistic(powersum 2) clear
assert reldif(exit_mass, (nu0 + exit) * (1 - (N0 / NC)^(1 / 2))) < 1e-10
assert reldif(entry_mass, (nu1 - entry) * ((N1 / NC)^(1 / 2) - 1)) < 1e-10
use `panel', clear
ddd x, statistic(average 0.2 0.7) clear
assert abs(exit_mass) < 1e-12 & abs(entry_mass) < 1e-12
use `panel', clear
rcof "ddd x, statistic(average 0.9 1) normalize clear" == 198
rcof "ddd x, statistic(var) normalize clear" == 198
di "Test 11 passed"

di _newline "All tests passed"
