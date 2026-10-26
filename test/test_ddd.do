/***************************************************************************************************
test_ddd.do

Purpose:
	Tests for ddd: the accounting identities, the closed forms of Table 1 for the
	built-in indices, the agreement of user-supplied indices (ado and Mata, numerical drift) with
	the built-ins (analytic drift), weights, by(), horizons, lagged weights, the log decomposition,
	and the CES index.

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
	assert abs(demography - (mass + selection)) < 1e-10 if !missing(drift)
end

/***************************************************************************************************
Test 1: Variance. Drift = 2Cov(x, dx) and displacement = Var(dx) among survivors
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
ddd x, index(var) clear
_check
assert abs(drift[1] - 2 * `cov') < 1e-10
assert abs(displacement[1] - `var') < 1e-10
di "Test 1 passed"

/***************************************************************************************************
Test 2: User-supplied indices (ado and Mata) give the built-in decomposition
***************************************************************************************************/
di _newline "=== Test 2: ado and Mata indices ==="
capture program drop myvar
program define myvar, rclass
	syntax varname, w(varname)
	tempvar d
	qui sum `varlist' [iw = `w'], meanonly
	qui gen double `d' = (`varlist' - r(mean))^2
	qui sum `d' [iw = `w'], meanonly
	return scalar nu = r(mean)
end
capture program drop mygini
program define mygini, rclass
	* Gini as the mean absolute difference over twice the mean (unweighted)
	syntax varname, w(varname)
	tempvar r
	sort `varlist'
	qui gen double `r' = (2 * _n - 1 - _N) * `varlist'
	qui sum `r', meanonly
	local num = r(sum)
	qui sum `varlist', meanonly
	return scalar nu = `num' / (_N * r(sum))
end
mata:
real scalar mytheil(real matrix X, real colvector w)
{
	real colvector xt
	xt = X / (quadsum(w :* X) / quadsum(w))
	return(quadsum(w :* xt :* ln(xt)) / quadsum(w))
}
end

foreach spec in "var|myvar" "ge 1|mata: mytheil()" "gini|mygini"{
	gettoken a b : spec, parse("|")
	local b = substr("`b'", 2, .)
	use `panel', clear
	ddd x, index(`a') saving(`"`c(tmpdir)'/dd_a.dta"', replace) normalize
	use `panel', clear
	ddd x, index(`b') clear normalize
	_check
	rename (drift displacement demography mass) b_=
	merge 1:1 year0 year1 using `"`c(tmpdir)'/dd_a.dta"', assert(match) nogen
	foreach v in drift displacement demography mass{
		assert reldif(`v', b_`v') < 1e-6
	}
	di "  `a' = `b': ok"
}
di "Test 2 passed"

/***************************************************************************************************
Test 3: Closed forms of Table 1 (normalized positions)
***************************************************************************************************/
di _newline "=== Test 3: closed forms ==="
* Survivors between 2000 and 2001, positions normalized by their mean at each date
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
gen double s1 = e1 / r(mean)
gen double kl = s0 * ln(s0 / s1)
qui sum kl, meanonly
local expm = r(mean) / 2
* Gini: displacement = 2 E[n1 (R1 - R0)], midrank percentile ranks among survivors
foreach t in 0 1{
	sort x`t'
	gen double R`t' = (_n - 0.5) / _N
}
gen double gk = 2 * n1 * (R1 - R0)
qui sum gk, meanonly
local gini = r(mean)

use `panel', clear
ddd x, index(theil) clear normalize
assert abs(displacement[1] - `theil') < 1e-10
use `panel', clear
gen double lx = ln(x)
ddd lx, index(expmean 2) clear
assert abs(displacement[1] - `expm') < 1e-10
use `panel', clear
ddd x, index(gini) clear normalize
assert abs(displacement[1] - `gini') < 1e-10
di "Test 3 passed"

/***************************************************************************************************
Test 4: Every built-in satisfies the identities; analytic drift = numerical drift of the same index
***************************************************************************************************/
di _newline "=== Test 4: built-ins ==="
* Mata functions defined in an ado-file are private to it; loading the file with run makes the
* built-ins visible to the program bi below
capture program drop ddd
mata: mata clear
run "`c(pwd)'/ddd.ado"
capture program drop bi
program define bi, rclass
	* A built-in index called through the ado interface, so that drift is computed numerically
	syntax varname, w(varname) name(string) [par(real 0)]
	tempname nu
	mata: st_numscalar("`nu'", _ddd_builtin("`name'", `par', st_data(., "`varlist'"), ///
		st_data(., "`w'"), 0, J(0, 1, .)))
	return scalar nu = `nu'
end
foreach spec in "mean|" "var|" "power|1.5" "expmean|0.5" "ces|0.5" "ge|2" "ge|0" "ge|1" "atkinson|0.5" ///
	"atkinson|1" "hhi|" "gini|" "top|0.1"{
	gettoken a p : spec, parse("|")
	local p = substr("`p'", 2, .)
	local q = cond("`p'" == "", "", "par(`p')")
	use `panel', clear
	ddd x [fw = 1 + mod(id, 3)], index(`a' `p') saving(`"`c(tmpdir)'/dd_a.dta"', replace)
	use `panel', clear
	gen w = 1 + mod(id, 3)
	ddd x [fw = w], index(bi name(`a') `q') clear
	_check
	rename (drift displacement) b_=
	merge 1:1 year0 year1 using `"`c(tmpdir)'/dd_a.dta"', assert(match) nogen
	_check
	* Rank indices are piecewise linear along the path, so numerical drift is exact unless units
	* cross within the step
	assert reldif(drift, b_drift) < 1e-8
	di "  `a' `p': ok"
}
di "Test 4 passed"

/***************************************************************************************************
Test 5: Frequency weights are equivalent to duplicated units
***************************************************************************************************/
di _newline "=== Test 5: weights ==="
foreach a in hhi gini "top 0.1" var{
	use `panel', clear
	gen w = 1 + mod(id, 3)
	ddd x [fw = w], index(`a') saving(`"`c(tmpdir)'/dd_a.dta"', replace)
	use `panel', clear
	gen w = 1 + mod(id, 3)
	expand w
	bys id year: replace id = id * 10 + _n
	tsset id year
	ddd x, index(`a') clear
	rename (drift displacement entry exit mass) b_=
	merge 1:1 year0 year1 using `"`c(tmpdir)'/dd_a.dta"', assert(match) nogen
	foreach v in drift displacement entry exit mass{
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
ddd x, index(gini) by(grp) horizon(2) saving(`"`c(tmpdir)'/dd_a.dta"', replace)
use `panel', clear
gen str1 grp = cond(mod(id, 2), "a", "b")
keep if grp == "b"
ddd x, index(gini) horizon(2) clear
gen grp = "b"
merge 1:1 grp year0 year1 using `"`c(tmpdir)'/dd_a.dta"', keep(match) nogen
assert _N == 1
_check
* Lagged weights on a balanced panel: drift_lag = 2Cov(x at t-1, dx)
_sim 300 3 2
bys id: keep if _N == 3
tsset id year
gen double dx = F.x - x
qui corr L.x dx if year == 2001, cov
local cov = r(cov_12) * (r(N) - 1) / r(N)
ddd x, index(var) lag(1) clear
assert abs(drift_lag[2] - 2 * `cov') < 1e-10
assert abs(drift + displacement - drift_lag - displacement_lag) < 1e-10 if !missing(drift_lag)
assert missing(drift_lag[1])
di "Test 6 passed"

/***************************************************************************************************
Test 7: Two variables per unit (sales-weighted average productivity)
***************************************************************************************************/
di _newline "=== Test 7: several variables ==="
mata:
real scalar wavg(real matrix X, real colvector w)
{
	return(quadsum(w :* X[., 2] :* X[., 1]) / quadsum(w :* X[., 2]))
}
end
_sim 300 3 3
gen double s = exp(rnormal())
ddd x s, index(mata: wavg()) clear
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
ddd x s, index(mata: wavg()) clear
assert reldif(drift[1], `d') < 1e-8
di "Test 7 passed"

/***************************************************************************************************
Test 8: Log decomposition. Logs of the split in levels
***************************************************************************************************/
di _newline "=== Test 8: log ==="
use `panel', clear
ddd x, index(gini) clear
gen double nuS0 = nu0 + exit
gen double nuS1 = nu1 - entry
keep year0 year1 nu0 nu1 nuS0 nuS1 drift
rename drift drift_lev
tempfile lev
save `lev'
use `panel', clear
ddd x, index(gini) log clear
_check
merge 1:1 year0 year1 using `lev', assert(match) nogen
assert abs(total - ln(nu1 / nu0)) < 1e-10
assert abs(drift - ln((nuS0 + drift_lev) / nuS0)) < 1e-10
assert abs(displacement - ln(nuS1 / (nuS0 + drift_lev))) < 1e-10
assert abs(entry - ln(nu1 / nuS1)) < 1e-10
assert abs(exit - ln(nuS0 / nu0)) < 1e-10
* The Gini is population invariant: no mass
assert abs(mass) < 1e-10
* The mean is linear: initial weights predict the final index, so no displacement
use `panel', clear
ddd x, index(mean) log clear
_check
assert abs(displacement) < 1e-10
* Lagged weights split the same log survivor change
_sim 300 3 2
ddd x, index(var) lag(1) log clear
_check
assert abs(drift + displacement - drift_lag - displacement_lag) < 1e-10 if !missing(drift_lag)
* A nonpositive index has missing log terms
use `panel', clear
replace x = x - 100
ddd x, index(mean) log clear
assert missing(total) & missing(drift)
di "Test 8 passed"

/***************************************************************************************************
Test 9: CES index. Expmean plus ln(N)/g: same survivor terms, mass is the taste for variety
***************************************************************************************************/
di _newline "=== Test 9: ces ==="
use `panel', clear
gen double lx = ln(x)
ddd lx, index(expmean 3) clear
rename (drift displacement demography mass selection) e_=
tempfile em
save `em'
use `panel', clear
gen double lx = ln(x)
ddd lx, index(ces 3) clear
_check
merge 1:1 year0 year1 using `em', assert(match) nogen
assert abs(drift - e_drift) < 1e-10
assert abs(displacement - e_displacement) < 1e-10
assert abs(mass - ln(N1 / N0) / 3) < 1e-10
assert abs(e_mass) < 1e-10
assert abs(selection - e_selection) < 1e-10
* Demography is ln of the survivors' share of exp(3 lx) at t over that at t + 1, divided by 3
use `panel', clear
gen double e = x^3
gen byte s0 = !missing(F.x)
gen byte s1 = !missing(L.x)
qui sum e if year == 2000, meanonly
local E0 = r(sum)
qui sum e if year == 2000 & s0, meanonly
local share0 = r(sum) / `E0'
qui sum e if year == 2001, meanonly
local E1 = r(sum)
qui sum e if year == 2001 & s1, meanonly
local share1 = r(sum) / `E1'
use `panel', clear
gen double lx = ln(x)
ddd lx, index(ces 3) clear
assert abs(demography[1] - ln(`share0' / `share1') / 3) < 1e-10
di "Test 9 passed"

di _newline "All tests passed"
