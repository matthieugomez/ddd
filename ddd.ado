* version 0.3.0, November 2026
* Decompose the change in a distributional statistic into drift, displacement, and demography
* (Gomez, "The Micro Origins of Aggregate Trends: Drift, Displacement, and Demography").
*
* Drift is the sum of continuing units' changes times their weights at the start, where a unit's weight is the
* derivative of the statistic with respect to its position; displacement is the rest of the change among
* continuing units; and demography is the change in the full population minus the change among continuing
* units. Built-in statistics use their analytic weights (_ddd_builtin). A user-supplied statistic (an ado
* program) only needs to compute the statistic: drift is then the derivative of the statistic along
* continuing units' changes, computed by central differences (_ddd_drift).
*
* All computations are by cell, a group at a date: the statistic of a set of units in each cell (_ddd_nu)
* and drift (_ddd_drift), with the masses of the units passed in mass().

program define ddd
	version 16
	syntax varlist(numeric) [if] [in] [fweight iweight pweight], STATistic(string asis) ///
		[saving(string asis) clear HORizon(integer 1) lag(integer 0) NORMalize step(real 1e-5) ///
		by(varlist)]

	/* 0: Check inputs */
	cap qui tsset
	if _rc{
		di as error "The dataset must be declared as a panel before applying the decomposition (e.g. tsset panelvariable timevariable)"
		exit 198
	}
	local id = r(panelvar)
	local time = r(timevar)
	local delta = r(tdelta)
	if "`id'" == "" | "`time'" == ""{
		di as error "The dataset must be declared as a panel before applying the decomposition (e.g. tsset panelvariable timevariable)"
		exit 198
	}

	if `"`saving'`clear'"' == ""{
		di as error "You need to specify either the option saving(filename) or the option clear. The first saves the output in an external file while the second replaces the existing dataset."
		exit 198
	}
	if `"`saving'"' != "" & "`clear'" != ""{
		di as error "Options saving(filename) and clear are mutually exclusive; specify only one."
		exit 198
	}
	if `"`saving'"' != ""{
		_ddd_saving `saving'
		local savefile `"`s(filename)'"'
		local replace `s(replace)'
	}
	if `horizon' <= 0{
		di as error "horizon() must be a positive number of periods"
		exit 198
	}
	if `lag' < 0{
		di as error "lag() must be a positive number of periods"
		exit 198
	}
	if `step' <= 0{
		di as error "step() must be positive"
		exit 198
	}

	local others `id' `time' `varlist'
	local overlap : list by & others
	if "`overlap'" != ""{
		di as error "by() must not contain the panel variable, the time variable, or the decomposed variables"
		exit 198
	}
	if `: list posof "_ddd_w" in varlist'{
		di as error "_ddd_w is a reserved name; rename this variable"
		exit 198
	}

	/* 1: Parse the statistic */
	* kind is builtin (analytic weights) or ado (an rclass program)
	local statistic `statistic'
	local par .
	local par2 .
	local cmd
	gettoken name rest : statistic
	local builtins mean var powermean powersum expmean expsum ge theil mld atkinson gini share average hhi
	if `: list name in builtins'{
		local kind builtin
		if inlist("`name'", "share", "average"){
			* share lo hi and average lo hi: the share held by, and the average of, the units whose percentile
			* rank lies between lo and hi
			local rest = trim(`"`rest'"')
			if `: word count `rest'' != 2{
				di as error "statistic(`name' # #) requires two ranks between 0 and 1, e.g. statistic(`name' 0.9 1) for the top decile"
				exit 198
			}
			local lo : word 1 of `rest'
			local hi : word 2 of `rest'
			cap confirm number `lo'
			local rc1 = _rc
			cap confirm number `hi'
			if `rc1' | _rc | !(`lo' >= 0 & `lo' < `hi' & `hi' <= 1){
				di as error "statistic(`name' # #) requires two ranks 0 <= lo < hi <= 1, e.g. statistic(`name' 0.9 1) for the top decile"
				exit 198
			}
			local par = `lo'
			local par2 = `hi'
		}
		else if inlist("`name'", "powermean", "powersum", "expmean", "expsum", "ge", "atkinson"){
			local rest = trim(`"`rest'"')
			cap confirm number `rest'
			if _rc{
				di as error "statistic(`name' #) requires a numeric parameter"
				exit 198
			}
			local par = `rest'
			if inlist("`name'", "powersum", "expmean", "expsum") & `par' == 0{
				di as error "statistic(`name' #) requires a nonzero parameter"
				exit 198
			}
		}
		else if trim(`"`rest'"') != ""{
			di as error "statistic(`name') takes no parameter"
			exit 198
		}
		if "`name'" == "theil"{
			local name ge
			local par 1
		}
		if "`name'" == "mld"{
			local name ge
			local par 0
		}
		if `: word count `varlist'' > 1{
			di as error "Built-in statistics take a single variable"
			exit 198
		}
		* normalize splits the change in positions divided by their mean, which leaves the statistic unchanged
		* only if it is scale invariant
		if "`normalize'" != "" & inlist("`name'", "mean", "var", "powermean", "powersum", "expmean", "expsum", ///
			"average"){
			di as error "normalize requires a scale-invariant statistic (ge, theil, mld, atkinson, gini, share, hhi)"
			exit 198
		}
	}
	else{
		local kind ado
		cap program list `name'
		if _rc{
			cap findfile `name'.ado
			if _rc{
				di as error "`name' is neither a built-in statistic nor a program (built-in statistics: `builtins')"
				exit 111
			}
		}
		* No comma without options: a program that declares none rejects a trailing comma
		local cmd `name' `varlist' [iw=_ddd_w]
		if trim(`"`rest'"') != "" local cmd `cmd', `rest'
	}

	/* 2: Sample and masses */
	marksample touse
	if "`by'" != "" markout `touse' `by', strok
	tempvar wgt
	if "`weight'" != "" qui gen double `wgt' `exp' if `touse'
	else qui gen double `wgt' = 1 if `touse'
	cap assert `wgt' > 0 & !missing(`wgt') if `touse'
	if _rc{
		di as error "Weights must be positive and nonmissing"
		exit 198
	}
	if "`kind'" == "builtin"{
		local positive = inlist("`name'", "ge", "atkinson", "powermean", "powersum")
		if `positive'{
			cap assert `varlist' > 0 if `touse'
			if _rc{
				di as error "statistic(`statistic') requires `varlist' > 0"
				exit 198
			}
		}
		if "`name'" == "hhi"{
			cap assert `varlist' >= 0 if `touse'
			if _rc{
				di as error "statistic(hhi) requires `varlist' >= 0"
				exit 198
			}
		}
	}
	if mod(`lag', `delta') != 0{
		di as error "lag() must be a multiple of the panel's time delta (`delta')"
		exit 198
	}

	local tfmt : format `time'
	local stat kind(`kind') name(`name') par(`par') par2(`par2') cmd(`cmd') vnames(`varlist')

	/* 3: The units, in a separate frame. A cell is a group at a date. */
	tempvar gid cell Nall nuall cont N1 NC nuC0 nuC1 drift disp m0 m1 mC0 mC1 three NA corr
	tempname work
	tempfile output
	frame put `by' `id' `time' `varlist' `wgt' if `touse', into(`work')
	frame `work'{
		if "`by'" != "" qui egen long `gid' = group(`by'), missing
		else qui gen long `gid' = 1
		qui egen long `cell' = group(`gid' `time')
		qui egen double `Nall' = total(`wgt'), by(`cell')

		* Continuing units are the units observed in the same group at the start and at the end of the
		* interval, horizon() periods later. They keep their mass at the start at both dates. x1 holds their
		* positions at the end, and N1 the mass of all units at the end.
		sort `id' `time'
		qui tsset `id' `time', delta(`delta')
		qui gen byte `cont' = F`horizon'.`gid' == `gid'
		local x1
		foreach v of local varlist{
			tempvar v1
			qui gen double `v1' = F`horizon'.`v' if `cont'
			local x1 `x1' `v1'
		}
		qui gen double `N1' = F`horizon'.`Nall' if `cont'

		* With lag(), xl holds the positions lag() time units before the start of the continuing units also
		* observed then in the same group (the units observed at the three dates)
		if `lag' > 0{
			local k = `lag' / `delta'
			qui gen byte `three' = `cont' & L`k'.`gid' == `gid'
			local xl
			foreach v of local varlist{
				tempvar vl
				qui gen double `vl' = L`k'.`v' if `three'
				local xl `xl' `vl'
			}
		}

		/* 4: The statistic of all units at each date, and of continuing units at the start and the end */
		* Sort once by cell: egen, by(cell) then does not sort the data again
		sort `cell'
		_ddd_nu `varlist', mass(`wgt') cell(`cell') gen(`nuall') `stat'
		qui egen double `NC' = total(`wgt' * `cont'), by(`cell')
		_ddd_nu `varlist' if `cont', mass(`wgt') cell(`cell') gen(`nuC0') `stat'
		_ddd_nu `x1' if `cont', mass(`wgt') cell(`cell') gen(`nuC1') `stat'

		* Drift: the continuing units' changes times their weights at the start. Displacement: the rest of the
		* change among continuing units.
		_ddd_drift if `cont', at(`varlist') from(`varlist') to(`x1') mass(`wgt') cell(`cell') gen(`drift') ///
			`stat' step(`step') `normalize'
		qui gen double `disp' = `nuC1' - `nuC0' - `drift' if `cont'

		* Mass parts of entry and exit: the statistic of the continuing units, scaled to the mass of all units
		* at the start and at the end, as if entrants and exiters were distributed like the continuing units.
		* They are zero for a statistic that does not depend on the number of units.
		if "`kind'" == "ado" | inlist("`name'", "hhi", "powersum", "expsum"){
			qui gen double `m0' = `wgt' * `Nall' / `NC' if `cont'
			qui gen double `m1' = `wgt' * `N1' / `NC' if `cont'
			_ddd_nu `varlist' if `cont', mass(`m0') cell(`cell') gen(`mC0') `stat'
			_ddd_nu `x1' if `cont', mass(`m1') cell(`cell') gen(`mC1') `stat'
		}
		else{
			qui gen double `mC0' = `nuC0'
			qui gen double `mC1' = `nuC1'
		}

		* Lagged weights: the correction to drift (drift with lagged weights minus drift). The displacement of a
		* rank band can be close to zero, so its correction combines those of the upper bands [lo, 1] and
		* [hi, 1]: share lo hi is share lo 1 minus share hi 1, and average lo hi is
		* ((1 - lo) average lo 1 - (1 - hi) average hi 1) / (hi - lo). An upper band that holds all units or
		* none has no correction.
		if `lag' > 0{
			qui egen double `NA' = total(`wgt' * `three'), by(`cell')
			local args x0(`varlist') x1(`x1') xl(`xl') cont(`cont') three(`three') mass(`wgt') cell(`cell') ///
				step(`step') `normalize'
			if inlist("`name'", "share", "average"){
				foreach j in 1 2{
					local p = cond(`j' == 1, `par', `par2')
					tempvar c`j'
					if inlist(`p', 0, 1) qui gen double `c`j'' = 0
					else _ddd_correction, `args' gen(`c`j'') kind(builtin) name(`name') par(`p') par2(1)
				}
				if "`name'" == "share" qui gen double `corr' = `c1' - `c2'
				else qui gen double `corr' = ((1 - `par') * `c1' - (1 - `par2') * `c2') / (`par2' - `par')
			}
			else _ddd_correction, `args' disp(`disp') gen(`corr') `stat'
		}

		/* 5: One observation per group and interval */
		local cellvars `Nall' `nuall' `NC' `nuC0' `nuC1' `drift' `mC0' `mC1'
		if `lag' > 0 local cellvars `cellvars' `NA' `corr'
		collapse (firstnm) `cellvars', by(`cell' `gid' `time' `by')
		* The interval starting at each date exists if the group is observed horizon() periods later
		qui tsset `gid' `time', delta(`delta')
		qui gen double nu1 = F`horizon'.`nuall'
		qui gen double N1 = F`horizon'.`Nall'
		qui keep if !missing(N1)
		if _N == 0{
			di as error "No pair of dates `horizon' periods apart in the data"
			exit 2000
		}
		rename (`time' `Nall' `nuall' `NC' `drift') (`time'0 N0 nu0 NC drift)
		qui gen double `time'1 = `time'0 + `horizon' * `delta'
		qui gen double total = nu1 - nu0
		qui gen double displacement = `nuC1' - `nuC0' - drift
		qui gen double entry = nu1 - `nuC1'
		qui gen double exit = `nuC0' - nu0
		qui gen double demography = entry + exit
		qui gen double entry_mass = `mC1' - `nuC1'
		qui gen double entry_selection = entry - entry_mass
		qui gen double exit_mass = `nuC0' - `mC0'
		qui gen double exit_selection = exit - exit_mass
		* With lag(), drift and displacement use lagged weights, and the intervals without continuing units
		* observed at the lagged date, where the correction is not defined, are dropped. chi_lag is the
		* correction as a share of displacement with weights at the start.
		local lagvars
		if `lag' > 0{
			rename `NA' NC_lag
			qui keep if NC_lag > 0
			if _N == 0{
				di as error "No interval with continuing units observed `lag' before its start"
				exit 2000
			}
			qui gen double chi_lag = `corr' / displacement
			qui replace drift = drift + `corr'
			qui replace displacement = displacement - `corr'
			local lagvars NC_lag chi_lag
		}
		keep `by' `time'0 `time'1 nu0 nu1 total drift displacement demography entry exit entry_mass ///
			entry_selection exit_mass exit_selection N0 N1 NC `lagvars'
		order `by' `time'0 `time'1 nu0 nu1 total drift displacement demography entry exit entry_mass ///
			entry_selection exit_mass exit_selection N0 N1 NC `lagvars'
		format `time'0 `time'1 `tfmt'
		label var `time'0 "Start of the interval"
		label var `time'1 "End of the interval"
		label var N0 "Mass of units at the start"
		label var N1 "Mass of units at the end"
		label var NC "Mass of continuing units"
		label var nu0 "Statistic at the start"
		label var nu1 "Statistic at the end"
		label var total "Change in the statistic (nu1 - nu0)"
		label var drift "Drift"
		label var displacement "Displacement"
		label var demography "Demography (= entry + exit)"
		label var entry "Entry (= entry_mass + entry_selection)"
		label var exit "Exit (= exit_mass + exit_selection)"
		label var entry_mass "Mass part of entry"
		label var entry_selection "Selection part of entry"
		label var exit_mass "Mass part of exit"
		label var exit_selection "Selection part of exit"
		if `lag' > 0{
			label var drift "Drift, with weights `lag' before the start"
			label var displacement "Displacement, with weights `lag' before the start"
			label var NC_lag "Mass of continuing units also observed `lag' before the start"
			label var chi_lag "Correction as a share of displacement with weights at the start"
		}
		sort `by' `time'0
		qui save `output'
	}
	if `"`savefile'"' != "" frame `work': save `"`savefile'"', `replace'
	else use `output', clear
end

/* Parse saving(filename[, replace]) */
program define _ddd_saving, sclass
	gettoken filename 0 : 0, parse(",")
	local filename = strtrim(`"`filename'"')
	if inlist(`"`filename'"', "", ","){
		di as error "saving() requires a filename"
		exit 198
	}
	syntax [, replace]
	sreturn local filename `"`filename'"'
	sreturn local replace `replace'
end

/* The statistic in each cell, among the units that satisfy if, with masses mass() */
program define _ddd_nu
	syntax varlist(numeric) [if], mass(varname) cell(varname) gen(name) kind(string) name(string) ///
		par(string) par2(string) [cmd(string asis) vnames(string)]
	marksample touse, novarlist
	if "`kind'" == "builtin"{
		_ddd_builtin `varlist' if `touse', mass(`mass') cell(`cell') nu(`gen') name(`name') par(`par') ///
			par2(`par2')
	}
	else _ddd_ado `varlist' if `touse', mass(`mass') cell(`cell') gen(`gen') cmd(`cmd') vnames(`vnames')
end

/* Drift in each cell: the sum, over the units that satisfy if, of each unit's change from positions from() to
   positions to(), times its weight, the derivative of the statistic with respect to its position, evaluated
   at positions at(). With normalize, each position is divided by its mean at its date. */
program define _ddd_drift
	syntax [if], at(varlist numeric) from(varlist numeric) to(varlist numeric) mass(varname) ///
		cell(varname) gen(name) kind(string) name(string) par(string) par2(string) ///
		[cmd(string asis) vnames(string) step(real 1e-5) NORMalize]
	marksample touse, novarlist
	if "`normalize'" != ""{
		tempvar N
		qui egen double `N' = total(`mass') if `touse', by(`cell')
		foreach list in at from to{
			local normalized
			foreach v of local `list'{
				tempvar mean n
				qui egen double `mean' = total(`mass' * `v' / `N') if `touse', by(`cell')
				qui gen double `n' = `v' / `mean' if `touse'
				local normalized `normalized' `n'
			}
			local `list' `normalized'
		}
	}
	* The change of each variable
	local K : word count `from'
	local D
	forvalues k = 1/`K'{
		tempvar d
		qui gen double `d' = `: word `k' of `to'' - `: word `k' of `from'' if `touse'
		local D `D' `d'
	}
	if "`kind'" == "builtin"{
		tempvar nu weight
		_ddd_builtin `at' if `touse', mass(`mass') cell(`cell') nu(`nu') deriv(`weight') name(`name') ///
			par(`par') par2(`par2')
		qui egen double `gen' = total(`weight' * `D') if `touse', by(`cell')
	}
	else{
		* A user-supplied statistic: drift is the derivative of the statistic along the changes,
		* d/ds nu(X + s D) at s = 0, from central differences with steps h and h/2 combined by Richardson
		* extrapolation
		local j 0
		foreach s in 1 -1 0.5 -0.5{
			local ++j
			local Xs
			forvalues k = 1/`K'{
				tempvar x
				qui gen double `x' = `: word `k' of `at'' + `s' * `step' * `: word `k' of `D'' if `touse'
				local Xs `Xs' `x'
			}
			tempvar nu`j'
			_ddd_ado `Xs' if `touse', mass(`mass') cell(`cell') gen(`nu`j'') cmd(`cmd') vnames(`vnames')
		}
		qui gen double `gen' = (4 * (`nu3' - `nu4') / `step' - (`nu1' - `nu2') / (2 * `step')) / 3 if `touse'
	}
end

/* Lagged-weights correction in each cell, on all its observations: drift with lagged weights minus drift,
   estimated on the units observed at the three dates (three) as a share chi of their displacement, times the
   displacement of all continuing units (cont), computed here unless given in disp(). x0, x1, and xl hold
   positions at the start, at the end, and at the lagged date. */
program define _ddd_correction
	syntax, x0(varlist) x1(varlist) xl(varlist) cont(varname) three(varname) mass(varname) cell(varname) ///
		gen(name) kind(string) name(string) par(string) par2(string) [disp(varname) cmd(string asis) ///
		vnames(string) step(real 1e-5) NORMalize]
	local stat kind(`kind') name(`name') par(`par') par2(`par2') cmd(`cmd') vnames(`vnames')
	local opts mass(`mass') cell(`cell') `stat'
	tempvar nuA0 nuA1 dA dL chi
	* On the units observed at the three dates: drift with weights at the start and with weights at the
	* lagged date
	_ddd_nu `x0' if `three', gen(`nuA0') `opts'
	_ddd_nu `x1' if `three', gen(`nuA1') `opts'
	_ddd_drift if `three', at(`x0') from(`x0') to(`x1') gen(`dA') `opts' step(`step') `normalize'
	_ddd_drift if `three', at(`xl') from(`x0') to(`x1') gen(`dL') `opts' step(`step') `normalize'
	qui gen double `chi' = (`dL' - `dA') / (`nuA1' - `nuA0' - `dA')
	* On all continuing units: displacement
	if "`disp'" == ""{
		tempvar nuC0 nuC1 dC disp
		_ddd_nu `x0' if `cont', gen(`nuC0') `opts'
		_ddd_nu `x1' if `cont', gen(`nuC1') `opts'
		_ddd_drift if `cont', at(`x0') from(`x0') to(`x1') gen(`dC') `opts' step(`step') `normalize'
		qui gen double `disp' = `nuC1' - `nuC0' - `dC'
	}
	* chi is defined on the units observed at the three dates, displacement on all continuing units
	tempvar c d
	qui egen double `c' = max(`chi'), by(`cell')
	qui egen double `d' = max(`disp'), by(`cell')
	qui gen double `gen' = `c' * `d'
end

/* Built-in statistics: the value nu() of the statistic in each cell, among the units with positions x that
   satisfy if, with masses m = mass(), and with deriv() each unit's weight, the derivative of nu with respect
   to the unit's position. N is the mass of the cell, so that a sum over units divided by N is an average. */
program define _ddd_builtin
	syntax varname(numeric) [if], mass(varname) cell(varname) nu(name) name(string) par(string) ///
		par2(string) [deriv(name)]
	marksample touse, novarlist
	local x `varlist'
	local m `mass'
	local a `par'
	tempvar N mu c S lM M X A om
	qui egen double `N' = total(`m') if `touse', by(`cell')
	if "`name'" == "mean"{
		qui egen double `nu' = total(`m' * `x' / `N') if `touse', by(`cell')
		local weight `m' / `N'
	}
	else if "`name'" == "var"{
		qui egen double `mu' = total(`m' * `x' / `N') if `touse', by(`cell')
		qui egen double `nu' = total(`m' * (`x' - `mu')^2 / `N') if `touse', by(`cell')
		local weight 2 * `m' * (`x' - `mu') / `N'
	}
	else if inlist("`name'", "powermean", "powersum"){
		* Computed in logs to avoid overflow: lM is the log of the statistic. powersum sums over units rather
		* than averaging across them.
		local n = cond("`name'" == "powermean", "`N'", "1")
		if `a' == 0 qui egen double `lM' = total(`m' * ln(`x') / `N') if `touse', by(`cell')
		else{
			qui egen double `c' = max(`a' * ln(`x')) if `touse', by(`cell')
			qui egen double `S' = total(`m' * exp(`a' * ln(`x') - `c')) if `touse', by(`cell')
			qui gen double `lM' = (`c' + ln(`S' / `n')) / `a' if `touse'
		}
		qui gen double `nu' = exp(`lM') if `touse'
		local weight `m' * exp((`a' - 1) * (ln(`x') - `lM')) / `n'
	}
	else if inlist("`name'", "expmean", "expsum"){
		* Computed relative to the largest exp(g x) to avoid overflow. expsum sums over units rather than
		* averaging across them.
		local n = cond("`name'" == "expmean", "`N'", "1")
		qui egen double `c' = max(`a' * `x') if `touse', by(`cell')
		qui egen double `S' = total(`m' * exp(`a' * `x' - `c')) if `touse', by(`cell')
		qui gen double `nu' = (`c' + ln(`S' / `n')) / `a' if `touse'
		local weight `m' * exp(`a' * `x' - `c') / `S'
	}
	else if "`name'" == "ge"{
		qui egen double `mu' = total(`m' * `x' / `N') if `touse', by(`cell')
		if `a' == 1{
			qui egen double `nu' = total(`m' * `x' / `mu' * ln(`x' / `mu') / `N') if `touse', by(`cell')
			local weight `m' / (`N' * `mu') * (ln(`x' / `mu') - `nu')
		}
		else if `a' == 0{
			qui egen double `nu' = total(-`m' * ln(`x' / `mu') / `N') if `touse', by(`cell')
			local weight `m' / `N' * (1 / `mu' - 1 / `x')
		}
		else{
			qui egen double `A' = total(`m' * (`x' / `mu')^(`a') / `N') if `touse', by(`cell')
			qui gen double `nu' = (`A' - 1) / (`a' * (`a' - 1)) if `touse'
			local weight `m' / (`N' * (`a' - 1) * `mu') * ((`x' / `mu')^(`a' - 1) - `A')
		}
	}
	else if "`name'" == "atkinson"{
		* M is the power mean of order 1 - e
		qui egen double `mu' = total(`m' * `x' / `N') if `touse', by(`cell')
		if `a' == 1{
			qui egen double `M' = total(`m' * ln(`x') / `N') if `touse', by(`cell')
			qui replace `M' = exp(`M')
		}
		else{
			qui egen double `M' = total(`m' * `x'^(1 - `a') / `N') if `touse', by(`cell')
			qui replace `M' = `M'^(1 / (1 - `a'))
		}
		qui gen double `nu' = 1 - `M' / `mu' if `touse'
		local weight `m' / (`N' * `mu') * (`M' / `mu' - `M'^(`a') * `x'^(-`a'))
	}
	else if inlist("`name'", "gini", "share", "average"){
		* Rank-weighted statistics: om is each unit's rank weight
		_ddd_rankweight `x' if `touse', mass(`m') cell(`cell') gen(`om') name(`name') lo(`a') hi(`par2')
		if "`name'" == "average"{
			qui egen double `nu' = total(`m' * `om' * `x' / (`N' * (`par2' - `a'))) if `touse', by(`cell')
			local weight `m' * `om' / (`N' * (`par2' - `a'))
		}
		else{
			qui egen double `X' = total(`m' * `x') if `touse', by(`cell')
			qui egen double `nu' = total(`m' * `om' * `x' / `X') if `touse', by(`cell')
			local weight `m' / `X' * (`om' - `nu')
		}
	}
	else if "`name'" == "hhi"{
		qui egen double `X' = total(`m' * `x') if `touse', by(`cell')
		qui egen double `nu' = total(`m' * (`x' / `X')^2) if `touse', by(`cell')
		local weight 2 * `m' / `X' * (`x' / `X' - `nu')
	}
	if "`deriv'" != "" qui gen double `deriv' = `weight' if `touse'
end

/* Rank weight of each unit: the average of the rank-weight function over the percentile interval [a, b]
   that the unit occupies, with tied units sharing their combined interval. For the Gini, the function is
   2R - 1; for share lo hi and average lo hi, it is one between the ranks lo and hi and zero elsewhere.
   Sorting by cell first keeps the data sorted by cell, so that egen does not sort them again. */
program define _ddd_rankweight
	syntax varname [if], mass(varname) cell(varname) gen(name) name(string) lo(string) hi(string)
	marksample touse, novarlist
	tempvar N cum a b
	qui egen double `N' = total(`mass') if `touse', by(`cell')
	sort `cell' `touse' `varlist'
	qui by `cell' `touse': gen double `cum' = sum(`mass') / `N'
	qui by `cell' `touse' `varlist': gen double `b' = `cum'[_N]
	qui by `cell' `touse' `varlist': gen double `a' = `cum'[1] - `mass'[1] / `N'
	if "`name'" == "gini" qui gen double `gen' = `a' + `b' - 1 if `touse'
	else qui gen double `gen' = max(0, min(`b', `hi') - max(`a', `lo')) / (`b' - `a') if `touse'
end

/* A user-supplied statistic in each cell: the program runs on a frame that holds exactly the units of the
   cell that satisfy if, with the user's variable names and their masses in _ddd_w */
program define _ddd_ado
	syntax varlist(numeric) [if], mass(varname) cell(varname) gen(name) cmd(string asis) vnames(string)
	marksample touse, novarlist
	tempvar obs range
	tempname f nu
	qui gen double `gen' = .
	* The units of each cell that satisfy if are in a block of observations, from obs[1] to obs[_N]
	sort `cell' `touse'
	qui gen long `obs' = _n
	qui by `cell' `touse': gen `range' = string(`obs'[1], "%12.0f") + "/" + string(`obs'[_N], "%12.0f") ///
		if _n == 1 & `touse'
	qui levelsof `range', local(ranges) clean
	foreach r of local ranges{
		frame put `varlist' `mass' in `r', into(`f')
		if "`varlist'" != "`vnames'" frame `f': rename (`varlist') (`vnames')
		frame `f': rename `mass' _ddd_w
		* Clear r() before calling the program (count is r-class)
		qui count in 1
		frame `f': cap noi `cmd'
		local rc = _rc
		frame drop `f'
		if `rc'{
			di as error "The program computing the statistic exited with an error: `cmd'"
			exit `rc'
		}
		if "`r(nu)'" == ""{
			di as error "The program must return the statistic in r(nu): `cmd'"
			exit 498
		}
		scalar `nu' = r(nu)
		qui replace `gen' = scalar(`nu') in `r'
	}
end
