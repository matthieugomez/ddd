* version 0.1.0, October 2026
* Decompose the change in a distributional index into drift, displacement, and demography
* (Gomez, "The Micro Origins of Aggregate Trends: Drift, Displacement, and Demography").
*
* Drift is the change in the index with each continuing unit's marginal weight held at its initial value,
* displacement is the rest of the change among continuing units, and demography is the difference between
* the change in the full population and the change among continuing units. Built-in indices use their
* analytic weights. A user-supplied index (an ado program) only needs to compute
* the index: drift is then the directional derivative of the index along continuing units' changes,
* computed by a central difference with Richardson extrapolation.

program define ddd
	version 16
	syntax varlist(numeric) [if] [in] [fweight iweight pweight], Index(string asis) ///
		[saving(string asis) clear by(varlist) HORizon(integer 1) lag(integer 0) ///
		NORMalize step(real 1e-5)]

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

	/* 1: Parse the index */
	* kind is builtin (analytic weights) or ado (an rclass program)
	local index `index'
	local par .
	local cmd
	gettoken name rest : index
	local builtins mean var power expmean ces ge theil mld atkinson gini top hhi
	if `: list name in builtins'{
		local kind builtin
		if inlist("`name'", "power", "expmean", "ces", "ge", "atkinson", "top"){
			local rest = trim(`"`rest'"')
			cap confirm number `rest'
			if _rc{
				di as error "index(`name' #) requires a numeric parameter"
				exit 198
			}
			local par = `rest'
			if "`name'" == "top" & !(`par' > 0 & `par' < 1){
				di as error "index(top #) requires a share between 0 and 1, e.g. index(top 0.01)"
				exit 198
			}
			if inlist("`name'", "expmean", "ces") & `par' == 0{
				di as error "index(`name' #) requires a nonzero parameter"
				exit 198
			}
		}
		else if trim(`"`rest'"') != ""{
			di as error "index(`name') takes no parameter"
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
			di as error "Built-in indices take a single variable"
			exit 198
		}
	}
	else{
		local kind ado
		cap program list `name'
		if _rc{
			cap findfile `name'.ado
			if _rc{
				di as error "`name' is neither a built-in index nor a program (built-in indices: `builtins')"
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
		local positive = inlist("`name'", "ge", "atkinson") | ("`name'" == "power" & `par' != round(`par'))
		if `positive'{
			cap assert `varlist' > 0 if `touse'
			if _rc{
				di as error "index(`index') requires `varlist' > 0"
				exit 198
			}
		}
		if "`name'" == "hhi"{
			cap assert `varlist' >= 0 if `touse'
			if _rc{
				di as error "index(hhi) requires `varlist' >= 0"
				exit 198
			}
		}
	}
	if mod(`lag', `delta') != 0{
		di as error "lag() must be a multiple of the panel's time delta (`delta')"
		exit 198
	}

	local tfmt : format `time'

	/* 3: Compute the decomposition on a copy of the variables used */
	tempvar gid iid
	tempname work ef out res
	tempfile map output
	frame put `by' `id' `time' `varlist' `wgt' if `touse', into(`work')
	frame create `ef'
	frame `work'{
		if "`by'" != "" qui egen long `gid' = group(`by'), missing
		else qui gen long `gid' = 1
		qui egen long `iid' = group(`id')
		sort `gid' `time' `iid'
		if "`by'" != ""{
			preserve
			qui keep `gid' `by'
			qui duplicates drop `gid', force
			qui save `map'
			restore
		}
		mata: _ddd_run("`varlist'", "`wgt'", "`gid'", "`iid'", "`time'", `delta', ///
			`horizon', `lag', "`normalize'" != "", "`kind'", "`name'", `par', ///
			`"`cmd'"', "`ef'", `step', "`res'")
	}

	/* 4: Output dataset */
	if rowsof(`res') == 0 | `res'[1, 1] == .{
		di as error "No pair of dates `horizon' periods apart in the data"
		exit 2000
	}
	frame create `out'
	frame `out'{
		qui svmat double `res'
		local names _ddd_gid `time'0 `time'1 N0 N1 NC nu0 nu1 total drift displacement ///
			demography entry exit entry_mass entry_selection exit_mass exit_selection
		if `lag' > 0 local names `names' NC_lag drift_lag displacement_lag
		local j 0
		foreach v of local names{
			local ++j
			rename `res'`j' `v'
		}
		order N0 N1 NC, after(exit_selection)
		if "`by'" != ""{
			rename _ddd_gid `gid'
			qui merge m:1 `gid' using `map', keep(master matched) nogen
			drop `gid'
			order `by'
		}
		else drop _ddd_gid
		format `time'0 `time'1 `tfmt'
		label var `time'0 "Start of the interval"
		label var `time'1 "End of the interval"
		label var N0 "Mass of units at the start"
		label var N1 "Mass of units at the end"
		label var NC "Mass of continuing units"
		label var nu0 "Index at the start"
		label var nu1 "Index at the end"
		label var total "Change in the index (nu1 - nu0)"
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
			label var NC_lag "Mass of continuing units also observed `lag' before the start"
			label var drift_lag "Drift with lagged weights"
			label var displacement_lag "Displacement with lagged weights"
		}
		sort `by' `time'0 `time'1
		qui save `output'
	}
	if `"`savefile'"' != ""{
		frame `out': save `"`savefile'"', `replace'
	}
	else{
		use `output', clear
	}
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


mata:
mata set matastrict on

struct _ddd_index{
	string scalar kind
	string scalar name
	real scalar par
	string scalar cmd
	string rowvector vnames
	string scalar evalframe
	real scalar step
}

/* Driver: loop over groups and start dates */
void _ddd_run(string scalar vars, string scalar wvar, string scalar gvar, string scalar ivar,
	string scalar tvar, real scalar delta, real scalar h, real scalar lag, real scalar norm,
	string scalar kind, string scalar name, real scalar par, string scalar cmd,
	string scalar ef, real scalar step, string scalar outmat)
{
	struct _ddd_index scalar ix
	real matrix X, ginfo, tinfo, R
	real colvector W, G, I, T, rg, Tg, times, r0, r1, rl, j1, jl
	real rowvector res
	real scalar g, j, M, t1

	ix.kind = kind
	ix.name = name
	ix.par = par
	ix.cmd = cmd
	ix.vnames = tokens(vars)
	ix.evalframe = ef
	ix.step = step

	X = st_data(., vars)
	W = st_data(., wvar)
	G = st_data(., gvar)
	I = st_data(., ivar)
	T = st_data(., tvar)
	M = max(I)

	R = J(0, 18 + 3 * (lag > 0), .)
	ginfo = panelsetup(G, 1)
	for (g = 1; g <= rows(ginfo); g++){
		rg = ginfo[g, 1]::ginfo[g, 2]
		Tg = T[rg]
		tinfo = panelsetup(Tg, 1)
		times = Tg[tinfo[., 1]]
		for (j = 1; j <= rows(times); j++){
			r0 = rg[tinfo[j, 1]::tinfo[j, 2]]
			rl = J(0, 1, .)
			if (lag > 0){
				jl = selectindex(times :== times[j] - lag)
				if (rows(jl) == 1) rl = rg[tinfo[jl, 1]::tinfo[jl, 2]]
			}
			t1 = times[j] + h * delta
			j1 = selectindex(times :== t1)
			if (rows(j1) != 1) continue
			r1 = rg[tinfo[j1, 1]::tinfo[j1, 2]]
			res = _ddd_pair(ix, X, W, I, r0, r1, rl, lag > 0, norm, M)
			R = R \ (G[rg[1]], times[j], t1, res)
		}
	}
	st_matrix(outmat, R)
}

/* Decomposition between two dates */
real rowvector _ddd_pair(struct _ddd_index scalar ix, real matrix X, real colvector W,
	real colvector I, real colvector r0, real colvector r1, real colvector rl,
	real scalar dolag, real scalar norm, real scalar M)
{
	real colvector pos, c1, keep, rc0, rc1, WC, cl, k3, WA
	real matrix C0, C1, A0, A1, AL
	real rowvector res
	real scalar N0, N1, NC, nu0, nu1, total, nuC0, nuC1, mC0, mC1, drift, disp, entry, exit,
		entrymass, exitmass, NA, dA, dispA, dL, share, driftlag, displag

	N0 = quadsum(W[r0])
	N1 = quadsum(W[r1])
	nu0 = _ddd_nu(ix, X[r0, .], W[r0])
	nu1 = _ddd_nu(ix, X[r1, .], W[r1])
	total = nu1 - nu0

	// Continuing units: units observed at both dates; they keep their mass at the start date
	pos = J(M, 1, .)
	pos[I[r1]] = r1
	c1 = pos[I[r0]]
	keep = c1 :< .
	if (!any(keep)){
		return((N0, N1, 0, nu0, nu1, total, ., ., ., ., ., ., ., ., .) , J(1, 3 * dolag, .))
	}
	rc0 = select(r0, keep)
	rc1 = select(c1, keep)
	WC = W[rc0]
	NC = quadsum(WC)
	C0 = X[rc0, .]
	C1 = X[rc1, .]
	nuC0 = _ddd_nu(ix, C0, WC)
	nuC1 = _ddd_nu(ix, C1, WC)

	drift = _ddd_drift(ix, C0, C0, C1, WC, norm)
	disp = nuC1 - nuC0 - drift
	entry = nu1 - nuC1
	exit = nuC0 - nu0
	// Mass parts of entry and exit: entrants and exiters that replicate the continuing units
	mC0 = _ddd_nu(ix, C0, WC * (N0 / NC))
	mC1 = _ddd_nu(ix, C1, WC * (N1 / NC))
	entrymass = mC1 - nuC1
	exitmass = nuC0 - mC0

	// Lagged weights: estimated on the continuing units also observed at the lagged date, as a share of
	// their displacement, and applied to all continuing units
	driftlag = .
	displag = .
	NA = 0
	if (dolag & rows(rl) > 0){
		pos = J(M, 1, .)
		pos[I[rl]] = rl
		cl = pos[I[rc0]]
		k3 = cl :< .
		if (any(k3)){
			A0 = select(C0, k3)
			A1 = select(C1, k3)
			AL = X[select(cl, k3), .]
			WA = select(WC, k3)
			NA = quadsum(WA)
			dA = _ddd_drift(ix, A0, A0, A1, WA, norm)
			dispA = _ddd_nu(ix, A1, WA) - _ddd_nu(ix, A0, WA) - dA
			dL = _ddd_drift(ix, AL, A0, A1, WA, norm)
			share = (dA - dL) / dispA
			driftlag = drift - share * disp
			displag = disp + share * disp
		}
	}

	res = (N0, N1, NC, nu0, nu1, total, drift, disp, entry + exit, entry, exit, entrymass,
		entry - entrymass, exitmass, exit - exitmass)
	if (dolag) res = res, (NA, driftlag, displag)
	return(res)
}

/* Drift with weights evaluated at Xw and changes Xb - Xa; with norm, positions are divided by
   their mean at each date */
real scalar _ddd_drift(struct _ddd_index scalar ix, real matrix Xw, real matrix Xa, real matrix Xb,
	real colvector w, real scalar norm)
{
	real scalar N
	N = quadsum(w)
	if (norm){
		return(_ddd_dd(ix, Xw :/ (quadcolsum(w :* Xw) / N), w,
			Xb :/ (quadcolsum(w :* Xb) / N) - Xa :/ (quadcolsum(w :* Xa) / N)))
	}
	return(_ddd_dd(ix, Xw, w, Xb - Xa))
}

/* Directional derivative of the index at X along D */
real scalar _ddd_dd(struct _ddd_index scalar ix, real matrix X, real colvector w, real matrix D)
{
	real colvector g
	real scalar h, d1, d2
	if (ix.kind == "builtin"){
		(void) _ddd_builtin(ix.name, ix.par, X, w, 1, g)
		return(quadsum(g :* D))
	}
	h = ix.step
	d1 = (_ddd_nu(ix, X + h * D, w) - _ddd_nu(ix, X - h * D, w)) / (2 * h)
	d2 = (_ddd_nu(ix, X + (h / 2) * D, w) - _ddd_nu(ix, X - (h / 2) * D, w)) / h
	return((4 * d2 - d1) / 3)
}

/* Value of the index */
real scalar _ddd_nu(struct _ddd_index scalar ix, real matrix X, real colvector w)
{
	real colvector g
	if (ix.kind == "builtin") return(_ddd_builtin(ix.name, ix.par, X, w, 0, g))
	return(_ddd_ado(ix, X, w))
}

/* Index computed by an rclass ado program, run on a frame that holds the units */
real scalar _ddd_ado(struct _ddd_index scalar ix, real matrix X, real colvector w)
{
	string scalar cur
	string rowvector names
	real scalar rc
	real matrix nu

	names = (ix.vnames, "_ddd_w")
	cur = st_framecurrent()
	st_framecurrent(ix.evalframe)
	stata("quietly clear")
	st_addobs(rows(X))
	(void) st_addvar("double", names)
	st_store(., names, (X, w))
	st_rclear()
	rc = _stata(ix.cmd, 1)
	nu = st_numscalar("r(nu)")
	st_framecurrent(cur)
	if (rc){
		errprintf("The index program exited with an error: %s\n", ix.cmd)
		exit(rc)
	}
	if (rows(nu) == 0){
		errprintf("The index program must return the index in r(nu): %s\n", ix.cmd)
		exit(498)
	}
	return(nu)
}

/* Built-in indices: value, and with wantg the derivative with respect to each unit's position */
real scalar _ddd_builtin(string scalar name, real scalar a, real colvector x, real colvector w,
	real scalar wantg, real colvector g)
{
	real scalar N, m, nu, A, c, S, Mn, X
	real colvector xt, e, om

	N = quadsum(w)
	if (name == "mean"){
		nu = quadsum(w :* x) / N
		if (wantg) g = w / N
	}
	else if (name == "var"){
		m = quadsum(w :* x) / N
		nu = quadsum(w :* (x :- m):^2) / N
		if (wantg) g = 2 * w :* (x :- m) / N
	}
	else if (name == "power"){
		nu = quadsum(w :* x:^a) / N
		if (wantg) g = a * w :* x:^(a - 1) / N
	}
	else if (name == "expmean" | name == "ces"){
		// ces sums over units rather than averaging across them
		c = max(a * x)
		e = exp(a * x :- c)
		S = quadsum(w :* e)
		nu = (c + ln(name == "ces" ? S : S / N)) / a
		if (wantg) g = w :* e / S
	}
	else if (name == "ge"){
		m = quadsum(w :* x) / N
		xt = x / m
		if (a == 1){
			nu = quadsum(w :* xt :* ln(xt)) / N
			if (wantg) g = w / (N * m) :* (ln(xt) :- nu)
		}
		else if (a == 0){
			nu = -quadsum(w :* ln(xt)) / N
			if (wantg) g = w / N :* (1 / m :- 1 :/ x)
		}
		else{
			A = quadsum(w :* xt:^a) / N
			nu = (A - 1) / (a * (a - 1))
			if (wantg) g = w / (N * (a - 1) * m) :* (xt:^(a - 1) :- A)
		}
	}
	else if (name == "atkinson"){
		m = quadsum(w :* x) / N
		if (a == 1) Mn = exp(quadsum(w :* ln(x)) / N)
		else Mn = (quadsum(w :* x:^(1 - a)) / N)^(1 / (1 - a))
		nu = 1 - Mn / m
		if (wantg) g = w / (N * m) :* (Mn / m :- Mn^a * x:^(-a))
	}
	else if (name == "gini" | name == "top"){
		om = _ddd_rankweight(x, w, name, a)
		X = quadsum(w :* x)
		nu = quadsum(w :* om :* x) / X
		if (wantg) g = w / X :* (om :- nu)
	}
	else if (name == "hhi"){
		X = quadsum(w :* x)
		nu = quadsum(w :* x:^2) / X^2
		if (wantg) g = 2 * w / X :* (x / X :- nu)
	}
	return(nu)
}

/* Rank weight of each unit: the average of the rank-weight function over the percentile interval
   that the unit occupies, with tied units sharing their combined interval */
real colvector _ddd_rankweight(real colvector x, real colvector w, string scalar name, real scalar p)
{
	real colvector o, xs, ws, cum, gid, a, b, om, res
	real matrix info
	real scalar N, n

	n = rows(x)
	N = quadsum(w)
	o = order(x, 1)
	xs = x[o]
	ws = w[o]
	cum = quadrunningsum(ws) / N
	if (n == 1) gid = 1
	else gid = runningsum(1 \ (xs[|2 \ n|] :!= xs[|1 \ n - 1|]))
	info = panelsetup(gid, 1)
	a = cum[info[., 1]] - ws[info[., 1]] / N
	b = cum[info[., 2]]
	if (name == "gini") om = a + b :- 1
	else om = rowmax((J(rows(a), 1, 0), b - rowmax((a, J(rows(a), 1, 1 - p))))) :/ (b - a)
	res = J(n, 1, .)
	res[o] = om[gid]
	return(res)
}
end
