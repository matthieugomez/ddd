{smcl}
{* *! version 0.1.0 October 2026}{...}
{viewerjumpto "Syntax" "ddd##syntax"}{...}
{viewerjumpto "Description" "ddd##description"}{...}
{viewerjumpto "Indices" "ddd##indices"}{...}
{viewerjumpto "User-defined indices" "ddd##user"}{...}
{viewerjumpto "Options" "ddd##options"}{...}
{viewerjumpto "Output" "ddd##output"}{...}
{viewerjumpto "Examples" "ddd##examples"}{...}
{viewerjumpto "Author" "ddd##contact"}{...}


{title:Title}

{p2colset 4 12 12 8}{...}
{p2col :{cmd:ddd} {hline 2}}Decompose the change in a distributional index into drift, displacement, and demography{p_end}
{p2colreset}{...}


{marker syntax}{...}
{title:Syntax}

{p 8 15 2} {cmd:ddd} {varlist} {ifin} [{it:{help weight}}] {cmd:,} {cmd:index(}{it:index}{cmd:)}
{{cmd:saving(}{it:filename}[{cmd:, replace}]{cmd:)} | {cmd:clear}}
[{cmd:by(}{it:varlist}{cmd:)} {cmd:horizon(}{it:#}{cmd:)} {cmd:lag(}{it:#}{cmd:)} {cmd:normalize} {cmd:log} {cmd:step(}{it:#}{cmd:)}]{p_end}

{pstd}
{it:index} is a built-in index ({help ddd##indices:below}), the name of an rclass program,
or {cmd:mata:} {it:fcn}{cmd:()}, a Mata function ({help ddd##user:below}).
{cmd:fweight}s, {cmd:iweight}s, and {cmd:pweight}s are allowed; all are treated as masses.
The dataset must be declared as a panel ({help tsset}).


{marker description}{...}
{title:Description}

{pstd}
Productivity, concentration, inequality, and other distributional indices summarize a distribution of units,
so their changes come from how units move, enter, and exit.
Between two dates {it:t} and {it:t}+{it:h}, {cmd:ddd} splits the change in the index into

{p2colset 8 24 26 2}{...}
{p2col :{cmd:drift}}the sum over survivors of each survivor's change times the derivative of the index with respect to its position, evaluated at the initial distribution of survivors: the change predicted by holding each unit's marginal weight at its initial value{p_end}
{p2col :{cmd:displacement}}the rest of the change in the index among survivors: the effect of the weights changing as units move{p_end}
{p2col :{cmd:demography}}the change in the full population minus the change among survivors, the sum of {cmd:entry} (adding the entrants to the survivors at {it:t}+{it:h}) and {cmd:exit} (removing the exiters from the population at {it:t}){p_end}
{p2colreset}{...}

{pstd}
Survivors are the units observed at both dates. The decomposition is exact:
{cmd:total} = {cmd:drift} + {cmd:displacement} + {cmd:demography}.
Demography splits in turn into {cmd:mass}, the effect of entry and exit that replicate the survivors,
and {cmd:selection}, the rest, which is zero when entrants and exiters are distributed like survivors.

{pstd}
A unit is in the population at a date when it has an observation with nonmissing values of {varlist}.
With weights, each unit counts with its weight at each date, and survivors keep their weight at the start date at both dates,
so that reweighting of survivors is part of demography.


{marker indices}{...}
{title:Built-in indices}

{pstd}
Built-in indices take one variable {it:x} and use the analytic derivative of the index.
E[.] is the average over units, weighted by their mass, and R is the percentile rank (ties share the average of their percentile interval).

{p2colset 8 26 28 2}{...}
{p2col :{cmd:mean}}E[x]{p_end}
{p2col :{cmd:var}}Var(x){p_end}
{p2col :{cmd:power} {it:a}}E[x^a]{p_end}
{p2col :{cmd:expmean} {it:g}}ln(E[exp(g x)]) / g{p_end}
{p2col :{cmd:ces} {it:g}}ln(N E[exp(g x)]) / g, with N the mass of units, e.g. industry productivity with log productivity x and g = elasticity - 1;
it is {cmd:expmean} plus ln(N) / g, so it has the same drift and displacement, and its mass term is the taste for variety{p_end}
{p2col :{cmd:ge} {it:a}}generalized entropy (E[(x/E[x])^a] - 1) / (a (a - 1)){p_end}
{p2col :{cmd:theil}}{cmd:ge 1}, E[(x/E[x]) ln(x/E[x])]{p_end}
{p2col :{cmd:mld}}{cmd:ge 0}, mean log deviation -E[ln(x/E[x])]{p_end}
{p2col :{cmd:atkinson} {it:e}}1 - E[x^(1-e)]^(1/(1-e)) / E[x]{p_end}
{p2col :{cmd:gini}}E[(2R - 1) x] / E[x]{p_end}
{p2col :{cmd:top} {it:p}}share of the top fraction {it:p} of units, E[1(R >= 1 - p) x] / E[x], e.g. {cmd:top 0.01}{p_end}
{p2col :{cmd:hhi}}Herfindahl-Hirschman index, the sum of squared shares of x; it depends on the number of units{p_end}
{p2colreset}{...}

{pstd}
{cmd:ge}, {cmd:theil}, {cmd:mld}, {cmd:atkinson}, and {cmd:power} with a non-integer exponent require x > 0; {cmd:hhi} requires x >= 0.


{marker user}{...}
{title:User-defined indices}

{pstd}
Any other index can be supplied as a program that computes it. Drift is then the derivative of the index along survivors' changes,
d/ds nu(x + s dx) at s = 0, computed by a central difference with step {cmd:step()} and Richardson extrapolation.
It takes four evaluations of the index, whatever the number of units.
Several variables per unit are allowed (e.g. productivity and sales); drift then sums the changes of all variables.

{pstd}
{bf:Ado program.} {cmd:index(}{it:myprog} [{it:options}]{cmd:)} calls

{p 8 8 2}{cmd:myprog} {varlist}{cmd:, w(_ddd_w)} [{it:options}]{p_end}

{pstd}
on a dataset that contains exactly the units at which the index is evaluated, with their mass in the variable {cmd:_ddd_w}.
The program must be rclass and return the index in {cmd:r(nu)}. Two rules:

{phang2}1. Use {cmd:double}s for every variable you create ({cmd:gen double}). The numerical derivative needs about ten significant digits, and a {cmd:float} has seven.{p_end}
{phang2}2. Treat {cmd:_ddd_w} as masses, not probabilities. {cmd:[iweight=_ddd_w]} does this in {cmd:summarize, meanonly};
{cmd:aweight}s rescale the weights and {cmd:fweight}s require integers.{p_end}

{pstd}
{bf:Mata function.} {cmd:index(mata: }{it:fcn}{cmd:())} calls {it:fcn}{cmd:(X, w)}, where {it:fcn} is declared as

{p 8 8 2}{cmd:real scalar} {it:fcn}{cmd:(real matrix X, real colvector w)}{p_end}

{pstd}
with one column of X per variable of {varlist} and w the masses. A Mata function is much faster than an ado program.

{pstd}
{bf:Indices with kinks.} For an index whose derivative jumps, such as a rank-based index, the numerical derivative is exact
as long as no two units cross within {cmd:step()} times their changes; the default step makes this rare, and the built-ins avoid the issue.
Indices whose influence function jumps, such as quantiles, have no drift in panel data.


{marker options}{...}
{title:Options}

{phang}
{cmd:index(}{it:index}{cmd:)} the index to decompose; required.

{phang}
{cmd:saving(}{it:filename}[{cmd:, replace}]{cmd:)} saves the decomposition to a dataset; {cmd:clear} replaces the dataset in memory with it.
One of the two is required.

{phang}
{cmd:by(}{it:varlist}{cmd:)} computes the decomposition separately within each group. A unit that changes group exits the first and enters the second.

{phang}
{cmd:horizon(}{it:#}{cmd:)} length of the intervals, in periods of the panel; default {cmd:horizon(1)}.
A term of first order in time grows in proportion to the interval: to compare horizons, run {cmd:ddd} once per horizon.

{phang}
{cmd:lag(}{it:#}{cmd:)} also reports drift and displacement with lagged weights: each survivor's weight is evaluated at its position {it:#} time units before the start of the interval,
which removes the bias from transitory measurement error and from flows (such as annual sales).
The correction is estimated on the survivors also observed at the lagged date, as a share of their displacement, and applied to all survivors.
{it:#} must be a multiple of the panel's delta.

{phang}
{cmd:normalize} divides positions by the survivors' mean at each date before splitting the survivor change into drift and displacement.
For a scale-invariant index (Gini, top shares, generalized entropy, Atkinson, HHI), the index, the survivor change, and demography are unchanged,
a common growth rate is not drift, and drift is the covariance between a weight and a change in normalized position.

{phang}
{cmd:log} decomposes the log change in a positive index, ln {cmd:nu1} - ln {cmd:nu0}, by taking logs of the split in levels.
Among survivors, the change passes through the index predicted by the initial weights, the survivors' initial index plus {cmd:drift}:
drift is the log growth from the survivors' initial index to this prediction, and displacement the log growth from the prediction to their final index.
{cmd:entry}, {cmd:exit}, {cmd:mass}, and {cmd:selection} are the corresponding log changes, so the identities hold in logs.
A common proportional change of positions is all drift for an index homogeneous of degree one.
To first order, each term is its counterpart in levels divided by the index.
Terms are missing when an index they involve is not positive.

{phang}
{cmd:step(}{it:#}{cmd:)} step of the numerical derivative for user-defined indices; default {cmd:step(1e-5)}.


{marker output}{...}
{title:Output}

{pstd}
One observation per group and interval:

{p2colset 8 30 32 2}{...}
{p2col :{it:byvars}}the {cmd:by()} variables{p_end}
{p2col :{it:timevar}{cmd:0}, {it:timevar}{cmd:1}}start and end of the interval{p_end}
{p2col :{cmd:N0}, {cmd:N1}, {cmd:NS}}mass of units at the start, at the end, and of survivors{p_end}
{p2col :{cmd:nu0}, {cmd:nu1}}index at the start and at the end{p_end}
{p2col :{cmd:total}}{cmd:nu1} - {cmd:nu0}{p_end}
{p2col :{cmd:drift}, {cmd:displacement}}survivor terms{p_end}
{p2col :{cmd:demography}}{cmd:entry} + {cmd:exit} = {cmd:mass} + {cmd:selection}{p_end}
{p2col :{cmd:NS_lag}, {cmd:drift_lag}, {cmd:displacement_lag}}with {cmd:lag()}: mass of survivors also observed at the lagged date, and the survivor terms with lagged weights{p_end}
{p2colreset}{...}

{pstd}
With {cmd:log}, {cmd:total} and the terms that follow it are in logs; {cmd:nu0} and {cmd:nu1} remain the index.


{marker examples}{...}
{title:Examples}

{pstd}Gini of wealth in a household panel, in normalized wealth, with survey weights and lagged weights{p_end}
{phang2}{cmd:. tsset hid year}{p_end}
{phang2}{cmd:. ddd wealth [pw = hw], index(gini) normalize lag(1) clear}{p_end}

{pstd}Concentration of sales within each market, over intervals of one to four years{p_end}
{phang2}{cmd:. tsset firm year}{p_end}
{phang2}{cmd:. forvalues h = 1/4 {c -(}}{p_end}
{phang2}{cmd:  2.     ddd sales, index(hhi) by(market) horizon(`h') saving(hhi_ddd`h', replace)}{p_end}
{phang2}{cmd:  3. {c )-}}{p_end}

{pstd}A user-defined index as an ado program: the coefficient of variation{p_end}
{phang2}{cmd:. program define cv, rclass}{p_end}
{phang2}{cmd:  1.     syntax varname, w(varname)}{p_end}
{phang2}{cmd:  2.     tempvar d}{p_end}
{phang2}{cmd:  3.     summarize `varlist' [iw = `w'], meanonly}{p_end}
{phang2}{cmd:  4.     local m = r(mean)}{p_end}
{phang2}{cmd:  5.     gen double `d' = (`varlist' - `m')^2}{p_end}
{phang2}{cmd:  6.     summarize `d' [iw = `w'], meanonly}{p_end}
{phang2}{cmd:  7.     return scalar nu = sqrt(r(mean)) / `m'}{p_end}
{phang2}{cmd:  8. end}{p_end}
{phang2}{cmd:. ddd wealth, index(cv) clear}{p_end}

{pstd}A user-defined index of two variables as a Mata function: average productivity weighted by sales{p_end}
{phang2}{cmd:. mata:}{p_end}
{phang2}{cmd:: real scalar wavg(real matrix X, real colvector w) return(quadsum(w :* X[., 2] :* X[., 1]) / quadsum(w :* X[., 2]))}{p_end}
{phang2}{cmd:: end}{p_end}
{phang2}{cmd:. ddd tfp sales, index(mata: wavg()) clear}{p_end}


{marker references}{...}
{title:References}

{pstd}
Gomez, Matthieu. "The Micro Origins of Aggregate Trends: Drift, Displacement, and Demography."
{browse "http://www.matthieugomez.com/files/ddd.pdf"}


{marker contact}{...}
{title:Author}

{pstd}
Matthieu Gomez{break}
Columbia University{break}
Please report issues on GitHub{break}
{browse "https://github.com/matthieugomez/ddd"}{p_end}
