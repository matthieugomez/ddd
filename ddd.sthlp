{smcl}
{* *! version 0.3.0 November 2026}{...}
{viewerjumpto "Syntax" "ddd##syntax"}{...}
{viewerjumpto "Description" "ddd##description"}{...}
{viewerjumpto "Statistics" "ddd##statistics"}{...}
{viewerjumpto "User-defined statistics" "ddd##user"}{...}
{viewerjumpto "Options" "ddd##options"}{...}
{viewerjumpto "Output" "ddd##output"}{...}
{viewerjumpto "Standard errors" "ddd##se"}{...}
{viewerjumpto "Examples" "ddd##examples"}{...}
{viewerjumpto "Author" "ddd##contact"}{...}


{title:Title}

{p2colset 4 12 12 8}{...}
{p2col :{cmd:ddd} {hline 2}}Decompose the change in a distributional statistic into drift, displacement, and demography{p_end}
{p2colreset}{...}


{marker syntax}{...}
{title:Syntax}

{p 8 15 2} {cmd:ddd} {varlist} {ifin} [{it:{help weight}}] {cmd:,} {cmd:statistic(}{it:statistic}{cmd:)}
{{cmd:saving(}{it:filename}[{cmd:, replace}]{cmd:)} | {cmd:clear}}
[{cmd:horizon(}{it:#}{cmd:)} {cmd:lag(}{it:#}{cmd:)} {cmd:normalize} {cmd:step(}{it:#}{cmd:)} {cmd:by(}{it:varlist}{cmd:)}]{p_end}

{pstd}
{it:statistic} is a built-in statistic ({help ddd##statistics:below}) or the name of an rclass program ({help ddd##user:below}).
{cmd:fweight}s, {cmd:iweight}s, and {cmd:pweight}s are allowed; all are treated as masses.
The dataset must be declared as a panel ({help tsset}).


{marker description}{...}
{title:Description}

{pstd}
Productivity, concentration, inequality, and other distributional statistics summarize a distribution of units,
so their changes come from how units move, enter, and exit.
Between two dates {it:t} and {it:t}+{it:h}, {cmd:ddd} splits the change in the statistic into

{p2colset 8 24 26 2}{...}
{p2col :{cmd:drift}}the sum over continuing units of each unit's change times the derivative of the statistic with respect to its position (the unit's weight), evaluated at the initial distribution of continuing units: the change predicted by holding each unit's weight at its initial value{p_end}
{p2col :{cmd:displacement}}the rest of the change in the statistic among continuing units: the effect of the weights changing as units move{p_end}
{p2col :{cmd:demography}}the change in the full population minus the change among continuing units, the sum of {cmd:entry} (adding the entrants to the continuing units at {it:t}+{it:h}) and {cmd:exit} (removing the exiters from the population at {it:t}){p_end}
{p2colreset}{...}

{pstd}
Continuing units are the units observed at both dates. The decomposition is exact:
{cmd:total} = {cmd:drift} + {cmd:displacement} + {cmd:demography}.
Entry and exit each split in turn into a mass part ({cmd:entry_mass}, {cmd:exit_mass}), the effect of adding entrants, or removing exiters, distributed like the continuing units,
which changes the size of the population but not the shape of its distribution,
and a selection part ({cmd:entry_selection}, {cmd:exit_selection}), the rest, the effect of where entrants and exiters actually are, which is zero when they are distributed like the continuing units.

{pstd}
A unit is in the population at a date when it has an observation with nonmissing values of {varlist}.
With weights, each unit counts with its weight at each date, and continuing units keep their weight at the start date at both dates,
so that reweighting of continuing units is part of demography.


{marker statistics}{...}
{title:Built-in statistics}

{pstd}
Built-in statistics take one variable {it:x} and use the analytic derivative of the statistic.
E[.] is the average over units, weighted by their mass, and R is the percentile rank (ties share the average of their percentile interval).

{p2colset 8 26 28 2}{...}
{p2col :{cmd:mean}}E[x]{p_end}
{p2col :{cmd:var}}Var(x){p_end}
{p2col :{cmd:powermean} {it:a}}power mean E[x^a]^(1/a); {cmd:powermean 0} is the geometric mean exp(E[ln x]) and {cmd:powermean 1} the mean{p_end}
{p2col :{cmd:powersum} {it:a}}(N E[x^a])^(1/a), the sum of x^a to the power 1/a, with N the mass of units: the CES aggregate, e.g. industry productivity with productivity x and a = elasticity of substitution - 1;
it is N^(1/a) times {cmd:powermean}, and its mass parts are the love of variety{p_end}
{p2col :{cmd:expmean} {it:g}}exponential mean ln(E[exp(g x)]) / g, the log of {cmd:powermean} {it:g} of exp(x){p_end}
{p2col :{cmd:expsum} {it:g}}ln(N E[exp(g x)]) / g, the log of {cmd:powersum} {it:g} of exp(x): the CES aggregate in logs, e.g. industry log productivity with log productivity x;
it is {cmd:expmean} plus ln(N) / g, so it has the same drift and displacement, and its mass parts are the love of variety{p_end}
{p2col :{cmd:ge} {it:a}}generalized entropy (E[(x/E[x])^a] - 1) / (a (a - 1)){p_end}
{p2col :{cmd:theil}}{cmd:ge 1}, E[(x/E[x]) ln(x/E[x])]{p_end}
{p2col :{cmd:mld}}{cmd:ge 0}, mean log deviation -E[ln(x/E[x])]{p_end}
{p2col :{cmd:atkinson} {it:e}}1 - E[x^(1-e)]^(1/(1-e)) / E[x]{p_end}
{p2col :{cmd:gini}}E[(2R - 1) x] / E[x]{p_end}
{p2col :{cmd:share} {it:lo hi}}share of the units with percentile rank between {it:lo} and {it:hi}, E[1(lo <= R < hi) x] / E[x], e.g. {cmd:share 0.99 1} for the top 1% and {cmd:share 0.4 0.5} for the fifth decile; for a band in the middle of the distribution, the weight is not monotone in rank, so displacement and demography can have either sign; decomposing every decile shows how the whole distribution changes{p_end}
{p2col :{cmd:average} {it:lo hi}}average of the units with percentile rank between {it:lo} and {it:hi}, E[x | lo <= R < hi], e.g. {cmd:average 0.9 1} for the top decile; {cmd:average 0 1} is the mean;
over percentile bins that cover all units, each component, weighted by the bins' widths {it:hi} - {it:lo}, adds up to that of the mean{p_end}
{p2col :{cmd:hhi}}Herfindahl-Hirschman index, the sum of squared shares of x; it depends on the number of units;
with {cmd:normalize}, drift is 2 sum(s ds) and displacement sum(ds^2), with s the shares of continuing units{p_end}
{p2colreset}{...}

{pstd}
{cmd:powermean}, {cmd:powersum}, {cmd:ge}, {cmd:theil}, {cmd:mld}, and {cmd:atkinson} require x > 0; {cmd:hhi} requires x >= 0.


{marker user}{...}
{title:User-defined statistics}

{pstd}
Any other statistic can be supplied as a program that computes it. Drift is then the derivative of the statistic along continuing units' changes,
d/ds nu(x + s dx) at s = 0, computed by a central difference with step {cmd:step()} and Richardson extrapolation.
It takes four evaluations of the statistic, whatever the number of units.
Several variables per unit are allowed (e.g. productivity and sales); drift then sums the changes of all variables.

{pstd}
{cmd:statistic(}{it:myprog} [{it:options}]{cmd:)} calls

{p 8 8 2}{cmd:myprog} {varlist} {cmd:[iw=_ddd_w]} [{cmd:,} {it:options}]{p_end}

{pstd}
on a dataset that contains exactly the units at which the statistic is evaluated, with their mass as an {cmd:iweight}.
The program must be rclass, declare {cmd:[iw]} in its {cmd:syntax}, and return the statistic in {cmd:r(nu)}. Two rules:

{phang2}1. Use {cmd:double}s for every variable you create ({cmd:gen double}). The numerical derivative needs about ten significant digits, and a {cmd:float} has seven.{p_end}
{phang2}2. Treat the weights as masses, not probabilities. Passing them on with {cmd:[`weight'`exp']} does this in {cmd:summarize, meanonly};
{cmd:aweight}s rescale the weights and {cmd:fweight}s require integers.
To use the mass in a formula, declare {cmd:[iw/]}: {cmd:`exp'} is then the name of the variable.{p_end}

{pstd}
{bf:Statistics with kinks.} For a statistic whose derivative jumps, such as a rank-based statistic, the numerical derivative is exact
as long as no two units cross within {cmd:step()} times their changes; the default step makes this rare, and the built-ins avoid the issue.
Statistics whose influence function jumps, such as quantiles, have no drift in panel data.


{marker options}{...}
{title:Options}

{phang}
{cmd:statistic(}{it:statistic}{cmd:)} the statistic to decompose; required.

{phang}
{cmd:saving(}{it:filename}[{cmd:, replace}]{cmd:)} saves the decomposition to a dataset; {cmd:clear} replaces the dataset in memory with it.
One of the two is required.

{phang}
{cmd:horizon(}{it:#}{cmd:)} length of the intervals, in periods of the panel; default {cmd:horizon(1)}.
A term of first order in time grows in proportion to the interval: to compare horizons, run {cmd:ddd} once per horizon.

{phang}
{cmd:lag(}{it:#}{cmd:)} computes {cmd:drift} and {cmd:displacement} with lagged weights, which remove the bias from transitory measurement error and from flows (such as annual sales).
Each continuing unit's weight is then the derivative of the statistic at the positions, {it:#} time units before the start of the interval, of the continuing units observed at that date.
The correction, drift with lagged weights minus drift with weights at the start, is estimated on these units as a share {cmd:chi_lag} of their displacement, and this share of displacement is moved into drift for all continuing units.
For {cmd:share} and {cmd:average}, whose displacement can be close to zero in the middle of the distribution, the correction is that of the upper bands [{it:lo}, 1] and [{it:hi}, 1], combined.
{cmd:chi_lag} is not capped; capping it at one keeps displacement of the same sign as displacement with weights at the start.
Intervals without continuing units observed at the lagged date, such as the first ones, are dropped.
Drift and displacement with weights at the start are {cmd:drift} - {cmd:chi_lag} {cmd:displacement} / (1 - {cmd:chi_lag}) and {cmd:displacement} / (1 - {cmd:chi_lag}), or the output of {cmd:ddd} without {cmd:lag()}.
{it:#} must be a multiple of the panel's delta.

{phang}
{cmd:normalize} divides positions by the continuing units' mean at each date when splitting the change among continuing units into drift and displacement.
It changes only this split: the statistic, its change among continuing units, and demography are unchanged.
It requires a scale-invariant statistic ({cmd:ge}, {cmd:theil}, {cmd:mld}, {cmd:atkinson}, {cmd:gini}, {cmd:share}, {cmd:hhi}).
A common growth rate is then not drift, and drift is the covariance between a weight and a change in normalized position.

{phang}
{cmd:step(}{it:#}{cmd:)} step of the numerical derivative for user-defined statistics; default {cmd:step(1e-5)}.

{phang}
{cmd:by(}{it:varlist}{cmd:)} computes the decomposition separately within each group. A unit that changes group exits the first and enters the second.


{marker output}{...}
{title:Output}

{pstd}
One observation per group and interval:

{p2colset 8 30 32 2}{...}
{p2col :{it:byvars}}the {cmd:by()} variables{p_end}
{p2col :{it:timevar}{cmd:0}, {it:timevar}{cmd:1}}start and end of the interval{p_end}
{p2col :{cmd:nu0}, {cmd:nu1}}statistic at the start and at the end{p_end}
{p2col :{cmd:total}}{cmd:nu1} - {cmd:nu0}{p_end}
{p2col :{cmd:drift}, {cmd:displacement}}terms among continuing units, with lagged weights with {cmd:lag()}{p_end}
{p2col :{cmd:demography}}{cmd:entry} + {cmd:exit}{p_end}
{p2col :{cmd:entry_mass}, {cmd:entry_selection}}mass and selection parts of {cmd:entry}{p_end}
{p2col :{cmd:exit_mass}, {cmd:exit_selection}}mass and selection parts of {cmd:exit}{p_end}
{p2col :{cmd:N0}, {cmd:N1}, {cmd:NC}}mass of units at the start, at the end, and of continuing units{p_end}
{p2col :{cmd:NC_lag}}with {cmd:lag()}: mass of continuing units also observed at the lagged date{p_end}
{p2col :{cmd:chi_lag}}with {cmd:lag()}: the correction for lagged weights as a share of displacement with weights at the start{p_end}
{p2colreset}{...}

{pstd}
{bf:Growth of the statistic.} For a positive statistic, decompose ln {cmd:nu1} - ln {cmd:nu0} by taking logs of the split in levels.
With the continuing units' indices nuC0 = {cmd:nu0} + {cmd:exit} and nuC1 = {cmd:nu1} - {cmd:entry},
log drift is ln((nuC0 + {cmd:drift}) / nuC0), log displacement ln(nuC1 / (nuC0 + {cmd:drift})),
log entry ln({cmd:nu1} / nuC1), and log exit ln(nuC0 / {cmd:nu0}).


{marker se}{...}
{title:Standard errors}

{pstd}
{cmd:ddd} reports no standard errors, since what to resample depends on the data.
To get them, rerun {cmd:ddd} on bootstrap samples and take the standard deviation, across draws, of the quantity you report:
a component, its average over intervals or groups, or a ratio of components.
The draws give the standard error of such a combination of observations of the output, which the standard errors of each observation would not.

{pstd}
Draw whole histories, so that each draw keeps every date of a unit, including the lagged positions used by {cmd:lag()}.
In a household survey, draw households.
In a near-census of firms, the uncertainty is rather over groups, such as markets or industries: draw groups, with all their units.
This also matters for the statistics that depend on the number of units ({cmd:hhi}, {cmd:powersum}, {cmd:expsum}), for which a unit drawn twice counts as two units.
With survey replicate weights, run {cmd:ddd} once per replicate weight instead, and use the survey's variance formula.

{pstd}Standard errors of the components of the Gini of wealth, averaged over intervals, from 200 draws of households{p_end}
{phang2}{cmd:. tsset hid year}{p_end}
{phang2}{cmd:. set seed 1}{p_end}
{phang2}{cmd:. local draws}{p_end}
{phang2}{cmd:. forvalues b = 1/200 {c -(}}{p_end}
{phang2}{cmd:  2.     preserve}{p_end}
{phang2}{cmd:  3.     bsample, cluster(hid) idcluster(newid)}{p_end}
{phang2}{cmd:  4.     tsset newid year}{p_end}
{phang2}{cmd:  5.     tempfile d`b'}{p_end}
{phang2}{cmd:  6.     ddd wealth [pw = hw], statistic(gini) normalize lag(1) saving(`d`b'')}{p_end}
{phang2}{cmd:  7.     restore}{p_end}
{phang2}{cmd:  8.     local draws `draws' `d`b''}{p_end}
{phang2}{cmd:  9. {c )-}}{p_end}
{phang2}{cmd:. clear}{p_end}
{phang2}{cmd:. append using `draws', generate(draw)}{p_end}
{phang2}{cmd:. collapse (mean) drift displacement demography, by(draw)}{p_end}
{phang2}{cmd:. summarize drift displacement demography}{p_end}

{pstd}
The standard deviations are the standard errors; the point estimates come from {cmd:ddd} on the original data.
To draw markets instead, give each drawn market, and each of its firms, its own identifier:{p_end}
{phang2}{cmd:  3.     bsample, cluster(market) idcluster(newmarket)}{p_end}
{phang2}{cmd:  4.     egen newfirm = group(newmarket firm)}{p_end}
{phang2}{cmd:  5.     tsset newfirm year}{p_end}
{phang2}{cmd:  6.     tempfile d`b'}{p_end}
{phang2}{cmd:  7.     ddd sales, statistic(hhi) normalize lag(1) by(newmarket) saving(`d`b'')}{p_end}


{marker examples}{...}
{title:Examples}

{pstd}Gini of wealth in a household panel, in normalized wealth, with survey weights and lagged weights{p_end}
{phang2}{cmd:. tsset hid year}{p_end}
{phang2}{cmd:. ddd wealth [pw = hw], statistic(gini) normalize lag(1) clear}{p_end}

{pstd}Average wealth in the top 1%, with survey weights and lagged weights{p_end}
{phang2}{cmd:. ddd wealth [pw = hw], statistic(average 0.99 1) lag(1) clear}{p_end}

{pstd}Concentration of sales within each market, in shares, over intervals of one to four years, with lagged weights since annual sales are flows{p_end}
{phang2}{cmd:. tsset firm year}{p_end}
{phang2}{cmd:. forvalues h = 1/4 {c -(}}{p_end}
{phang2}{cmd:  2.     ddd sales, statistic(hhi) normalize lag(1) horizon(`h') by(market) saving(hhi_ddd`h', replace)}{p_end}
{phang2}{cmd:  3. {c )-}}{p_end}

{pstd}A user-defined statistic: the coefficient of variation{p_end}
{phang2}{cmd:. program define cv, rclass}{p_end}
{phang2}{cmd:  1.     syntax varname [iw]}{p_end}
{phang2}{cmd:  2.     tempvar d}{p_end}
{phang2}{cmd:  3.     summarize `varlist' [`weight'`exp'], meanonly}{p_end}
{phang2}{cmd:  4.     local m = r(mean)}{p_end}
{phang2}{cmd:  5.     gen double `d' = (`varlist' - `m')^2}{p_end}
{phang2}{cmd:  6.     summarize `d' [`weight'`exp'], meanonly}{p_end}
{phang2}{cmd:  7.     return scalar nu = sqrt(r(mean)) / `m'}{p_end}
{phang2}{cmd:  8. end}{p_end}
{phang2}{cmd:. ddd wealth, statistic(cv) clear}{p_end}

{pstd}A user-defined statistic of two variables: average productivity weighted by sales{p_end}
{phang2}{cmd:. program define wavg, rclass}{p_end}
{phang2}{cmd:  1.     syntax varlist(min=2 max=2) [iw]}{p_end}
{phang2}{cmd:  2.     gettoken x s : varlist}{p_end}
{phang2}{cmd:  3.     tempvar sx}{p_end}
{phang2}{cmd:  4.     gen double `sx' = `s' * `x'}{p_end}
{phang2}{cmd:  5.     summarize `sx' [`weight'`exp'], meanonly}{p_end}
{phang2}{cmd:  6.     local num = r(mean)}{p_end}
{phang2}{cmd:  7.     summarize `s' [`weight'`exp'], meanonly}{p_end}
{phang2}{cmd:  8.     return scalar nu = `num' / r(mean)}{p_end}
{phang2}{cmd:  9. end}{p_end}
{phang2}{cmd:. ddd tfp sales, statistic(wavg) clear}{p_end}


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
