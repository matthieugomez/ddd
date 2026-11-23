# ddd

This Stata command implements the decomposition of Gomez, ["The Micro Origins of Aggregate Trends: Drift, Displacement, and Demography"](http://www.matthieugomez.com/files/ddd.pdf). Between two dates, it splits the change in any distributional statistic (a Gini, a top share, a Herfindahl index, an aggregate productivity) computed on a panel of units into three terms:

- **Drift**: the sum over continuing units of each unit's change times the derivative of the statistic with respect to its position, at the initial distribution
- **Displacement**: the rest of the change among continuing units, the effect of the weights changing as units move
- **Demography**: the change in the full population minus the change among continuing units (= entry + exit)

The three terms sum up to the total change of the statistic.

## Installation

```stata
net install ddd, from("https://raw.githubusercontent.com/matthieugomez/ddd/main/")
```

## Syntax

```stata
ddd varlist [if] [in] [weight], statistic(statistic) {saving(filename[, replace]) | clear} [horizon(#) lag(#) normalize step(#) by(varlist)]
```

The dataset must be `tsset`. Weights (`fweight`, `iweight`, `pweight`) are treated as masses.

| Option | Meaning |
|---|---|
| `statistic()` | mean, var, powermean a, powersum a (CES aggregate), expmean g, expsum g (CES aggregate in logs), ge a, theil, mld, atkinson e, gini, share lo hi (share of the units with percentile rank between lo and hi, such as the top 1% or a decile), average lo hi (their average), hhi, or [your own statistic](#define-your-own-statistic) |
| `horizon(#)` | length of the intervals in periods (default 1) |
| `lag(#)` | drift and displacement with lagged weights (measurement error, flows); drops the intervals without continuing units observed at the lagged date |
| `normalize` | divide positions by the continuing units' mean at each date when splitting drift and displacement (scale-invariant statistics) |
| `by()` | separate decomposition within each group |

## Output

One observation per group and interval:

- with `by()`: the group variables
- the start and end of the interval: `{timevar}0`, `{timevar}1`
- the statistic at the start and at the end: `nu0`, `nu1`
- the components: `total` = `drift` + `displacement` + `demography`
- demography = `entry` + `exit`, with `entry` = `entry_mass` + `entry_selection` and `exit` = `exit_mass` + `exit_selection`: mass is the effect of adding entrants, or removing exiters, distributed like the continuing units, and selection is the rest
- the masses at the start, at the end, and of continuing units: `N0`, `N1`, `NC`
- with `lag()`: `drift` and `displacement` use lagged weights; `NC_lag` is the mass of continuing units also observed at the lagged date, and `chi_lag` the correction as a share of displacement with weights at the start, so that drift and displacement with weights at the start are `drift - chi_lag * displacement / (1 - chi_lag)` and `displacement / (1 - chi_lag)`

## Standard errors

`ddd` reports no standard errors, since what to resample depends on the data. Rerun it on bootstrap samples and take the standard deviation, across draws, of the quantity you report: a component, its average over intervals or groups, or a ratio of components. Draw whole histories, so that each draw keeps every date of a unit, including the lagged positions used by `lag()`: households in a survey; groups such as markets or industries, with all their units, in a near-census of firms (this also matters for `hhi`, `powersum`, and `expsum`, which depend on the number of units, so that a unit drawn twice counts as two units). With survey replicate weights, run `ddd` once per replicate weight instead, and use the survey's variance formula.

```stata
tsset hid year
set seed 1
local draws
forvalues b = 1/200 {
    preserve
    bsample, cluster(hid) idcluster(newid)
    tsset newid year
    tempfile d`b'
    ddd wealth [pw = hw], statistic(gini) normalize lag(1) saving(`d`b'')
    restore
    local draws `draws' `d`b''
}
clear
append using `draws', generate(draw)
collapse (mean) drift displacement demography, by(draw)
summarize drift displacement demography
```

The standard deviations are the standard errors; the point estimates come from `ddd` on the original data. To draw markets instead, give each drawn market, and each of its firms, its own identifier:

```stata
    bsample, cluster(market) idcluster(newmarket)
    egen newfirm = group(newmarket firm)
    tsset newfirm year
    tempfile d`b'
    ddd sales, statistic(hhi) normalize lag(1) by(newmarket) saving(`d`b'')
```

## Define your own statistic

You can supply the statistic as an ado program. The program must return the statistic in `r(nu)`. Use `double`s throughout:

```stata
program define cv, rclass
    syntax varname [iw]
    tempvar d
    summarize `varlist' [`weight'`exp'], meanonly
    local m = r(mean)
    gen double `d' = (`varlist' - `m')^2
    summarize `d' [`weight'`exp'], meanonly
    return scalar nu = sqrt(r(mean)) / `m'
end
ddd wealth, statistic(cv) clear
```

## Relation to `topdecompose`

[`topdecompose`](https://github.com/matthieugomez/topdecompose) is an earlier command, written for Gomez, ["Decomposing the Growth of Top Wealth Shares"](https://doi.org/10.3982/ECTA21396) (*Econometrica*, 2024). It applies this decomposition to one statistic, the average wealth in a top percentile (`statistic(average 0.99 1)` for the top 1%): `within` is drift, `between` is displacement, and `demography` is demography. It reports growth rates, so on a balanced panel its terms equal those of `ddd` divided by `nu0`.

`topdecompose` was designed for a setup where the wealth of deaths and births is only observed in the top percentile, not in the whole population, as in the Forbes 400. The mean among continuing units is then unknown, so for top shares `topdecompose` is applied to wealth divided by its mean in the whole population, whereas `ddd` with `normalize` divides by the mean among continuing units. Unless entrants and exiters are on average as rich as continuing units, part of demography then moves into `within`, even in continuous time.

Another smaller difference is that `between` measures not only the displacement among continuing units but also the distance to the threshold of the marginal units that cross it because of demography (this is second order in time).
