# ddd

This Stata command implements the decomposition of Gomez, ["The Micro Origins of Aggregate Trends: Drift, Displacement, and Demography"](http://www.matthieugomez.com/files/ddd.pdf). Between two dates, it splits the change in any distributional index (a Gini, a top share, a Herfindahl index, an aggregate productivity) computed on a panel of units into three terms:

- **Drift**: the sum over continuing units of each unit's change times the derivative of the index with respect to its position, at the initial distribution
- **Displacement**: the rest of the change among continuing units, the effect of the weights changing as units move
- **Demography**: the change in the full population minus the change among continuing units (= entry + exit)

The three terms sum up to the total change of the index.

## Installation

```stata
net install ddd, from("https://raw.githubusercontent.com/matthieugomez/ddd/main/")
```

## Syntax

```stata
ddd varlist [if] [in] [weight], index(index) {saving(filename[, replace]) | clear} [by(varlist) horizon(#) lag(#) normalize step(#)]
```

The dataset must be `tsset`. Weights (`fweight`, `iweight`, `pweight`) are treated as masses.

| Option | Meaning |
|---|---|
| `index()` | `mean`, `var`, `power a`, `expmean g`, `ces g`, `ge a`, `theil`, `mld`, `atkinson e`, `gini`, `top p`, `hhi`, or [your own index](#define-your-own-index) |
| `by()` | separate decomposition within each group |
| `horizon(#)` | length of the intervals in periods (default 1) |
| `lag(#)` | also report drift and displacement with lagged weights (measurement error, flows) |
| `normalize` | divide positions by the continuing units' mean at each date (scale-invariant indices) |

## Output

One observation per group and interval:

- with `by()`: the group variables
- the start and end of the interval: `{timevar}0`, `{timevar}1`
- the index at the start and at the end: `nu0`, `nu1`
- the components: `total` = `drift` + `displacement` + `demography`
- demography = `entry` + `exit`, with `entry` = `entry_mass` + `entry_selection` and `exit` = `exit_mass` + `exit_selection`
- the masses at the start, at the end, and of continuing units: `N0`, `N1`, `NC`
- with `lag()`: the mass of continuing units also observed at the lagged date, and drift and displacement with lagged weights: `NC_lag`, `drift_lag`, `displacement_lag`

## Define your own index

An index can be supplied as a program that only *computes* the index. Drift is the directional derivative of the index along continuing units' changes, `d/ds nu(x + s dx)` at `s = 0`, computed numerically with four evaluations of the index whatever the number of units. Several variables per unit are allowed.

As an ado program, returning the index in `r(nu)`, with masses in `_ddd_w` (use `double`s throughout):

```stata
program define cv, rclass
    syntax varname, w(varname)
    tempvar d
    summarize `varlist' [iw = `w'], meanonly
    local m = r(mean)
    gen double `d' = (`varlist' - `m')^2
    summarize `d' [iw = `w'], meanonly
    return scalar nu = sqrt(r(mean)) / `m'
end
ddd wealth, index(cv) clear
```

As a Mata function (faster), here of two variables:

```stata
mata:
real scalar wavg(real matrix X, real colvector w) return(quadsum(w :* X[., 2] :* X[., 1]) / quadsum(w :* X[., 2]))
end
ddd tfp sales, index(mata: wavg()) clear
```

## Tests

`do test/test_ddd.do` from the package folder checks the accounting identities, the closed forms of the paper's Table 1, that numerical drift equals analytic drift for every built-in index, weights, `by()`, `horizon()`, and `lag()`.

## Relation to `topdecompose`

[`topdecompose`](https://github.com/matthieugomez/topdecompose) is an earlier command, written for Gomez, ["Decomposing the Growth of Top Wealth Shares"](https://doi.org/10.3982/ECTA21396) (*Econometrica*, 2024). It should be seen as an application of this decomposition to one index, the average wealth in a top percentile: `within` is drift, `between` is displacement, and `demography` is demography. The built-in `top p` is the share of the top; the average can be supplied as [your own index](#define-your-own-index). `topdecompose` reports growth rates, so compare its terms with those of `ddd` divided by `nu0`. On a balanced panel, the two then coincide.

`topdecompose` was designed for a setup where the wealth of deaths and births is only observed in the top percentile, not in the whole population, as in the Forbes 400, whereas `ddd` needs the whole distribution. So, for top shares, `topdecompose` divides wealth by its mean in the whole population, whereas `ddd` with `normalize` divides by the mean among continuing units. Unless entrants and exiters are on average as rich as continuing units, part of demography then moves into `within`, even in continuous time.

Another smaller difference, second order in time, is that `between` does not purely measure displacement among continuing units: deaths and births move the threshold of the top, so continuing units also cross it without any re-ranking among them.
