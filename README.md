# ddd

This Stata command implements the decomposition of Gomez, ["The Micro Origins of Aggregate Trends: Drift, Displacement, and Demography"](http://www.matthieugomez.com/files/ddd.pdf). Between two dates, it splits the change in any distributional index (a Gini, a top share, a Herfindahl index, an aggregate productivity) computed on a panel of units into three terms:

- **Drift**: the sum over survivors of each survivor's change times the derivative of the index with respect to its position, at the initial distribution
- **Displacement**: the rest of the change among survivors, the effect of the weights changing as units move
- **Demography**: the change in the full population minus the change among survivors (= entry + exit = mass + selection)

The decomposition is exact: **total = drift + displacement + demography**.

## Installation

```stata
net install ddd, from("https://raw.githubusercontent.com/matthieugomez/ddd/main/")
```

## Syntax

```stata
ddd varlist [if] [in] [weight], index(index) {save(filename) [replace] | clear} [by(varlist) horizon(numlist) lag(#) normalize step(#)]
```

The dataset must be `tsset`. Weights (`fweight`, `iweight`, `pweight`) are treated as masses.

| Option | Meaning |
|---|---|
| `index()` | a built-in index, an rclass program, or `mata: fcn()` |
| `by()` | separate decomposition within each group |
| `horizon(1 2 4)` | intervals of 1, 2, and 4 periods |
| `lag(#)` | also report drift and displacement with lagged weights (measurement error, flows) |
| `normalize` | divide positions by the survivors' mean at each date (scale-invariant indices) |

## Indices

Built-in indices use analytic derivatives: `mean`, `var`, `power a`, `expmean g`, `ge a`, `theil`, `mld`, `atkinson e`, `gini`, `top p`, `hhi`.

You can also define your own index (see [below](#define-your-own-index)).

## Output

One observation per group and interval: `time0`, `time1`, the masses `N0`, `N1`, `NS`, the index `nu0`, `nu1`, `total`, `drift`, `displacement`, `demography`, `entry`, `exit`, `mass`, `selection`, and with `lag()`, `NS_lag`, `drift_lag`, `displacement_lag`.

## Define your own index

An index can be supplied as a program that only *computes* the index. Drift is the directional derivative of the index along survivors' changes, `d/ds nu(x + s dx)` at `s = 0`, computed numerically with four evaluations of the index whatever the number of units. Several variables per unit are allowed.

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
