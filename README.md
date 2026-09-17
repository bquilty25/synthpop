# synthpop

R package for generating synthetic versions of sensitive microdata for
statistical disclosure control.

Synthetic data replaces sensitive original values with values drawn from
conditional distributions fitted to the original data, using parametric
models or classification and regression trees. Variables are synthesised
sequentially, each conditioned on those already generated. The result
preserves the statistical properties of the original while reducing
disclosure risk.

For a description of the method, see Nowok, Raab and Dibben (2016),
doi:[10.18637/jss.v074.i11](https://doi.org/10.18637/jss.v074.i11).

## Installation

Install the released version from CRAN:

```r
install.packages("synthpop")
```

Or install the development version from this repository:

```r
# install.packages("remotes")
remotes::install_github("bquilty25/synthpop")
```

## Usage

```r
library(synthpop)

syn_obj <- syn(mtcars, method = "cart", seed = 2024)
compare(syn_obj, mtcars)
```

`syn()` returns a `synds` object containing the synthetic data in `$syn`.
`compare()` produces per-variable distribution plots of the real and
synthetic data.

## Key functions

| Function | Purpose |
|---|---|
| `syn()` | Generate synthetic data |
| `compare()` | Per-variable distribution comparison |
| `lm.synds()` / `glm.synds()` | Fit models to synthetic data |
| `compare.fit.synds()` | Compare model coefficients (real vs synthetic) |
| `utility.gen()` | Propensity score utility measure |
| `utility.tables()` | Table-based utility measures |
| `replicated.uniques()` | Check for records matching unique real individuals |
| `disclosure()` | Identity and attribute disclosure risk |
| `sdc()` | Remove or replace high-risk synthetic records |

## Vignettes

The package includes vignettes on synthesis, utility, inference,
disclosure, and using synthetic data with AI coding agents:

```r
browseVignettes("synthpop")
```

## Authors

Beata Nowok (maintainer), Gillian M. Raab, Chris Dibben, Joshua Snoke,
Caspar van Lissa, Lotte Pater.

## Licence

GPL-2 | GPL-3

## Links

- Package website: <https://www.synthpop.org.uk/>
- CRAN: <https://cran.r-project.org/package=synthpop>
