# synthpop (fork)

This is a fork of [synthpop](https://github.com/bnowok/synthpop), an R
package for generating synthetic versions of sensitive microdata for
statistical disclosure control. The upstream package is maintained by
Beata Nowok and is available on
[CRAN](https://cran.r-project.org/package=synthpop).

This fork adds a vignette on using `synthpop` to generate synthetic data
for AI coding agents, and modernises the package's plotting code
(updated ggplot2 API usage, consistent theming, accessible colour
palette).

## Installation

```r
# install.packages("remotes")
remotes::install_github("bquilty25/synthpop")
```

## Changes from upstream

- **New vignette** (`vignettes/agentic-development.qmd`): using
  synthetic data to develop analysis pipelines with AI coding agents
  without exposing patient-identifiable data.
- **Plotting updates**: replaced deprecated `aes_string()` and
  `eval(parse())` with `.data` pronoun, switched `lwd` to `linewidth`,
  applied `theme_minimal()` consistently, updated default colour palette.

## Upstream

- Repository: <https://github.com/bnowok/synthpop>
- Package website: <https://www.synthpop.org.uk/>
- Reference: Nowok, Raab and Dibben (2016),
  doi:[10.18637/jss.v074.i11](https://doi.org/10.18637/jss.v074.i11)

## Licence

GPL-2 | GPL-3
