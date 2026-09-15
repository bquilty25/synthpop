# Handoff: synthpop vignette for agentic development

## What was done

Investigated using synthetic data to let AI coding agents develop R analysis pipelines without seeing real patient data. Evaluated `synthpop`, `fabricatr`, and `charlatan`. `synthpop` (parametric method) was the clear winner for preserving correlations when the real data is available. `fabricatr` generates independent marginals only, so all correlations collapse to near zero.

## Deliverables

### Vignette (`vignettes/agentic-development.Rmd`)

On branch `vignette-agentic-development`. Covers:

- `syn()` with parametric and CART methods
- Preparing categorical variables (factor conversion before synthesis)
- `compare()` with consistent colour scheme
- GGally pairwise comparison plots (numeric and mixed factor/numeric)
- `utility.gen()` with S_pMSE and SPECKS metrics
- Correlation preservation comparison
- Model fit comparison with `lm.synds()` and `compare()`
- Disclosure control: `replicated.uniques()`, `disclosure()`, `sdc()`
- Data governance section on population-level information leakage
- Integration workflow: env var toggle (`SYNTHETIC_DATA`), agent permissions
- All output formatted as clean data frames

Renders cleanly with `rmarkdown::render("vignettes/agentic-development.Rmd")`.

### Related repos (can be cleaned up)

- `bquilty25/scramble`: thin wrapper package around synthpop. We concluded the vignette approach is better. Can be archived or deleted.
- `bquilty25/synthetic`: predecessor to scramble. Should be deleted (`gh repo delete bquilty25/synthetic --yes`).

## Key findings

- **Parametric vs CART**: parametric preserves pairwise correlations better for most datasets. CART is better when categorical variables define subpopulations with distinct patterns (e.g. iris Species). Default to parametric.
- **Factor handling matters**: numeric columns with few levels (cyl, vs, am, gear in mtcars) should be converted to factors before synthesis so synthpop uses logistic/multinomial models rather than linear.
- **Disclosure risk**: CART produces more replicated uniques than parametric. Use `sdc(result, real, rm.replicated.uniques = TRUE)` to remove them.
- **Population-level risk**: synthetic data preserves the statistical structure of the original. For small or rare-disease cohorts, this structure itself is sensitive, carrying the same information as detailed cross-tabulations. Synthesis does not solve the governance problem for these datasets.
- **The practical test**: if you would be comfortable sharing summary statistics and cross-tabulations with the API, synthetic data is comparable. If not, synthesis is not enough.

## What remains

- PR the vignette to upstream `gillian-raab/synthpop`, or keep it as a standalone guide
- Test on a real epidemiological dataset with dates, mixed types, and more rows
- Get DPO input on whether synthetic data derived from specific cohorts qualifies as non-personal data under GDPR
- Consider adding a section on `fabricatr` for schema-only generation (when real data is unavailable)
- Consider adding a correlation heatmap comparison (we built one earlier but it didn't make it into the vignette)
