# popgen-r-utils

**`popgenutils`** — reusable R helpers for SSR (microsatellite) population-genetic
analysis.

This package grew out of the one-off scripts in
[genetic_diversity_analysis](https://github.com/seema6203/genetic_diversity_analysis),
[NJ_CLUSTER](https://github.com/seema6203/NJ_CLUSTER) and
[PCOA](https://github.com/seema6203/PCOA), written for a study of Indian
*Nymphaea* (water lily) populations. Those scripts work, but each one hardcodes
its own paths, its own colour vector and its own axis labels, so none of them
can be reused on a second dataset without editing. This package is the same
analyses, turned into documented functions with arguments.

## What it fixes

| In the original scripts | Here |
| --- | --- |
| `setwd("/Users/Seema/Documents/pcoa/")` | Paths are function arguments. Nothing is hardcoded. |
| `colr <- c(rep("purple", 2), rep("red", 39), ...)` — breaks silently if a sample is added or reordered | `cluster_colours(pop(gen))` looks colours up by group. Correct by construction, and colour-blind-safe by default. |
| `labs(x = "PCo1 (18.74%)")` — typed in by hand, wrong as soon as the data changes | `plot_pcoa()` computes the percentage from `pco$eig`. |
| `q$data[[3]]$colour <- ...` — reaching into `ggplot_build()` internals, tied to layer order | A plain `ggplot` object you extend with `+`. |
| Nei distance fed straight to `dudi.pco()`, producing negative eigenvalues | `run_pcoa()` applies a Cailliez or Lingoes correction by default. |
| Stats spread over `basic.stats()`, `poppr()` and `Ho()` printouts | `diversity_summary()` returns one tidy data frame. |
| `set.seed()` sometimes present, sometimes not | `run_amova()` takes and records a `seed`. |

## Installation

```r
# install.packages("remotes")
remotes::install_github("seema6203/popgen-r-utils")
```

Dependencies: `adegenet`, `ade4`, `ape`, `ggplot2`, `hierfstat`, `poppr`.

```r
install.packages(c("adegenet", "ade4", "ape", "ggplot2", "hierfstat", "poppr"))
```

## Quick start

```r
library(popgenutils)
library(adegenet)

gen    <- read_genalex_csv("data/genotypes.csv")
groups <- pop(gen)

# Diversity table: n, MLG, Ho, Hs, Fis per population plus an overall row
diversity_summary(gen)

# Pairwise differentiation
pairwise_fst(gen, method = "WC")

# AMOVA with a reproducible permutation test
run_amova(gen, ~Pop, nrepet = 999, seed = 1999)

# Bootstrapped NJ tree
fit <- nj_bootstrap(gen, bootstrap = 1000, support_cutoff = 70)
save_figure("nj_tree.pdf", plot_nj(fit, groups = groups))

# PCoA, with axis labels computed from the eigenvalues
pc <- run_pcoa(gen, correction = "cailliez")
save_figure("pcoa.pdf", plot_pcoa(pc, groups = groups))
```

A full worked script is in
[`inst/examples/run_analysis.R`](inst/examples/run_analysis.R).

## Function reference

### Reading and writing

| Function | Purpose |
| --- | --- |
| `read_genalex_csv(path)` | Read a GenAlEx CSV, with a readable error if the format is wrong and a warning if no populations are assigned. |
| `save_figure(filename, plot)` | Write a base or ggplot figure to PDF, PNG, TIFF or SVG. Closes the device even if plotting fails. |

### Diversity and differentiation

| Function | Purpose |
| --- | --- |
| `diversity_summary(gen)` | Data frame of `n`, `mlg`, `Ho`, `Hs`, `Fis` per population, plus an overall row. |
| `pairwise_fst(gen, method)` | Pairwise F<sub>ST</sub>, Weir–Cockerham (`"WC"`) or Nei (`"nei"`). |
| `genetic_distance(gen, method)` | Distance between individuals: `"nei"`, `"euclidean"` or `"bruvo"`. |
| `run_amova(gen, hier, nrepet, seed)` | AMOVA plus its permutation test, returned together with the seed used. |

### Figures

| Function | Purpose |
| --- | --- |
| `nj_bootstrap(gen, bootstrap, support_cutoff)` | NJ tree with bootstrap support; support below the cutoff is set to `NA`. |
| `plot_nj(fit, groups)` | Draw the tree with coloured tips, node support, legend and scale bar. |
| `run_pcoa(gen, nf, correction)` | PCoA, returning coordinates **and** variance explained per axis. |
| `plot_pcoa(pc, groups, axes)` | ggplot ordination with data-derived axis labels and optional ellipses. |
| `cluster_colours(groups, palette)` | One colour per individual, looked up by group. Okabe–Ito palette by default. |

## Input format

Genotypes are read with `poppr::read.genalex()`, so the input must be a
**GenAlEx-formatted CSV**: three header rows (number of loci, number of samples,
population sizes) followed by one row per individual, two columns per locus for
a diploid.

## A note on bootstrapping

`nj_bootstrap()` defaults to 1000 replicates, which is the right number for a
final figure and slow enough to be annoying while you are still adjusting
colours and layout. Pass `bootstrap = 0` (or a small number) while iterating,
then raise it for the run that goes in the paper.

## Status

Version 0.1.0 — early. The functions are documented and unit-tested where they
can be tested without genotype data, but the package has not yet been run
end-to-end against a full dataset. Treat the first run on your own data as a
check: compare the output against the results from the original scripts before
relying on it.

If something is wrong or missing, please open an issue.

## License

MIT — see [LICENSE](LICENSE).
