# rRMSAnalyzer 2.0.1

## New feature

-   `get_annotation()`: get a dataframe with all current annotation in a RiboClass.

-   `remove_annotation()`: Added a new parameter to select a subset of annotation to remove.

# rRMSAnalyzer 2.0.0

-   Added a `NEWS.md` file to track changes to the package.

## Breaking change

-   `plot_PCA()`has been renamed to `plot_pca()`.

## New features

### Plots

-   Correspondence analysis for samples' counts : `plot_coa()`
-   Heatmaps: `plot_heatmap()` and `plot_heatmap_corr()`
-   "Simple" boxplots : `boxplot_count()` (by sample) and `boxplot_cscores()` (by site)
-   Boxplot comparing c-score by condition for sites considered as differential `plot_diff_sites()`
-   And others: `plot_counts_env()`, `plot_counts_fraction()`...

### Error messages

-   Error messages are now using Cli package for better clarity

-   More error messages to help the users

## Bugfixes & improvements

-   `compute_cscore()` has drastically been improved in terms of performance

# rRMSAnalyzer 1.0.0

First public version of rRMSAnalyzer 🥳 !

## New features

Only main features will be shown here :

-   C-score computation for all positions (`compute_cscore()`)

-   Technical biases adjustment with CombatSeq (`adjust_bias()`)

-   Principal Component Analysis (`plot_PCA()`)
