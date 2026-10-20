# Statistical review of the methods (what is appropriate, what is only descriptive, what was changed)

| Component | Method | Assessment | Caveat / action taken |
|---|---|---|---|
| Proportions | Wilson interval; cluster bootstrap by collector | Wilson is correct for independent samples; samples are clustered by collector | Both reported; the cluster-bootstrap interval is the one to quote |
| Concentration | Gini, Lorenz, Hill N0–N2, bootstrap CI | Standard, descriptive | Bootstrap resamples collectors' counts |
| pH digit preference | Observed/expected whole-unit ratio by interpolation of adjacent half units; chi-square with Monte-Carlo p among collectors with ≥10 samples | Transparent but ad hoc; descriptive indices (Whipple/Myers are alternatives) | Reported as ratio with cluster-bootstrap CI; only 14 collectors enter the heterogeneity test |
| Earthworm counts | Poisson dispersion test, NB2 fit, expected zeros | Standard | Heaping tests against a geometric interpolation are exploratory (p-values indicative only) |
| Spatial autocorrelation | Moran's I on mid-ranks, binary distance bands, row-standardised; permutation tests | Appropriate for ordinal/count data; exact under exchangeability | Global null mixes regional gradient and local structure → restricted nulls (50/100-km blocks) added; closest band is biased towards zero by the ±1 km coordinate displacement |
| Nominal variables | Join-count type agreement, chance-corrected ("spatial kappa") | Sound but not a textbook statistic; interpret as spatial agreement beyond chance | Null from the same permutation schemes |
| Multiplicity | Benjamini–Hochberg within each null model | Valid under independence/positive dependence | Tests across bands are positively dependent; accepted |
| Neighbour consistency | Spearman ρ between a sample and the median of its k nearest eligible neighbours | Intuitive but non-standard; neighbourhoods overlap, so naive p-values and CIs are invalid | Inference by permutation (global, within 100-km blocks); CIs by collector-cluster and by spatial-block bootstrap (the latter reported in figures) |
| Collector vs place | One-way ICC(1) on ranks | Descriptive; collector and place are confounded | Replaced as the main estimate by a crossed random-effects model (collector + 25-km cell + 100-km block) on rank-based inverse-normal scores, REML, LRT with 50:50 χ² mixture. Powell optimiser used (L-BFGS did not converge for reduced models) |
| Collector effect, local | Pairwise dissimilarity same vs different collector within 25-km cells; collector labels permuted within cells | Good design: removes place by restriction | Informative only where several collectors share a cell |
| Kruskal–Wallis among collectors | Rank test | Confounded with place; descriptive | Only 14 collectors with ≥10 samples |
| Missing OM | Nearest-3 neighbour mean, snapped to class | Neighbour imputation builds spatial dependence | Primary analyses use complete cases; imputed version used as sensitivity analysis (no change in any conclusion) |

Not corrected (limitations to state): preferential and opportunistic sampling (72% of samples in one month), possible temporal confounding, unknown form of the coordinate displacement beyond "fixed random value within 1 km", and the use of Gaussian scores for ordinal variables with few classes (many ties).
