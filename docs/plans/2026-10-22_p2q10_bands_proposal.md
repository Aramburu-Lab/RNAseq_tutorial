# P2Q.10 bands: top_gene_share and rrna_fraction, a proposal for Rytis

**Status:** open for Rytis (D-052: bands are his call). Nothing in `conf/qc_metrics.toml`
has a band from this file; both metrics stay `bands = "dataset_only"` until he answers.

**Data.** Step-1 (E-MTAB-11190), 210 runs, every one Salmon libType U (210 of 210
`aux_info/meta_info.json`). Measured 2026-10-22 by `ctx qc` from this branch into a scratch
root that links `~/ctx_tier1/data` and `resources` read-only; GENCODE 50 primary-assembly
annotation. The judge's strata are the sheet's cell type groups (LDN, NDN, Monocyte, BM_1 to
BM_4). Quantiles are numpy's linear interpolation over the 210 values.

## The definitions the numbers are on

| metric | definition | changed by P2Q.10 on this branch |
|---|---|---|
| top_gene_share | largest count among genes other than the annotation's Mt_rRNA genes and the registry's `exclude_gene_ids` (ENSG00000280614, ENSG00000281383) / library size | the two IDs moved from code into the registry with their source; with no annotation the metric is now left out instead of measured over every gene |
| rrna_fraction | counts on gene_type rRNA and rRNA_pseudogene off chrM, plus the same two lncRNAs / library size | not changed; see below |

Every Step-1 value, flag, band and status row is equal to the bit before and after the
change (210 of 210 values of each metric; 282 of 282 flag rows).

## top_gene_share

| | min | q05 | q25 | median | q75 | q90 | q95 | q99 | max |
|---|---|---|---|---|---|---|---|---|---|
| over every gene, pre-P2Q.10 | 0.0165 | | | 0.0534 | | | | | 0.1788 |
| corrected | 0.0104 | 0.0156 | 0.0223 | 0.0322 | 0.0496 | 0.0755 | 0.0872 | 0.1439 | 0.1788 |

Pre-P2Q.10, ENSG00000280614 was the top gene in 102 of 210 runs and MT-RNR2 in 16. Corrected,
the top gene is S100A9 in 63, ACTB, B2M and MALAT1 in 19 each, DEFA1B in 18, LYZ in 17, and 14
other genes in the rest; the largest share is DEFA1B in ERX6994400, 0.1788.

| stratum | n | median | 3 scaled MAD upper limit |
|---|---|---|---|
| NDN | 73 | 0.0332 | 0.0675 |
| Monocyte | 68 | 0.0202 | 0.0347 |
| LDN | 38 | 0.0585 | 0.1488 |
| BM_1 | 10 | 0.1078 | 0.2691 |
| BM_3 | 10 | 0.0557 | 0.1382 |
| BM_4 | 8 | 0.0393 | 0.0708 |
| BM_2 | 3 | 0.0550 | 0.1351 |

The stratum rule flags 7 runs: ERX6994436, ERX6994439, ERX6994404, ERX6994443, ERX6994460,
ERX6994508, ERX6994355.

| candidate cut-off | runs above it |
|---|---|
| 0.02 | 176 |
| 0.05 | 52 |
| 0.10 | 7 |
| 0.15 | 2 |

**Sources.** None found. The 2026-10-01 proposal searched and found no sourced cut-off for this
metric (its row: "none sourced"); this proposal made no new search, as it ran with no network.

**Options.**

1. Keep `dataset_only`: the stratum rule alone, 7 flags. The tab says why there is no band.
2. A course convention that shades only, `plot_bands = [0.10, 0.20]` (good to 0.10, usable to
   0.20): 7 runs in usable, none in bad. Method: round numbers above Step-1's q95 (0.087) and
   its max (0.179), so the shading marks the tail without judging. It is drawn from the data it
   shades, so it describes this cohort and is labelled a course convention.
3. A judging band at the same values: 7 runs would read usable on top of the stratum flags.
   Neutrophil and bone-marrow libraries are dominated by S100A9 and DEFA1B by biology, and
   BM_1's median is 0.108, so a fixed band judges cell type, not library quality.

**Recommended: option 2**, shading only; judging stays with the stratum rule, as for globin.

## rrna_fraction

**The item's definition, 18S and 28S reads, cannot be computed from what Step-1's run has.**

- GENCODE 50 (primary assembly, 78,941 genes) has no 18S, 28S or 45S gene: no gene_name
  RNA18S*, RNA28S* or RNA45S*. Its rRNA features are 53 rRNA genes, all 5S or 5.8S (17 RNA5S,
  18 RNA5SP typed rRNA, 4 5S_rRNA, 5 RNA5-8SN, 9 5_8S_rRNA), 497 rRNA_pseudogene and 2 Mt_rRNA.
- Of these the counts hold 27 nuclear rRNA genes, 487 rRNA_pseudogene and both Mt_rRNA, of
  78,470 genes; 3 of the 27 sit on unplaced contigs (GL000220.1, KI270442.1), at most
  6.05e-05 of a library.
- 18S reaches the counts only through the two 18S-antisense lncRNAs, and only because Step-1 is
  unstranded. In a stranded library the measure is a lower bound.
- 28S reaches no counted feature. Salmon's index decoys are the genome's 194 contigs, so a 28S
  read lands among the decoy fragments (3,437,107 of 38,268,783 in ERX6994339), with intronic
  and intergenic reads, and no file of the run says which contig it hit.

What is missing: a 45S rDNA reference sequence quantified with the reads, or alignments to
the genome counted over its rDNA copies. Both need the reads quantified again. A 45S accession
is not named here: none was looked up, and no network was available.

| rrna_fraction | min | q05 | q25 | median | q75 | q90 | q95 | q99 | max |
|---|---|---|---|---|---|---|---|---|---|
| rRNA and rRNA_pseudogene alone, pre-P2Q.10 | 5.94e-07 | | | 6.48e-06 | | | | | 6.79e-05 |
| the two lncRNAs alone | 0.00180 | | | 0.0410 | | | | | 0.1673 |
| as computed now, both | 0.00181 | 0.0062 | 0.0211 | 0.0410 | 0.0652 | 0.0867 | 0.0992 | 0.1228 | 0.1673 |

As computed now the stratum rule flags 2 runs: ERX6994528 (0.167, Monocyte) and ERX6994345
(0.0898, BM_2, 3 runs). Above 0.05: 81 runs; above 0.10: 11; above 0.15: 1.

**Options.**

1. No band, as now, until the run has an rRNA reference: the measure counts 18S in unstranded
   libraries and never 28S.
2. Accept the 18S proxy as the definition, renamed for what it measures, then band it.
3. Quantify Step-1 against a 45S reference on Kebnekaise, then propose a band on 18S + 28S.

**Recommended: option 1.** No band is proposed, since the definition the item asks for does not
exist in this run.

## Not done

- No band is written into `conf/qc_metrics.toml`.
- The registry has `action = "flag"` on rrna_fraction, as since b1fe960; D-052 decided keep until
  the definition is corrected, and the registry never took it. Not changed here.
- No source was searched for either metric.
