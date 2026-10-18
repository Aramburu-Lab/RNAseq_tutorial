# R3 proof: the external check, registered before any protein value of PXD038437 is compared

Registered 2026-10-18, before the check script exists or runs. Committed alone; the check's
commit follows it.

## The study

Ramberger E. et al. 2024, "The proteogenomic landscape of multiple myeloma reveals insights into
disease biology and therapeutic opportunities", Nature Cancer, doi:10.1038/s43018-024-00784-3,
PMID 38942927, PMC11358022 (open access). Protein: PRIDE PXD038437, `proteinGroups.txt`,
MaxQuant 2.0.3.0, TMTpro 16-plex, 11 plexes, each with an internal standard channel. RNA: GEO
GSE222727 (123 per-sample RSEM `genes.results`, not fetched), and the paper's own normalised RNA
in Supplementary Table S2, keyed by TMT channel.

The PMC full text was read in full on 2026-10-18 (59,830 characters; figure legends are not in
that text). The supplement was read through Europe PMC's supplementaryFiles zip (209,347,597 B).

## The statement

Results, "Proteomic landscape of newly diagnosed MM": RNA-to-protein correlation was moderate,
with a median Pearson correlation coefficient of 0.29. Methods, "RNA-protein correlation":
RNA samples filtered for at least 80 % plasma cells and over 20 million mapped reads; protein
collapsed to gene by median and matched to RNA by gene name; per-gene Pearson across MM samples
on genes with at least ten paired values. Protein values: corrected reporter intensities, log2,
minus each plex's internal standard, then median-MAD normalised; protein groups filtered for
reverse hits, contaminants, only-by-site, and at least two peptides.

## The check, R3-1

The protein table read by ctx's generic reader through the guide's route (MaxQuant adapter,
`Reporter intensity corrected`, the internal standard as each plex's reference), against the
paper's Table S2 RNA on the columns its `use_for_correlation` row marks TRUE, matched by TMT
channel; per-gene Pearson on genes with at least ten paired values; the median.

- Recorded: our median beside the paper's 0.29.
- Gate (a course convention, not the paper's): the median within 0.05 of 0.29.

What differs from the paper, stated before the run: no median-MAD normalisation; one row per
gene (the best-quantified) instead of the median of a gene's groups; no two-peptide filter; a
group naming several genes left out. The RNA is the paper's own table, not an RNA route of
ctx, so this checks the protein reader only. The paired M7 run waits for an RNA route (see the
R3 report).
