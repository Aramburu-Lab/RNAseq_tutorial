# lessons.md — Cancer Transcriptomics Teaching Suite

Rules written after corrections, so the same mistake is not repeated. Reviewed at the start of every session on this project.

## 2026-09-26

- **Cohort choice: immune cells first.** Do not default to TCGA solid tumours when proposing cancer datasets. Order of preference here is innate immune cells (neutrophils, monocytes) → other immune cells → solid tumours with large cancer and healthy groups. TCGA is tier 3, processed counts only.
- **Integration is in scope by default.** Transcriptomics plans for this project include the RNA→protein question (matched proteomics, translation as the missing layer) and non-coding RNA nuances (library type, lncRNA filters, miRNA-seq) from the first draft, not as later additions.
- **Verify accessions before proposing them.** Dataset candidates go into the plan only with a "verify" column; accessions asserted from memory were replaced by searched ones (E-MTAB-11190, Cassetta 2019). ENA/BioStudies pages are not searchable from the sandbox, so FASTQ availability checks are handed to Rytis with an exact query.
- **Document edits are targeted replacements with grep verification**, applied through one scripted pass with an exactly-one-match assertion per edit; no whole-file rewrites.
- **Check-ins use tappable decisions**, at most three per turn, recommended default marked; everything else is stated as an assumption in the plan.
