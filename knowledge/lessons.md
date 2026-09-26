# lessons.md — Cancer Transcriptomics Teaching Suite

Rules written after corrections, so the same mistake is not repeated. Reviewed at the start of every session on this project.

## 2026-09-26

- **Cohort choice: immune cells first.** Do not default to TCGA solid tumours when proposing cancer datasets. Order of preference here is innate immune cells (neutrophils, monocytes) → other immune cells → solid tumours with large cancer and healthy groups. TCGA is tier 3, processed counts only.
- **Integration is in scope by default.** Transcriptomics plans for this project include the RNA→protein question (matched proteomics, translation as the missing layer) and non-coding RNA nuances (library type, lncRNA filters, miRNA-seq) from the first draft, not as later additions.
- **Verify accessions before proposing them.** Dataset candidates go into the plan only with a "verify" column; accessions asserted from memory were replaced by searched ones (E-MTAB-11190, Cassetta 2019). ENA/BioStudies pages are not searchable from the sandbox, so FASTQ availability checks are handed to Rytis with an exact query.
- **Document edits are targeted replacements with grep verification**, applied through one scripted pass with an exactly-one-match assertion per edit; no whole-file rewrites.
- **Check-ins use tappable decisions**, at most three per turn, recommended default marked; everything else is stated as an assumption in the plan.
- **Prove a template update with a real one.** An update to the same template version is a no-op and cannot fail. Tag a synthetic release in a scratch clone that edits next to every project addition, and update a scratch project clone to it. Lines appended at the end of CLAUDE.md and AGENTS.md conflicted with the template's own appends; at the top of each file they merge.
- **Measure the finished figure file, not the code.** On macOS, pyplot's GUI backend made an 85 mm page 84.84 mm, and seaborn desaturated every box fill to 75 %. Both were invisible in the code and caught only by reading the PDF. Build figures with `matplotlib.figure.Figure`, pass `saturation=1`, and assert page size and painted colours.
- **Files built from the same coordinates cannot check each other.** A wrong rebase offset kept slice.fa, slice.gtf and transcripts.fa mutually consistent and passed every offline test. Anchor a fixture to external truth: its start codons must read ATG.
- **A pre-registered rule needs its content checked, not only its commit order.** A rule rewritten between the rule commit and the record commit passed an ancestry test. Compare the rule as first committed with the rule in force.
- **Eye-check a figure at its real size before accepting a style.** 14 pt text in an 85 mm column overlapped for 8 genes. The sizes stay, the layout adapts, and a test measures every text box for overlap.
- **Raw-read download is nf-core/fetchngs 1.12.0, not custom code.** The first plan proposed ENA API / `fastq-dl` download logic; Rytis replaced it with fetchngs. The suite consumes the fetchngs output layout through `ctx gather` alone and teaches accession discovery around it.
- **Read a pipeline's output layout from its code, not its docs.** The fetchngs docs list `*.runinfo.tsv` and untimestamped `pipeline_info` names that 1.12.0 does not publish, and md5 files exist only for FTP and Aspera downloads. Pin the layout with a real run of the pinned release, committed as a fixture.
- **Assert byte identity only on what a tool fixes.** A fetchngs rebuild reproduced every file except the row order of one `runinfo_ftp.tsv`; fetchngs does not fix that order. Check what the tool guarantees, compare the rest by content, and retry a network proof once before calling it failed: the first download rebuild stalled on NCBI and the retry passed.
