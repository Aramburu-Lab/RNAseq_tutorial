# Dataset survey: sizes and terms (work-order item 13 answer 6)

Read 2026-10-09. No data file was downloaded. Bytes are `Content-Length` from HEAD requests (GEO FTP over HTTPS, Xena S3) or `file_size` / `size` from documented metadata APIs (PDC GraphQL `filesPerStudy`, GDC `/files` at Data Release 46.0, GCS JSON listing). Kit cap: 64 MB (commit ecf03fc). Nothing in the repository was edited; pixi was not run.

## Table

| Dataset | Course use | File | Source | Bytes | Terms (URL) | Kit | Recommendation |
|---|---|---|---|---|---|---|---|
| SCAN-B | M10.12 breast OS validation (TODO.md:281) | GSE96058_gene_expression_3273_samples_and_136_replicates_transformed.csv.gz | https://ftp.ncbi.nlm.nih.gov/geo/series/GSE96nnn/GSE96058/suppl/ | 591,676,211 | NCBI: "NCBI itself places no restrictions on the use or distribution" of the data; submitters may claim IP rights. https://www.ncbi.nlm.nih.gov/home/about/policies/ . GSE96058 record states no licence. | Derived tables only (591 MB exceeds the cap) | Later: when M10.12 is scheduled |
| SCAN-B | OS, follow-up, clinical per sample | GSE96058-GPL11154_series_matrix.txt.gz | .../GSE96058/matrix/ | 218,719 | as above | Yes (small), or derived | Later, with the expression file |
| SCAN-B | as above, second platform | GSE96058-GPL18573_series_matrix.txt.gz | .../GSE96058/matrix/ | 42,331 | as above | Yes | Later, with the expression file |
| SCAN-B | not needed | GSE96058_transcript_expression_3273_samples_and_136_replicates.csv.gz | .../suppl/ | 860,717,799 | as above | No | Not |
| SCAN-B | not needed | GSE96058_UCSC_hg38_knownGenes_22sep2014.gtf.gz | .../suppl/ | 20,808,999 | as above | No | Not, unless gene IDs need its model |
| TOIL | Step-3 TCGA versus GTEx (TODO.md:98, 701) | TcgaTargetGtex_rsem_gene_tpm.gz | https://toil-xena-hub.s3.us-east-1.amazonaws.com/download/ | 1,323,254,426 | Xena: "there is no Xena Data Use Agreement"; cite Xena, check the providers. https://ucsc-xena.gitbook.io/project/data-use-agreement . TCGA: GDC open access, no authorisation needed, redistribution not addressed, https://gdc.cancer.gov/access-data/data-access-processes-and-tools . GTEx: licence below. | Derived tables only | Later: Phase 4 Step-3 |
| TOIL | counts, if DE runs on TOIL | TcgaTargetGtex_gene_expected_count.gz | same | 1,289,124,490 | same | Derived only | Later, only if DE on counts is chosen |
| TOIL | sample sheet (tissue, study, type) | TcgaTargetGTEX_phenotype.txt.gz | same | 135,753 | same | Yes | Later, with the matrix |
| TOIL | gene IDs, GENCODE v23 | probeMap/gencode.v23.annotation.gene.probemap | same | 3,244,244 | same | Yes | Later, with the matrix |
| TOIL | not needed | TcgaTargetGtex_rsem_isoform_tpm.gz | same | 4,467,218,883 | same | No | Not |
| TOIL | not needed | TcgaTargetGtex_RSEM_Hugo_norm_count.gz | same | 1,656,112,846 | same | No | Not |
| CPTAC ccRCC | M7 tumour versus NAT concordance, next after breast (TODO.md:706) | Supplementary_Data_Proteome_TMT.tar.gz (paper's tables; contents not opened) | PDC000127 (also listed under PDC000128, PDC000200) | 15,052,030, md5 0fa5f5b83969d0b12a177c912a206a77 | CPTAC DUA: acceptance required, publication guidelines, an acknowledgement statement; redistribution not addressed. https://cptac-data-portal.georgetown.edu/cptac/study/disclaimer/73217 | Derived only | Download now, on approval |
| CPTAC ccRCC | protein, CDAP gene report | CPTAC3_Clear_Cell_Renal_Cell_Carcinoma_Proteome.tmt10.tsv | PDC000127 | 79,198,141 | CPTAC DUA | Derived only | Later, only if the paper's tables do not suffice |
| CPTAC ccRCC | paper's RNA (RPKM) | CPTAC_CCRCC_Transcriptome_rpkm.tar.gz | PDC000127 | 26,807,229 | CPTAC DUA | Derived only | Later: cross-check only |
| CPTAC ccRCC | RNA counts | GDC STAR-Counts, CPTAC-3, PDC000127's cases: 202 tumour files (110 cases), 75 normal (75 cases) | GDC `/files` | 855,228,071 + 317,742,586 = 1,172,970,657 | GDC open access | Derived only | Download now, on approval |
| CPTAC ccRCC | design, clinical | S044_TMT10_Label_to_Sample_Mapping_File_CCRCC_r2_Jan2019.xlsx; S050_S044_CPTAC_ccRCC_Discovery_Cohort_Clinical_Data_r4_Sept2019.xlsx | PDC000127 | 23,811; 139,294 | CPTAC DUA | Derived only | Download now, on approval |
| CPTAC ccRCC | not needed | CPTAC_CCRCC_Methylation.tar.gz | PDC000127 | 8,258,306,942 | CPTAC DUA | No | Not |
| CPTAC LUAD | M7, after ccRCC | CPTAC3_Lung_Adeno_Carcinoma_Proteome.tmt10.tsv | PDC000153 | 95,159,608 | CPTAC DUA | Derived only | Later |
| CPTAC LUAD | design, clinical | S046_S056_BI_CPTAC3_TMT10_Label_to_Sample_Mapping_File_LUAD_v2.1-sample-annotation_plexes_01-25_20191219.xlsx; S046_S056_BI_CPTAC3_LUAD_Discovery_Cohort_Clinical_Data_r2_July2020.xlsx | PDC000153 | 57,866; 158,024 | CPTAC DUA | Derived only | Later |
| CPTAC LUAD | RNA counts | GDC STAR-Counts: 111 tumour, 102 normal files | GDC | 470,047,599 + 431,990,692 = 902,038,291 | GDC open access | Derived only | Later |
| CPTAC LUAD | paper's freeze | none listed at PDC000153, PDC000149, PDC000224 or PDC000489 | PDC | not found | n/a | n/a | Open question: where Gillette 2020's tables are |
| CPTAC LSCC | M7, after ccRCC | Satpathy_Supplementary_Data.tar.gz (contents not opened) | PDC000234 (also PDC000232, 233, 237) | 67,044,939, md5 936ceae4a5973d3f20ff60159ce2f4ca | CPTAC DUA | Derived only | Later |
| CPTAC LSCC | protein, CDAP gene report | CPTAC3_Lung_Squamous_Cell_Carcinoma_Proteome.tmt11.tsv | PDC000234 | 98,078,567 | CPTAC DUA | Derived only | Later |
| CPTAC LSCC | design, clinical | S063_BI_CPTAC3_TMT11_Label_to_Sample_Mapping_File_LSCC_Satpathy_Cell2021.xlsx; S063_S058_CPTAC_LSCC_Discovery_Cohort_Clinical_Data_r2_July2021.xlsx | PDC000234 | 42,529; 127,185 | CPTAC DUA | Derived only | Later |
| CPTAC LSCC | RNA counts | GDC STAR-Counts: 109 tumour files (108 cases), 95 normal | GDC | 461,540,020 + 402,313,794 = 863,853,814 | GDC open access | Derived only | Later |
| GTEx V11 | M9 junction panel (TODO.md:442) | GTEx_Analysis_2025-08-22_v11_STARv2.7.11b_junctions.gct.gz (per-sample counts; no median junction file exists in V11) | https://storage.googleapis.com/storage/v1/b/adult-gtex/o?prefix=bulk-gex/v11/ | 4,329,808,666, md5 (base64) awD/7xuai4SC9qcdsFfI5g== | GTEx licence: acknowledge with a dated phrase; imply no endorsement; a redistributor keeps data current or states it does not. https://gtexportal.org/home/license | Derived medians only | Not on the laptop; later, if V11 junctions are required, as a cluster job that writes per-tissue medians |

## Facts behind the table

- PDC studies (CPTAC3 Discovery), GraphQL `getPaginatedUIStudy` and `biospecimenPerStudy`: ccRCC PDC000127, TMT10, 124 cases, 208 aliquots (110 Primary Tumor, 84 Solid Tissue Normal, 9 Not Reported, 5 Cell Lines); LUAD PDC000153, TMT10, 115 cases, 217 aliquots (111 tumour, 102 NAT); LSCC PDC000234, TMT11, 115 cases, 217 aliquots (110 tumour, 102 NAT). Phospho: PDC000128, PDC000149, PDC000232. A LUAD Confirmatory proteome, PDC000489 (TMT11, 131 cases), and a ccRCC DIA proteome, PDC000200, also exist.
- PDC-to-GDC matching by case submitter ID: ccRCC 110 of 110 tumour cases and 75 of 84 NAT cases have GDC STAR-Counts; LUAD 111 of 111 and 102 of 102; LSCC 108 of 110 and 95 of 102. ccRCC has 202 tumour files for 110 cases, so several aliquots per case: choose one per case before M5.
- Check of the method: PDC000120's freeze tarball reads 34,878,978 B here, equal to the M7 plan.
- CPTAC freedom to publish: met once a global analysis paper is out; Clark 2019, Gillette 2020 and Satpathy 2021 are those papers. Requested statement: "Data used in this publication were generated by the Clinical Proteomic Tumor Analysis Consortium (NCI/NIH)." The OCCPR DUA URL (proteomics.cancer.gov/data-portal/about/data-use-agreement) now redirects to the OCCPR home page.
- TOIL: Xena hub API metadata gives version 2016-09-03, unit log2(tpm+0.001) for TPM and log2(expected_count+1) for counts, GENCODE v23 probemap; `rows` 19,131 for TPM and phenotype, 19,109 for expected counts (Xena's field; read as sample counts, UNVERIFIED). Files last modified 2021-04-09. The values are already log-transformed: undo before any count model.
- GTEx API v2 refuses `datasetId=gtex_v11` for `medianJunctionExpression` (accepts gtex_v8, gtex_snrnaseq_pilot, gtex_v10). A per-gene V11 junction panel is impossible through the API today; the alternatives are the 4.33 GB file or V10 medians labelled as V10.
- GTEx V11 bucket also holds a 2026-05-19 re-release of gene reads and TPM with LCM samples excluded (gene_tpm.parquet 2026-05-19, gene_reads.gct.gz 1,008,045,909 B); the course's on-disk copy is the 2025-08-22 release. Worth a decision on which release M9 cites.

## Not done

- NCBI GEO's own disclaimer and FAQ pages (www.ncbi.nlm.nih.gov/geo/info/...) returned a reCAPTCHA challenge; not bypassed. The NCBI-wide policy page above was read instead.
- No archive was opened, so the contents of the ccRCC and LSCC supplementary tarballs are UNVERIFIED.
- Gillette 2020's LUAD data freeze was not located; Cell supplementary terms were not read.
- SCAN-B: no study-specific licence was found beyond GEO; Brueffer et al.'s article licences were not read.
