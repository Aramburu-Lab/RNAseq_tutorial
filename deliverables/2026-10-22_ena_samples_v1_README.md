# 2026-10-22_ena_samples_v1

ENA's record of every sample in six studies, metadata only, for sorting samples by tissue,
stage and sex. Written by `scripts/python/ena_samples.py`; do not edit, rerun instead.

| Field | Value |
| --- | --- |
| Created | 2026-10-22, 14:56 local time |
| Command | `python -I scripts/python/ena_samples.py PRJNA955664 PRJNA393443 PRJNA483877 PRJEB49000 PRJNA792551 PRJNA344450 --out-dir deliverables --version 1` |
| Code | the commit that adds this file, on f19d576 |
| Renamed | 2026-10-22, by `git mv`, to carry each study's label (Rytis 2026-10-22); no byte changed, every SHA-256 below is the one measured at the first write |
| Python | 3.12, standard library only |
| Sources | ENA portal API `filereport?result=read_run` and ENA browser API `xml/<sample accessions>`, read 2026-10-22; ENA states no release number |
| Wall time, peak RAM | 27.4 s, 41.5 MB maximum resident set (`/usr/bin/time -l`) |

## Files

Each name carries the study's label, then its accession, lowercased: the deliverable
rule allows lowercase words joined by underscores only. A future run names files this way
with `--label <study>=<label>`, e.g. `--label PRJNA955664=lung_tissue_luac`. PRJNA344450's
two files keep their first names, untouched: kept, not used by the course.

| File | Label | Project | Rows | Bytes | SHA-256 |
| --- | --- | --- | --- | --- | --- |
| `2026-10-22_ena_samples_lung_tissue_luac_prjna955664_attributes_long_v1.tsv` | Lung tissue (LUAC) | PRJNA955664 | 3444 | 135753 | `ca0c7c41fa7e2757983c5521a516d9fa874c480f7e2f725014aa200cc17a622b` |
| `2026-10-22_ena_samples_lung_tissue_luac_prjna955664_v1.tsv` | Lung tissue (LUAC) | PRJNA955664 | 246 | 73504 | `fac249b99b226bdb20e6ba247c1596e721f78796b1226fa419f1803680e5c673` |
| `2026-10-22_ena_samples_breast_tissue_brca_prjna393443_attributes_long_v1.tsv` | Breast tissue (BRCA) | PRJNA393443 | 1648 | 56898 | `aedfcb49762caf884fac36a683614c239c19e2618ac0f18aa411b366ea70fa2b` |
| `2026-10-22_ena_samples_breast_tissue_brca_prjna393443_v1.tsv` | Breast tissue (BRCA) | PRJNA393443 | 103 | 24811 | `00aaaf576efc6185bf6b799b2aad74700fe70d0d27f5727ee4bd0b9175ef9fad` |
| `2026-10-22_ena_samples_brca_immune_cells_prjna483877_attributes_long_v1.tsv` | BRCA immune cells | PRJNA483877 | 1351 | 50663 | `d665f5d5c1c552516bb1ad8dd91857921469c21cbaedee4ff32504e04aca9905` |
| `2026-10-22_ena_samples_brca_immune_cells_prjna483877_v1.tsv` | BRCA immune cells | PRJNA483877 | 116 | 35856 | `db285cfd8a5d3ad3e9b67191fd3aebb0b84bc825888c678cf79cfd312138c05c` |
| `2026-10-22_ena_samples_noncancerous_blood_prjeb49000_attributes_long_v1.tsv` | Non-cancerous blood | PRJEB49000 | 4410 | 185728 | `085112e1c3675251361afe3e0f799004fc2d588f17c73de25200cc04ecbf8a97` |
| `2026-10-22_ena_samples_noncancerous_blood_prjeb49000_v1.tsv` | Non-cancerous blood | PRJEB49000 | 210 | 565774 | `c7852ec501763ea7769c7e7e29d797c7731f8fcda14595db009e4ca0eb6341d7` |
| `2026-10-22_ena_samples_pancreas_tissue_paad_prjna792551_attributes_long_v1.tsv` | Pancreas tissue (PAAD) | PRJNA792551 | 252 | 10583 | `9db302285143c9785ac98eaa6554cb6b8a21a1b1a0a9a6bea5069489ad7226ce` |
| `2026-10-22_ena_samples_pancreas_tissue_paad_prjna792551_v1.tsv` | Pancreas tissue (PAAD) | PRJNA792551 | 28 | 6837 | `4af6fd848468ff40c33bfb1ead356096f3d8b1fb4b74178f5e776e0719401970` |
| `2026-10-22_ena_samples_PRJNA344450_attributes_long_v1.tsv` | kept, not used by the course | PRJNA344450 | 918 | 30735 | `36cdd47884214c0ceac246d9d1587cae8dd84731edb044e95894e78ff2a97c56` |
| `2026-10-22_ena_samples_PRJNA344450_v1.tsv` | kept, not used by the course | PRJNA344450 | 54 | 13384 | `fe1ce36fe14f5f6468dfcb93fea89c359dc91273259b86c501c29255bf731c4d` |

## Data dictionary: `<date>_ena_samples_<label>_<study>_v1.tsv`, one row per sample

| Column | Type | Unit | Meaning |
| --- | --- | --- | --- |
| sample_accession | string | | BioSample accession, the XML's PRIMARY_ID |
| secondary_accession | string | | SRS or ERS accession, the XML's SECONDARY_ID |
| sample_alias | string | | submitter's alias, the XML's SUBMITTER_ID |
| sample_title | string | | the XML's TITLE |
| scientific_name | string | | SAMPLE_NAME/SCIENTIFIC_NAME |
| center_name | string | | the SAMPLE element's center_name attribute |
| broker_name | string | | the SAMPLE element's broker_name attribute |
| runs | string | | the sample's run accessions from the file report, joined by `;` |
| layout | string | | its distinct library layouts, joined by `;` |
| description | string | | the XML's DESCRIPTION, empty where absent |
| one column per tag | string | as submitted | the tag's value, verbatim; empty where the sample lacks the tag; a tag a sample carries twice is joined by `; ` |

Tag columns follow in the order the study first shows each tag. Tag names are ENA's raw XML
tags; ENA's browser title-cases them for display. A tag named as a fixed column heads its
column as `attribute_<tag>`: PRJEB49000's `scientific_name` is `attribute_scientific_name`.
PRJEB49000 samples each carry `organism`, `ENA-FIRST-PUBLIC` and `ENA-LAST-UPDATE` twice, and
their two date forms are joined in one cell.

## Data dictionary: `<date>_ena_samples_<label>_<study>_attributes_long_v1.tsv`

| Column | Type | Unit | Meaning |
| --- | --- | --- | --- |
| sample_accession | string | | BioSample accession |
| tag | string | | the attribute's TAG, as ENA writes it |
| value | string | as submitted | the attribute's VALUE, verbatim |

Tab-separated, UTF-8. The csv module quotes a value only when it holds a tab, a newline or a
double quote; no value in these twelve files holds one, so none is quoted.
