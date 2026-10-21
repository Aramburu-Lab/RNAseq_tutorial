# M9 gene card fixtures: NPTX2's real responses

Bodies exactly as received by `ctx genecard` on 2026-10-08, gzipped; SHA-256 of each body
before gzip, as `resources/cache/manifest.json` recorded it.

| File | URL | Retrieved (UTC) | SHA-256 of the body |
|---|---|---|---|
| `hpa_v25_ENSG00000106236.xml.gz` | https://v25.proteinatlas.org/ENSG00000106236.xml | 2026-10-08T12:25:43 | e6e9db0992d9ecdad81e406b62fce7f5747263dfcd7fd94194873e8519b669fb |
| `hpa_v25_ENSG00000106236.json.gz` | https://v25.proteinatlas.org/ENSG00000106236.json | 2026-10-08T12:25:43 | d4d20a7e8dd4237d62ebb3fdecc0e63212eb2fcacec3207d80dd479a2d6743c6 |
| `ensembl_116_ENSG00000106236_lookup.json.gz` | https://rest.ensembl.org/lookup/id/ENSG00000106236?expand=1;mane=1;content-type=application/json | 2026-10-08T12:25:46 | fc0e8c0a16dc28f39e0620b43d6ecf3097d00722fa7eb43d6fa0a2fd54fd42ba |
| `ensembl_116_info_data.json.gz` | https://rest.ensembl.org/info/data?content-type=application/json | 2026-10-08T12:25:42 | 86e87c46ae1a424702676493d024d0508ce60f6c951b0b2803ce10269e09c20c |

Two HPA search responses, fetched with curl under `ctx genecard`'s user agent on 2026-10-09
(v25 host, redirected to www.proteinatlas.org; the search names no release), at the URLs
`hpa.search_url` builds:

| File | URL | Retrieved (UTC) | SHA-256 of the body |
|---|---|---|---|
| `hpa_v25_search_ENSG00000106236_cell_RNA_leukemia.json.gz` | https://v25.proteinatlas.org/api/search_download.php?search=ENSG00000106236&format=json&columns=g,eg,cell_RNA_leukemia&compress=no | 2026-10-09T07:11:26 | 9344fd3a0e73352c3931cccaadd51e5bc553934b5ba87bea1ff47c0106b13668 |
| `hpa_v25_search_name_NDUFA4L2.json.gz` | https://v25.proteinatlas.org/api/search_download.php?search=NDUFA4L2&format=json&columns=g,gs,eg&compress=no | 2026-10-09T07:11:36 | 630c72b6a010a5ad426e5d2b4fe66c629eecbeee37da0e48101253e81173a2d8 |

GTEx API v2 responses for NPTX2, fetched by `ctx genecard` on 2026-10-09 and copied from its
cache unchanged (a gzipped body, mtime 0):

| File | URL | Retrieved (UTC) | SHA-256 of the body |
|---|---|---|---|
| `gtex_v10_ENSG00000106236_reference_v39.json.gz` | https://gtexportal.org/api/v2/reference/gene?geneId=ENSG00000106236&gencodeVersion=v39&genomeBuild=GRCh38%2Fhg38 | 2026-10-09T07:28:49 | b5b33c3c2a0feb9053ed1a562e1cf27456ae4058fb9788533ec70782b92da76f |
| `gtex_v10_ENSG00000106236_junctions.json.gz` | https://gtexportal.org/api/v2/expression/medianJunctionExpression?gencodeId=ENSG00000106236.4&datasetId=gtex_v10&itemsPerPage=100000 | 2026-10-09T07:28:49 | 28ebdea9ff452c8ee8f5275ab923391a019ef11dc805400b37af59eb30eef9f4 |

One GTEx file here is made up, not a response: `gtex_v10_ENSG00000106236_median_tpm_madeup.json.gz`
(E.1c, D-106), NPTX2's V10 median TPM per tissue in the shape the GTEx API v2 reference gives
for `expression/medianGeneExpression` (no live response was read). Its 54 tissues are the
junction response's, in its order. Pituitary 240, Uterus 3 and Whole_Blood 0.05 are the V11
fixture's all-donor medians (`tests/genecard_root.py`, `NPTX2_TPM`); EBV-transformed
lymphocytes is set to 0.08; every other tissue is 0.6 times its mean of the four junction
medians, to 2 decimals. Seven tissues fall below 1 TPM. SHA-256 of the body before gzip:
e49dd73b656b990c48ea291366b53268ec4fc0336f22b78884f0740ccaf1504b. Written by a script
kept outside the repository; the values carry no meaning beyond the tests.

The two HPA search bodies above are byte-identical (same SHA-256) to the ones `ctx genecard`
cached at 2026-10-09T07:28:37 and 07:27:03.

Licences: the Human Protein Atlas data are CC BY 4.0 (Uhlén et al., Science 2015,
doi:10.1126/science.1260419; https://v25.proteinatlas.org/ENSG00000106236). Ensembl REST
under the EMBL-EBI terms of use (Ensembl 2026, doi:10.1093/nar/gkaf1239).

GTEx V11's files are not here: the tests build made-up GCT, parquet and annotation files of
the same shape in `tests/genecard_root.py`. The two GTEx V10 API responses above are the
GTEx Portal's, obtained on 2026-10-09, used under its licence with that acknowledgement.
