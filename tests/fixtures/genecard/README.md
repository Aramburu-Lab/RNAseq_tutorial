# M9 gene card fixtures: NPTX2's real responses

Bodies exactly as received by `ctx genecard` on 2026-10-08, gzipped; SHA-256 of each body
before gzip, as `resources/cache/manifest.json` recorded it.

| File | URL | Retrieved (UTC) | SHA-256 of the body |
|---|---|---|---|
| `hpa_v25_ENSG00000106236.xml.gz` | https://v25.proteinatlas.org/ENSG00000106236.xml | 2026-10-08T12:25:43 | e6e9db0992d9ecdad81e406b62fce7f5747263dfcd7fd94194873e8519b669fb |
| `hpa_v25_ENSG00000106236.json.gz` | https://v25.proteinatlas.org/ENSG00000106236.json | 2026-10-08T12:25:43 | d4d20a7e8dd4237d62ebb3fdecc0e63212eb2fcacec3207d80dd479a2d6743c6 |
| `ensembl_116_ENSG00000106236_lookup.json.gz` | https://rest.ensembl.org/lookup/id/ENSG00000106236?expand=1;mane=1;content-type=application/json | 2026-10-08T12:25:46 | fc0e8c0a16dc28f39e0620b43d6ecf3097d00722fa7eb43d6fa0a2fd54fd42ba |
| `ensembl_116_info_data.json.gz` | https://rest.ensembl.org/info/data?content-type=application/json | 2026-10-08T12:25:42 | 86e87c46ae1a424702676493d024d0508ce60f6c951b0b2803ce10269e09c20c |

Licences: the Human Protein Atlas data are CC BY 4.0 (Uhlén et al., Science 2015,
doi:10.1126/science.1260419; https://v25.proteinatlas.org/ENSG00000106236). Ensembl REST
under the EMBL-EBI terms of use (Ensembl 2026, doi:10.1093/nar/gkaf1239).

GTEx is not here: the tests build made-up GCT, parquet and annotation files of the same
shape in `tests/genecard_root.py`.
