# ENA sample fixtures for `scripts/python/ena_samples.py`

FIXTURES, not real data. Shaped on two real responses read on 2026-10-22: ENA's portal API
`filereport?result=read_run&format=tsv` with `scripts/python/ena_samples.py`'s seven fields, and
ENA's browser API `xml/<BioSample accession>` for one human tumour sample of an NCBI-brokered
study. Every accession, alias, lab name and value is invented; the `000000` accessions resolve
to nothing at ENA.

| File | Shape of | Planted for the tests |
|---|---|---|
| `read_run.tsv` | read_run file report | sample `SAMN00000001` has two runs |
| `samples.xml` | browser-API `SAMPLE_SET` of three samples | `SAMN00000002` lacks `recurrence` and has no `stage`; `SAMN00000003` lacks `Sex`, and its `stage` value holds a tab |
