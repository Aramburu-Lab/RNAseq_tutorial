# NF.5: the six datasets renamed everywhere to their labels

Draft for Rytis, 2026-10-23. Decided: "literally everything" (Rytis 2026-10-23), the labels
without spaces and colons, and parentheses dropped from every file, folder and job name (Rytis
2026-10-23). Pages keep the labels as they read, "Breast tissue (BRCA)".

## The names

| Slug today | Name from NF.5 | Label, unchanged |
|---|---|---|
| `tier1` | `Non-cancerous-blood` | Non-cancerous blood |
| `tier1_teach` | `Non-cancerous-blood-teaching-subset` | Non-cancerous blood, teaching subset |
| `four_runs` | `Non-cancerous-blood-four-runs` | Non-cancerous blood, four runs |
| `breast_tissue` | `Breast-tissue-BRCA` | Breast tissue (BRCA) |
| `breast_mac` | `BRCA-immune-cells-macrophages` | BRCA immune cells: macrophages |
| `breast_mono` | `BRCA-immune-cells-monocytes` | BRCA immune cells: monocytes |
| `lung_gse229705` | `Lung-tissue-LUAD` | Lung tissue (LUAD) |
| `pancreas_gse192647` | `Pancreas-tissue-PAAD` | Pancreas tissue (PAAD) |

The two subsets are named from their own labels, since `tier1` is in both of their slugs.
Spaces become `-`; a comma or a colon is dropped. One rule makes the name from the label, the
same one `job_name.sh` already applies to job names, so a ninth dataset follows it unasked.

## The file-name rule changes

An output's name is `YYYY-MM-DD_<module>_<description>_vN`, and the description today is
lower-case words joined by `_` (`ctx.io.naming.DESCRIPTION`), the dataset's slug first. The new
names carry capitals and `-`, so a word becomes `[A-Za-z0-9]+(-[A-Za-z0-9]+)*`, words still
joined by `_`: `2026-10-23_de_Breast-tissue-BRCA_overview_v1.html`. A D-entry records the change;
CONVENTIONS' naming line and the naming lint follow it.

## What is renamed

- `conf/datasets/<slug>.toml`, `.tsv`, `.csv` by `git mv`, and every reference to them:
  `conf/cohorts/*.toml`, `conf/editions/core.toml`, the guide, the kit's dataset list.
- Every slug in `src/`, `scripts/`, `tests/` and their fixtures, and in the live docs:
  CONVENTIONS, AGENTS and TODO's open items.
- Step-1's frozen dataset file and its sanity record's path: a second change to a frozen file
  after the record, its own D-entry as D-112 was, with `test_frozen_since_registration`
  updated to the new names and nothing else.
- On this Mac, the ctx roots and their pages: `~/ctx_tier1` becomes `~/ctx_Non-cancerous-blood`
  and so on, by `mv`. Each root's pages are rebuilt under the new names from M1, and the pages
  with the old names move to the root's `superseded_2026-10-23/`, never deleted.
- On Kebnekaise, once NF.2's six runs end: each run folder `mv`'d to its new name, with a
  symlink left at the old path, because nf-core's and fetchngs' own outputs (the samplesheet's
  FASTQ paths, `pipeline_info`) hold absolute paths, and another step's output is not edited.
- The memory file of the names, and the release kit, rebuilt.

## What is not renamed, and why

Dated records keep the names they were written with: the 73 plans under `docs/plans`, the
D-entries in `knowledge/decisions.md`, `knowledge/lessons.md`, TODO's done items and the
deliverables already dated (`deliverables/2026-09-30_de_tier1_sanity_v1_README.md`). Rewriting
them would make a record say something that was not true when it was written. The D-entry holds
the old-to-new table, so each old name is found. If Rytis wants those rewritten too, it is one
more pass of the same substitution.

Accessions (PRJEB49000 and the rest), sample and run IDs, and GEO series stay as they are.

## How

1. One agent in a worktree does the repository part, cut from pushed main: the naming rule and
   its tests first, then the `git mv`s and the substitution, then a grep that finds no old slug
   outside the dated records, then whole-repo ruff and the full test suite. I review and merge;
   preflight before the push.
2. On the Mac, the roots: `mv`, rebuild, the old pages moved aside; the lite review set rebuilt.
3. On Kebnekaise, after NF.2: the `mv` and symlink commands, for Rytis to run.
4. The kit rebuilt from the renamed roots.

Checks: the naming lint passes on every output; a test plants an old slug in a page file name
and the lint refuses it; Step-1's sanity rules still pass on the renamed record; every page of
the review set opens in Chromium with no console error.

## Step 3's block (written 2026-10-25; run only once no run of the six is queued or running)

A `set -euo pipefail` script, run with `bash`. It refuses while any `RNAseq-Tutorial_` job is in
the queue, since a running head holds the old path open. Each folder moves to its new name and
the old path becomes a symlink to it. A run resumed later from either path reads the same files,
because every path nf-core and fetchngs recorded is absolute and still resolves. A name that
exists already is left alone, so a second run of the script changes nothing.

```
S=/proj/nobackup/hpc2n2024-150/Rytis/Tasks/RNAseq_tutorial
cat > "$S/rename_runs.sh" <<'SH'
set -euo pipefail
S=/proj/nobackup/hpc2n2024-150/Rytis/Tasks/RNAseq_tutorial
if squeue -u "$USER" -h -o '%j' | grep -q '^RNAseq-Tutorial_'; then echo "a run is queued; nothing renamed"; exit 1; fi
while read -r old new; do
    if [[ -L "$S/$old" || -e "$S/$new" ]]; then echo "$old: done already"; continue; fi
    [[ -d "$S/$old" ]] || { echo "$old: no such folder"; continue; }
    mv "$S/$old" "$S/$new" && ln -s "$new" "$S/$old" && echo "$old -> $new"
done <<'MAP'
tier1 Non-cancerous-blood
breast_tissue Breast-tissue-BRCA
breast_mac BRCA-immune-cells-macrophages
breast_mono BRCA-immune-cells-monocytes
lung_gse229705 Lung-tissue-LUAD
pancreas_gse192647 Pancreas-tissue-PAAD
MAP
ls -l "$S" | grep -- '->'
SH
bash "$S/rename_runs.sh"
```
