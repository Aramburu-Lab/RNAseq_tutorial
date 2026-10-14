# Page text

Every sentence a reader sees on a page is in one of these files. Edit the text here; the
code only names a key and fills in the numbers. Rebuild a page by running its command again.

| To change… | Edit | How |
|---|---|---|
| Headings, buttons, labels and table heads on every page | `ui.toml` | edit the strings; `[sections]` holds the eight skeleton headings |
| What a command-line option means, and what moving it changes | `cards.toml` | `[options."--name"]`; `[effects.<command>]` overrides the effect on one command |
| The M0 Setup page | `setup.toml` | `purpose`, `input`, `method`, `interpretation`, `pitfalls`, `[[terms]]`, `[[questions]]`, then one table per results block |
| The M1 Gather page | `gather.toml` | the same skeleton, then the checks, curation, counts and the two figures |
| The custom-data guide (`ctx guide`) | `guide.toml` | the same skeleton, then the key numbers, `[find]` (the accession levels, routes and class exercise before Part 1), each part's heading, lead and block titles; the commands are not here: they come from the launchers' Usage blocks and the CLI |
| The teacher guide (`python -m ctx.report.teacher`) | `teacher.toml` | the same skeleton, its two argument cards, then the time table's words; every time is read from a module page |
| The M1 Gather page of a GDC project (`--gdc-dir`) | `gather_gdc.toml` | the same skeleton, then the checks, the one-file-per-patient flow, the annotations and each table's columns |
| The M3 Trim, align, count page | `process.toml` | the same skeleton, then the notes, groups, key numbers, outcomes, settings, comparisons, each figure's words and each table's columns; `[preview]` holds what the course preview alone says |
| The M5 Differential expression page | `de.toml` | the same skeleton, `[questions.*]` asked by what the design holds, how an unlabelled contrast reads, then the groups, the contrasts table, the key numbers and each figure's words; the volcano's reading lines are filled live in the browser |
| The M6 Function page | `function.toml` | the same skeleton, `[questions.*]` by the background's direction, a gene set in words (`[gene_sets]`), then the groups, the key numbers and each figure's words |
| The course hub (`ctx report`) | `hub.toml` | its title and lead, each module's name, the key numbers |
| The statistics every figure's caption shares | `annotate.toml` | the test names, the box figure's caption and labels, PERMANOVA's sentence |

A contrast reads by the `label` its dataset file gives it; without one, by its two levels as
the SDRF spells them.

## Rules the tests hold

- **Placeholders.** `{name}` is filled by the code. Keep the name; `{name:,}` and
  `{name:.1f}` format a number. A placeholder the code does not fill stops the page with an
  error naming the file and key, so no page ships with a hole in it.
- **Variants.** A sentence that changes with the data (one run or several, one batch or
  more) is two keys, and the code picks one: `lead_one` and `lead_many`.
- **No dead text.** Every key must be read somewhere; a key no code names fails
  `tests/test_text.py`.
- **Pitfalls.** In `pitfalls`, the words before the first `": "` print in bold.
- **Prose rules.** CONVENTIONS.md, Report prose: word budgets, sentence length, the banned list, numbers from `{placeholders}`. `pixi run lint-prose` lists every violation.
- **Glossary marks.** `[[term]]` marks a glossary term, `[[term|shown words]]` shows other words, and `[[term!]]` gives a later mention its icon too. The term is an id, term or alias in `../glossary.yaml`, in any case. Marks work in the sections, questions, pitfalls, figure captions and cards; the first mention per section gets the icon. A mark inside backticks is code, shown as written. `python -m ctx.report.glossary_lint` refuses an unknown term or a mark anywhere else, and lists glossary terms left unmarked or marked only after a bare first mention (D-048).
- **No em dashes or arrows** anywhere, as on the rest of the pages.

## An edited copy, without touching the package

Copy this folder, edit the copy, and point `CTX_CONTENT` at it:

```bash
CTX_CONTENT=/path/to/my_content ctx setup
```

With `CTX_CONTENT` set, every file is read from that folder and none from the package, so
the copy must hold all of them.
