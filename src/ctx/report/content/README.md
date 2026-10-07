# Page text

Every sentence a reader sees on a page is in one of these files. Edit the text here; the
code only names a key and fills in the numbers. Rebuild a page by running its command again.

| To change… | Edit | How |
|---|---|---|
| Headings, buttons, labels and table heads on every page | `ui.toml` | edit the strings; `[sections]` holds the eight skeleton headings |
| What a command-line option means, and what moving it changes | `cards.toml` | `[options."--name"]`; `[effects.<command>]` overrides the effect on one command |
| The M0 Setup page | `setup.toml` | `purpose`, `input`, `method`, `interpretation`, `pitfalls`, `[[terms]]`, `[[questions]]`, then one table per results block |
| The M1 Gather page | `gather.toml` | the same skeleton, then the checks, curation, counts and the two figures |
| The M3 Trim, align, count page | `process.toml` | the same skeleton, then the notes, groups, key numbers, outcomes, settings, comparisons, each figure's words and each table's columns; `[preview]` holds what the course preview alone says |
| The M5 Differential expression page, the sentences moved so far | `de.toml` | the volcano's adjustment sentence, and what a thinned page calls the points it may leave out |
| The M6 Function page, the sentences moved so far | `function.toml` | the ORA lead, caption and passing line; the GSEA highlight line; the background lesson's words in each direction, and its self-checks |

Not here yet: the rest of the M5 and M6 pages, the hub and the other figure captions. Their
text moves with P2S.4 (TODO.md).

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
