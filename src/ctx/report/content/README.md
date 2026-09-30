# Page text

Every sentence a reader sees on a page is in one of these files. Edit the text here; the
code only names a key and fills in the numbers. Rebuild a page by running its command again.

| To change… | Edit | How |
|---|---|---|
| Headings, buttons, labels and table heads on every page | `ui.toml` | edit the strings; `[sections]` holds the eight skeleton headings |
| What a command-line option means, and what moving it changes | `cards.toml` | `[options."--name"]`; `[effects.<command>]` overrides the effect on one command |
| The M0 Setup page | `setup.toml` | `purpose`, `input`, `method`, `interpretation`, `pitfalls`, `[[terms]]`, `[[questions]]`, then one table per results block |
| The M1 Gather page | `gather.toml` | the same skeleton, then the checks, curation, counts and the two figures |

Not here yet: the M5 and M6 pages, the hub and the figure captions. Their code is frozen
until the Step-1 record, and their text moves after it (TODO.md, P2S.4).

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
- **No em dashes or arrows** anywhere, as on the rest of the pages.

## An edited copy, without touching the package

Copy this folder, edit the copy, and point `CTX_CONTENT` at it:

```bash
CTX_CONTENT=/path/to/my_content ctx setup
```

With `CTX_CONTENT` set, every file is read from that folder and none from the package, so
the copy must hold all of them.
