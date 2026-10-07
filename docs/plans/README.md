# Plans approved by Rytis

Each file is the plan or proposal as it stood when Rytis approved it. TODO.md and
knowledge/decisions.md record what was decided; these hold the detail the records cite.

| File | Approved | Record |
|---|---|---|
| [2026-09-30_m3_plan.md](2026-09-30_m3_plan.md) | M3 plan, answered 2026-09-30 and 2026-10-01 | TODO.md M3.0 to M3.14 |
| [2026-09-30_m3_facts.md](2026-09-30_m3_facts.md) | the measurements behind it | the same |
| [2026-09-30_gene_selection_plan.md](2026-09-30_gene_selection_plan.md) | 2026-09-30 | D-049, TODO.md GS.0 to GS.10 |
| [2026-10-01_qc_bands_proposal_approved.md](2026-10-01_qc_bands_proposal_approved.md) | 2026-10-01, all 16 recommended answers | D-052, TODO.md P2Q.9 and P2Q.10 |
| [2026-10-01_qc_bands_proposal.json](2026-10-01_qc_bands_proposal.json) | the registry fields of that proposal, one object per metric | the same |
| [2026-10-01_figure_controls_plan_approved.md](2026-10-01_figure_controls_plan_approved.md) | 2026-10-01, all 12 recommended answers | D-051, TODO.md P2S.6 and the work order |
| [2026-10-01_figure_controls_critique.md](2026-10-01_figure_controls_critique.md) | the critique that plan was written from | the same |

## How these were recovered

They were first kept in the session's scratch folder under `/private/tmp`, which macOS's
periodic cleanup emptied of files unused for three days; on 2026-10-07 all seven were gone.
They were rebuilt from the agent transcripts under `~/.claude/projects`:

- the bands proposal from its write, with the `&gt;` to `>` replacement the agent made
  after it, matching the agent's own printout of the final file;
- `proposal.json` by running the agent's own `build_json.py`, recovered from the command that
  wrote it and never edited after;
- the controls plan from its write plus the four-phrase edit script the agent ran after it,
  giving the original 49,842 bytes;
- the critique and the gene-selection plan from their writes, with one later edit replayed;
- the M3 plan and facts from an agent's full read of each on 2026-10-01 12:38 UTC.

One gap: the M3.6 review corrected one sentence of the M3 plan's section 4 on 2026-10-02 by
a command the transcripts do not hold. The sentence was rewritten here to the correction that
review recorded in decisions.md D-039: Salmon's chr22 share is context, not rule 7's ceiling.
