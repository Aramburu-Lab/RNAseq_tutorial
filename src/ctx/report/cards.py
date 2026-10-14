"""Parameter cards: what each CLI option means and what moving it changes, from cards.toml.

Kept apart from ``ctx.cli`` so the command file holds commands; imports nothing heavy, as
``ctx.cli`` must not (tests/test_runtime.py).
"""

from typing import Any

from ctx.report import text
from ctx.report.builder import ParameterCard
from ctx.report.marks import plain

# (meaning, what changes if you move it), per CLI option, from content/cards.toml; the default
# is read from click. The meaning doubles as the --help text, so it describes the option on
# every command; --help prints a glossary mark as its display text (D-048).
CARDS: dict[str, tuple[str, str]] = {
    option: (card["meaning"], card["effect"])
    for option, card in text.raw("cards", "options").items()
}
HELP = {option: plain(meaning) for option, (meaning, _) in CARDS.items()}
# Where one option does different things on different commands, the effect is per command.
CARD_EFFECTS: dict[tuple[str, str], str] = {
    (command, option): effect
    for command, effects in text.raw("cards", "effects").items()
    for option, effect in effects.items()
}
_UNCARDED = {"--help", "--version"}


def parameter_cards(*commands: Any) -> tuple[ParameterCard, ...]:
    """One card per CLI option of the given commands, with the default read from click.

    Args:
        commands: click commands (typer bundles its own click, so typed loosely): the group
            and the subcommand whose options a page documents.

    Returns:
        The cards, in the order click lists the options.

    Raises:
        ValueError: If an option has no entry in ``CARDS``; a parameter without a card
            is a teaching page that is silently incomplete.
    """
    options = [p for c in commands for p in c.params if p.opts and p.opts[0] not in _UNCARDED]
    if missing := [p.opts[0] for p in options if p.opts[0] not in CARDS]:
        raise ValueError(f"no parameter card for {missing}")
    command = commands[-1].name if commands else ""
    cards = []
    for option in options:
        meaning, effect = CARDS[option.opts[0]]
        effect = CARD_EFFECTS.get((command, option.opts[0]), effect)
        if option.required:
            default = "required"
        elif option.secondary_opts:  # an on/off switch, such as --lfc-shrink/--no-lfc-shrink
            default = "on" if option.default else "off"
        else:
            default = "unset" if option.default is None else str(option.default)
        cards.append(ParameterCard(option.opts[0], meaning, default, effect))
    return tuple(cards)
