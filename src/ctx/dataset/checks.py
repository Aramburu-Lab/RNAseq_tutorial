"""What every table parser of a dataset file shares: the error, the typed read, filters."""

from collections.abc import Mapping
from typing import Any

# An own table's endings: comma-separated, or tab-separated as an SDRF is (ctx.gather.curation).
TABLE_SUFFIXES = (".csv", ".tsv", ".txt")


class DatasetError(ValueError):
    """A dataset file that cannot be used; the message lists every problem."""


def _need(table: Mapping[str, Any], key: str, kind: type, where: str, problems: list[str]) -> Any:
    """Read one required key of a given type, recording a problem instead of raising."""
    value = table.get(key)
    # bool is an int to isinstance; a count written as true is still a mistake.
    if not isinstance(value, kind) or (isinstance(value, bool) and kind is not bool):
        problems.append(f"{where}.{key}: expected {kind.__name__}, found {value!r}")
        return None
    return value


def _table(value: Any, where: str, problems: list[str]) -> Mapping[str, Any]:
    """A value that must be a TOML table; anything else is a problem and reads as empty."""
    if isinstance(value, Mapping):
        return value
    problems.append(f"{where}: expected a table, found {value!r}")
    return {}


def _list(value: Any, where: str, problems: list[str]) -> list[Any]:
    """A value that must be a TOML array; anything else is a problem and reads as empty."""
    if isinstance(value, list):
        return value
    problems.append(f"{where}: expected a list, found {value!r}")
    return []


def _values(raw: Any, where: str, problems: list[str]) -> dict[str, tuple[str, ...]]:
    """A filter table: column to allowed values, each a list."""
    table = _table(raw, where, problems)
    return {k: tuple(str(x) for x in _list(v, f"{where}.{k}", problems)) for k, v in table.items()}
