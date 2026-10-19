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


def _values(raw: Mapping[str, Any]) -> dict[str, tuple[str, ...]]:
    """A filter table: column to allowed values."""
    return {k: tuple(str(x) for x in v) for k, v in raw.items()}
