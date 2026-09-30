"""M2: read QC and the checkpoint on counts (`ctx qc`)."""

from dataclasses import dataclass


class QcError(ValueError):
    """A QC input that cannot be used; the message names every problem found."""


@dataclass(frozen=True)
class Fastp:
    """fastp's read filter; the defaults are fastp 1.1.0's own (`fastp --help`)."""

    q: int = 15
    u: int = 40
    l: int = 15  # noqa: E741 - fastp's own name for the length filter
    cut_tail: bool = False
    window: int = 4
    mean: int = 20

    def args(self) -> list[str]:
        """The settings as fastp options."""
        args = ["-q", str(self.q), "-u", str(self.u), "-l", str(self.l)]
        return args + (["--cut_tail", "-W", str(self.window), "-M", str(self.mean)]
                       if self.cut_tail else [])  # fmt: skip
