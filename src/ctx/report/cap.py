"""The page cap: one definition of the limit, and of the order a page's optional parts leave.

CONVENTIONS.md, HTML reports: a page is self-contained and under 15 MB, and ``builder.write``
refuses one over ``MAX_BYTES``. A page with optional parts (M2's flag view, M9's cards after the
first) keeps as many of them, in list order, as its page as written holds: each candidate is
rendered, Reproduce panels and all, and measured; parts leave from the end of the list. Y.1
found two pages over the cap with every test green: M9 had fitted its cards before their panels
existed, and M2's view was never dropped, as its test named a stem no real map has.

Every reader takes ``MAX_BYTES`` from this module when it runs, so a test that lowers it here
lowers it for the writer and for every fit alike.
"""

from collections.abc import Callable, Sequence

MAX_BYTES = 15 * 10**6


def over(html: str) -> int:
    """How many bytes a page as written passes the cap; zero or less where it fits."""
    return len(html.encode()) - MAX_BYTES


def kept(pages: Callable[[int], tuple[str, Sequence[int]]], parts: int, least: int = 0) -> int:
    """How many optional parts, in list order, the page as written holds under the cap.

    ``pages(k)`` renders the page as it would be written with the first ``k`` parts and gives
    its HTML and the bytes of each part from ``least`` on, as far as it measures them. Parts
    leave from the end, as many at once as their bytes take to cover the excess, one at a time
    where their bytes are not given, and the page is measured again, since what replaces them (a
    note, the names of the cards left out) has bytes too. The first ``least`` parts always
    stay; a page still over the cap with them is left for ``builder.write`` to refuse.
    """
    count = parts
    while count > least:
        html, sizes = pages(count)
        excess = over(html)
        if excess <= 0:
            break
        freed = 0
        while count > least and freed < excess:
            count -= 1
            freed += sizes[count - least] if count - least < len(sizes) else excess
    return count
